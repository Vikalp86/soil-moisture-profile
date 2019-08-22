"""
Licensed under the MIT License (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

https://opensource.org/licenses/MIT

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Authors:  Kel Markert, Vikalp Mishra, W. Lee Ellenburg
Contact: kel.markert@uah.edu
Copyright (c) 2019, Authors

"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pylab as plt
import matplotlib.dates as mdates
import seaborn as sns

def nse(df):
    """
    Nash-sutcliffe model efficiency coefficient function

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Returns: out (float): calculated NSE value
    """

    df = df.dropna()
    sim = df['simulation']
    obs = df['observation']
    out = 1 - (np.nansum((sim-obs)**2)/np.nansum((obs-np.nanmean(obs))**2))
    return out


def rmse(df):
    """
    Root mean square error function

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Returns: out (float): calculated RMSE value
    """

    df = df.dropna()
    sim = df['simulation']
    obs = df['observation']
    out = np.sqrt(np.nanmean((sim-obs)**2))
    return out

def r(df):
    """
    Correlation coefficient function

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Returns: out (float): calculated correlation coefficient value
    """

    df = df.dropna()
    sim = df['simulation']
    obs = df['observation']
    out = stats.pearsonr(sim,obs)[0]
    return out

def bias(df):
    """
    Bias function

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Returns: out (float): calculated bias value
    """

    df = df.dropna()
    sim = df['simulation']
    obs = df['observation']
    out = np.nanmean(sim-obs)
    return out

def ubRmse(df):
    """
    Un-biased RMSE function

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Returns: out (float): calculated un-biased rmse value
    """

    rms = rmse(df)
    b = bias(df)
    out = np.round(np.sqrt(rms**2-b**2),3)
    return out


def calc_Topt(sur,obs,Tvals,objfunc='nse'):
    """
    Function to calibrate the T parameter using a brute-force method

    Args: sur (pandas.Series): pandas series of the surface soil moisture
          obs (pandas.Series): pandas series of the soil moisture at layer x to calibrate
          Tvals (list,tuple,set,np.array): sequence of values to test for optimal value

    Kwargs: objfuc (string): objective function used to search for optimal value;
                             options: "nse","rmse","bias",and "r"; default: "nse"

    Returns: out (dict): dictionary with the optimal T value keyed at 'T' and the
                          objective function value keyed at 'objval'
    """

    objOpts = dict(nse=nse,rmse=rmse,bias=bias,r=r,ubrmse=ubRmse)
    objectiveFunc = objOpts[objfunc]

    df = pd.concat([sur,obs],axis=1)
    df.columns = ('surface','depth')
    df.dropna(inplace = True)
    # new_df = new_df[~new_df.index.duplicated(keep='first')]

    results = []
    for T in Tvals:
        Ttest = expFilter(df['surface'],T=T)
        tempDf = pd.concat([Ttest,df['depth']],axis=1)
        tempDf.columns = ('simulation','observation')
        N = objectiveFunc(tempDf)
        results.append(N)

    # check to either find the min or max depending on objectivFunc
    if objfunc in ('nse','r'):
        best = np.array(results).argmax()
        objVal = np.nanmax(results)
    else:
        best = np.array(results).argmin()
        objVal = np.nanmin(results)

    out = dict(T=Tvals[best],objval=objVal)
    return out

def expFilter(series,T=1):
    """
    Function to calculate the exponential filter function for a time series
    Expects the soil moisture to be effective soil moisture

    Args: series (pandas.Series): pandas series of the surface observations

    Kwargs: T (int): integer value for T parameter used in the exponential filter

    Returns: out (pandas.Series): series of simulated soil moisture
    """

    sim = series.copy()

    K = 1
    SWIn = series.iloc[0]

    sim = []
    for i in range(series.size):
        if i == 0:
            gap = 1
        else:
            dd = series.index[i] - series.index[i-1]
            gap = dd // np.timedelta64(1, 'D')

        if gap <= 12 :
            p = series.iloc[i]
            recurSWI = SWIn + K * (p-SWIn)
            K = K / (K+np.exp(-1*gap/T))

        if gap > 12:
            p = series.iloc[i]
            K = 1
            SWIn = series.iloc[i-1]
            recurSWI = SWIn + K * (p-SWIn)

        SWIn = recurSWI
        sim.append(recurSWI)

    out = pd.Series(sim,index=series.index,name='simulation')
    return out

def normalize(series,lower=None,upper=None):
    """
    Function to normalize values between a lower and upper bounds

    Args: series (np.array | pd.Series | pd.DataFrame): continuous values to normalize

    Kwargs: lower (float): lower bound to normalize to; default: None and will use
                           the series minimum value
            upper (float): upper bound to normalize to; default: None and will use
                           the series maximum value

    Returns: normalized series with values 0-1
    """

    if lower == None:
        lower = np.nanmin(series)

    if upper == None:
        upper = np.nanmax(series)

    return (series-lower)/(upper-lower)


def make_plot(df,plotType='scatter',outFile=None,**kwargs):
    """
    Function to create a plot of observed and simulated time series

    Args: df (pandas.DataFrame): pandas dataframe with 'simulation' and 'observation' columns

    Kwargs: plotType (string): type of plot to create; options: "scatter", "series";
                                default: scatter
            outFile (string): output file path to save figure to; if outFile is
                              specified then the plot will not show; default: None
            **kwargs (dict): keyword arguments for saving matplotlib figure

    Returns: None
    """

    # Compute Time series statistics
    ns = np.round(nse(df),2)   # NS efficiency
    rms = np.round(rmse(df),3) # RMSE
    r2 = np.round(r(df)**2,3)  # R^2
    b = np.round(bias(df),3)   # Bias
    ub = np.round(np.sqrt(rms**2-b**2),3)

    x,y = df['observation'].values,df['simulation'].values

    if plotType == "scatter":
        sns.regplot(x,y, ci = None, truncate = False)
        plt.ylim(0,1)
        plt.xlim(0,1)
        plt.xlabel('SCAN')
        plt.ylabel('Exponential Filter')
        xTxt = 0.82
        yTxt = 0.05

    elif plotType == 'series':
        plt.plot(y,ls = '--', color = 'b',label = 'simulation')
        plt.plot(x, color = 'r',label = 'observation')
        plt.ylim(0,1)

        plt.ylabel('SM (Effective)')
        plt.xlabel('Year')
        legend = plt.legend()
        xTxt = 0.82
        yTxt = 0.80

    else:
        raise NotImplementedError()

    plt.text(xTxt,yTxt+(0.05*3),'R: '+str(np.sqrt(r2)))
    plt.text(xTxt,yTxt+(0.05*2),'RMSE: '+str(rms))
    plt.text(xTxt,yTxt+(0.05),'ubRMSE: '+str(ub))
    plt.text(xTxt,yTxt,'Bias: '+str(b))

    # plt.title('Layer ' + str(lyr+2))
    if outFile != None:
        plt.savefig(outFile,**kwargs)
        plt.close()
    else:
        plt.show()

    return
