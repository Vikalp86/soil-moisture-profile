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

Authors:  Vikalp Mishra, W. Lee Ellenburg, Kel Markert
Contact: vmishra@nsstc.uah.edu
Copyright (c) 2019, Authors

::::::::::::::::::::::::::::::::::::::::::::::::::::
Created on Mon Feb 26 10:12:12 2018
-----------------------
A code to develop dry case SM profile
using the Principle of Maximum Entropy Model - POME

The model developed by Al-Hamdan & Cruise, 2010
Original Code by - Al_Hamdan,UAH
Citiation: http://ascelibrary.org/doi/abs/10.1061/(ASCE)HE.1943-5584.0000196
----------
Modified by- Vikalp Mishra
NASA-SERVIR

Change: instead of while loop for computing the min of
mass balance, this code uses the for loop to select the
min of mass balance error -  the best case scenario
starting from the field capacity.[Jul 28, 2018]
-----------------------------------------------------
:::::::::::::::::::::::::::::::::::::::::::::::::::::
"""


import numpy as np
import pandas as pd


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------


def rmse(sim,obs):
    sim = np.array(sim)
    obs = np.array(obs)
    return np.sqrt(((sim-obs)**2).mean())


def bias(sim,obs):
    sim = np.array(sim)
    obs = np.array(obs)
    return np.mean(sim-obs)

def solver(sur,avg,bott):
    temp = 0.
    temp2 = 0.

    if (sur+bott)> 2*avg:
        dat1 = -50
        dat2 = -0.01
    else:
        dat1 = 0.01
        dat2 = 50

    while abs(dat2-dat1)>0.0005:
        f1 = (math.exp(dat1*bott)-math.exp(dat1*sur))/(bott*math.exp(dat1*bott)-sur*math.exp(dat1*sur)-avg*(math.exp(dat1*bott)-math.exp(dat1*sur))) - dat1
        f2 = (math.exp(dat2*bott)-math.exp(dat2*sur))/(bott*math.exp(dat2*bott)-sur*math.exp(dat2*sur)-avg*(math.exp(dat2*bott)-math.exp(dat2*sur))) - dat2

        if (f1*f2)<=0:
            temp = dat1
            dat1 = (dat1+dat2)/2.
        else:
            temp2 = dat2
            if temp == 0:
                dat2 = dat1
            else:
                dat2 = temp
            dat1 = (dat1+temp2)/2.

    return dat1

def wet_case(sur,avg,bott,D):
    depths = np.arange(5,int(D)+5,5)

    # check to see if avg is within limits -- in not then provide fix values
    if avg < sur or avg > bott:
        avg = sur*0.65 + bott*0.35

    dat = solver(sur,avg,bott)

    try:
        x = np.log(dat/(np.exp(dat*sur)-np.exp(dat*bott)))+1
    except Exception:
        x = None

    theta = np.log(np.exp(dat*sur)-(dat*np.exp(1.-x)*(depths-5.)/(D-5.)))/dat

    return theta

def dry_case(sur,avg,bott,D):
    depths = np.arange(5,int(D)+5,5)

    # check to see if avg is within limits -- in not then provide fix values
    if avg < sur or avg > bott:
        avg = sur*0.65 + bott*0.35

    dat = solver(sur,avg,bott)

    try:
        x = np.log(dat/(np.exp(dat*bott)-np.exp(dat*sur)))+1
    except Exception:
        x = None

    theta = np.log((np.exp(dat*sur)+(dat*np.exp(1.-x)*(depths-5.)/(D-5.))))/dat

    return theta


def run_case(sur,avg,bott,inf_val,D,inflection,avg_out=False):

    # assumes the average SM is uniformly distributed throughout the column
    upp_avg = avg * (D/inflection)
    low_avg = avg * (1- (D/inflection))

    low_sur = upp_bott = inf_val

    if sur < avg and avg > bott:
        # case I
        upp = dry_case(sur,upp_avg,upp_bott,inflection)
        low = wet_case(low_sur,low_avg,bott,D-inflection+5)

    elif sur > avg and avg < bott:
        # case II
        upp = wet_case(sur,upp_avg,upp_bott,inflection)
        low = dry_case(low_sur,low_avg,bott,D-inflection+5)
    else:
        raise ValueError()

    col = np.concatenate([upp,low[1:]])

    if avg_out:
        out = np.nanmean(col)
    else:
        out = col

    return out

# vectorize the function to run dynamic case
vcases= np.vectorize(run_case)


def dynamic_case(sur,avg,bott,D,inflection,fc,resid,sat):
    Inf_val = (fc - resid) / (sat - resid)
    if avg > Inf_val: Inf_val = 1
    if sur < avg and avg > bott:
        jj = np.arange(np.max([sur,bott])+0.001,Inf_val+0.001,0.001)
    elif sur > avg and avg < bott:
        jj = np.arange(Inf_val,np.min([sur,bott]),0.001)

    trys = vcases(sur,avg,bott,jj,D,inflection,avg_out=True)
    errors = (np.abs(trys-avg)*100) / avg
    minimalInf = jj[np.nanargmin(errors)]

    data = run_case(sur,minimalInf,bott,minimalInf,D,inflection)

    return data

def lookup_soilType(soilType):
    usda = {
    1:'blah',
    }
    return usda[soilType]


def pome(sur,avg,bott,D,inflection=None,fc=None,resid=None,sat=None,soilType=None):
    """
    Function to

    soilType = USDA soiltype
    
    sur, avg & bott are SM values in effective soil moisture
    D           - total layer depth in cm (~100)
    inflection  - depth in cm at which inflection point is observed 
    fc          - field capacity
    resid       - residual soil moisture
    sat         - satuarated SM 
    soilType    - USDA soil classification (to calculate soil characteristics if fc, resid & sat values are not know)

    """

    if sur < avg and avg < bott:
        # dry case
        profile = dry_case(sur,avg,bott,D)

    elif sur > avg and avg > bott:
        # wet case
        profile = wet_case(sur,avg,bott,D)


    elif ((sur < avg) and (avg > bott)) or ((sur > avg) and (avg < bott)):
        if np.any([fc,resid,sat]==None):
            if soilType == None:
                fc,resid,sat = lookup_soilType(soilType)
            else:
                raise ValueError('When')

        profile = dynamic_case(sur,avg,bott,D,inflection,fc,resid,sat)

    else:
        raise RuntimeError('Counld not determine soil moisture profile please check input data')

    return profile

# vectorize the pome function to accept inputs as numpy arrays
vpome = np.vectorize(pome,otypes=[np.ndarray])
