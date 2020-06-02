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

Authors:  Vikalp Mishra
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

-----------------------------------------------------
:::::::::::::::::::::::::::::::::::::::::::::::::::::
"""


import math, timeit, os, sys
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
st1 = timeit.default_timer()



def eqs_wet(sur,avg,bott,D):
    data = []
    temp = 0.
    temp2 = 0.
    
    if avg<((sur+bott)*.5):
        dat1 = -50
        dat2 = -0.01
    else:
        dat1 = 0.01
        dat2 = 50
        
    while abs(dat2-dat1)>0.0005:
        f1 = (math.exp(dat1*sur)-math.exp(dat1*bott))/(sur*math.exp(dat1*sur)-bott*math.exp(dat1*bott)-avg*math.exp(dat1*sur)+avg*math.exp(dat1*bott)) - dat1
        f2 = (math.exp(dat2*sur)-math.exp(dat2*bott))/(sur*math.exp(dat2*sur)-bott*math.exp(dat2*bott)-avg*math.exp(dat2*sur)+avg*math.exp(dat2*bott)) - dat2

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
    try:
        dao1 = math.log(dat1/(math.exp(dat1*sur)-math.exp(dat1*bott)))+1
    except Exception:
        pass
                    
    for z in range(5,int(D)+5,5):
        try:
            theta = math.log(math.exp(dat1*sur)-(dat1*math.exp(1.-dao1)*(z-5.)/(D-5.)))/dat1
            data.append(theta)
        except Exception:
            pass      
    
    return data
#-----------------------------------------------------------------

print('-------------------------------------------------')
print('POME Model - Wet Case')
print('---------------------')

sur = input('Enter surface effective soil moisture: ')
avg = input('Enter profile mean effective soil moisture: ')
bott = input('Enter bottom effective soil moisture: ')
D = input('Enter total profile depth in cm: ')


sur = float(sur)
avg = float(avg)
bott = float(bott)
D = int(D)

if sur > avg and avg > bott:
    # Wet Case
    profile = eqs_wet(sur,avg,bott,D)
    
    df = pd.DataFrame()
    df['Depth'] = np.arange(5,D+5,5)
    df['SM'] = profile
    
    # saving and plotting the profile result
    df.to_csv('DryCase.csv')
    plt.plot(df.SM, df.Depth, color = 'blue', marker = '*')
    plt.ylim(D,5)
    plt.xlim(0,1)
    plt.ylabel('Depth (cm)')
    plt.xlabel('Effective Soil Moisture')

else:    
    print('**********************************')
    print('ERROR :: Incorrect inputs for Wet case, please re-enter values!')
    print('**********************************')
    sys.exit()
    

