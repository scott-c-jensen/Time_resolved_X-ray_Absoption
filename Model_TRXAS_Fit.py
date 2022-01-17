#!/usr/bin/env python2
"""
Library of functions for fitting the series of time resolved x-ray absorption of photosystem II.
"""
import numpy as np
import lmfit

def singleRate(t, amp1, tau1, amp2=None, tau2=None):
    """This model treats the absortion changes as a simple chemical process:
    population A goes changes into B and the x-ray absorption changes by amp1
    t: time array in microseconds
    amp1: x-ray absorption difference bewteen the two states
    tau1: the decay time constant for the lifetime after being illuminated by a laser
    amp2,tau2: dummy variables"""
    value = amp1*(1-np.exp(-t/tau1))
    return value

def sequentialRate(t, amp1, amp2, tau1, tau2):
    """This model treats the absortion changes as a simple chemical process:
    population A goes changes into B then B changes in to C
    t: time array in microseconds
    amp1: x-ray absorption difference bewteen the initial and intermediate state
    amp2: x-ray absorption difference bewteen the intermediate and final state
    tau1: the decay time constant for the lifetime of the initial state after being illuminated by a laser
    tau2: the decay time constant for the lifetime of the intermediate state
    """
    k1=tau1**(-1)
    k2=tau2**(-1)
    value = amp2*(1+(k1*np.exp(-t*k2)-k2*np.exp(-t*k1))/(k2-k1))+amp1*(k1/(k2-k1)*(np.exp(-k1*t)-np.exp(-k2*t)))
    return value
     
def simultaneousRates(t, amp1, amp2, tau1,tau2):
    """This model treats the absortion changes as a simple chemical process:
    population A goes changes into either B or C at different rates
    t: time array in microseconds
    amp1: x-ray absorption difference bewteen the initial and state B
    amp2: x-ray absorption difference bewteen the initial and state C
    tau1: the decay time constant for the initial state to state B
    tau2: the decay time constant for the initial state to state c
    """
    k1=tau1**(-1)
    k2=tau2**(-1)
    value = amp1*(1-np.exp(-k1*t)) + amp2*(1-np.exp(-k2*t))
    return value    
    
def sequentialExp (t, amp1, tau1, tau2, amp2=None):
    """This model treats the absortion changes as a simple chemical process:
    population A goes changes into B then B changes in to C but state B doesn't affect absorption
    t: time array in microseconds
    amp1: x-ray absorption difference bewteen the initial and intermediate state
    amp2: x-ray absorption difference bewteen the intermediate and final state
    tau1: the decay time constant for the lifetime of the initial state after being illuminated by a laser
    tau2: the decay time constant for the lifetime of the intermediate state
    """
    k1=tau1**(-1)
    k2=tau2**(-1)
    value = amp1*(1+(k1*np.exp(-t*k2)-k2*np.exp(-t*k1))/(k2-k1))
    return value

    
def sStatePopulation(initS1, pAdvance, numFlashes):
    """
    Takes floats for all variables except flash num which is an int, outputs numpy array
    Finds the population of the S-states based on the inital conditions and advancement probability
    Technically the final population is unecessary but I left it in for reference/validation
    """
    if pAdvance>1 or initS1<0 or initS1>1:
        raise Exception('PSII Advancement and initial populations cannot be over 1')
        
    sStatePop = np.zeros((numFlashes+1,4))
    
    sStatePop[0,0]=1-initS1 #Initial population
    sStatePop[0,1]=initS1 #Initial population
    
    for f in range(1,numFlashes+1):
        for s in range(0,4):
            sStatePop[f,s] = sStatePop[f-1,s-1]*pAdvance+sStatePop[f-1,s]*(1-pAdvance)
    return (sStatePop)
    
def sStateAdvancement(sStatePop, pAdvance):
    """
    Finds the percentage of centers advancing in each sState advancing given the current sState Population
    sStatePop: numpy array of the population in each state 
    pAdvance: the percentage of advancement as a fraction 
    """
    return(sStatePop*pAdvance)
    
def getAdvance(initS1, pAdvance, numFlashes=5):
    """Given the initial s-states and the number of flashes to consider
    initS1: Initial population in the S1 state
    pAdvance: Percent of population that advances 
    numFlashes: Number of laser flashes to advance the states
    """
    return(sStateAdvancement(sStatePopulation(initS1, pAdvance, numFlashes),pAdvance))
    
def get_flashModel(s3Model):    
    """This function changes which model is used for fitting the S3 model
    S3Model: function to model the S3-S0 Transition"""
    def flashModel(t,  tauS01, tauS12, tauS23, tau1S30, tau2S30, ampS01, ampS12, ampS23=None,  amp2S30=None, amp1S30=None):
        """
        This function calculates all s-state advancements which are used for fitting flash transitions elsewhere
        t in this case is just the range of x-values used in the fit (ie not a list)
        This assumes the S0 is the final state in the transition S3 to S0
        """
        valuesS = np.zeros((np.shape(t)[0],4))
        #The value contributions based on the S-state contributions
        #S0 to S1 uses the index 0
        valuesS[0] = singleRate(t, ampS01, tauS01)
        valuesS[1] = singleRate(t, ampS12, tauS12)
        valuesS[2] = singleRate(t, ampS23, tauS23)
        valuesS[3] = s3Model(t, amp1=amp1S30, amp2=amp2S30, tau1=tau1S30, tau2=tau2S30) #For different S3-S0 models this will have to change
        return valuesS
    return(flashModel)

def makeTotalModel(s3Model):
    """Creates the full function to fit, the advancement probability, the amplitudes and the time constants for each state
    s3Model: model for the S3-S0 transition"""
    flashModel = get_flashModel(s3Model)
    def totalModel(tList, initS1, pAdvance, tauS01, tauS12, tauS23, tau1S30, tau2S30, ampS01, ampS12, ampS23, amp2S30, amp1S30):
        """The model that will be fit, contains all variables for the states and advancement probability.
        All states are single rate transitions except S3
        tList: List of times 
        initS1: the initial population of the S1 state
        pAdvance: the fractional probability of each state advancing
        tauSxx: the time constant for the transition Sxx (ie S30 is S3-S0) and used in the models
        ampSxx: the amplitude for the transition between Sxx"""
        advM = getAdvance(initS1, pAdvance)
        #get the theoretical kinetics using the longest time duration in the data
        tLongest = tList[np.argmax(tList)]
        traceFits = flashModel(tLongest,  tauS01, tauS12, tauS23, tau1S30, tau2S30,ampS01, ampS12, ampS23,  amp2S30, amp1S30)
        
        values = np.zeros(np.shape(tList))
        
        for i in range(len(tList)):
            values[i] = advM[i,1]*traceFits[:,:len(tList[i])]
        return values
    return(totalModel)
    
def getParams(yData,paramDict):
    """Gets the parameters for initial guesses
    This combines passed params with the amplitudes estimated from data
    yData: time resolved data
    paramDict: dictionary of parameters for fitting
    """
    params = lmfit.Parameters()
    
    #get amplitude guess based on the actual kinetics
    amp = ['ampS12', 'ampS23',  'amp2S30', 'ampS01']
    for i in range(len(amp)):
        ampGuess = np.average(yData[i][-1000:]) 
        params.add(amp[i],value = ampGuess, min=None, max=None)
        if i ==2:
            params.add('amp1S30', value=-ampGuess/3., min=None, max=None)      

    #Get dictionary of parameters
    for i in paramDict:
        params.add(i, value = paramDict[i][0], min = paramDict[i][1], max = paramDict[i][2])

    return params

def fitData(xData, yData, weightList, s3Model, paramD):
    """Fits x, Y data with weights on the initial data points according to the weightList
    xData, yData: x and y data for the absorption
    weightList: weighted list for weighting data in yData during fitting
    s3Model: model for S3-S0 transition
    paramD: parameters for fitting"""
    modelUse = makeTotalModel(s3Model)
    fitModel = lmfit.Model(modelUse)
    paramsIn = getParams(yData,paramD)
    fit = fitModel.fit(yData, t=xData, params=paramsIn, method='leastsq', weights=weightList)
    return(fit)
    
    
       