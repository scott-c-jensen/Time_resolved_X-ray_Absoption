"""
@author: Scott C. Jensen
2018

This program takes the time resolved x-ray absorption data, recorded with 1 µs time resolution.
Fits the data based on which kinetic model is of interest. The manually chosen for the s3_model.

There are different oxidation states (S-states) in photosystem II. They cycle from S0-S4. 
The S4-S0 transition doesn't require light to transition, returning to the most reduced state.
Not only are there multiple changes in the system at each state, but only the S1 (the dark adapted state)
is a pure state, others are mixed. This program therefore, has global fitting of all s-state transitions
for the advancement and individual state transitions with the percent of each state.

So for example:
If we have a 80% advancement then after each flash (F) we have
F0 = 100% S1
F1 = 20% S1, 80% S2

But what we are looking at is the change in absorption... so we only care about the changes
So if we make a matrix for the S-state transitions for each flash, such as the first transition S1 to S2 (or S1-S2)
This can be written for each flash as an array [S0-S1,S1-S2,S2-S3,S3-S0]:

1F = [   0, 0.80,   0.,   0.,]
2F = [   0, 0.16, 0.64,   0.,]
3F = [   0, 0.03, 0.26, 0.51,]
4F = [0.41, 0.01, 0.08, 0.31,]
5F = [0.33, 0.33, 0.02, 0.12,]

Numbers/names of amplidues
S3 model to use
"""

import Model_TRXAS_Data as md
import Model_TRXAS_Fit as mf
import functools

def main(d, s3_fn):
    """
    Loads in the data, combines it, and fits it to one of a number of kinetic models.
    Reports on fits and plots.
    Uses a few hard-coded parameters as the file locations, and initial guesses don't change.
    """
    inDir = 'C:\\Users\\Pushkar\\Desktop\\Beamtime_Data\\TR_XAS\\Output\\'
    fileNames = ['April_2017_2det.h5','April_2017_1Det.h5','July2017_1Det.h5','July2017_2det.h5','feb_19_BothDet.h5']
    filesIn = [inDir+fn for fn in fileNames] 

    #dictionary with lists for [inital guess, min value, max value]
    #Amplitudes are taken based on the data
    paramsDict = {'initS1':[.5,0,1], 
                    'pAdvance':[.8,0,1], 
                    'tauS01':[50,0,None], 
                    'tauS12':[90,0,None], 
                    'tauS23':[400,None,None], 
                    'tau1S30':[50,0,None], 
                    'tau2S30':[1300,0,None]} #example tauS01 is the decay constant for S0 to S1 transition
                    
    S3Model = functools.partial(s3_fn, **d)
    xData, yData, weightList, sd = md.getData(filesIn)
    fit = mf.fitData(xData, yData, weightList, mf.totalModel, s3_fn paramsDict)
    mf.getAll(fit)
    #fit.report()
    fit.plot()

s3_function = mf.sequentialRate # Model to test for the S3-S0 transition
    
if __name__ == "__main__":
    main(d, s3_function)