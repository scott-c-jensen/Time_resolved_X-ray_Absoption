#!/usr/bin/env python
"""
Library to read in the x-ray absorption data, process it, background subtract it, align it etc.
Basically takes the 5 laser flash sequence data and conditions it to have:
    single kinetics in each list
    Background subtracted
    Aligned by the laser
    rebinned

Scott Jensen
"""
import numpy as np
import matplotlib.pyplot as plt  
import scipy.signal as sci
#from scipy.interpolate import UnivariateSpline
import h5py
import scottLib as sl
import lmfit

def AlignData(dataList, scatterList, flashList):
    """
    Takes sets of data and shifts them so they align by flashes
    They are aligned to the minimum flash index
    The initial points are removed and zeros are added to the end
    A new data list and the flash times are output
    dataList: List of data Lists
    scatterList: List of scatter lists (scatter off the sample)
    flashList: List of times where the laser was illuminated
    """
    alignedDataList =[]
    alignedScatterList = []
    alignIdx = np.min(flashList) #lowest flash idx in all lists
    minListIdx = np.argmin(flashList)//len(flashList[0]) #gets the list of flashTimes that has the min value
    print(alignIdx)
    print(minListIdx)
    for i in range(len(flashList)):
        offset = np.min(flashList[i])-alignIdx
        print(offset)
        print(i)
        if offset == 0:
            print('zero Offset')
            alignedDataList.append(dataList[i])
            alignedScatterList.append(scatterList[i])
        else:
            alignedDataList.append(np.concatenate((dataList[i][offset:],np.zeros(offset))))
            alignedScatterList.append(np.concatenate((scatterList[i][offset:],np.zeros(offset))))
    return(alignedDataList, alignedScatterList, flashList[minListIdx])
    
def sumData(array):
    return(np.array(array).sum(axis = 0))
    
def bkgSubtract(Data, flashTimes):
    """
    Takes data and subtracts the background after each flash (15000 before 25000 after)
    bkgRange is the range of time before and after the flash that is used for establishing baseline bkg removal
    Data: list of data
    flashTimes: times where the laser illumination occured
    """
    #bkgRange=[[8000+15000,15000+15000]]
    #bkgRange2 = [[3000+15000,8000+15000]]
    bkgRange = [[-6000+15000,-1000+15000]]
    bkgRange2 = [[-6000+15000,-1000+15000]]

    bkgSubData = np.array(Data)
    for i,flashT in enumerate(flashTimes):
        idxStart, idxEnd = flashT-15000, flashT+25000
        if i==2:
            bkgSubData[idxStart:idxEnd] = sl.bkgLinear(np.arange(len(Data[idxStart:idxEnd])), Data[idxStart:idxEnd], bkgRange)
        else:
            bkgSubData[idxStart:idxEnd] = sl.bkgLinear(np.arange(len(Data[idxStart:idxEnd])), Data[idxStart:idxEnd], bkgRange2)

    return(bkgSubData)

def getLaserPos(laserTrace, clockRate=1000000,laserFrequency=10, numFlashes=5):
    """Finds the laser position for the laser trace
    laserTrace: laser intensity over time
    clockRate: recording clock frequency
    laserFrequency: repetition rate of the laser
    numFlashes: number of laser flashes on each sample"""
    #find max peak of laser flashes
    flashMax=[]
    usPerFlash = clockRate/laserFrequency
    
    #Find the max num of flashes in the range
    for flashNum in range(numFlashes):
        rangeStart=usPerFlash*flashNum #I shift all data to be at 100000 for the first flash so I put that in the center of the range... (ie start 50000)
        if flashNum == numFlashes - 1 and len(laserTrace)<rangeStart+usPerFlash:
            rangeStop = len(laserTrace)
        else:
            rangeStop = rangeStart+usPerFlash
        flashMax.append(np.argmax(laserTrace[rangeStart:rangeStop])+rangeStart)
    return(flashMax)

def truncate(xIn, yIn, xStart, xStop):
    """Removes data from xIn and yIn that are outside xStart and xStop
    xIn, yIn: x and y arrays for data
    xStart and xStop: x Values that define the range included in the output arrays"""
    xx = np.asarray(xIn)
    yy = np.asarray(yIn)
    xtemp = xx[xx>=xStart]
    ytemp = yy[xx>=xStart]
    xout = xtemp[xtemp<=xStop]
    yout = ytemp[xtemp<=xStop]
    return xout, yout
      
def getLinearBkg(xIn, yIn, idxStart, idxStop):
    """Finds the linear background from data based on the range between idxStart and idxStop
    XIn, yIn: the data input
    idxStart,idxStop: index values into xIn that determine the background range"""
    xtemp, ytemp = truncate(xIn, yIn, idxStart, idxStop)
       
    linFit = np.polyfit(xtemp, ytemp, 1)
    return (linFit)
    
def getKineticData (xIn, yIn, start, stop, step):
    """
    This function gets the range of data and offsets the data based on 1ms before the flash
    xIn,yIn: absorption data as arrays
    start, stop, step: microsecond positon in array to extract (also equivalently the index in the array)
    """
    
    zeroDataAverage = 1000 #Number of points to average before laser flash to set as zero
    #Rebin, Set the average of previous ~50 us to be the initial zero, then manually set first point as zero
    if step ==1:
        xData, yData =xIn[start:stop],np.array(yIn[start:stop])
        print ('This is yOffset before setting to zero: {}').format(yData[0])
    else:
        start = start+1 #offset so the us where the sample was hit is not included
        xData, yData = sl.rebin (xIn, yIn, start, stop, step)
    xData = xData-start
    yData = yData -np.mean(yIn[start-zeroDataAverage:start-1])
    xData[0] = 0
    yData[0] = 0
    
    return(xData,yData)
    
    
def addWeights(xTrace, usWeightRange, var, weightFactor):
    """Gets weights for fitting the data
    xTrace: list of x data
    useWeightRange: the range of values to weight
    var: variance of the data
    weightFactor: how much to weight the data by"""
    lowRange = (xTrace<=usWeightRange)*1/np.sqrt(var)*weightFactor
    highRange = (xTrace>usWeightRange)*1/np.sqrt(var)
    weightList = lowRange+highRange
    return weightList

    
def getVar(oldX, inputData, varStart, varEnd, step):
    """gets the variance for the data
    oldX: range of x values
    inputData: y values
    varStart: where to start looking at the variance
    varEnd: 
    step: x step size to rebin to"""
    if step ==1:
        var = np.var(inputData[varStart:varEnd])
    else:
        _, binY = sl.rebin(oldX, inputData, varStart, varEnd, step)
        var = np.var(binY)
    return var

    
def plotList(listIn):
    """Plots the data in listIn
    listIn: data as a list
    """
    for plotData in listIn:
        plt.figure()
        plt.clf()
        
        plt.plot(plotData, linewidth=1)
        
        plt.legend()
        plt.legend(bbox_to_anchor=(0, 1), loc=2, borderaxespad=0., frameon=False,  fontsize = '16')
        plt.xlabel('Time bins (20us)', fontsize='18')
        plt.ylabel('Normalized Intensity', fontsize= '18')
        #plt.figtext(x=.1, y=.05, fontdict=None,  fontsize = '16') 
        plt.rc('xtick', labelsize=13.5) 
        plt.rc('ytick', labelsize=13.5)
        plt.rcParams.update({'font.size': 13.5})
        plt.show()
    return

def loadData(fileIn):
    """ Loads data traces (time resolved x-ray absorption spectra)
    fileIn: Input file path (str)
    outputs (Mn fluorescence, Scattering, flashTimes)
    """
    with h5py.File(fileIn, 'r') as f:
        traceData=np.squeeze(np.asarray(f['trace']))
    flashTimes = getLaserPos(traceData[2,:])
    #MnFluorescence, Scatter, LaserFlash times (bin numbers)
    return(traceData[0,:],traceData[1,:], flashTimes)

def loadList(fileList):
    """
    Takes the list of file names for h5 files
    Outputs the DataTrace and the flashTimestra
    """
    traceDList, flashTList, traceScatterList = [],[],[]
    for fileIn in fileList:
        traceD, traceScatter, flashT = loadData(fileIn)
        traceDList.append(traceD)
        traceScatterList.append(traceScatter)
        flashTList.append(flashT)
    return(traceDList,traceScatterList,flashTList)

def getData(filesIn,
            timeBin = 1, 
            kineticTime=[3000, 3000, 6000, 3000, 3000],
            finalBkg = [5000,5000,5000,5000,5000],
            usWeightRange = [500, 300, 1500, 300, 500], 
            smoothNum=0,
            weightFactor=1 ):
    """Gets the original time sequence, and laser sequence
    Takes the data around the laser flash to isolate the kinetic changes
    removes the background based baseline data after and before laser flash
    determines the var and standard deviation of the data
    filesIn: file path to files
    timeBin: x bin size in microseconds
    kineticTime: how much time to ensure all kinetics are caught
    finalBkg: where to start taking background data
    useWeightRange: the range to weight data
    smoothNum: window length for a linear interpolation/smoothing (default 0 meaning no filter)
    weightFactor: how much to weight the data by at the begining of the transition (default 1 meaning no weighting)
    returns data, weights, and standard deviation of the baseline data
    """
    inputDataList, scatterList, flashTimesList = loadList(filesIn) #load Data
    
    alignedData, scatterTrace, flashTimes = AlignData(inputDataList, scatterList, flashTimesList) #align Data
    #alignedData = [alignedData[i]/scatterTrace[i] for i in range(len(alignedData))]
    inputData=sumData(alignedData) # combine data sets

    Std = np.std(inputData[flashTimes[1]-5000:flashTimes[1]])
    print('Std: {}'.format(Std))
    
    kineticList, xList, weights = [],[], []

    for i, flash in enumerate(flashTimes):
    
        oldX = range(len(inputData))
    
        # Here I get data that is offset to zero, rebinned and has a range for bkg
        xTraceBkg, yTraceBkg = getKineticData (oldX, inputData, flash, flash+kineticTime[i]+finalBkg[i], timeBin)
        
        #remove Data before the laser flash
        xTrace = xTraceBkg[xTraceBkg<=kineticTime[i]]
        yTrace = yTraceBkg[xTraceBkg<=kineticTime[i]]
        
        #BkgRemoval
        BkgSlope, BkgConst = getLinearBkg(xTraceBkg, yTraceBkg, kineticTime[i],kineticTime[i]+finalBkg[i])
        yNoBkg = yTrace - BkgSlope*xTrace #Note that we want yBkg[0] to remain at zero so the constant offset is not used
    
        #smooth if desired
        if smoothNum != 0:
            yTrace = sci.savgol_filter(yNoBkg, smoothNum, 1)
    
        #Making the weights with the variance
        if i == 0:
            var = getVar(oldX, inputData, flash-1000, flash, timeBin)
        weightList = addWeights(xTrace, usWeightRange[i], var, weightFactor)
        #weightList = [var**0.5]*len(xTrace)
    
        weights.append(weightList)
        kineticList.append(yTrace)
        xList.append(xTrace)
    return(xList, kineticList, weights, Std)
