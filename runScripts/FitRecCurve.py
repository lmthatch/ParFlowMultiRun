


from .RecCurveModels import exponentialMshift, boussinesqshift, coutagneshift
from scipy.optimize import curve_fit
import numpy as np

def subsetRec(allSto):
    ''' 
    Subset the recession period from storage data array/series 
    Recession period is defined as all storage values after max storage 
    
    Args:  
        allSto: pandas Series with storage data over time
    
    Returns:
        pandas Series with recession storage data
    '''

    # subset recession
    maxStoIdx = allSto.idxmax() # find index of max storage
    recSto = allSto[maxStoIdx:]
    return recSto

def fitModel(model,xvalues,yvalues,p0s):
    '''
    Fits Data (xvalues, yvalues) to desired model
    using Scipy's curve_fit function

    Args:
        model: function of desired model to fit
        xvalues: x values to fit (time values)
        yvalues: y values to fit (storage values)
        p0s: model parameter initial 'guess' values
    
    Returns: 
        fitted model parameters (or empty array if fit fails)
        model fit statistics from 'getfitstats' (or empty array if fit fails)
    '''

    # try to fit model
    try:
        popt = curve_fit(model,xvalues,yvalues,p0=p0s)[0]
        model_stats = getFitStats(model(xvalues, *popt), yvalues)
    
    # return empty array if model fit fails
    except:
        popt = np.empty(len(p0s))
        model_stats = np.empty(3)

    return popt, model_stats

def as_csv(array):
    '''
    Converts numpy array to csv seperated string

    Args:
        array: numpy array 
    
    Returns:
        string of array values seperated by commas
    '''

    return ','.join([str(i) for i in array]) + '\n'

def getFitStats(yest,yData):
    '''
    Calculates model fit statistics:
        Mean Squared Error (MSE)
        Root Mean Squared error (RMSE)
        Rsquared 

    Args:
        yest: estimated y value from model fit (model est storage)
        yData: actual y value (parflow storage)

    Returns:
        Tuple with MSE, RMSE, and Rsquared
    '''

    absError = yest - yData
    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    return MSE,RMSE,Rsquared

def fitRecCurve(stoDat, outFile):
    '''Fits Recession Curves to Parflow Storage Recession Data
    Fits Multiple Recession Curve Models
        Mal
        Boussinesq
        Coutainge
    Writes out model fit statistics to given file

    Args:
        stoDat: pandas series with storage over time
        outFile: string with outfile full path location to save storage statistics
    
    No Returns
    '''

    # subset the 'rain-rec period' from the full storage data
    recSto = subsetRec(stoDat.squeeze()).to_numpy() # and convert to numpy array

    # get min and max storage
    maxSto = stoDat.max()
    minSto = stoDat.min()

    # get model optimal parameters
    recLen = len(recSto)
    xvalues = np.arange(0,recLen)
    popt_mals, stats_mals = fitModel(exponentialMshift,xvalues,recSto, [maxSto,minSto, 0.047])
    popt_boss, stats_boss = fitModel(boussinesqshift,xvalues,recSto, [maxSto,minSto, 0.036])
    popt_cous, stats_cous = fitModel(coutagneshift, xvalues, recSto, [maxSto,minSto,0.045, 1.012])

    # merge data together to write out
    mergeStats = np.concatenate([popt_mals,popt_boss,popt_cous,stats_mals,stats_boss, stats_cous])
        
    # save the data w/ Header
    filehead=np.asarray(['Mal_Qo','Mal_Qf','Mal_a',
    'Bous_Qo','Bous_Qf','Bous_a',
    'Cou_Qo','Cou_Qf','Cou_a','Cou_b',
    'Mal_MSE','Mal_RMSE','Mal_Rsq',
    'Bous_MSE','Bous_RMSE','Bous_Rsq',
    'Cou_MSE','Cou_RMSE','Cou_Rsq'])

    with open(outFile,'w') as f:
        f.write(as_csv(filehead))
        f.write(as_csv(mergeStats))
