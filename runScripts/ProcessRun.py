# Run Processing
# Gets parameters from parameter file generated by SC_GenerateParamSet.py

import pandas as pd
import numpy as np
import random
import sys
import os
import math
import pfio #Hoang's Parflow Read/Write Module
import time
import datetime

from FitRecCurve import *

def pullTimeSeries(runLen,nz,nclm=0,pullVars=[('sat','test.out.satur'),('press','test.out.press')]):
    ''' 
    Pulls all Parflow output data for all variables including:
        pressure
        saturation 
        CLM/LSM (if applicable)

    Args:
        runLen: length, in hours, of parflow run
        nz: number of vertical layers in the parflow model
        nclm: number of soil layers shared between parflow & clm (if applicable)
        pullVars: list of tuples with variable name and file prefix
        clmColnames: column names for clm data (if applicable)

    Returns: 
        Pandas Dataframe with all Parflow output data
    '''

    print('Pulling Time Series data')

    # list of all soil layers (based on nz)
    allLayers = list(np.arange(1,nz+1)) 

    # pull data for all variables
    firstVar = True
    for var, varpre in pullVars:
        print('current variable: ' + var)

        # create column names (w/ more complexity for CLM variable)
        if var == 'clm':
            clmLayers = list(np.arange(1,nclm+1)) # list of all CLM soil layers (based on nclm)
            varColnames =  ['eflx_lh_tot','eflx_lwrad_out','eflx_sh_tot',
                            'eflx_soil_grnd','qflx_evap_tot','qflx_evap_grnd',
                            'qflx_evap_soi','qflx_evap_veg','qflx_trans_veg',
                            'qflx_infl','swe_out','t_grnd','qflx_qirr']
            varColnames.extend(['tsoil_' + str(l) for l in clmLayers]) # adds a soil temperature list
            varpost = 'C.pfb'
        else:
            varColnames = [var + '_' + str(l) for l in allLayers]
            varpost='pfb'

        # pull variable data  
        varData = pullSingleVar(varpre,runLen,varColnames,post=varpost)

        # merge into data frame with all variables 
        if firstVar:
            allVarData = varData
            firstVar = False
        else:
            allVarData = pd.concat([allVarData,varData],axis=1)

    return allVarData


def pullSingleVar(prefix,runLen,colnames,post='pfb'):
    '''
    Pulls full time series for a single variable

    Args:
        prefix: parflow output file prefix
        runLen: length, in hours, of parflow run (and output data) 
        colnames: column names for output data/variable
        post: file extension (default 'pfb', clmfiles 'C.pfb')
    
    Returns:
        Pandas Dataframe with all variable data
    '''

    # Loop through all files by 'hour' and merge
    for k in range(1,(runLen+1)):
        
        # get hourly data
        fin = prefix + ".{:05d}.".format(k) + post # set file name
        #print('pulling: ',fin)
        
        incompleteRun = False
        try:
            datINpf = pfio.pfread(fin) # read in parflow file
        except:
            print('PROBLEM PULLING FILE')
            incompleteRun = True

        if incompleteRun:
            break

        datDF = pd.DataFrame(np.transpose(datINpf[:,:,0]),columns=colnames) 

        # merge hourly data into single dataframe
        if k==1:
            datAll = datDF
        else:
            datAll = datAll.append(datDF, ignore_index=True) 

    datAll = datAll.reset_index() # reset index for easier future merging
    
    return datAll

def addStorage(allData,nz,dz,ss,por):
    '''
    Calculates storage data for all parflow timesteps, at all layers
    Adds to Pandas DataFrame with all parflow output data

    Args:
        allData: all Parflow output data (Press, Saturation, and CLM if applicable)
        nz: number of vertical layers in the parflow model
        dz: vertical layer depth (currently only constant layer depths is configued)
        ss: specific storage
        por: porosity
    
    Returns:
        tuple:
            Panda Dataframe: updated parflow output and storage data
            List of Storage column names (for easier access later)
    '''

    allLayers = list(np.arange(1,nz+1)) # list of all soil layers
    sat_colnames = ['sat_' + str(l) for l in allLayers]
    pres_colnames = ['press_' + str(l) for l in allLayers]
    sto_colnames = ['sto_' + str(l) for l in allLayers]

    # loop through each layer and calculate storage for that layer
    for i in range(nz):
        sat = allData[sat_colnames[i]]
        pres = allData[pres_colnames[i]]
        sto =  pres.multiply(sat) * dz * ss + sat * por * dz
        allData[sto_colnames[i]] = sto

    # total storage over time
    allData['totalSto'] = allData[sto_colnames].sum(axis='columns')
    sto_colnames.append('totalSto') # add to storage column names 

    return allData, sto_colnames

def processDataSC(rpars,parDict): #,saveAllPFData,saveTotStoSL,saveRecCurve_Total, saveRecCurve_Layers, saveCLMSL, saveStoStats):
    '''
    Processes Parflow Run Output. 
    Data is processed based on processing selections in the SCInput.txt file (read in through ParDict)
    
    Args:
        rpars: Pandas Dataframe with Parflow Run Parameters
        parDict: Dictionary w/ ParflowMultiRun Parameters, including processing requirements
    
    Returns:
        None, all data written out to file.
    
    '''

    # Key Parflow Run Parameters
    n = rpars['n']
    nz = rpars['ComputationalGrid.NZ']
    runLen = rpars['TimingInfo.StopTime']

    # evaluate which parflow output to collect (varies depending on whether CLM/LSM was used)
    try:
        LSM = rpars['Solver.LSM'] == 'CLM' # check if clm is being used
        nclm = rpars['Solver.CLM.RootZoneNZ']
        pullFiles = [('sat','test.out.satur'),('press','test.out.press'),('clm','test.out.clm_output')]     
    except:
        LSM = False
        nclm = 0
        pullFiles = [('sat','test.out.satur'),('press','test.out.press')]

    # pull all Parflow Run Output Data
    allData = pullTimeSeries(runLen,nz,nclm=nclm,pullVars=pullFiles)

    # save all press/sat/clm data (if applicable)
    if parDict['saveAllPFData']:
        fileOut = '../FullRunData/FullRunData_run' + str(n) + ".csv"
        allData.to_csv(fileOut,index=False) 

    # add storage data (as needed)
    if any([parDict['saveTotStoSL'], parDict['saveRecCurve_Total'], parDict['saveRecCurve_Layers'], parDict['saveStoStats']]):
        # calculate storage
        # storage is dx*dy*dz*press*Ss*Sat + Sat*porosity
        dz = rpars['ComputationalGrid.DZ']
        por = rpars['Geom.domain.Porosity.Value']
        ss = rpars['Geom.domain.SpecificStorage.Value']
        allData, sto_colnames = addStorage(allData,nz,dz,ss,por)

    # save storage single line data (if applicable)
    if parDict['saveTotStoSL']:
        # save single line total storage
        fileOut = '../SingleLineOutput/SL_TotStoAllHr_run' + str(n) + ".csv"
        allData.totalSto.to_csv(fileOut,index=False)

    # save Total Storage Recession Curve Data as single line (if applicable)
    if parDict['saveRecCurve_Total']:
        fileOut = '../SingleLineOutput/SL_TotStoRecCurveFit_run' + str(n) + '.csv'
        fitRecCurve(allData.totalSto, fileOut)
    
    # save Layer Recession Curve Data as single line file (if applicable)
    if parDict['saveRecCurve_Layers']:
        # save storage curves for [top 0-0.1m (if dz <= 0.1), 0-1m, 0-2m, 1m-2m, 2m-10m]
        
        # calculate each layer's 'centroid' depth
        layCent = np.arange(0,nz) * dz + dz/2

        ## top 1m
        layNums = np.where(layCent > 9)[0] + 1 #layer numbers are not zero-indexed
        layNames = ['sto_' + str(l) for l in layNums]
        laySto = allData[layNames].mean(axis=1).to_numpy()
        # rec1
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1m_Rec1_run' + str(n) + '.csv'
        fitRecCurve(laySto[17:90], fileOut)
        # rec2
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1m_Rec2_run' + str(n) + '.csv'
        fitRecCurve(laySto[95:168], fileOut)
        # rec3
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1m_Rec3_run' + str(n) + '.csv'
        fitRecCurve(laySto[173:246], fileOut)

        # top 2m
        layNums = np.where(layCent > 8)[0] +1 
        layNames = ['sto_' + str(l) for l in layNums]
        laySto = allData[layNames].mean(axis=1).to_numpy()
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_2m_Rec1_run' + str(n) + '.csv'
        fitRecCurve(laySto[17:90], fileOut)
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_2m_Rec2_run' + str(n) + '.csv'
        fitRecCurve(laySto[95:168], fileOut)
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_2m_Rec3_run' + str(n) + '.csv'
        fitRecCurve(laySto[173:246], fileOut)

        # 1-2m
        layNums = np.where((layCent < 9) & (layCent > 8))[0] + 1
        layNames = ['sto_' + str(l) for l in layNums]
        laySto = allData[layNames].mean(axis=1).to_numpy()
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1-2m_Rec1_run' + str(n) + '.csv'
        fitRecCurve(laySto[17:90], fileOut)
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1-2m_Rec2_run' + str(n) + '.csv'
        fitRecCurve(laySto[95:168], fileOut)
        fileOut = '../SingleLineOutput/SL_StoRecCurveFit_1-2m_Rec3_run' + str(n) + '.csv'
        fitRecCurve(laySto[173:246], fileOut)

        # bottom 8m
        #layNums = np.where(layCent < 8)[0] +1 
        #layNames = ['sto_' + str(l) for l in layNums]
        #laySto = allData[layNames].mean(axis=1)
        #fileOut = '../SingleLineOutput/SL_StoRecCurveFit_2-10m_run' + str(n) + '.csv'
        #fitRecCurve(laySto, fileOut)

        # top 0.1 m
        #if dz <= 0.1:
        #    layNums = np.where(layCent < 9.9)[0] + 1
        #    layNames = ['sto_' + str(l) for l in layNums]
        #    laySto = allData[layNames].mean(axis=1)
        #    fileOut = '../SingleLineOutput/SL_StoRecCurveFit_01m_run' + str(n) + '.csv'
        #    fitRecCurve(laySto, fileOut)

    # save all CLM hourly data
    if parDict['saveCLMAll']: 
        clmVars = ['qflx_evap_tot','qflx_evap_grnd',
            'qflx_evap_soi','qflx_evap_veg','qflx_trans_veg',
            'qflx_infl','qflx_qirr','swe_out','t_grnd'] # ignoring soil temperatures for now
        fileOut = '../FullRunData/FullCLMData_run' + str(n) + ".csv"
        allData[clmVars].to_csv(fileOut,index=False)


    # calculate/save all CLM single line data (if applicable)
    if parDict['saveCLMSL']:

        if LSM: # double check, can't save CLM data if CLM didn't run...
            fluxVars = ['qflx_evap_tot','qflx_evap_grnd',
                        'qflx_evap_soi','qflx_evap_veg','qflx_trans_veg',
                        'qflx_infl','qflx_qirr']
            fluxData = allData[fluxVars] #%>% summarise_all(sum)
            fluxTots = fluxData.agg(['sum'])
            fluxTots.columns = ['Total_m_' + str(s) for s in fluxTots.columns]

            sweData = allData[['swe_out']]
            sweTots = sweData.agg(['mean'])
            sweTots.columns = ['Mean_m_' + str(s) for s in sweTots.columns]

            tempVars = ['t_grnd']
            clmLayers = list(np.arange(1,nclm+1))
            tempVars.extend(['tsoil_' + str(l) for l in clmLayers])
            tempData = allData[tempVars]
            tempTots = tempData.agg(['mean'])
            tempTots.columns = ['Mean_' + s for s in tempTots.columns]
        
            fluxTots=fluxTots.reset_index(drop=True)
            sweTots=sweTots.reset_index(drop=True)    
            tempTots=tempTots.reset_index(drop=True)
            singleLineOut = pd.concat([fluxTots,sweTots,tempTots], axis=1) # merge the datasets
            fileOut = '../SingleLineOutput/SL_CLM_run' + str(n) + ".csv"
            singleLineOut.to_csv(fileOut,index=False)
        
        else:
            print('WARNING: saveCLMSL set to true but CLM did not run...')

    # calculate and save storage statistics (if applicable)
    if parDict['saveStoStats']:
        
        # calculate storage statistics
        # initial, max, and final stroage for all layers

        # subset just storage data
        stoData = allData[sto_colnames]
        hrsPreRain = rpars['Cycle.prerainrec.pre.Length'] - 1 # get hour index before rain starts
        allSto_init = stoData.iloc[hrsPreRain,] # get intial storage before rain
        allSto_final = stoData.iloc[-1,] # get final storage at end of PF run
        allSto_max = stoData.max(axis=0) # get maximum storage
        allSto_hourMax = stoData.idxmax(axis=0) # get hour index of maximum storage

        # make column names for new variables
        initCols = pd.Series([s + '_init' for s in sto_colnames]) 
        finalCols = pd.Series([s + '_final' for s in sto_colnames])
        maxCols = pd.Series([s + '_max' for s in sto_colnames])
        hourMaxCols = pd.Series([s + '_hrMax' for s in sto_colnames])
        allCols = pd.concat([initCols,finalCols,maxCols,hourMaxCols]) 

        # pull data together and write it out for Total Storage Data
        allStoData_out = pd.concat([allSto_init,allSto_final,allSto_max,allSto_hourMax])
        allStoDataDF = pd.DataFrame(allStoData_out).transpose()
        allStoDataDF.columns=allCols # update Pandas DataFrame w/ column names
        fileOut = '../SingleLineOutput/SL_Sto_run' + str(n) + ".csv"
        allStoDataDF.to_csv(fileOut,index=False)

    if parDict['saveSM']:

        por = rpars['Geom.domain.Porosity.Value']

        # calculate each layer's 'centroid' depth
        layCent = np.arange(0,nz) * dz + dz/2

        # top 1m
        layNums = np.where(layCent > 9)[0] + 1 #layer numbers are not zero-indexed
        layNames = ['sat_' + str(l) for l in layNums]
        laySM = allData[layNames].mean(axis=1) * por
        fileOut = '../SingleLineOutput/SL_SM_1m_run' + str(n) + '.csv'
        laySM.to_csv(fileOut,index=False)

        # top 2m
        layNums = np.where(layCent > 8)[0] + 1 
        layNames = ['sat_' + str(l) for l in layNums]
        laySM = allData[layNames].mean(axis=1) * por
        fileOut = '../SingleLineOutput/SL_SM_2m_run' + str(n) + '.csv'
        laySM.to_csv(fileOut,index=False)
    
    if parDict['processFlow']:
        fileOut = '../SingleLineOutput/TopPress_run' + str(n) + '.csv'
        allData['press_' + str(nz)].to_csv(fileOut,index=False)
        #mannings = rpars['Mannings.Geom.domain.Value']
        #slope = rpars['TopoSlopesX.Geom.domain.Value']
        #dx = rpars['ComputationalGrid.DX']
        #q = [(p**(5/3))*(abs(slope)**(1/2))*dx/mannings for p in topPress]
        #q = pd.Series(q)