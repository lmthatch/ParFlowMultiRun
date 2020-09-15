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


def pullTimeSeries(runLen,nz,nclm=0,pullVars=[('sat','test.out.satur'),('press','test.out.press'),('clm','test.out.clm_output')]):
    ''' get full time series of pressure, saturation, and LSM data '''
    ''' input data:
        runLen is the number of hour of run files
        nz is the number of verticle layers
        pull vars is a list of tuples with the name for the column header output, and the file name prefixes
    '''

    allLayers = list(np.arange(1,nz+1)) # list of all soil layers

    # pull data for all variables
    firstVar = True
    for var, varpre in pullVars:

        # create column names (w/ more complexity for CLM variable)
        if var == 'clm':
            clmLayers = list(np.arange(1,nclm+1)) # list of all soil layers
            varColnames =  ['eflx_lh_tot','eflx_lwrad_out','eflx_sh_tot',
                            'eflx_soil_grnd','qflx_evap_tot','qflx_evap_grnd',
                            'qflx_evap_soi','qflx_evap_veg','qflx_trans_veg',
                            'qflx_infl','swe_out','t_grnd','qflx_qirr']
            varColnames.extend(['tsoil_' + str(l) for l in clmLayers]) # adds a soil temperature list
            varpost = 'C.pfb'
        else:
            varColnames = [var + '_' + str(l) for l in allLayers]
            varpost='.pfb'
            
        varData = pullSingleVar(varpre,runLen,varColnames,post=varpost)

        if firstVar:
            allData = varData
            firstVar = False
        else:
            allData = pd.concat([allData,varData],axis=1)

    return allData


def pullSingleVar(prefix,runLen,colnames,post='pfb'):
    # Loop through all files and merge
    for k in range(1,(runLen+1)):
        fin = prefix + ".{:05d}.".format(k) + post
        datINpf = pfio.pfread(fin)
        datDF = pd.DataFrame(np.transpose(datINpf[:,:,0]),columns=colnames)

        if k==1:
            datAll = datDF
        else:
            datAll = datAll.append(datDF, ignore_index=True) 

    datAll = datAll.reset_index()
    
    return datAll

def addStorage(allData,nz,dz,ss,por):

    allLayers = list(np.arange(1,nz+1)) # list of all soil layers
    sat_colnames = ['sat_' + str(l) for l in allLayers]
    pres_colnames = ['press_' + str(l) for l in allLayers]
    sto_colnames = ['sto_' + str(l) for l in allLayers]

    for i in range(nz):
        sat = allData[sat_colnames[i]]
        pres = allData[pres_colnames[i]]
        sto =  pres.multiply(sat) * dz * ss + sat * por * dz
        allData[sto_colnames[i]] = sto

    # total storage over time
    allData['totalSto'] = allData[sto_colnames].sum(axis='columns')

    return allData

def processDataSC(rpars, saveAllData=True, saveAllSL=True, calcSto = True, saveAllSto = True, calcStoStats=True, saveStoSL=True):
    '''process parflow output data'''
    n = rpars['n']
    nz = rpars['ComputationalGrid.NZ']
    runLen = rpars['TimingInfo.StopTime']

    # check if clm is being used
    try:
        LSM = rpars['Solver.LSM'] == 'CLM'
        nclm = rpars['Solver.CLM.RootZoneNZ']
    except:
        LSM = False
        nclm = 0

    pullFiles = [('sat','test.out.satur'),('press','test.out.press'),('clm','test.out.clm_output')]
    allData = pullTimeSeries(runLen,nz,nclm=nclm,pullVars=pullFiles)

    # add storage data
    if calcSto:
        # calculate storage
        # storage is dx*dy*dz*press*Ss*Sat + Sat*porosity
        dz = rpars['ComputationalGrid.DZ']
        por = rpars['Geom.domain.Porosity.Value']
        ss = rpars['Geom.domain.SpecificStorage.Value']
        allData = addStorage(allData,nz,dz,ss,por)

    # save all press/sat/clm data
    if saveAllData:
        fileOut = '../FullRunData/FullRunData_test' + str(n) + ".csv"
        allData.to_csv(fileOut,index=False) 

    # save storage single line data
    if saveStoSL:
        # save single line total storage
        fileOut = '../TotalStorageData/TotSto_test' + str(n) + ".csv"
        allData.totalSto.to_csv(fileOut,index=False)

    # calculate/save all slingle line data
    if saveAllSL:

        if LSM:
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
            tempVars.extend(['tsoil_' + str(l) for l in clmLayers])
            tempData = allData[tempVars]
            tempTots = tempData.agg(['mean'])
            tempTots.columns = ['Mean_' + s for s in tempTots.columns]
        
        chngStoData = allData[sto_colnames]
        initSto = chngStoData.iloc[0]
        finalSto = chngStoData.iloc[-1]
        diffSto = finalSto.subtract(initSto)
        chngStoTots = pd.DataFrame(diffSto).transpose()
        chngStoTots.columns = ['Chng_' + s for s in chngStoTots.columns] 

        # merge data frames for Single Line Output
        chngStoTots=chngStoTots.reset_index(drop=True)

        if LSM:
            fluxTots=fluxTots.reset_index(drop=True)
            sweTots=sweTots.reset_index(drop=True)    
            tempTots=tempTots.reset_index(drop=True)
            singleLineOut = pd.concat([fluxTots,sweTots,tempTots,chngStoTots], axis=1) # merge the datasets
        else:
            singleLineOut = chngStoTots
        #chngPTots=chngPTots.reset_index(drop=True)
        #chngSTots=chngSTots.reset_index(drop=True)    
        #merge1 = pd.concat([fluxTots,sweTots],axis=1,ignore_index = True)
        #merge2 = pd.concat([merge1,tempTots],axis=1,ignore_index = True)
        #merge3 = pd.concat([merge2,chngPTots],axis=1,ignore_index = True)
        #singleLineOut = pd.concat([merge2,chngSTots],axis=1,ignore_index = True)
        
        fileOut = '../SingleLineOutput/SingleLineOutput_test' + str(n) + ".csv"
        singleLineOut.to_csv(fileOut,index=False)


        # calculate and save storage statistics
        # add data for 'peak' storage, for all layers and for total storage
        # 'output' variables
        # initial storage - gets time right before rain
        if calcStoStats:
            hrsPreRain = rpars['Cycle.prerainrec.pre.Length'] - 1
            allLayers = list(np.arange(1,nz+1)) # list of all soil layers
            sto_colnames = ['sto_' + str(l) for l in allLayers]
            sto_colnames.append('totalSto')
            stoData = allData[sto_colnames]
            allSto_init = stoData.iloc[hrsPreRain,]
            allSto_final = stoData.iloc[-1,]
            allSto_max = stoData.max(axis=0)
            allSto_hourMax = stoData.idxmax(axis=0)

            # make column names
            initCols = pd.Series([s + '_init' for s in sto_colnames])
            finalCols = pd.Series([s + '_final' for s in sto_colnames])
            maxCols = pd.Series([s + '_max' for s in sto_colnames])
            hourMaxCols = pd.Series([s + '_hrMax' for s in sto_colnames])
            allCols = pd.concat([initCols,finalCols,maxCols,hourMaxCols])

            # pull data together and write it out for Total Storage Data
            allStoData_out = pd.concat([allSto_init,allSto_final,allSto_max,allSto_hourMax])
            allStoDataDF = pd.DataFrame(allStoData_out).transpose()
            allStoDataDF.columns=allCols
            fileOut = '../SingleLineOutput/SL_Storage_test' + str(n) + ".csv"
            allStoDataDF.to_csv(fileOut,index=False)
