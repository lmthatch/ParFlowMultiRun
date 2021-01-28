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

def waterTableDepth(satCol,pressCol,dz):
    '''
    get watertable depth in a saturation column
    asumes wtd is the depth of all saturated cells + pressure at the next partially saturated cell

    Args:
        satCol: array of the cell saturations in the column
        pressCol: array of the cell pressure in the column
        dz: cell depth
    '''

    # reshape the columns
    satCol.shape=(1,len(satCol))
    satCol = satCol.flatten()
    #print('satCol')
    #print(satCol.shape)
    #print(satCol)
    pressCol.shape=(1,len(pressCol))
    pressCol=pressCol.flatten()
    zindex = np.arange(len(satCol))
    if len(zindex[satCol>= 1]) > 0:
        lastSatCell = max(zindex[satCol>= 1])
        satDepth = (lastSatCell + 1)*dz
    else:
        lastSatCell = 0
        satDepth = 0
    nextPress = pressCol[lastSatCell]

    return nextPress + satDepth # total depth


def wtdCrossSection(satCS,pressCS,dz):
    '''
    gets full wtd cross section from tilted 

    Args:
        satCS: array of saturation for the 2D cross section
        pressCS: array of pressure for the 2d cross section
        dz: vertical layer depths
    Returns:
        wtd: 1D numpy array of wtd over the cross section
    '''
    #print('Saturated Cross Section Size is....')
    #print(satCS.shape)
    nx = satCS.shape[2]
    wtd = []
    for x in range(nx):
        wtd.append(waterTableDepth(satCS[:,:,x],pressCS[:,:,x],dz))
    wtdArray = np.array(wtd)
    wtdArray.shape = (1,len(wtdArray))
    #print(wtdArray.shape)
    #print(wtdArray)
    return wtdArray


def processDataTV(rpars): #,saveAllPFData,saveTotStoSL,saveRecCurve_Total, saveRecCurve_Layers, saveCLMSL, saveStoStats):
    '''
    Processes Tilted V Parflow Run Output. 
    Data is processed based on processing selections in the SCInput.txt file (read in through ParDict)
    
    Args:
        rpars: Pandas Dataframe with Parflow Run Parameters
        parDict: Dictionary w/ ParflowMultiRun Parameters, including processing requirements
    
    Returns:
        None, all data written out to files.
    
    '''

    # Output Variables:
    #   Streamflow
    #   SM:
    #       Upper Hill: 1, 2m depths
    #       Lower Hill: 1, 2m depths
    #   Total Storage
    #       Full Domain
    #       Upper Hill
    #       Lower Hill

    dz = rpars['ComputationalGrid.DZ']
    dx = rpars['ComputationalGrid.DX']
    dy = rpars['ComputationalGrid.DY']
    nx = rpars['ComputationalGrid.NX']
    ny = rpars['ComputationalGrid.NY']
    por = rpars['Geom.domain.Porosity.Value']
    ss = rpars['Geom.domain.SpecificStorage.Value']
    runLen = rpars['TimingInfo.StopTime']
    mannings = rpars['Mannings.Geom.domain.Value']
    slope = rpars['TopoSlopesY.Geom.channel.Value']


    # get geometry information for upper and lower hillslopes and center channel
    centerX = int((nx - 1)/2)
    upHillXs = np.concatenate((np.arange(0,int(centerX/2)), np.arange(int(centerX/2 + centerX + 1), int(nx))))
    lowHillXs = np.concatenate((np.arange(int(centerX/2), centerX),np.arange(centerX + 1,int(centerX + 1 + centerX/2))))
    zpoints = np.arange(10 - dz/2, 0, -dz)
    #xpoints = np.arange(dx/2,nx*dx,dx)
    ypoints = np.arange(dy/2,ny*dy,dy)
    #zelev = np.array([(abs(x - centerX) - dx/2)*slope for x in xpoints]) # elevation for each grid cell along the x
    #print(ypoints)
    
    # Loop through all files by 'hour' and merge
    firstFile = True
    incompleteRun = False
    for k in range(1,(runLen+1)):
        
        # get hourly data
        pressFN = 'test.out.press.'+ "{:05d}".format(k) + '.pfb'
        satFN = 'test.out.satur.'+ "{:05d}".format(k) + '.pfb'
        
        try:
            pressDat = pfio.pfread(pressFN) # read in parflow file
            satDat = pfio.pfread(satFN) # read in parflow file
        except:
            print('PROBLEM PULLING FILE')
            incompleteRun = True

        if incompleteRun:
            break

        # calculate storages
        stoDat =  np.multiply(pressDat,satDat) * dz * ss + satDat * por * dz
        satDat2m = satDat[zpoints<2,:,:]
        satDat1m = satDat[zpoints<1,:,:]

        # cross section wtd
        if 500 in ypoints:
            wtd500 = wtdCrossSection(satDat[:,ypoints==500,:],pressDat[:,ypoints==500,:],dz)
            wtd1500 = wtdCrossSection(satDat[:,ypoints==1500,:],pressDat[:,ypoints==1500,:],dz)
            wtd2500 = wtdCrossSection(satDat[:,ypoints==2500,:],pressDat[:,ypoints==2500,:],dz)
            wtd3500 = wtdCrossSection(satDat[:,ypoints==3500,:],pressDat[:,ypoints==3500,:],dz)
            wtd4500 = wtdCrossSection(satDat[:,ypoints==4500,:],pressDat[:,ypoints==4500,:],dz)
        else: # get average of two cells to the sides
            wtd500a = wtdCrossSection(satDat[:,ypoints==(500-dx/2),:],pressDat[:,ypoints==(500-dx/2),:],dz)
            wtd500b = wtdCrossSection(satDat[:,ypoints==(500+dx/2),:],pressDat[:,ypoints==(500+dx/2),:],dz)
            wtd500 = np.array([np.mean([i,j]) for i,j in zip(wtd500a.flatten(),wtd500b.flatten())])
            wtd500.shape = wtd500a.shape
            wtd1500a = wtdCrossSection(satDat[:,ypoints==(1500-dx/2),:],pressDat[:,ypoints==(1500-dx/2),:],dz)
            wtd1500b = wtdCrossSection(satDat[:,ypoints==(1500+dx/2),:],pressDat[:,ypoints==(1500+dx/2),:],dz)
            wtd1500 = np.array([np.mean([i,j]) for i,j in zip(wtd1500a.flatten(),wtd1500b.flatten())])
            wtd1500.shape = wtd1500a.shape
            wtd2500a = wtdCrossSection(satDat[:,ypoints==(2500-dx/2),:],pressDat[:,ypoints==(2500-dx/2),:],dz)
            wtd2500b = wtdCrossSection(satDat[:,ypoints==(2500+dx/2),:],pressDat[:,ypoints==(2500+dx/2),:],dz)
            wtd2500 = np.array([np.mean([i,j]) for i,j in zip(wtd2500a.flatten(),wtd2500b.flatten())])
            wtd2500.shape = wtd2500a.shape
            wtd3500a = wtdCrossSection(satDat[:,ypoints==(3500-dx/2),:],pressDat[:,ypoints==(3500-dx/2),:],dz)
            wtd3500b = wtdCrossSection(satDat[:,ypoints==(3500+dx/2),:],pressDat[:,ypoints==(3500+dx/2),:],dz)
            wtd3500 = np.array([np.mean([i,j]) for i,j in zip(wtd3500a.flatten(),wtd3500b.flatten())])
            wtd3500.shape = wtd3500a.shape
            wtd4500a = wtdCrossSection(satDat[:,ypoints==(4500-dx/2),:],pressDat[:,ypoints==(4500-dx/2),:],dz)
            wtd4500b = wtdCrossSection(satDat[:,ypoints==(4500+dx/2),:],pressDat[:,ypoints==(4500+dx/2),:],dz)
            wtd4500 = np.array([np.mean([i,j]) for i,j in zip(wtd4500a.flatten(),wtd4500b.flatten())])
            wtd4500.shape = wtd4500a.shape

        if firstFile == True:
            #storage
            allSto = np.array(np.sum(stoDat))
            nchanSto = np.array(np.sum(stoDat[:,:,np.arange(stoDat.shape[2]) != centerX]))
            upSto = np.array(np.sum(stoDat[:,:,upHillXs]))
            lowSto = np.array(np.sum(stoDat[:,:,lowHillXs]))
            #soil moisture
            low1m = np.array(np.mean(satDat1m[:,:,lowHillXs]))
            up1m = np.array(np.mean(satDat1m[:,:,upHillXs]))
            all1m = np.array(np.mean(satDat1m[:,:,np.arange(satDat.shape[2]) != centerX]))
            low2m = np.array(np.mean(satDat2m[:,:,lowHillXs]))
            up2m = np.array(np.mean(satDat2m[:,:,upHillXs]))
            all2m = np.array(np.mean(satDat2m[:,:,np.arange(satDat.shape[2]) != centerX]))
            #flow/channel pressure
            topChanPress = np.array(pressDat[pressDat.shape[0]-1,pressDat.shape[1]-1,centerX])
            #wtd
            wtd500All = wtd500
            wtd1500All = wtd1500
            wtd2500All = wtd2500
            wtd3500All = wtd3500
            wtd4500All = wtd4500
            firstFile=False
        else:
            allSto = np.append(allSto,np.sum(stoDat))
            nchanSto = np.append(nchanSto,np.sum(stoDat[:,:,np.arange(stoDat.shape[2]) != centerX]))
            upSto = np.append(upSto,np.sum(stoDat[:,:,upHillXs]))
            lowSto = np.append(lowSto,np.sum(stoDat[:,:,lowHillXs]))
            low1m = np.append(low1m,np.mean(satDat1m[:,:,lowHillXs]))
            up1m = np.append(up1m,np.mean(satDat1m[:,:,upHillXs]))
            all1m = np.append(all1m,np.mean(satDat1m[:,:,np.arange(satDat.shape[2]) != centerX]))
            low2m = np.append(low2m,np.mean(satDat2m[:,:,lowHillXs]))
            up2m = np.append(up2m,np.mean(satDat2m[:,:,upHillXs]))
            all2m = np.append(all2m,np.mean(satDat2m[:,:,np.arange(satDat.shape[2]) != centerX]))
            topChanPress = np.append(topChanPress,pressDat[pressDat.shape[0]-1,pressDat.shape[1]-1,centerX])
            #wtd
            wtd500All = np.concatenate((wtd500All,wtd500),axis=0)
            wtd1500All = np.concatenate((wtd1500All,wtd1500))
            wtd2500All = np.concatenate((wtd2500All,wtd2500))
            wtd3500All = np.concatenate((wtd3500All,wtd3500))
            wtd4500All = np.concatenate((wtd4500All,wtd4500))

    # do final multiplying at the end to avoid unnecearry multiplication
    allSto = allSto * dx * dy
    nchanSto = nchanSto * dx * dy
    upSto = upSto*dx*dy
    lowSto = lowSto*dx*dy
    low1m = low1m*por
    up1m = up1m*por
    all1m = all1m*por
    low2m = low2m*por
    up2m = up2m*por
    all2m = all2m*por

    # flow 
    topChanPress[topChanPress < 0] = 0
    q = [(p**(5/3))*(abs(slope)**(1/2))*dx/mannings for p in topChanPress]

    # write out files
    outPre = '../SingleLineOutput/'
    np.savetxt(outPre + 'TotSto' + '_run' + str(rpars['n']) + '.csv',allSto,delimiter=',')
    np.savetxt(outPre + 'TotStowoChan' + '_run' + str(rpars['n']) + '.csv',nchanSto,delimiter=',')
    np.savetxt(outPre + 'SM_1m' + '_run' + str(rpars['n']) + '.csv',all1m,delimiter=',')
    np.savetxt(outPre + 'SM_2m' + '_run' + str(rpars['n']) + '.csv',all2m,delimiter=',')
    np.savetxt(outPre + 'LowerHillSM_1m' + '_run' + str(rpars['n']) + '.csv',up1m,delimiter=',')
    np.savetxt(outPre + 'LowerHillSM_2m' + '_run' + str(rpars['n']) + '.csv',up2m,delimiter=',')
    np.savetxt(outPre + 'UpperHillSM_1m' + '_run' + str(rpars['n']) + '.csv',low1m,delimiter=',')
    np.savetxt(outPre + 'UpperHillSM_2m' + '_run' + str(rpars['n']) + '.csv',low2m,delimiter=',')
    np.savetxt(outPre + 'Outflow' + '_run' + str(rpars['n']) + '.csv',q,delimiter=',')
    np.savetxt(outPre + 'WTD500m' + '_run' + str(rpars['n']) + '.csv',wtd500All,delimiter=',')
    np.savetxt(outPre + 'WTD1500m' + '_run' + str(rpars['n']) + '.csv',wtd1500All,delimiter=',')
    np.savetxt(outPre + 'WTD2500m' + '_run' + str(rpars['n']) + '.csv',wtd2500All,delimiter=',')
    np.savetxt(outPre + 'WTD3500m' + '_run' + str(rpars['n']) + '.csv',wtd3500All,delimiter=',')
    np.savetxt(outPre + 'WTD4500m' + '_run' + str(rpars['n']) + '.csv',wtd4500All,delimiter=',')

    # add clm processing
    # if needed...

    # get all factors at upper and lower hillslopes
    if rpars['Solver.LSM'] == 'CLM':

        nclm = 10 # rpars['Solver.CLM.RootZoneNZ']
        clmLayers = list(np.arange(1,nclm+1)) # list of all CLM soil layers (based on nclm)

        clmColnames =  ['eflx_lh_tot','eflx_lwrad_out','eflx_sh_tot',
                'eflx_soil_grnd','qflx_evap_tot','qflx_evap_grnd',
                'qflx_evap_soi','qflx_evap_veg','qflx_trans_veg',
                'qflx_infl','swe_out','t_grnd','qflx_qirr']
        clmColnames.extend(['tsoil_' + str(l) for l in clmLayers]) # adds a soil temperature list
        
        clmOutFile =  '../FullRunData/allCLMDat' + '_run' + str(rpars['n']) + '.csv'

        # loop through all hours
        firstFile = True
        incompleteRun = False
        for k in range(1,(runLen+1)):
            
            # get hourly data
            clmFN = 'test.out.clm_output.' + "{:05d}".format(k) + '.C.pfb'
            
            try:
                clmDat = pfio.pfread(clmFN) # read in parflow file
            except:
                print('PROBLEM PULLING CLM FILE')
                incompleteRun = True

            if incompleteRun:
                break
            
            firstVar = True
            for i,var in enumerate(clmColnames):

                varDat = clmDat[i,:,:]
                vartot = np.mean(varDat)
                varLowHill = np.mean(varDat[:,lowHillXs])
                varUpHill = np.mean(varDat[:,upHillXs])

                # create DF
                varDF = pd.DataFrame(data=[[vartot,varLowHill,varUpHill]])
                varDF.columns = [var + '_' + area for area in ['total','lower','upper']]

                if firstVar:
                    hourDF = varDF
                    firstVar = False
                else:
                    hourDF = pd.concat([hourDF,varDF],axis=1)
                    
            # write out data
            if firstFile: # write out header
                hourDF.to_csv(clmOutFile, sep=',',header=True,index=False)
                firstFile = False
            else:
                hourDF.to_csv(clmOutFile, sep=',',mode='a', header=False,index=False)      
        