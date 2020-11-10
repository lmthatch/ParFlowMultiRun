# Compile ParFlow MultiRun SC Error Data
# L. Thatch
# 10.26.2020

# %% [markdown]
# ## Pull all SC PF Run Parameters


#%%
 
import pandas as pd
import numpy as np

# import all run parameter data
parDir = '/home/lmthatch/Documents/01_ScalingTests/01_SC/01_ParFlowOnly/FailureTesting/Run2_decreaseResTol/'

# pull parameter files together
import glob
pattern = parDir + 'ParameterSets_AutoGenPY_*.csv'

firstFile = True
for file in glob.glob(pattern):
    fileDat = pd.read_csv(file)
    if firstFile:
        parDat = fileDat
        firstFile = False
        print('FirstFile')
    else:
        parDat = parDat.append(fileDat)



#%% [markdown]

# ## Subset 'Variable' Run Parameters
# And save single file with 'Constant' Parameters

#%%

nUnique = parDat.nunique(axis=0) 
uniqueVals = nUnique[nUnique > 1].index.values.tolist()
nUniqueVals = nUnique[nUnique == 1].index.values.tolist()

#%%
# output 'constant' values
constOutFile = parDir + 'ErrorAnalysis/ConstantVariables.csv'
constantVals = parDat[nUniqueVals].iloc[0]
#constantVals.to_csv(constOutFile,index=True,header=True,sep=',')


#%%

# reduce parameter dataframe for Variable Parameters Only

# Clearer/Shorter Column Names - (This part has to be done manually :( ))
selectKeys = ['n',
    'Geom.domain.ICPressure.Value',
    'Mannings.Geom.domain.Value',
    'Geom.domain.Perm.TensorValZ',
    'Geom.domain.Perm.Value',
    'Patch.z-upper.BCPressure.rain.Value',
    'TopoSlopesX.Geom.domain.Value',
    'ComputationalGrid.DX',
    'ComputationalGrid.NZ',
    'ComputationalGrid.DZ']

newColnames = ['n','GW','Mannings',
    'PermTZ','Perm','Rain','SlopeX',
    'DX','NZ','DZ']

colDict = dict(zip(selectKeys,newColnames))
parDat = parDat[selectKeys]
parDat = parDat.rename(columns=colDict)

# set n as index
parDat = parDat.set_index('n')



#%%

# Extra Clean Up - convert to categorical variables
parDat.DX = parDat.DX.astype('category')
parDat.DY = parDat.DZ.astype('category')
parDat.DZ = parDat.DZ.astype('category')
parDat.NZ = parDat.DZ.astype('category')



#%% [markdown]

# ## Add Success/Failure Indicator for Each Parameter Run

#%%
import glob

# failed runs
errorDir = parDir + 'Errors/'
errorPost =  '*parflow.test.log'
errorList = glob.glob(errorDir + errorPost) # grab all error files
nFail = [f.split('/')[-1] for f in errorList] # extract the run number from the file name
nFail = [int(f.split('_')[0]) for f in nFail] 

# successfull runs
fileDir = parDir + 'SingleLineOutput/'
filePrefix = 'SL_TotStoAllHr_run*'
allFiles = glob.glob(fileDir + filePrefix) # grab all output files
nSuccess = [s.split('run')[1] for s in allFiles] # extract the run number from the file name
nSuccess = [int(s.split('.csv')[0]) for s in nSuccess] 

# couple double checks to make sure failed and successfull are mutually exclusive
print('Runs in Fail and Success? ')
print(any(item in nSuccess for item in nFail))

# check 'success rate' 
print('\nSuccess Rate: ' + str(len(nSuccess)/(len(nFail) + len(nSuccess))))

# update parameter data frame
parDat['Success'] = np.NaN
parDat.Success.loc[nSuccess] = 1
parDat.Success.loc[nFail] = 0

#%%
# export the data
varOutFile = parDir + 'ErrorAnalysis/VariableVariables.csv'
#parDat.to_csv(varOutFile,sep=',')


#%% [markdown]

# ## Get KINSOL Error Data


#%%

parFail = parDat[parDat.Success == 0].copy()

# add new error value columns
parFail['error'] = np.nan
parFail['lastHour'] = np.nan
parFail['nni'] = np.nan
parFail['nform'] = np.nan
parFail['nfe'] = np.nan

# loop through each failed runs and add Data
for n in nFail:

    errorFn = errorDir + str(n) + '_test.out.kinsol.log'

    with open(errorFn, 'r') as ef:
        Lines = ef.readlines()
        for line in Lines: 
            if 'KINSOL starting step' in line:
                lastHourLine = line.strip()
            elif 'KINSol nni=' in line:
                lastSolverLine = line.strip()
            elif '---KINSOL' in line:
                lastError = line.strip()
            elif 'KINSol return' in line:
                lastError = line.strip()

    # strip the time from the file
    parFail.lastHour.loc[n] = float(lastHourLine.split('time')[1])

    # strip solver numbers
    parFail.nni.loc[n] = int(lastSolverLine.split("fnorm=")[0].split('nni=')[1])
    parFail.nform.loc[n] = float(lastSolverLine.split('fnorm=')[1].split('nfe=',)[0])
    parFail.nfe.loc[n] = int(lastSolverLine.split('nfe=')[1])

    # save error
    parFail.error.loc[n] = lastError

print(parFail.error.value_counts())

#%%
# save data
errorOut = parFail.drop(inplace=False,columns = ['GW', 'Mannings', 'PermTZ', 'Perm', 'Rain', 'SlopeX', 'DX', 'NZ', 'DZ',
       'Success'])
errorOutFN = parDir + 'MergedData/KinsolErrors.csv'
errorOut.to_csv(errorOutFN,index=True,sep=',')

#%% [markdown]

# ## Get ParFlow Run Data (Pressure/Saturtion) for Failed Runs

#%% 

frDir = parDir + 'FullRunData/'

# porosity
por = 0.25
ss = 9.80E-05

outMat = np.zeros((len(nFail),302)) * np.NaN

outMatTopP = np.zeros((len(nFail),302)) * np.NaN
outMatBottomP = np.zeros((len(nFail),302)) * np.NaN

for n,i in zip(nFail,range(len(nFail))):
    
    dz = parFail.DZ.loc[n]

    runFN = frDir + 'FullRunData_run' + str(n) + '.csv'
    
    try:
        runDat = pd.read_csv(runFN)
    except:
        print('PROBLEM: ' + str(n))
        next

    # calculate and Save Storage
    pcols = runDat.columns[runDat.columns.str.contains(pat='press')]
    scols = runDat.columns[runDat.columns.str.contains(pat='sat')]

    pdat = runDat[pcols]
    sdat = runDat[scols]

    stoTot = pd.DataFrame(pdat.values*sdat.values * dz * ss)  + pd.DataFrame(sdat.values * por * dz)
    stoTotAvg = stoTot.sum(axis=1) #row sum for each hour

    outMat[i,0:len(stoTotAvg)] = stoTotAvg.to_numpy()

    # get top and bottom pressure values
    bottomP = runDat[pcols[0]]
    topP = runDat[pcols[-1]]

    outMatTopP[i,0:len(topP)] = topP.to_numpy()
    outMatBottomP[i,0:len(bottomP)] = bottomP.to_numpy()


#%%
# Save Total Storage
allStoDat = pd.DataFrame(outMat) # convert to DataFrame
allStoDat['n'] = nFail
allStoDat = allStoDat.set_index('n') # set n as index
outFile = parDir + 'MergedData/TotalStorageData.csv'
allStoDat.to_csv(outFile,index=True) # save data

# Save Top Pressure
topPDat = pd.DataFrame(outMatTopP) # convert to DataFrame
topPDat['n'] = nFail
topPDat = topPDat.set_index('n') # set n as index
outFile = parDir + 'MergedData/TopPressData.csv'
topPDat.to_csv(outFile,index=True) # save data

# Save Bottom Pressure
bottomPDat = pd.DataFrame(outMatBottomP) # convert to DataFrame
bottomPDat['n'] = nFail
bottomPDat = bottomPDat.set_index('n') # set n as index
outFile = parDir + 'MergedData/BottomPressData.csv'
bottomPDat.to_csv(outFile,index=True) # save data



# %%
