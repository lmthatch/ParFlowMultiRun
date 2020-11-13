
#%% [markdown]

# #Subset PFIDB Files


# %%



ns = [1,501,1001,1501,2001,2501,3001,3501,4001,4501,5001,5501,
    6001,6501,7001,7501,8001,8501,9001,9501,10001,10501,
    11001,11501,12001,12501,13001,13501,14001,14501]

outDir = '/home/lmthatch/Documents/88_TEMP/pfidbsForReed/'

#%% [markdown]

# ### pull in all parameter data

# %%

import glob
import pandas as pd

# file directory
#runDir = '/home/lmthatch/Documents/01_ScalingTests/01_SC/01_ParFlowOnly/FailureTesting/Run6_MultiRec/'
runDir = '/home/lmthatch/Documents/01_ScalingTests/01_SC/02_PFCLM/Run1/'

# get parameter data
parPat = runDir + 'ParameterSets_AutoGenPY_*.csv'

firstFile = True
for file in glob.glob(parPat):
    fileDat = pd.read_csv(file)
    if firstFile:
        parDat = fileDat
        firstFile = False
    else:
        parDat = parDat.append(fileDat)


#%% [markdown]

# generate pfidb
parDat = parDat.set_index('n')
parDat = parDat.loc[ns]

# reset index
parDat.reset_index(inplace=True)


# %% 

# fix some clm data

if 'Solver.CLM.MetFilePath' in parDat.columns:
    parDat['Solver.CLM.MetFilePath'] = ''

# %%
from RunParSet import pfidbGen
import os

os.chdir(outDir)

for i in parDat.index:
    parDatSub = parDat.loc[i]
    fn = 'clmtest_dx' + str(parDatSub['ComputationalGrid.DX']) + '_nz' + str(parDatSub['ComputationalGrid.NZ']) + '.pfidb'
    pfidbGen(parDat.loc[i])

    os.system('mv test.pfidb ' + fn)



# %%
