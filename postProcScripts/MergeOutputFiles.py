
# %% [markdown]

# # Merge Curve Fits
# L. Thatch
#
# created: 11.4.2020
#
# updated: 11.7.2020
#   
# %%
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% [markdown]

# # Merge PF Run Parameter Data

# ### Import Parameter Data

# %%

# file directory
runDir = '/home/lmthatch/CheyMount/'

# merged file out directory
outDir = runDir + 'MergedData/'

# get parameter data
parPat = runDir + 'origPars/ParameterSets_AutoGenPY_*.csv'

firstFile = True
for file in glob.glob(parPat):
    fileDat = pd.read_csv(file)
    if firstFile:
        parDat = fileDat
        firstFile = False
    else:
        parDat = parDat.append(fileDat)

#%% [markdown]

# ### Clean up Parameter Dataframe
# %%

# get constant parameters
nUnique = parDat.nunique(axis=0)
constantVals = nUnique[nUnique > 1].index.values.tolist()
constDat = parDat[constantVals]

# save constant values
outFN = outDir + 'ConstantParameters.csv'
constDat.to_csv(outFN)

#%% 

# get/cleanup variable parameters
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

parDat = parDat.set_index('n')

# %% [markdown]

# ### Add parameter sets references 
# (aka same porosity/perm/rain) but different dz/dx

# %%
uniquePars = parDat.drop_duplicates(subset=['GW','Mannings',
    'PermTZ','Perm','Rain','SlopeX'])
uniquePars['nsub'] = np.arange(len(uniquePars))

parDat['nsub'] = np.NaN

for i in range(len(parDat)):
    gw = parDat.GW.iloc[i]
    perm = parDat.Perm.iloc[i]
    parDat.nsub.iloc[i] = uniquePars[(uniquePars['GW'] == gw) & (uniquePars['Perm'] == perm)].nsub.iloc[0]

#%% save variable parameter data
outFN = outDir + 'VariableParameters.csv'
parDat.to_csv(outFN)


#%% [markdown]

# # Merge Curve Fit Data

#%%

fileDir = runDir + 'SingleLineOutput/'

curveFitPrefixs = [('01m','SL_StoRecCurveFit_01m_run'),
    ('1-2m','SL_StoRecCurveFit_1-2m_run'),
    ('1m','SL_StoRecCurveFit_1m_run'),
    ('2-10m','SL_StoRecCurveFit_2-10m_run'),
    ('2m','SL_StoRecCurveFit_2m_run'),
    ('TotalSto','SL_TotStoRecCurveFit_run')
    ]

# loop through all runs
firstFile = True
for n in parDat.index:
    print(n)
    for ctype,prefix in curveFitPrefixs:
        fn = fileDir + prefix + str(n) + '.csv'
    
        try:
            fileDat = pd.read_csv(fn)
        except:
            print('no curve for: ' + str(n))
            continue
    
        fileDat['n'] = n
        fileDat['CurveType'] = ctype

        if firstFile:
            curveDat = fileDat
            firstFile = False
        else:
            curveDat = curveDat.append(fileDat)
        
# %%

# save data
curveDat = curveDat.set_index('n')
outFN = runDir + 'MergedData/RecCurveFitMerge.csv'
curveDat.to_csv(outFN)


# %% [markdown]

# # Get Top 1M and 2M SM changes

#%%

depths = (1,2)
fns = (outDir + 'Top1mSM_Merged.csv',outDir + 'Top2mSM_Merged.csv')

for depth,fn in zip(depths,fns):
    
    with open(fn,'w') as f:
        
        firstFile = True

        for n in range(5):

            # domain parameters
            por = 0.25
            dz = parDat.DZ.loc[n]

            # evaluate which layers to pull
            nlay = 10/dz
            nlayperM = 10/dz/10
            layerNs = np.arange(nlay-nlayperM*depth,nlay) + 1
            layers = ['sat_' + str(int(n)) for n in layerNs]

            # get run data
            fn = runDir + 'FullRunData/FullRunData_run' + str(n) + '.csv'
            frDat = pd.read_csv(fn)

            # calculate soil moisture
            avgSM = frDat[top1MLays].mean(axis=1) * por

            # write out data
            if firstFile:
                f.write(str(n)+',')
                firstFile=False
            else:
                f.write("\n")
                f.write(str(n)+',')

            f.write(".".join([str(i) for i in avgSM]))

            
 