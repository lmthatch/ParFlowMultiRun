
# Converted from 'SC_GenerateParameterSets.R'

def genParSet(inputData):
   '''
   Generates parameter sets based on Input File
   Files with each of the parameter sets are written into the test directory (testDir)

   Args:
      InputData (dict): dictionary with ParFlowMultiRun keys, including
         n: number of parameter sets to generate
         inputFile: file name w/ parflow parameter values
         testDir: folder where ParFlow MultiRun will run
         numFold: number of folders ParFlow MultiRun Running in (number of files to split parameters in)
   Returns:
      none

   '''

   import pandas as pd
   import numpy as np
   import random
   import sys
   import os

   print('Inside genParSet Function')
   # needed dictionary keys
   n=inputData['totaln']
   inputFile=inputData['parFN']
   testDir=inputData['runFolder']
   numFold=inputData['nfold']
   randseed=inputData['randomseed']

   print('got all my data sets')

   # create parameter data sets
   parDataSets={'n': list(range(0,n))}
   parDF=pd.DataFrame(parDataSets)

   # load csv file of keys and values
   print('Name of input file: ' + inputFile)
   inputVarDF=pd.read_csv(inputFile, delimiter=',')
   inputVarDF.columns = ['KeyName', 'SCValue', 'inputType', 'set', 'MinRange', 'MaxRange']

   # evaluate variable parameters
   varVars = inputVarDF[inputVarDF['set'] == 'Variable']
   varNames = np.unique(varVars['SCValue']) # some variables may be set together, so read unique names

   random.seed(randseed) # set random seed for reproducibility
   
   # loop through each of the variables and create a set (size n) of random values
   for var in varNames: 
      varIdx = varVars[varVars['SCValue']==var]
      if varIdx['inputType'].iat[0] =='double': 
         varmin = int(float(varIdx['MinRange'].iat[0])*1000*n)
         varmax = int(float(varIdx['MaxRange'].iat[0])*1000*n)
         varValues = np.array(random.sample(range(varmin,varmax),n))/(1000*n)
      elif varIdx['inputType'].iat[0] =='integer':
         varmin = int(varIdx['MinRange'].iat[0])
         varmax = int(varIdx['MaxRange'].iat[0])
         intRange = [format(x, '01d') for x in range(varmin,varmax)]
         varValues = np.array(random.choices(intRange,k=n))
         varValues = varValues.astype(int)
      else:
         varPos = varIdx['MinRange'].iat[0].split(sep=" ")
         varValues = np.array(random.choices(varPos,k=n))
      for name in varIdx['KeyName']:
         parDF[name]=varValues

   # add constant variables, these won't need to be changed, can be read straight in
   constantVars = inputVarDF[inputVarDF['set'] == 'Constant']
   constantKeys = constantVars['KeyName']

   # add them to data frame with variable parameters
   for key in constantKeys:
      keyRow = constantVars[constantVars['KeyName']==key]
      keyValue = keyRow['SCValue'].iat[0]
      parDF[key] = [keyValue] * n

   # calculate variables
   calcVars = inputVarDF[inputVarDF['set'] == 'Calculate']

   ##### Irrigation Variables #### *** NEEDS WORK, MESSY.
   if 'Solver.CLM.IrrigationStartTime' in inputVarDF.KeyName.values:
      
      # randomly generate volumes, durrations, and timings
      irrVols = [2.5,5,10]
      irrDurrs = [1,2.5,5,10]
      irrTiming = ['morn','day','eve','night']

      irrVolValues = np.array(random.choices(irrVols,k=n))
      irrDurrValues = np.array(random.choices(irrDurrs,k=n))
      irrTimingValues = np.array(random.choices(irrTiming,k=n))

      # convert to irrigation start and stop times
      irrTimingDF = pd.read_csv('SC_IrrigationTiming_20200813.csv')
      irrStart = []
      irrEnd = []

      for i in range(len(irrDurrValues)):
         dfsub = irrTimingDF[(irrTimingDF['durration'] == irrDurrValues[i]) & (irrTimingDF['timing'] == irrTimingValues[i])] 
         irrStart.append(dfsub['start'].iloc[0])
         irrEnd.append(dfsub['end'].iloc[0])

      parDF['Solver.CLM.IrrigationStartTime'] = irrStart
      parDF['Solver.CLM.IrrigationStopTime'] = irrEnd
      irrRate = irrVolValues/irrDurrValues/3600
      parDF['Solver.CLM.IrrigationRate'] = irrRate
   
   # Quick Fidelity Check
   # make sure clm soil layer not greater than the number of layers
   # first check if these variables exist
   if (('Solver.CLM.SoiLayer' in inputVarDF.KeyName.values) & ('Solver.CLM.RootZoneNZ' in inputVarDF.KeyName.values)):
      clmLayDif = parDF['Solver.CLM.SoiLayer'] - parDF['Solver.CLM.RootZoneNZ'] #check difference between soilLayer and RootZoneNZ 
      #parDF['Solver.CLM.SoiLayer'][clmLayDif > 0] = parDF['Solver.CLM.RootZoneNZ'][clmLayDif > 0] # if soilLayer > RootZoneNZ then replace with the RootZoneNZ
      parDF['Solver.CLM.SoiLayer'] = [parDF['Solver.CLM.SoiLayer'][i] if clmLayDif[i] < 0 else parDF['Solver.CLM.RootZoneNZ'][i] for i in range(len(parDF))]

   # calculate variables
   # calculate variables should be in the same format [operator, var1, var2]
   # operator can be [divide, multiply, add, or subtract]
   
   # loop through each calc Var
   for i in range(len(calcVars)):
      currRow = calcVars.iloc[i,:]
      allargs = currRow.SCValue.split(" ")
      op = allargs[0]
      firstArg = allargs[1]

      # check if first argument is a string or a number
      try:
         firstValue = int(firstArg)
         parValue = pd.Series(np.repeat(firstValue,n))
      except:
         parValue = parDF[firstArg]

         # check variable type & convert
         vartype = inputVarDF[inputVarDF['KeyName']==firstArg].inputType
         if vartype.iat[0] == 'double':
            parValue = parValue.astype('float')

         else:
            parValue = parValue.astype('int')

      for j in range(len(allargs)-2):
          currarg = allargs[j + 2]
          currValue = parDF[currarg]

          # check variable type & convert
          vartype = inputVarDF[inputVarDF['KeyName']==currarg].inputType
          if vartype.iat[0] == 'double':
              currValue = currValue.astype('float')
          else:
              currValue = currValue.astype('int')

          if op == 'divide':
              parValue = parValue.divide(currValue)
          elif op == 'multiply':
              parValue = parValue.multiply(currValue)
          elif op == 'add':
              parValue = parValue.add(currValue)
          elif op == 'subtract':
              parValue = parValue.subtract(currValue)
      parDF[currRow.KeyName] = parValue 

   # create new test directory
   os.system('mkdir -p ' + testDir)
   os.system('mkdir -p ' + testDir + '/FullRunData')
   os.system('mkdir -p ' + testDir + '/SingleLineOutput')
   os.system('mkdir -p ' + testDir + '/RunTimeData')
   os.system('mkdir -p ' + testDir + '/Errors')

   # check if we need CLM files, copy if necessary
   if parDF['Solver.LSM'].iloc[0] == 'CLM':
      try:
         os.system('cp -r clm_input/ ' + testDir)
      except:
         print('Missing CLM Files')

   # write out parameter data set variables to number of files
   numFold = int(numFold)

   if numFold == 1:
      outfn = testDir + '/ParameterSets_AutoGenPY_0.csv'
      parDF.to_csv(outfn,index=False)
   else:
      filesPerOutF = n/numFold
      for nf in range(numFold):

         # cut the out parameter set
         fileStart = nf * filesPerOutF
         fileEnd = (nf+1) * filesPerOutF - 1
         parDF_sub = parDF.loc[fileStart:fileEnd,:]

         # save the data
         outfn = testDir + '/ParameterSets_AutoGenPY_' + str(nf) + '.csv'
         parDF_sub.to_csv(outfn,index=False)

def main(args):
   keyNames = ['totaln','parFN','runFolder','nfold','randomseed']
   keyValues = [int(args[0]),args[1],args[2],args[3],args[4]]

   genParSet(dict(zip(keyNames,keyValues)))

if __name__ == "__main__":

    import sys
    main(sys.argv[1:]) # main will take in command line input for genParSet parameters

