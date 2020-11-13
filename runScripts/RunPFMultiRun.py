
from SC_GenerateParamSet import genParSet
from RunParSet import runSingleFolder
import os

def readInputFile():
    '''
    Reads ParflowMultiRun Input Parameters from SCInput.txt

    Args: 
        none
    Returns:
        dictionary with ParflowMultiRun Input Parameters
    '''

    fn = 'SCInput.txt' # input file (this file name should not change and always be in folder where ParFlowMultiRun is called)
    keys=['runFolder','inputDir','clmDir','parFN',
        'irrFN','parfDir','totaln','nfold','randomseed',
        'saveAllPFData','saveTotStoSL','saveRecCurve_Total', 
        'saveRecCurve_Layers', 'saveCLMSL', 'saveStoStats',
        'saveCLMAll','saveSM','processFlow']
    values = []


    # get key values from inputfile
    with open(fn,'r') as f:
        for line in f:

            currValue = line.split('#')[0]
            currValue = currValue.strip()

            if currValue.isnumeric():
                currValue = int(currValue)
            elif currValue == 'True':
                currValue = True
            elif currValue == 'False':
                currValue = False

            values.append(currValue)

    return dict(zip(keys,values)) #dictionary of key names and values

def createNewRun():

    '''
    Generates new ParFlowMultiRun
    Creates Run Files, and Directory
    Based on Input File, Launches multiple parrallel (in progress) folder runs
    
    Args:
        none
    Returns:
        none
    '''

    # get run 'meta' parameters
    inputData = readInputFile()
    #print(inputData)


    print(inputData)

    # generate run simulation parameters
    genParSet(inputData)
    
    # move into run directory
    os.chdir(inputData['runFolder'])

    # print out a log file with metaparameters
    outfn = 'runMetaPars.txt'
    with open(outfn,'w') as f:
        f.write(str(inputData))

    # enter parallel code here VVVVVVV
    #runSingleFolder(1, inputData) 

def main():
    createNewRun()

if __name__ == "__main__":

    main() # main is passed the line in the parameter file that you want to pull