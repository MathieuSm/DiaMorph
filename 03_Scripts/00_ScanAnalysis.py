#%% #!/usr/bin/env python3

import time
import argparse
import numpy as np
from pathlib import Path
from Utils import ISQReader, Utils

Version = '01'

# Define the script description
Description = """
    This script runs a simple analysis of an ISQ scan for the DiaMorph study
    
    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel
            University of Bern
    
    Date: November 2022
    """

#%%
# For testing purpose
class ArgumentsClass:

    def __init__(self):
        self.File = 'C0019219_1.ISQ'
        self.id = str(Path.cwd() / '../02_Data/')
        self.od = str(Path.cwd() / '../04_Results/')
        return

Arguments = ArgumentsClass()

#%%
def Main(Arguments):

    DataDirectory = Arguments.id
    RestultsDirectory = Arguments.od

    Tic = time.time()
    print('\nRead ISQ scan')
    Scan, AdditionalData = ISQReader.Main(Arguments)
    Toc = time.time()
    Utils.PrintTime(Tic, Toc)
    Utils.ShowSlice(Scan, Axis='X')

    return

#%% Tests to reshape

import matplotlib.pyplot as plt

Step = 1024*1024
Build = np.zeros((816,1024,1024))

for i in range(Build.shape[0]):

    VM = VoxelModel[i*Step:(i+1)*Step].reshape((1024,1024))
    Build[i] = np.roll(VM, shift=-4*i, axis=0)

Figure, Axis = plt.subplots(1,1)
Axis.imshow(Build[400])
plt.show()

#%% Tests to rebuild

# Build = VM.copy()

for i in range(Build.shape[0]):
    for j in range(Build.shape[1]):
        Build[i,j:] = np.roll(Build[i,j:],shift=-4, axis=1)

Array = np.roll(Build, shift=(-750,50), axis=(2,1))

#%%
Figure, Axis = plt.subplots(1,1)
Axis.imshow(Array[800])
plt.show()

#%%
ISQReader.WriteMHD(Array,AdditionalData,Arguments.od,Arguments.File[:-4],'float')

#%%
if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add long and short argument
    ScriptVersion = Parser.prog + ' version ' + Version
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('File', help='File name of the ISQ scan (required)', type=str)

    # Define paths
    InputDirectory = str(Path.cwd() / '../02_Data/')
    OutputDirectory = str(Path.cwd() / '../04_Results/')
    Parser.add_argument('-id', help='Set input directory', type=str, default=InputDirectory)
    Parser.add_argument('-od', help='Set output directory', type=str, default=OutputDirectory)


    # Read arguments from the command line
    Arguments = Parser.parse_args()

# %%
