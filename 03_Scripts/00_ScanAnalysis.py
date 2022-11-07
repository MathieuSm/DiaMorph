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
    Utils.ShowSlice(Scan.reshape((816,1024,1023), order='C'))

    return

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
