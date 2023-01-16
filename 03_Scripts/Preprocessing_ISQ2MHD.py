#%% #!/usr/bin/env python3

import argparse
from Utils import *
Read.Echo = False
Write.Echo = False

Version = '01'

# Define the script description
Description = """
    This script aims to write MHD images of scans in order to transfer
    them via Microsoft Onedrive as the original ISQ format becomes
    corrupted somehow.
    
    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel
            University of Bern
    
    Date: January 2023
    """

#%%
# For testing purpose
class ArgumentsClass:

    def __init__(self):
        self.File = 'Example.ISQ'
        return

Arguments = ArgumentsClass()

#%%
def Main(Arguments):

    ProcessTiming(1, 'Convert ISQs to MHDs')

    # Set project directories
    WD, DD, SD, RD = SetDirectories('DiaMorph')
    MHD_Dir = RD / '00_Preprocessing'

    # If provided analyse 1 specific file only
    # Otherwise list all ISQ files present in Data directory (DD)
    try:
        Files = [Arguments.File]
    except:
        Files = [File for File in os.listdir(str(DD)) if File.endswith('.ISQ')]

    # For each listed file
    for i, File in enumerate(Files):

        # Read ISQ scan
        ISQ, AdditionalData = Read.ISQ(str(DD / File))

        # Create corresponding sample folder and write MHD
        os.makedirs(str(MHD_Dir / File[:-4]), exist_ok=True)

        FileName = str(MHD_Dir / File[:-4] / File)
        Write.MHD(ISQ, FileName)

        # Update progress
        ProgressNext((i + 1) / len(Files) * 20)

    ProcessTiming(0)

    return

#%%
if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + Version
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--File', help='File name of the ISQ scan (if only one specific)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)

# %%
