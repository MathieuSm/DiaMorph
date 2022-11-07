import numpy as np
import matplotlib.pyplot as plt

def PrintTime(Tic, Toc):

    """
    Print elapsed time in seconds to time in HH:MM:SS format
    :param Tic: Actual time at the beginning of the process
    :param Toc: Actual time at the end of the process
    """

    Delta = Toc - Tic

    Hours = np.floor(Delta / 60 / 60)
    Minutes = np.floor(Delta / 60) - 60 * Hours
    Seconds = Delta - 60 * Minutes - 60 * 60 * Hours

    print('Process executed in %02i:%02i:%02i (HH:MM:SS)' % (Hours, Minutes, Seconds))

    return

def ShowSlice(Array, Slice=None, Title=None, Axis='Z'):
    
    if Axis == 'Z':
        if Slice:
            Array = Array[Slice,:,:]
        else:
            Array = Array[Array.shape[0]//2,:,:]
    if Axis == 'Y':
        if Slice:
            Array = Array[:,Slice,:]
        else:
            Array = Array[:,Array.shape[1]//2,:]
    if Axis == 'X':
        if Slice:
            Array = Array[:,:,Slice]
        else:
            Array = Array[:,:,Array.shape[2]//2]

    Figure, Axis = plt.subplots()
    Axis.imshow(Array,interpolation=None, cmap='binary_r')
    Axis.axis('Off')
    
    if(Title):
        Axis.set_title(Title)

    plt.show(Figure)

    return