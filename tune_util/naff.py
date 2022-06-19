# Implements the NAFF algorithm in python [1]. Based off the MML/AT version
#
# 1. J. Laskar, "Frequency analysis for multi-dimensional systems. Global dynamics and diffusion", Physica D 67, 257281, 1993.
#
# Written by Levon D.
# University of Maryland, Department of Physics
# Oct 2018
#
# Notes. 
# The implementation uses an fft as an initial guess and solves NAFF around (+/- 1 point) the fft guess.
#
#
import numpy as np
from scipy.optimize import fminbound

# Predefined functions
HANNING = lambda t: np.sin(np.pi*np.arange(t)/(t-1))**2 # hanning window filter
FREQ = lambda t: np.arange(t*1.0)/t # norm. freq
EXPO = lambda f,n: np.exp(-1j*2*np.pi*f*n) # complex exponential
NAFFfreq = lambda f,X,W: -1*np.abs(np.sum(EXPO(f,np.arange(1,len(X)+1))*X*W)) #NAFF freq algorithm
NAFFamp = lambda f,X: np.sum(np.real(EXPO(f, np.arange(len(X)))*X))/len(X) + 1j*np.sum(np.imag(EXPO(f, np.arange(len(X)))*X))/len(X) # NAFF amp algorithm

def naff(X, window=[], skipFirst=1):
    '''
        INPUT:
            X - an MxN numpy array. Will calculate tunes across columns, so N tune values.
                NOTE: input data should be a np array
            window - a 2x1 tuple used as the search range [min,max]
            skipFirst - ignore the first f values in the fft
        OUTPUT:
            freq,amp - 1xN array for each (frequency and amplitude)
    '''
    
    # quick check to fix input data being a single array
    if X.ndim == 1:
        X = np.array([X]).T

    # predefined variables
    nturns, ndata = X.shape
    F = FREQ(nturns)
    W = HANNING(nturns)
    
    # ouput variables
    freq = []
    amp = []

    # fft
    Xfft = np.abs(np.fft.fft(X,axis=0))
    
    # choose a starting guess and +/- range to scan
    if window:
        guess = np.mean(window)
        delta = (window[1] - window[0])/2.

        idxrange = np.arange(np.floor(nturns*(guess-delta))-1,np.ceil(nturns*(guess+delta)),dtype=int)
        idxs = (np.argmax(Xfft[idxrange,:],axis=0) + np.floor(nturns*(guess-delta))-1).astype(dtype=int)
    else:
        idxs = np.argmax(Xfft[skipFirst:int(np.floor(nturns/2.)),:],axis=0)+skipFirst
        
    # iterate through each data set and find the freq + amp
    for data,idx in zip(X.T,idxs):
        freq.append(fminbound(NAFFfreq,F[idx-1],F[idx+1],args=(data,W),xtol=1/nturns**4.0))
        amp.append(NAFFamp(freq[-1],data))
        
    return freq,amp

