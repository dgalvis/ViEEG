"""
All code was created by Daniel Galvis except where otherwise noted (MEGmat/BNI_find.m, MEGmat/BNI_single.m MEGmat/NI_model.m, MEGmat/theta_model.m, MEGmat/bonf_holm.m, MEGpy/MutualInformationICA/*, and MEGpy/connectivity.py->out=MIhigherdim(d1,d2))

Copyright (C) 2019 Daniel Galvis
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

preprocessing of time series data

data, time = downsample(data_in, time_in, f_or, f_new)
    input:
    data_in : time_points x sources
    time_in : time_points
    f_or    : original sampling frequency in Hz
    f_new   : new sampling frequency

    output:
    data_in : downsample time series x sources
    time_in : downsample time points

data = butterworth_filter(data_in, f_new, f_min, f_max, btype)
    input:
    data_in : time_points x sources
    f_new   : sampling frequency (this is usually done after resample hence "new")
    f_min   : low frequency
    f_max   : high frequency
    btype= 'band' (bandpass) or 'bandstop'
    output:
    4th order butterworth filtered data
    data : time_points x source

data_out, time_out = time_segment(data_in, time_in, low_range, high_range)
    input:
    data_in    : array time_points x sources
    time_in    : time points
    low_range  : list of lower bounds for epochs
    high_range : list of upper bounds for epochs

    output:
    data_out   : list of arrays, each array is time_points in [ low_range[idx], high_range[idx] ] x sources
    time_out   : list of arrays, each array is time_points in [ low_range[idx], high_range[idx] ] 

out = generate_iAAFT(data_in, nsurr, max_workers)
    Univariate Surrogates
    inputs:
    data_in : time_points x sources
    nsurr   : number of surrogates to create
    max_workers

    output:
    out   : array time_points x sources x nsurr (

sn = generate_iAAFT_start(data_in, perms)
    single univariate surrogate for a single time series
    inputs:
    data_in : array (time_points, )
    perms   : array (time_points, ) - a random permuatation of data_in

    output:
    sn      : surrogate (exact amplitude spectrum as original data, close signal distribution, mixed up phase spectrum

out = generate_multivariate_iAAFT(data_in, nsurr, max_workers)
    Multivariate Surrogates
    inputs:
    data_in : time_points x sources
    nsurr   : number of surrogates to create
    max_workers

    output:
    out   : array time_points x sources x nsurr

sn = generate_multivariate_iAAFT_start(X, perm)
    Generate surrogate data with matching amplitude spectrum and
    amplitude distribution (Schreiber and Schmitz, 1996).
    Multivariate extension as described in Schreiber and Schmitz (2000)

    input:
    X : array (time_points x sources)
    perm : array (time_points x source) random permutations along of the time points

    output:
    sn : array (time_points x source) - exact amplitude spectrum, close signal distribution,
"""
import numpy as np
import scipy.signal
import scipy.fftpack
import itertools
from concurrent.futures import ProcessPoolExecutor
from numba import jit

def downsample(data_in, time_in, f_or, f_new):
    """
    data, time = downsample(data_in, time_in, f_or, f_new= 512)
    input:
    data_in : time_points x sources
    time_in : time_points
    f_or    : original sampling frequency in Hz
    f_new = 512 : new sampling frequency

    output:
    data_in : downsample time series x sources
    time_in : downsample time points
    """
    data = data_in.copy()
    time = time_in.copy()

    count= 0
    # Time series which are a multiple of 10 run faster
    while np.mod(data.shape[0], 10) != 0:
        data   = data[:-1,:]
        time   = time[:-1]
        count += 1
    print('Removed ', str(count), ' points so that data has factor of 10 points (for speed!)')

    # resample
    if(f_or != f_new):
        old_samp = data.shape[0]
        new_samp = int(old_samp*f_new/f_or)
        data = scipy.signal.resample(data,new_samp,axis=0)
        time = np.linspace(time[0],time[-1],data.shape[0])

    print('Downsample complete: ', str(f_new), ' Hz')

    return data, time
# 4th order butterworth filter
def butterworth_filter(data_in, f_new, f_min, f_max, btype):
    """
    data = butterworth_filter(data_in, f_new, f_min= 0.5, f_max= 150, btype= 'band')
    input:
    data_in : time_points x sources
    f_new   : sampling frequency (this is usually done after resample hence "new")
    f_min   : low frequency
    f_max   : high frequency
    btype= 'band' (bandpass) or 'bandstop'
    output:
    4th order butterworth filtered data
    data : time_points x source
    """
    data = data_in.copy()
    if(f_min != f_max):
        Norder = 4
        Wn = np.array([f_min,f_max])/(f_new/2.0)
        b,a = scipy.signal.butter(Norder,Wn,btype=btype)
        data = scipy.signal.filtfilt(b,a,data,axis=0)
    print('Butterworth Filter: [', f_min, ' ', f_max,'] Hz of type',btype)

    return data

def time_segment(data_in, time_in, low_range, high_range):
    """
    data_out, time_out = time_segment(data_in, time_in, low_range, high_range)
    input:
    data_in    : array time_points x sources
    time_in    : time points
    low_range  : list of lower bounds for epochs
    high_range : list of upper bounds for epochs

    output:
    data_out   : list of arrays, each array is time_points in [ low_range[idx], high_range[idx] ] x sources
    time_out   : list of arrays, each array is time_points in [ low_range[idx], high_range[idx] ] 

    """

    data_out = []
    time_out = []
    # iterate over the list of epochs
    for idx, _ in enumerate(low_range):
        lrange = low_range[idx]
        hrange = high_range[idx]

        # actual time points closest to low_range[idx] and high_range[idx]
        lidx = np.argmin(np.abs(time_in - lrange))
        hidx = np.argmin(np.abs(time_in - hrange))


        # Time series which are a multiple of 10 run faster
        while np.mod(data_in[lidx:hidx,:].shape[0], 10) != 0:
            hidx = hidx - 1

        data_out.append(data_in[lidx:hidx,:])
        time_out.append(time_in[lidx:hidx])

    # Print the time intervals from the config file (just checking it works properly)
    for tidx, _ in enumerate(time_out):
        print('Segment: [', "%.2f" %(time_out[tidx][0]),' ', "%.2f" %( time_out[tidx][-1]),'] s')

    return data_out, time_out

def generate_iAAFT(data_in, nsurr, max_workers):
    """
    out = generate_iAAFT(data_in, nsurr, max_workers)
    Univariate Surrogates
    inputs:
    data_in : time_points x sources
    nsurr   : number of surrogates to create
    max_workers

    output:
    out   : array time_points x sources x nsurr (

    """
    print('Generating ', nsurr, ' univariate surrogates')
    data = data_in.copy()

    # number of epochs
    nepoch = len(data)

    # Define point permutations for all sensors, epochs, and surrogates
    # Do this outside parallel processing to get unique values
    # This is random permutations of the time series for each source for each surrogate
    perms  = []
    out    = []
    for epc in range(nepoch):
        npt  = data[epc].shape[0]
        nsen = data[epc].shape[1]
        perm_aux = np.zeros([npt, nsen, nsurr], dtype = int)
        out_aux  = np.zeros([npt, nsen, nsurr])
        for sen in range(nsen):
            for sur in range(nsurr):
                perm_aux[:, sen, sur] = np.random.permutation(npt)
        perms.append(perm_aux)
        out.append(out_aux)


    # generate each surrogate
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(generate_iAAFT_start,data[epc][:,sen],perms[epc][:,sen,sur])
                   for epc, sen, sur in itertools.product(range(nepoch), range(nsen), range(nsurr))]


    # extract surrogate idxs
    idx = np.zeros([int(nsen*nepoch*nsurr),3],dtype=int)
    count = 0
    for epc, sen, sur in itertools.product(range(nepoch), range(nsen), range(nsurr)):
        idx[count,0] = epc
        idx[count,1] = sen
        idx[count,2] = sur
        count = count + 1

    # extract surrogates
    count = 0
    for resulter in results:
        out[idx[count,0]][:,idx[count,1],idx[count,2]] = resulter.result()
        count += 1


    return out

def generate_multivariate_iAAFT(data_in, nsurr, max_workers):
    """
    out = generate_multivariate_iAAFT(data_in, nsurr, max_workers)
    Multivariate Surrogates
    inputs:
    data_in : time_points x sources
    nsurr   : number of surrogates to create
    max_workers

    output:
    out   : array time_points x sources x nsurr

    """

    print('Generating ', nsurr, ' multivariate surrogates')
    data = data_in.copy()

    # number of epochs
    nepoch = len(data)

    # Define point permutations for all sensors, epochs, and surrogates
    # Do this outside parallel processing to get unique values
    # This is random permutations of the time series for each source for each surrogate
    perms  = []
    out    = []
    for epc in range(nepoch):
        npt  = data[epc].shape[0]
        nsen = data[epc].shape[1]
        perm_aux = np.zeros([npt, nsen, nsurr], dtype = int)
        out_aux  = np.zeros([npt, nsen, nsurr])
        for sen in range(nsen):
            for sur in range(nsurr):
                perm_aux[:, sen, sur] = np.random.permutation(npt)
        perms.append(perm_aux)
        out.append(out_aux)


    # calculate multivariate surrogates
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(generate_multivariate_iAAFT_start,data[epc],perms[epc][:,:,sur])
                   for epc, sur in itertools.product(range(nepoch), range(nsurr))]

    # get idxs
    idx = np.zeros([int(nepoch*nsurr),2],dtype=np.int_)
    count = 0
    for epc, sur in itertools.product(range(nepoch),range(nsurr)):
        idx[count,0] = epc
        idx[count,1] = sur
        count = count + 1

    # retrieve the surrogate results
    count = 0
    for resulter in results:
        out[idx[count,0]][:,:,idx[count,1]] = resulter.result()
        count = count+1


    return out


def generate_iAAFT_start(data_in,perms):
    """
    sn = generate_iAAFT_start(data_in, perms)
    single univariate surrogate for a single time series
    inputs:
    data_in : array (time_points, )
    perms   : array (time_points, ) - a random permuatation of data_in

    output:
    sn      : surrogate (exact amplitude spectrum as original data, close signal distribution, mixed up phase spectrum

    """
    data = data_in.copy()

    max_it = 10000
    sn = np.zeros([data.shape[0]])


    # Initial conditions
    rn = data[perms]
    Xsorted = np.sort(rn) # desired signal spectrum
    Yamp = np.abs(scipy.fftpack.fft(data)) # desired amplitude spectrum

    c = 1
    prev_err = 1e10
    err = prev_err - 1
    while ((prev_err>err) and (c<max_it)):
        # match amplitude spectrum
        Yrn = scipy.fftpack.fft(rn)
        Yang = np.angle(Yrn)
        sn = np.real(scipy.fftpack.ifft(Yamp*np.exp(1j*Yang)))

        # scale to original signal distribution
        idx = np.argsort(sn)
        rn[idx] = Xsorted

        # Eval convergence
        prev_err = err
        A2 = np.abs(Yrn)
        err = np.mean(np.abs(A2-Yamp))

        c = c+1


    return sn


def generate_multivariate_iAAFT_start(X, perm):
    """

    sn = generate_multivariate_iAAFT_start(X, perm)
    Generate surrogate data with matching amplitude spectrum and
    amplitude distribution (Schreiber and Schmitz, 1996).
    Multivariate extension as described in Schreiber and Schmitz (2000)

    input:
    X : array (time_points x sources)
    perm : array (time_points x source) random permutations along of the time points

    output:
    sn : array (time_points x source) - exact amplitude spectrum, close signal distribution,

    """
    max_it = 10000
    pp    = X.shape[0]
    dim   = X.shape[1]
    Y     = scipy.fftpack.fft(X, axis=0)
    Yamp  = np.abs(Y)
    Porig = np.angle(Y)

    rn = np.zeros(X.shape)

    for k in range(dim):
        rn[:, k] = X[perm[:,k], k]

    Xsorted = np.sort(X, axis= 0)

    prev_err = 10e10
    err      = prev_err -1

    c= 1
    Pcurr = Porig
    while (prev_err > err) and (c<max_it):
        # match amplitude spectrum
        Prn   =  np.angle(scipy.fftpack.fft(rn, axis= 0))
        goal  =  Prn - Porig
        AUX1  =  np.expand_dims(np.sum(np.cos(goal), axis= 1), axis=1)
        AUX2  =  np.expand_dims(np.sum(np.sin(goal), axis= 1), axis=1)
        alpha =  np.tile((AUX1 != 0)*np.arctan(AUX2/(AUX1 + (AUX1 == 0))), [1, dim])
        AUX3  =  np.expand_dims(np.pi*(np.sum(np.cos(alpha-goal), axis= 1)<0), axis= 1)
        alpha =  alpha + np.tile(AUX3, [1,dim])
        Pcurr =  Porig + alpha

        # match signal distribution
        sn    = np.real(scipy.fftpack.ifft(Yamp*(np.cos(Pcurr) + 1j*np.sin(Pcurr)), axis= 0))
        INDS  = np.argsort(sn, axis= 0)

        for k in range(dim):
            rn[INDS[:,k], k] = Xsorted[:,k]

        # evaluate convergence
        prev_err = err
        A2       = np.abs(scipy.fftpack.fft(rn, axis= 0))
        err      = np.mean(np.abs(A2 - Yamp))
        c = c+1

    return sn

def chop_segment(data_in, nseg, sub_size):
    """
    out = chop_segment(data_in, nseg= 10, sub_size= 0.25)
    input:
    data_in : list of arrays for each epoch, arrays: (time_points x sources x real_data + surrogates)
    nseg= int: cut into nseg minimally overlapping segments
    sub_size= a fraction: of sub_size*time_points size

    output
    out   : list of arrays for each epoch, arrays: (sub_size*time_points x sources x nseg x real_data + surrogates)

    """
    print('chopping segment into ', nseg,' segments of fraction ', sub_size)
    data   = data_in.copy()
    nepoch = len(data)

    out = []
    # iterate over epochs
    for epc in range(nepoch):
        data_aux = data[epc]
        # points per subsegment
        npt      = int(data_aux.shape[0]*sub_size)
        # slide for minimally overlapping
        slide    = int(np.floor((data_aux.shape[0] - npt)/(nseg - 1)))

        nsen     = data_aux.shape[1]
        nsurr    = data_aux.shape[2]
        out_aux  = np.zeros((npt, nsen, nseg, nsurr))

        for sg in range(nseg):
            out_aux[:,:,sg,:] = data_aux[(slide*sg) : (slide*sg + npt), :, :]

        out.append(out_aux)

    return out

