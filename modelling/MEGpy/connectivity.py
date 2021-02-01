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

out = Run(data_uni, data_mul, methods, TH, BH, aux_dir, max_workers)
    Calculate functional connectivity matrices
    inputs:
    data_uni : array (time_points x sources x real data + all univariate surrogates)
    data_mul : array (time_points x sources x real data + all multivariate surrogates)
    methods  : functional connectivity methods MI, MI5, MI5_UNI, CC, oCC, MCC, PLV, iPLV
    TH       : use surrogate correction? True/False
    BH       : use bonferroni correction? True/False
    aux_dir  : res_dur/aux (for intermediate results)
    max_workers

    output
    out      : dictionary out[meth] = array (sources x sources x epochs)

    surrogate correction
    if TH: if real connection i,j is larger than surrogate connections i,j by ranksum test keep connection
        This is why we need to chop the data, so that real connection has 10 values and surrogates have 10*number of surrogates values
    if BH: use Bonferroni-Holms correction of the p-values from the ranksum test

out = CONN(func, data, aux_dir, max_workers)
    input:
    func : functional connectivity function name (see Run())
    data : array (time_points*sub_size x sources x subsegments x real data + surrogates)
    aux_dir  : res_dur/aux (for intermediate results)
    max_workers

    output
    out : array (sources x sources x subsegments x real data + surrogates)

out = MI5_EXEC(data, seed)
    input:
    data : array (time_points x sources)
    seed : random numbers needed for MIhigherdim, but should be chosen outside parallelization
    aux_dir  : res_dur/aux (for intermediate results)
    info   : list [epc, seg, sur] - to save intermediate results

    output:
    out  : MILCA MIhigherdim mutual information

out = MIhigherdim(d1, d2)
    Copyright 2009 Alexander Kraskov, Harald Stoegbauer, Peter Grassberger
    Please Reference:
    A. Kraskov, H. Stogbauer, and P. Grassberger,
    Estimating mutual information.
    Phys. Rev. E 69 (6) 066138, 2004

    Modified by Marinho Lopes @ 2018
    Written for Python Daniel Galvis 2019

    out = MIhigherdim(d1, d2)
    input:
    d1, d2 : two time series
    kneig: k nearest neighbor for MI algorithm (default is 3)
    emb_dim: embedding dimension (default is 1, no embedding)
    emb_tau: time-delay, only relevant for emb_dim > 1 (default is 1)

    out = MILCA MIhigherdim mutual information
    Call the C executable with default parameters and two time series
    Returns mutual information

out = MI_EXEC(data, seed)
    input:
    data : array (time_points x sources)
    seed : random numbers needed for MIhigherdim, but should be chosen outside parallelization

    output:
    out  : lazy proof of concept mutual information (don't use this for production runs!!)

MI = calc_MI(dX, dY, X, Y)
    input:
    dX = first time series
    dY = second time series
    X = histogram for a time series
    Y = histogram for another time series

    output:
    MI : lazy proof of concept mutual information (don't use this for production runs!!)

H = shan_entropy(c)
    input:
    c - a histogram
    output
    H - shannon entropy

out = MCC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: maximum-cross correlation (directed)

out = PLV_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: phase locking value

out = iPLV_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: imaginary phase locking value

out = CC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: amplitude based correlation coefficient

out = oCC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: orthogonalized amplitude based correlation coefficient

xT, yT = orthogonalize(x,y,xP,yP)
    x,y signals (HILBERT TRANSFORMED ALREADY)
    xP,yP pseudoinvers of x,y respectively
    OUTPUT
    xT,yT orthogonalized x,y that is x dot yT = 0
                                     y dot xT = 0

out = correlate(x,y n_pts)
    x, y signals with n_pts points (HILBERT TRANSFORMED ALREADY!!)
    Pearson zero lag correlation of the amplitude of the signals

lag = find_lag(nets)
    This code considers the seg_num true data networks contained in nets[:,:,:,0]
    For a given connection between i,j, this algorithm determines if the connection should be from i->j or j->i

    inputs:
    nets: original network with surrogates
    output:
    lag: a network such that lag[i,j] or lag[j,i] = 1 (i = j), determines if i->j or j->i for directed networks.

network, p = create_network(nets, TH, ALPHA= 0.05)
    inputs:
    nets: original network with surrogates
    TH: True use ranksum test to compare network to surrogates
        False just normalize the network by surrogates
    ALPHA: significance threshold
    outputs:
    network: uncorrected functional connectivity network
    p: list of p values for all considered connections (used if Bonferroni-Holms correction is added)

network = correct_network(nets, network, h)
    inputs:
    nets: the original network with surrogates
    network: the uncorrected network derived from comparing the original network to the surrogates
    h : the significance values after Bonferroni-Holms correction

    outputs:
    network - corrected network

network = wilcoxon_ranksum_directed(nets, TH, BH, ALPHA= 0.05)
    inputs:
    nets : directed network array (sources x sources x subsegments x data + surrogates)
         : for directed networks nets[i,j]> 0 then  nets[j,i]=0  always!!
    TH   : use wilcoxon ranksum test
    BH   : use Bonferroni-Holms correction
    ALPHA: 95 percent or some other percent

    output:
    network : directed and surrogate corrected network (corrected with ranksum test

network = wilcoxon_ranksum_undirected(nets, TH, BH, ALPHA= 0.05)
    inputs:
    nets : undirected network array (sources x sources x subsegments xdata + surrogates)
         : for undirected networks nets[i,j] = nets[j,i]  always!!
    TH   : use wilcoxon ranksum test
    BH   : use Bonferroni-Holms correction
    ALPHA: 95 percent or some other percent

    output:
    network : directed and surrogate corrected network (corrected with ranksum test

h = Bonferroni_Holms(p, ALPHA= 0.05)
    Reference:
    Holm, S. (1979) A simple sequentially rejective multiple test procedure.
    Scandinavian Journal of Statistics. 6, 65-70.
    INPUTS:
    p - list of p-values
    ALPHA - initial significance threshold
    OUPUTS:
    h - list of significances in the same order as the p list
      - 0 not significant, 1 significant
"""

import numpy as np
import os
import subprocess
import scipy
import scipy.fftpack
import scipy.stats
import itertools
from concurrent.futures import ProcessPoolExecutor
import scipy.signal
from scipy.signal import hilbert
from numba import jit
from scipy.linalg import pinv2 as pseudoinverse
from scipy.io import savemat

def Run(data_uni, data_mul, methods, TH, BH, aux_dir, max_workers):
    """
    out = Run(data_uni, data_mul, methods, TH, max_workers)
    Calculate functional connectivity matrices
    inputs:
    data_uni : array (time_points x sources x real data + all univariate surrogates)
    data_mul : array (time_points x sources x real data + all multivariate surrogates)
    methods  : functional connectivity methods MI, MI5, MI5_UNI, CC, oCC, MCC, PLV, iPLV
    TH       : use surrogate correction? True/False
    BH       : use bonferroni correction? True/False
    aux_dir  : res_dur/aux (for intermediate results)
    max_workers

    output
    out      : dictionary out[meth] = array (sources x sources x epochs)

    surrogate correction
    if TH: if real connection i,j is larger than surrogate connections i,j by ranksum test keep connection
        This is why we need to chop the data, so that real connection has 10 values and surrogates have 10*number of surrogates values
    if BH: use Bonferroni-Holms correction of the p-values from the ranksum test

    """

    out = {}
    # iterate over functional connectivity methods
    for meth in methods:
        # undirected methods
        if meth in ['MI','MI5','MI5_UNI','CC','oCC']:
            if meth in ['MI']:
                aux = CONN(MI_EXEC, data_mul, aux_dir, max_workers)
                print('MI complete')
            elif meth in ['MI5']:
                aux = CONN(MI5_EXEC, data_mul, aux_dir, max_workers)
                print('MI5 complete')
            elif meth in ['MI5_UNI']:
                aux = CONN(MI5_EXEC, data_uni, aux_dir, max_workers)
                print('MI5_UNI complete')
            elif meth in ['CC']:
                aux = CONN(CC_EXEC, data_uni, aux_dir, max_workers)
                print('CC complete')
            elif meth in ['oCC']:
                aux = CONN(oCC_EXEC, data_uni, aux_dir, max_workers)
                print('oCC complete')
            # comparing real data to surrogates
            aux2 = np.zeros((aux.shape[0],aux.shape[1],aux.shape[2]))
            for epc in range(aux.shape[2]):
                aux2[:,:,epc] = wilcoxon_ranksum_undirected(aux[:,:,epc,:,:], TH, BH)
        # directed methods
        elif meth in ['MCC', 'PLV', 'iPLV']:
            if meth in ['MCC']:
                aux = CONN(MCC_EXEC, data_uni, aux_dir, max_workers)
                print('MCC complete')

            elif meth in ['PLV']:
                aux = CONN(PLV_EXEC, data_uni, aux_dir, max_workers)
                print('PLV complete')

            elif meth in ['iPLV']:
                aux = CONN(iPLV_EXEC, data_uni, aux_dir, max_workers)
                print('iPLV complete')

            # comparing real data to surrogates
            aux2 = np.zeros((aux.shape[0],aux.shape[1],aux.shape[2]))
            for epc in range(aux.shape[2]):
                aux2[:,:,epc] = wilcoxon_ranksum_directed(aux[:,:,epc,:,:], TH, BH)

        # one of the functional connectivity methods doesn't exist
        else:
            print('one of your connectivity methods doesnt exist returning []')
            aux2 = []


        out[meth] = aux2

    return out

def CONN(func, data, aux_dir, max_workers):
    """
    out = CONN(func, data, max_workers)
    input:
    func : functional connectivity function name (see Run())
    data : array (time_points*sub_size x sources x subsegments x real data + surrogates)
    aux_dir  : res_dur/aux (for intermediate results)
    max_workers

    output
    out : array (sources x sources x subsegments x real data + surrogates)
    """

    nepoch  = len(data)
    nsen    = data[0].shape[1]
    nseg    = data[0].shape[2]
    nsurr   = data[0].shape[3]
    seed    = np.random.randint(1, 2**32-1, (nepoch, nseg, nsurr))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(func, data[epc][:,:,seg,sur], seed[epc, seg, sur], aux_dir, [epc, seg, sur]) for epc, seg, sur in itertools.product(range(nepoch),range(nseg),range(nsurr))]

    idx = np.zeros([int(nepoch*nseg*nsurr),3],dtype=np.int_)
    count = 0
    for epc, seg, sur in itertools.product(range(nepoch),range(nseg),range(nsurr)):
        idx[count,0] = epc
        idx[count,1] = seg
        idx[count,2] = sur
        count = count + 1

    out = np.zeros([nsen, nsen, nepoch, nseg, nsurr])

    count = 0
    for resulter in results:
        out[:,:,idx[count,0],idx[count,1],idx[count,2]] = resulter.result()
        count = count+1

    return out

def MI5_EXEC(data, seed, aux_dir, info):
    """
    out = MI5_EXEC(data, seed)
    input:
    data : array (time_points x sources)
    seed : random numbers needed for MIhigherdim, but should be chosen outside parallelization
    aux_dir  : res_dur/aux (for intermediate results)
    info   : list [epc, seg, sur] - to save intermediate results

    output:
    out  : MILCA MIhigherdim mutual information

    """
    np.random.seed(seed)
    n_sens = data.shape[1]
    out = np.zeros([n_sens,n_sens])


    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            out[ch1,ch2] = np.sqrt(1-np.exp(-2*MIhigherdim(data[:,ch1],data[:,ch2])))
            out[ch2,ch1] = out[ch1,ch2]

    # save intermediate results
    savemat(os.path.join(aux_dir,'aux_'+str(info[0])+'_'+str(info[1])+'_'+str(info[2])+'.mat'),{'out':out})

    return out
def MIhigherdim(d1, d2):
    """
    Copyright 2009 Alexander Kraskov, Harald Stoegbauer, Peter Grassberger
    Please Reference:
    A. Kraskov, H. Stogbauer, and P. Grassberger,
    Estimating mutual information.
    Phys. Rev. E 69 (6) 066138, 2004

    Modified by Marinho Lopes @ 2018
    Written for Python Daniel Galvis 2019

    out = MIhigherdim(d1, d2)
    input:
    d1, d2 : two time series
    kneig: k nearest neighbor for MI algorithm (default is 3)
    emb_dim: embedding dimension (default is 1, no embedding)
    emb_tau: time-delay, only relevant for emb_dim > 1 (default is 1)

    out = MILCA MIhigherdim mutual information
    Call the C executable with default parameters and two time series
    Returns mutual information
    """

    kneig = 3
    emb_dim = 1
    emb_tau = 1

    r = np.random.randint(10e10)
    name = './MEGpy/zwspMIhigh' + str(r) + '.txt'

    f = open(name, 'w')
    text = [ str(d1[i]) + ' ' + str(d2[i])+'\n' for i in range(d1.shape[0]) ]
    f.writelines(text)
    f.close()

    args = ['./MEGpy/MIhigherdim',name, str(2),str(emb_dim),str(emb_tau),str(d1.shape[0]),str(kneig)]
    popen = subprocess.Popen(args, stdout= subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    out = float(output)
    os.remove(name)

    return out



@jit
def MI_EXEC(data, *args):
    """
    out = MI_EXEC(data, seed)
    input:
    data : array (time_points x sources)
    seed : random numbers needed for MIhigherdim, but should be chosen outside parallelization

    output:
    out  : lazy proof of concept mutual information (don't use this for production runs!!)

    """

    n_sens = data.shape[1]
    out = np.zeros([n_sens,n_sens])

    hists = []
    for ch1 in range(n_sens):
        hists.append(np.histogram(data[:,ch1],bins = 'auto'))

    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            out[ch1,ch2] = calc_MI(data[:,ch1],data[:,ch2],hists[ch1],hists[ch2])
            out[ch2,ch1] = out[ch1,ch2]

    return out

def calc_MI(dX,dY,X,Y):
    """
    MI = calc_MI(dX, dY, X, Y)
    input:
    dX = first time series
    dY = second time series
    X = histogram for a time series
    Y = histogram for another time series

    output:
    MI : lazy proof of concept mutual information (don't use this for production runs!!)
    """


    c_X = X[0]
    c_Y = Y[0]

    c_XY = np.histogram2d(dX,dY, bins = (X[1],Y[1]))[0]

    H_XY = shan_entropy(c_XY)
    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)

    MI = np.sqrt( 1 - np.exp(-2*(H_X + H_Y - H_XY)))

    return MI

def shan_entropy(c):
    """
    H = shan_entropy(c)
    input:
    c - a histogram
    output
    H - shannon entropy
    """
    c_normalized = c / float(np.sum(c))

    c_normalized = c_normalized[np.nonzero(c_normalized)]

    H = -np.sum(c_normalized* np.log2(c_normalized))
    return H

def MCC_EXEC(data, *args):
    """
    out = MCC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: maximum-cross correlation (directed)
    """

    n_sens = data.shape[1]
    out = np.zeros([n_sens,n_sens])
    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            xcors = np.correlate(data[:,ch2],data[:,ch1],mode = 'full')
            xcors = np.abs(xcors)/np.sqrt(np.sum(data[:,ch2]**2)*np.sum(data[:,ch1]**2))
            lag = range(xcors.shape[0]) - np.median(range(xcors.shape[0]))
            mx = np.max(xcors)
            tau =  lag[mx == xcors]
            if tau>0:
                out[ch1,ch2] = mx
            elif tau<0:
                out[ch2,ch1] = mx

    return out


@jit
def PLV_EXEC(data, *args):
    """
    out = PLV_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: phase locking value
    """

    n_sens = data.shape[1]

    yh = hilbert(data,axis=0)
    yh_ang = np.angle(yh)
    out = np.zeros([n_sens,n_sens])
    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            delPhi = yh_ang[:,ch1]-yh_ang[:,ch2]
            aux = np.mean(np.exp(1j*delPhi))
            if np.angle(aux) > 0:
                out[ch1,ch2] = np.abs(aux)
            elif np.angle(aux) <0:
                out[ch2,ch1] = np.abs(aux)

    return out

@jit
def iPLV_EXEC(data, *args):
    """
    out = iPLV_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: imaginary phase locking value
    """

    n_sens = data.shape[1]

    yh = hilbert(data,axis=0)
    yh_ang = np.angle(yh)
    out = np.zeros([n_sens,n_sens])
    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            delPhi = yh_ang[:,ch1]-yh_ang[:,ch2]
            aux_angle = np.mean(np.exp(1j*delPhi))
            aux = np.imag(aux_angle)

            if np.angle(aux_angle) > 0:
                out[ch1,ch2] = np.abs(aux)
            elif np.angle(aux_angle) <0:
                out[ch2,ch1] = np.abs(aux)

    return out

@jit
def CC_EXEC(data, *args):
    """
    out = CC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: amplitude based correlation coefficient
    """

    n_pts  = data.shape[0]
    n_sens = data.shape[1]
    data = hilbert(data, axis= 0)

    out = np.zeros([n_sens,n_sens])
    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            x = data[:,ch1]
            y = data[:,ch2]
            aux = correlate(x,y, n_pts)
            out[ch1,ch2] = aux
            out[ch2,ch1] = aux

    return out

def oCC_EXEC(data, *args):
    """
    out = oCC_EXEC(data)
    input:
    data: array (time_points x sources)
    output
    out: orthogonalized amplitude based correlation coefficient
    """

    n_pts  = data.shape[0]
    n_sens = data.shape[1]
    data = hilbert(data, axis= 0)

    pinv   = np.zeros(data.T.shape, dtype= complex)
    for ch in range(n_sens):
        pinv[ch,:] = pseudoinverse(data[:,ch:ch+1])

    out = np.zeros([n_sens,n_sens])
    for ch1 in range(n_sens-1):
        for ch2 in range(ch1+1,n_sens):
            x = data[:,ch1]
            y = data[:,ch2]

            xT, yT = orthogonalize(x,y, pinv[ch1,:], pinv[ch2,:])
            aux2a  = correlate(x,yT, n_pts)
            aux2b  = correlate(xT,y, n_pts)
            aux2   = (aux2a + aux2b) /2.0
            out[ch1,ch2] = aux2
            out[ch2,ch1] = aux2

    return out


def orthogonalize(x, y, xP, yP):
    """
    xT, yT = orthogonalize(x,y,xP,yP)
    x,y signals (HILBERT TRANSFORMED ALREADY)
    xP,yP pseudoinvers of x,y respectively
    OUTPUT
    xT,yT orthogonalized x,y that is x dot yT = 0
                                     y dot xT = 0
    """
    y_pj2x = x*np.dot(xP,y)
    x_pj2y = y*np.dot(yP,x)

    xT = x - x_pj2y
    yT = y - y_pj2x

    xT = np.squeeze(xT)
    yT = np.squeeze(yT)

    return xT, yT


def correlate(x,y, n_pts):
    """
    out = correlate(x,y n_pts)
    x, y signals with n_pts points (HILBERT TRANSFORMED ALREADY!!)
    Pearson zero lag correlation of the amplitude of the signals
    """
    Ax = np.abs(x)
    Ay = np.abs(y)

    meanx = np.mean(Ax)
    meany = np.mean(Ay)

    stdx  = np.std(Ax)
    stdy  = np.std(Ay)

    Ax = (Ax - meanx)
    Ay = (Ay - meany)

    out = scipy.signal.correlate(Ax, Ay, mode='valid')
    out = out/(n_pts*stdx*stdy)

    return np.abs(out)

def find_lag(nets):
    """
    lag = find_lag(nets)
    This code considers the seg_num true data networks contained in nets[:,:,:,0]
    For a given connection between i,j, this algorithm determines if the connection should be from i->j or j->i

    inputs:
    nets: original network with surrogates
    output:
    lag: a network such that lag[i,j] or lag[j,i] = 1 (i = j), determines if i->j or j->i for directed networks.

    """
    # we need to decide whether conn[i,j] or conn[j,i] will be 0 (assuming it is significantly greater than the surrogates)
    lag = np.zeros([nets.shape[0],nets.shape[1]])
    for ch1 in range(nets.shape[0]-1):
        for ch2 in range(ch1+1,nets.shape[0]):
            if np.sum(nets[ch1,ch2,:,0]) > np.sum(nets[ch2,ch1,:,0]):
                lag[ch1,ch2] = 1
            else:
                lag[ch2,ch1] = 1
    return lag
def create_network(nets,TH, ALPHA= 0.05):
    """
    network, p = create_network(nets, TH, ALPHA= 0.05)
    inputs:
    nets: original network with surrogates
    TH: True use ranksum test to compare network to surrogates
        False just normalize the network by surrogates
    ALPHA: significance threshold
    outputs:
    network: uncorrected functional connectivity network
    p: list of p values for all considered connections (used if Bonferroni-Holms correction is added)
    """


    network = np.zeros([nets.shape[0],nets.shape[1]])
    nets = nets.copy()
    p= []
    check_nans = []
    check_surr_nans = []
    # Calculate p-values using the mann-whitney U test
    for ch1 in range(nets.shape[0]-1):
        for ch2 in range(ch1+1,nets.shape[1]):

            # network values over all epochs for true signal
            vals = nets[ch1,ch2,:,0]
            if np.sum(np.isnan(vals)) > 0:
                check_nans.append([np.sum(np.isnan(vals)),ch1,ch2])
                vals = vals[~np.isnan(vals)]

            # network values over all epochs for surrogate signal
            vals_surr = np.ravel(nets[ch1,ch2,:,1:])
            if np.sum(np.isnan(vals_surr)) > 0:
                check_surr_nans.append([np.sum(np.isnan(vals_surr)),ch1,ch2])
                vals_surr = vals_surr[~np.isnan(vals_surr)]

            # p value
            _, p_aux= scipy.stats.mannwhitneyu(vals, vals_surr, alternative= 'greater')
            p.append(p_aux)

            # Initial uncorrected network
            if p_aux < ALPHA:
                network[ch1,ch2] = (np.median(vals)-np.median(vals_surr))/(1-np.median(vals_surr))
                network[ch2,ch1] = network[ch1,ch2]
            # if TH false ignore ranksum test
            if (p_aux >= ALPHA) and (not TH):
                network[ch1,ch2] = np.max((0,(np.median(vals)-np.median(vals_surr))/(1-np.median(vals_surr))))
                network[ch2,ch1] = network[ch1,ch2]

    print('number of nans in data: ', check_nans)
    print('number of nans in surrogates: ', check_surr_nans)
    return network, p

def correct_network(network, h):
    """
    network = correct_network(nets, network, h)
    inputs:
    network: the uncorrected network derived from comparing the original network to the surrogates
    h : the significance values after Bonferroni-Holms correction

    outputs:
    network - corrected network
    """

    network = network.copy()
    # Correct the network
    count = 0
    for ch1 in range(network.shape[0]-1):
        for ch2 in range(ch1+1,network.shape[1]):
            network[ch1,ch2] = network[ch1,ch2]*h[count]
            network[ch2,ch1] = network[ch1,ch2]
            count += 1
    return network

def wilcoxon_ranksum_undirected(nets, TH, BH, ALPHA= 0.05):
    """
    network = wilcoxon_ranksum_undirected(nets, TH, BH, ALPHA= 0.05)
    inputs:
    nets : undirected network array (sources x sources x subsegments xdata + surrogates)
         : for undirected networks nets[i,j] = nets[j,i]  always!!
    TH   : use wilcoxon ranksum test
    BH   : use Bonferroni-Holms correction
    ALPHA: 95 percent or some other percent

    output:
    network : directed and surrogate corrected network (corrected with ranksum test
    """
    print('TH: ', TH)
    print('BH: ', BH)

    # Create a single network from the network values and surrogate network values
    network, p = create_network(nets, TH, ALPHA= ALPHA)

    # bonferroni holms only makes sense if threshold is true
    # return if either is false
    if (not BH) or (not TH):
        return network
    # Re-evaluate significance using multiple comparisons correction
    h = Bonferroni_Holms(p, ALPHA= ALPHA)

    # Bonferroni-Holms correct the network
    network = correct_network(network, h)

    return network

def wilcoxon_ranksum_directed(nets, TH, BH, ALPHA= 0.05):
    """
    network = wilcoxon_ranksum_directed(nets, TH, BH, ALPHA= 0.05)
    inputs:
    nets : directed network array (sources x sources x subsegments x data + surrogates)
         : for directed networks nets[i,j]> 0 then  nets[j,i]=0  always!!
    TH   : use wilcoxon ranksum test
    BH   : use Bonferroni-Holms correction
    ALPHA: 95 percent or some other percent

    output:
    network : directed and surrogate corrected network (corrected with ranksum test
    """

    # we need to decide whether conn[i,j] or conn[j,i] will be 0 (assuming it is significantly greater than the surrogates)
    lag = find_lag(nets)

    # now make the networks undirected (so that we can compare all values for nodes i,j to all surrogates i,j)
    nets = nets + np.transpose(nets,(1,0,2,3))

    network = wilcoxon_ranksum_undirected(nets, TH, BH, ALPHA= ALPHA)

    # make the network directed
    network = network*lag
    return network



def Bonferroni_Holms(p, ALPHA= 0.05):
    """
    h = Bonferroni_Holms(p, ALPHA= 0.05)
    Reference:
    Holm, S. (1979) A simple sequentially rejective multiple test procedure.
    Scandinavian Journal of Statistics. 6, 65-70.
    INPUTS:
    p - list of p-values
    ALPHA - initial significance threshold
    OUPUTS:
    h - list of significances in the same order as the p list
      - 0 not significant, 1 significant

    """
    p = np.array(p.copy())
    idx = np.argsort(p)
    psort = p[idx]
    check = ALPHA/(p.shape[0] + 1 - np.linspace(1,p.shape[0],p.shape[0]))
    try:
        threshold= np.min(np.nonzero(psort>check)[0])

    except Exception:
        threshold = p.shape[0]

    h = np.zeros(p.shape[0])
    if threshold >= 0:
        h[0:threshold] = 1.0

    h[idx] = h.copy()


    return h

