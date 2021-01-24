#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt

# LSQ algorithm for calculating AR parameters
def LSQ(samples,p):
    N = samples.shape[0] - p
    H = np.zeros((N,p))
    for m in range(0,N):
        if m==0:
            H[m] = -samples[m+p-1::-1,0].T
        else:
            H[m] = -samples[m+p-1:m-1:-1,0].T
    
    theta0 = np.linalg.inv( H.T.dot(H) ).dot(H.T).dot(samples[p:N+p,0])
    return theta0.T

# DUTTER algorithm for calculating AR parameters
def DUTTER(samples,p,theta0):
    N = samples.shape[0] - p
    d0 = np.median( np.abs( samples-np.median(samples,axis = 0) ), axis = 0 )/0.6745
    H = np.zeros((N,p))
    theta0 = theta0.reshape(8,1)
    for m in range(0,N):
        if m==0:
            H[m] = -samples[m+p-1::-1,0].T
        else:
            H[m] = -samples[m+p-1:m-1:-1,0].T
    
    n = np.zeros((N,1))
    delta = np.zeros((N,1))
    k = 1.5
    a = 5.6
    a2 = 3
    q = 0.5175 # q = min(1/(2*erf(k)),1.9)
    epsilon = 1e-4
    
    continue_func = 1
    num_iter = 0
    
    while (continue_func and num_iter<10):
        num_iter +=1
        
        # Step 1
        for m in range(N):
            n[m] = samples[m+p] - H[m].dot(theta0).T
        
        # Step 2
        suma = 0
        for m in range(N):
            eps = n[m]/d0
            if np.abs(eps)<=k:
                ksi = eps
            else:
                ksi = k * np.sign(eps)
            suma += ksi**2
        d1 = (d0**2*suma/N/0.7785)**0.5
        
        
        # Step 3
        for m in range(N):
            eps = n[m]/d1
            if eps>k:
                delta[m] = k*d1
            elif eps<-k:
                delta[m] = -k*d1
            else:
                delta[m] = n[m]
        
        # Step 4
        d_theta = np.linalg.inv( H.T.dot(H) ).dot(H.T).dot(delta)
        # Step 5
        theta = theta0 + q*d_theta
        
        # Step 6
        continue_func = 0
        
        for m in range(p):
            if not( (np.abs( theta[m]-theta0[m] ) < np.abs( epsilon*theta0[m] ) ).any() or np.abs( d1-d0 ) < epsilon*np.abs(d0).any() ):
                continue_func=1
        
        d0=d1
        theta0=theta
        
    d = d0
    return theta.reshape(8,), d

# WLSQ algorithm for calculating AR parameters
def WLSQ(samples, p, theta0, d0):
    N = samples.shape[0] - p
    H = np.zeros((N,p))
    
    for m in range(0,N):
        if m==0:
            H[m] = -samples[m+p-1::-1,0].T
        else:
            H[m] = -samples[m+p-1:m-1:-1,0].T
    
    W = np.zeros((N,N))
    a4 = 0.5
    epsilon = 1e-6
    
    continue_func = 1
    num_iter = 0
    
    while (continue_func and num_iter<4):
        num_iter +=1
        suma = 0
        
        for m in range(N):
            eps = (samples[m+p] - H[m].dot(theta0).T)/d0
            
            if samples[m+p] != H[m].dot(theta0).T:
                if np.abs(eps) <= a4*math.pi:
                    ksi = math.sin(eps/a4)
                else:
                    ksi = 0
                suma += d0**2 * ksi**2
                W[m,m] = ksi/eps
            else:
                W[m,m] = 1
               
        d1 = (suma/N/0.561722)**0.5
        theta = np.linalg.inv( H.T.dot(W).dot(H) ).dot(H.T).dot(W).dot(samples[p:N+p,0])
        
        
        # Step 6
        continue_func = 0
        
        for m in range(p):
            if not( (np.abs( theta[m]-theta0[m] ) < np.abs( epsilon*theta0[m] ) ).any() ) :
                continue_func=1
        
        d0=d1
        theta0=theta
        
    return theta.reshape(8,)

# Code for calculating AR coefficients using LSQ, DUTTER and WLSQ algorithm
def CalculateCoefficients(s):
    length_signal = 256
    A_LSQ = np.zeros((s.shape[0],8))    
    A_DUTTER = np.zeros((s.shape[0],8))
    A_WLSQ = np.zeros((s.shape[0],8))
    for m in range(length_signal,s.shape[0]):
        A_LSQ[m]       = LSQ   ( s[m-length_signal:m], 8)
        A_DUTTER[m], d = DUTTER( s[m-length_signal:m], 8, A_LSQ[m].T )
        A_WLSQ[m]      = WLSQ  ( s[m-length_signal:m], 8, A_DUTTER[m].T, d )
        
    return A_LSQ[length_signal:], A_DUTTER[length_signal:], A_WLSQ[length_signal:]

# Code for synthetizing the signal when the excitation is Strube's glottal wave
# both case of noised and non-noised signal is considered
def SynthetizeSignalStrube(A):
    # Parameters of the strube wave excitation
    Tp = 8e-3
    Ts = 3.2e-3
    Tn = 1.2e-3
    Tp_odb = 80
    T_sample = 1e-4
    Glotal = np.zeros((Tp_odb,1))
    for m in range(Tp_odb):
        if (m+1)*T_sample<Ts:
            Glotal[m] = math.sin( math.pi*(m+1)*T_sample/2/Ts )**2
        elif (m+1)*T_sample<Ts+Tn:
            Glotal[m] = math.cos( math.pi*((m+1)*T_sample-Ts)/2/Tn )
        else:
            Glotal[m] = 0
    G_pom = Glotal
    for m in range(20):
        Glotal = np.append(Glotal, G_pom, axis=0)
    Glotal = Glotal[:1100]
    ug_prim = np.append(Glotal, np.array([[0]]) , axis = 0 ) - np.append(np.array([[0]]), Glotal , axis = 0 )
    ug_sek  = np.append(ug_prim, np.array([[0]]), axis = 0 ) - np.append(np.array([[0]]), ug_prim, axis = 0 )
    
    Strube = ug_sek[45:1045]
    Strube = Strube/np.max(Strube)
    s_strube = signal.lfilter([1],A,Strube,axis=0)
    s_pom = s_strube
    
    
    s_strube_orig = s_strube.copy() # Synthetized signal wihout noise
    eps = 0.05 # Probability of outlier occurence
    sig1 = 2e-4
    sig2 = 100*sig1
    pom1 = np.random.rand(1000)
    pom2 = np.random.randn(1000)
    for m in range(s_strube.shape[0]):
        if pom1[m]<eps:
            temp = sig2 * pom2[m]
        else:
            temp = sig1 * pom2[m]
        s_pom[m] = s_strube[m] + temp
    s_strube = s_pom # Synthetized signal with noise
    plt.plot(s_strube_orig)
    plt.xlabel('samples')
    plt.title('Synthetized signal when the excitation is Strube\'s glottal wave')
    plt.show()
    return s_strube_orig, s_strube

# Code for synthetizing the signal when the excitation is a train of Dirac pulses
# both case of noised and non-noised signal is considered
def SynthetizeSignalDirac(A):
    # Parameters of the strube wave excitation
    Ts = 1e-4    
    Tp = 8e-3
    Tp_odb = Tp/Ts
    
    Dirac = np.zeros((1000,1))
    height = 1
    Dirac[1] = height
    for m in range(Dirac.shape[0]):
        if (m+1)%Tp_odb == 0:
            Dirac[m] = height
    
    s_dirac = signal.lfilter([1],A,Dirac,axis=0)
    
    s_pom = s_dirac
    s_dirac_orig = s_dirac.copy()    
    
    eps = 0.05 # Probability of outlier occurence
    sig1 = 2e-4
    sig2 = 100*sig1
    pom1 = np.random.rand(1000)
    pom2 = np.random.randn(1000)
    for m in range(s_dirac.shape[0]):
        if pom1[m]<eps:
            temp = sig2 * pom2[m]
        else:
            temp = sig1 * pom2[m]
        s_pom[m] = s_dirac[m] + temp
        
    s_dirac = s_pom # Synthetized signal with noise
    plt.plot(s_dirac_orig)
    plt.xlabel('samples')
    plt.title('Synthetized signal when the excitation is train of Dirac pulses')
    plt.show()
    return s_dirac_orig, s_dirac