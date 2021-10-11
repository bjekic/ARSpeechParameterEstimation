import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt
from tqdm import tqdm


def LSQ_algorithm(samples: np.ndarray, p: int) -> np.ndarray:
    """
    LSQ algorithm for calculating AR parameters
    :param samples: the samples of the signal
    :param p: AR model order
    :return: calculated theta parameters using LSQ method
    """

    N = samples.shape[0] - p
    H = np.zeros((N, p))
    for m in range(N):
        if m == 0:
            H[m] = -samples[m + p - 1::-1, 0].T
        else:
            H[m] = -samples[m + p - 1:m - 1:-1, 0].T
    
    theta0 = np.linalg.inv(H.T.dot(H)).dot(H.T).dot(samples[p:N + p, 0])

    return theta0.T


def DUTTER_algorithm(samples: np.ndarray, p: int, theta0: np.ndarray) -> tuple:
    """
    DUTTER algorithm for calculating AR parameters
    :param samples: the samples of the signal
    :param p: AR model order
    :param theta0: starting theta parameters
    :return: tuple of calculated theta parameters using Dutter method and parameter d
    """

    N = samples.shape[0] - p
    H = np.zeros((N, p))
    d0 = np.median(np.abs(samples - np.median(samples, axis=0)), axis=0) / 0.6745
    theta0 = theta0.reshape(8, 1)

    for m in range(N):
        if m == 0:
            H[m] = -samples[m + p - 1::-1, 0].T
        else:
            H[m] = -samples[m + p - 1:m - 1:-1, 0].T
    
    n = np.zeros((N, 1))
    delta = np.zeros((N, 1))

    k = 1.5
    q = 0.5175  # q = min(1 / (2 * erf(k)), 1.9)
    epsilon = 1e-4
    
    continue_func = 1
    num_iter = 0
    
    while continue_func and num_iter < 10:
        num_iter += 1
        
        # Step 1
        for m in range(N):
            n[m] = samples[m + p] - H[m].dot(theta0).T
        
        # Step 2
        suma = 0
        for m in range(N):
            eps = n[m] / d0
            if np.abs(eps) <= k:
                ksi = eps
            else:
                ksi = k * np.sign(eps)
            suma += ksi ** 2
        d1 = (d0 ** 2 * suma / N / 0.7785) ** 0.5

        # Step 3
        for m in range(N):
            eps = n[m] / d1
            if eps > k:
                delta[m] = k * d1
            elif eps < -k:
                delta[m] = -k * d1
            else:
                delta[m] = n[m]
        
        # Step 4
        d_theta = np.linalg.inv(H.T.dot(H)).dot(H.T).dot(delta)

        # Step 5
        theta = theta0 + q * d_theta
        
        # Step 6
        continue_func = 0
        
        for m in range(p):
            if not((np.abs(theta[m] - theta0[m]) < np.abs(epsilon * theta0[m])).any() or
                   np.abs(d1 - d0) < epsilon * np.abs(d0).any()):
                continue_func = 1
        
        d0 = d1
        theta0 = theta
        
    d = d0
    return theta.reshape(8, ), d


def WLSQ_algorithm(samples: np.ndarray, p: int, theta0: np.ndarray, d0: int) -> np.ndarray:
    """
    WLSQ algorithm for calculating AR parameters
    :param samples: the samples of the signal
    :param p: AR model order
    :param theta0: starting theta parameters
    :param d0:
    :return: calculated theta parameters using WLSQ method
    """

    N = samples.shape[0] - p
    H = np.zeros((N, p))

    for m in range(0, N):
        if m == 0:
            H[m] = -samples[m + p - 1::-1, 0].T
        else:
            H[m] = -samples[m + p - 1:m - 1:-1, 0].T
    
    W = np.zeros((N, N))
    a4 = 0.5
    epsilon = 1e-6
    
    continue_func = 1
    num_iter = 0
    
    while continue_func and num_iter < 4:
        num_iter += 1
        suma = 0
        
        for m in range(N):
            eps = (samples[m + p] - H[m].dot(theta0).T) / d0
            
            if samples[m + p] != H[m].dot(theta0).T:
                if np.abs(eps) <= a4 * math.pi:
                    ksi = math.sin(eps / a4)
                else:
                    ksi = 0
                suma += d0 ** 2 * ksi ** 2
                W[m, m] = ksi / eps
            else:
                W[m, m] = 1
               
        d1 = (suma / N / 0.561722) ** 0.5
        theta = np.linalg.inv(H.T.dot(W).dot(H)).dot(H.T).dot(W).dot(samples[p:N+p, 0])

        # Step 6
        continue_func = 0
        
        for m in range(p):
            if not((np.abs(theta[m] - theta0[m]) < np.abs(epsilon * theta0[m])).any()):
                continue_func = 1
        
        d0 = d1
        theta0 = theta
        
    return theta.reshape(8, )


def calculate_coefficients(s: np.ndarray) -> tuple:
    """
    Code for calculating AR coefficients using LSQ, DUTTER and WLSQ algorithm
    :param s: input signal
    :return: a tuple of calculated AR parameters using LSQ, Dutter and WLSQ methods
    """
    length_signal = 256
    ar_parameters_LSQ = np.zeros((s.shape[0], 8))
    ar_parameters_DUTTER = np.zeros((s.shape[0], 8))
    ar_parameters_WLSQ = np.zeros((s.shape[0], 8))
    for m in tqdm(range(length_signal, s.shape[0])):
        ar_parameters_LSQ[m] = LSQ_algorithm(s[m - length_signal:m], 8)
        ar_parameters_DUTTER[m], d = DUTTER_algorithm(s[m - length_signal:m], 8, ar_parameters_LSQ[m].T)
        ar_parameters_WLSQ[m] = WLSQ_algorithm(s[m - length_signal:m], 8, ar_parameters_DUTTER[m].T, d)
        
    return ar_parameters_LSQ[length_signal:], ar_parameters_DUTTER[length_signal:], ar_parameters_WLSQ[length_signal:]


def synthesize_signal_strube(ar_parameters: list) -> tuple:
    """
    Code for synthesizing the signal when the excitation is Strube's glottal wave
    both case of noised and non-noised signal is considered
    :param ar_parameters: AR parameters
    :return: A tuple of synthesized strube's signals, first is without noise, while the second is with noise
    """
    # Parameters of the Strube wave excitation

    Ts = 3.2e-3
    Tn = 1.2e-3
    Tp_odb = 80  # Tp = 8e-3
    T_sample = 1e-4
    glotal_wave = np.zeros((Tp_odb, 1))
    for m in range(Tp_odb):
        if (m + 1) * T_sample < Ts:
            glotal_wave[m] = math.sin(math.pi * (m + 1) * T_sample / 2 / Ts) ** 2
        elif (m + 1) * T_sample < Ts + Tn:
            glotal_wave[m] = math.cos(math.pi * ((m + 1) * T_sample - Ts) / 2 / Tn)
        else:
            glotal_wave[m] = 0
    G_pom = glotal_wave
    for m in range(20):
        glotal_wave = np.append(glotal_wave, G_pom, axis=0)
    glotal_wave = glotal_wave[:1100]
    ug_prim = np.append(glotal_wave, np.array([[0]]), axis=0) - np.append(np.array([[0]]), glotal_wave, axis=0)
    ug_sek = np.append(ug_prim, np.array([[0]]), axis=0) - np.append(np.array([[0]]), ug_prim, axis=0)
    
    strube_wave = ug_sek[45:1045]
    strube_wave = strube_wave / np.max(strube_wave)
    signal_strube = signal.lfilter([1], ar_parameters, strube_wave, axis=0)
    signal_tmp = signal_strube
    
    signal_strube_without_noise = signal_strube.copy()  # Synthesized signal without noise
    eps = 0.05  # Probability of outlier occurrence
    sig1 = 2e-4
    sig2 = 100 * sig1
    pom1 = np.random.rand(1000)
    pom2 = np.random.randn(1000)
    for m in range(signal_strube.shape[0]):
        if pom1[m] < eps:
            temp = sig2 * pom2[m]
        else:
            temp = sig1 * pom2[m]
        signal_tmp[m] = signal_strube[m] + temp
    signal_strube_with_noise = signal_tmp  # Synthesized signal with noise
    plt.plot(signal_strube_without_noise)
    plt.xlabel('samples')
    plt.title('Synthesized signal when the excitation is Strube\'s glottal wave')
    plt.show()

    return signal_strube_without_noise, signal_strube_with_noise


def synthesize_signal_dirac(ar_parameters: list) -> tuple:
    """
    Code for synthesizing the signal when the excitation is a train of Dirac pulses both case of noised and
    non-noised signal is considered
    :param ar_parameters: AR parameters
    :return: A tuple of synthesized strube's signals, first is without noise, while the second is with noise
    """

    # Parameters of the Strube wave excitation
    Ts = 1e-4    
    Tp = 8e-3
    Tp_odb = Tp / Ts
    
    Dirac = np.zeros((1000, 1))
    height = 1
    Dirac[1] = height
    for m in range(Dirac.shape[0]):
        if (m + 1) % Tp_odb == 0:
            Dirac[m] = height
    
    signal_dirac = signal.lfilter([1], ar_parameters, Dirac, axis=0)
    
    signal_tmp = signal_dirac
    s_dirac_orig = signal_dirac.copy()
    
    eps = 0.05  # Probability of outlier occurrence
    sig1 = 2e-4
    sig2 = 100 * sig1
    pom1 = np.random.rand(1000)
    pom2 = np.random.randn(1000)
    for m in range(signal_dirac.shape[0]):
        if pom1[m] < eps:
            temp = sig2 * pom2[m]
        else:
            temp = sig1 * pom2[m]
        signal_tmp[m] = signal_dirac[m] + temp
        
    signal_dirac = signal_tmp  # Synthesized signal with noise
    plt.plot(s_dirac_orig)
    plt.xlabel('samples')
    plt.title('Synthesized signal when the excitation is train of Dirac pulses')
    plt.show()

    return s_dirac_orig, signal_dirac
