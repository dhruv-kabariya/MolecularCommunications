import numpy as np
import scipy.special as sp
# from scipy.stats import poisson
from numpy.random import poisson as ran

lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 100  # Channel length

SNR = np.arange(-5, 50, 1)
sj = np.random.randint(2, size=L)


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**0.5)
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**0.5)
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i) for i in range(1, len(L)+1)])

c0 = 2*lambda_0*T*(10**(SNR/10))

Ntx = c0/pj[0]


def Cj(ntx):
    return ntx*pj


cj = np.array(
    [[Cj(ntx, j) for ntx in Ntx] for j in range(1, len(sj)+1)]
)

cj = np.transpose(cj)


def writeData():

    for i in len(SNR):

        t = sj*Cj(Ntx[i])

        summation1 = t.sum()

        Co(ii) = snr(ii)*2*lambda*T
        average1 = (lambda*T + summation1)
        average2 = average1 + data1.*Co(ii)
