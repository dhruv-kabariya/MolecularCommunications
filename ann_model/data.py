import numpy as np
import scipy.special as sp
# from scipy.stats import poisson
from numpy.random import poisson
import pandas as pd

lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 100  # Channel length

SNR = np.arange(0, 55, 1)
sj = np.random.randint(2, size=L)


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**0.5)
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**0.5)
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i) for i in range(1, L+1)])

c0 = 2*lambda_0*T*(10**(SNR/10))

Ntx = c0/pj[0]


def Cj(ntx):
    return ntx*pj


ri = np.zeros((55, 100))


def writeData():

    for i in range(len(SNR)):

        t = sj*Cj(Ntx[i])

        summation1 = t.sum()

        co = SNR[i]*2*lambda_0*T
        average1 = (lambda_0*T + summation1)
        average2 = average1 + sj*co
        for j in range(len(average2)):

            global ri
            ri[i][j] = poisson(average2[j])
        # ri[i] = poisson()

    # writing

    ri = np.reshape(ri, -1)
    y = sj*np.ones((55, 100))
    y = np.reshape(y, -1)
    A = np.stack((ri, y), axis=1)

    data = pd.DataFrame(
        {

            "X": ri,
            "Y": y
        }
    )
    data.to_csv("x_y.csv")


writeData()
