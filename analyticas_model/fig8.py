import math
from typing import cast
import scipy.special as sp
from scipy.stats import poisson
import random
import numpy as np
import matplotlib.pyplot as plt

lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 5  # Channel length

ri = np.arange(0, 140, 1)
# print(ri)
SNR = 25
sij = np.array([random.randint(0, 1) for _ in range(140)])


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**(1/2))
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**(1/2))
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i) for i in range(1, 8)])

c0 = 2*lambda_0*T*(10**(SNR/10))

Ntx = c0/pj[0]


def Cj(Ntx, j):
    return Ntx*probablity(j+1)


cj = np.array(
    [[Cj(Ntx, j) for j in range(1, 6)] for _ in range(140)]
)
# thresold when P(ri|si=0) = P(ri|si=1)
# def get_tau(Ntx):

# cj = np.array([Cj(Ntx,j) for j in range()])


def p_BER(i, si):

    lam = (c0*si) + (cj[i].sum()/2) + (lambda_0*T)

    prob = (math.exp(-lam)*(lam**i))/math.factorial(i)

    return prob


sub_tau = c0 / \
    math.log(
        1+(c0/((cj[1:5].sum()/2) + lambda_0*T)))


def plot_graph():
    s0 = [p_BER(i, 0) for i in ri]
    s1 = [p_BER(i, 1) for i in ri]

    # plt.yscale('log')

    plt.xlabel('Recivied Particle')
    plt.ylabel("Probability")

    plt.plot(ri, s0, '.-', linewidth=0.4, label="si=0")
    plt.plot(ri, s1, '.-', linewidth=0.4, label="si=1")
    plt.plot(sub_tau, 0.9, 'o-', linewidth=0.4, label="theory thresold")
    # plt.plot(theory_thresold, b2, 'o-',
    #          linewidth=0.4, label="optimal thresold")
    # plt.legend(loc=1, fontsize='x-small')

    plt.show()


plot_graph()
