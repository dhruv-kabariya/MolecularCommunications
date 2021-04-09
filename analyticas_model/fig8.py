import math
from typing import cast
import scipy.special as sp
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

Optimal_tau = 10**5

SNR = 25
sij = np.array([random.randint(0, 1) for _ in range(31)])


tau = np.arange(20, 141)
ri = np.arange(0, 141)


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


# thresold when P(ri|si=0) = P(ri|si=1)
# def get_tau(Ntx):

cj = np.array([Cj(Ntx, j) for j in range(1, 6)])

theory_thresold = c0 / math.log(
    1+(c0/((cj.sum()/2) + lambda_0*T)))


def p_BER(i, tau):

    y = np.flipud(cj)*sij[i-L-1:i-1]
    x = lambda_0*T + y.sum()
    x1 = x + c0

    ans = (0.5)*(sp.gammaincc(x, math.ceil(tau)) +
                 1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER():

    min_ber = 10**5
    for t in tau:

        ans = (1/(2**L))*(sum([p_BER(i, t)
                               for i in range(6, len(sij))]))
        if(min_ber >= ans):
            min_ber = ans
            Optimal_tau = t
    # print(min_ber)
    # print(Optimal_tau)
    return min_ber, Optimal_tau


# min of BER form all tau

# distribution

def p_ri_si(si, ri):

    lambda_p = lambda_0*T + c0*si + \
        (sum([Cj(Ntx, j) for j in range(1, L+1)])/2)
    return (((math.e)**(-lambda_p)) * (lambda_p**ri))/math.factorial(ri)


def plot_graph():
    b2 = (1/(2**L))*(sum([p_BER(i, theory_thresold)
                          for i in range(6, len(sij))]))
    b1, Optimal_tau = BER()
    s0 = [p_ri_si(0, i) for i in ri]
    s1 = [p_ri_si(1, i) for i in ri]
    # plt.yscale('log')
    plt.xlabel('Recivied Particle')
    plt.ylabel("BER")

    plt.plot(ri, s0, '.-', linewidth=0.4, label="si=0")
    plt.plot(ri, s1, '.-', linewidth=0.4, label="si=1")
    plt.plot(Optimal_tau, b1, 'o-', linewidth=0.4, label="theory thresold")
    plt.plot(theory_thresold, b2, 'o-',
             linewidth=0.4, label="optimal thresold")
    plt.legend(loc=1, fontsize='x-small')

    plt.show()


plot_graph()
