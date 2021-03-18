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
sij = np.array([random.randint(0, 1) for _ in range(11)])


tau = np.array([i for i in range(0, 141)])
ri = np.array([i for i in range(0, 141)])


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**(1/2))
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**(1/2))
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


P0 = probablity(1)

Ntx = (2*lambda_0*T*(10**(SNR/10)))/P0


def Cj(Ntx, j):
    return Ntx*probablity(j+1)


# thresold when P(ri|si=0) = P(ri|si=1)
# def get_tau(Ntx):

c0 = Cj(Ntx, 0)
theory_thresold = c0 / math.log(
    1+(c0/(sum([Cj(Ntx, i)/2 for i in range(1, L+1)]) + lambda_0*T)))


def p_BER(i, tau, Ntx):

    y = [Cj(Ntx, j)*sij[i-j] for j in range(1, L+1)]
    x = lambda_0*T + sum(y)
    x1 = x + Cj(Ntx, 0)

    ans = (1/2)*(sp.gammaincc(x, math.ceil(tau)) +
                 1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER(Ntx):

    min_ber = 10**5
    for t in tau:

        ans = (1/(2**L))*(sum([p_BER(i, t, Ntx)
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
    b2 = (1/(2**L))*(sum([p_BER(i, theory_thresold, Ntx)
                          for i in range(0, len(sij))]))
    b1, Optimal_tau = BER(Ntx)
    s0 = [p_ri_si(0, i) for i in ri]
    s1 = [p_ri_si(1, i) for i in ri]

    plt.plot(ri, s0, '.-')
    plt.plot(ri, s1, '.-', linewidth=0.4)
    plt.plot(Optimal_tau, b1, 'o')
    plt.plot(theory_thresold, b2, 'o')
    plt.show()


plot_graph()
