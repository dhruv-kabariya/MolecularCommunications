import math
from typing import cast
import scipy.special as sp
import random
import matplotlib.pyplot as plt

lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 5  # Channel length

SNR = [i for i in range(-5, 51)]
tau = [i for i in range(61)]
sij = [random.randint(0, 1) for i in range(50)]


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**(1/2))
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**(1/2))
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


P01 = probablity(1)  # q

Ntx = [(2*lambda_0*T*(10**(i/10)))/P01 for i in SNR]


def Cj(Ntx, j):
    return Ntx*probablity(j+1)


# def get_tau(Ntx):

#     c0 = Cj(Ntx, 0)
#     ans = c0 / \
#         math.log(
#             1+(c0/(sum([Cj(Ntx, i)/2 for i in range(1, L+1)]) + lambda_0*T)))
#     return ans


def p_BER(i, tau, Ntx):

    y = [Cj(Ntx, j)*sij[i-j] for j in range(1, L+1)]
    x = lambda_0*T + sum(y)
    x1 = x + Cj(Ntx, 0)

    ans = 1/2*(sp.gammaincc(x, math.ceil(tau)) +
               1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER(Ntx):

    ans = [(1/(2**L))*(sum([p_BER(i, t, Ntx)
                            for i in range(5, len(sij)-5)])) for t in tau]

    return min(ans)


def plot_graph():

    b1 = [BER(i) for i in Ntx]
    plt.plot(SNR, b1)
    plt.show()


plot_graph()
