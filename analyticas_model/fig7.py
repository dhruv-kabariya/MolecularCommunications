import math
from typing import cast
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')
plt.grid(b=True, linewidth=0.5)


lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 5  # Channel length

SNR = np.arange(-5, 56, 1)
tau = np.arange(0, 121, 0.5)
sij = np.random.randint(2, size=50)


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**0.5)
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**0.5)
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i) for i in range(1, 8)])

c0 = 2*lambda_0*T*(10**(SNR/10))

Ntx = c0/pj[0]


def Cj(Ntx, j):
    return Ntx*pj[j+1]


cj = np.array(
    [[Cj(ntx, j) for ntx in Ntx] for j in range(1, 6)]
)


def get_tau(co, cj):

    ans = co / \
        math.log(
            1+(co/((cj.sum()/2) + lambda_0*T)))
    return ans


def p_BER(i, tau, co, cj):

    y = np.flipud(cj)*sij[i-L-1:i-1]
    x = lambda_0*T + y.sum()
    x1 = x + co

    ans = 0.5*(sp.gammaincc(x, math.ceil(tau)) +
               1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER(cj, co):

    ans = [(1/(2**L))*(sum([p_BER(i, t, co, cj)
                            for i in range(6, len(sij)-5)])) for t in tau]

    return min(ans)


def BER2(cj, co):
    ans = (1/(2**L))*(sum([p_BER(i,  get_tau(co, cj), co, cj)
                           for i in range(6, len(sij)-6)]))

    return ans


def plot_graph(cj):

    plt.xlabel("SNR (db)")
    plt.ylabel("BER")
    plt.title("BER of actual vs theory reciver at T=30 \u0394")

    cj = np.transpose(cj)
    b1 = [BER(cj[i], c0[i]) for i in range(len(c0))]
    b2 = [BER2(cj[i], c0[i]) for i in range(len(c0))]
    # print(b1)
    plt.plot(SNR, b1, '.-', label="optimal thresold value", linewidth=0.5)
    plt.plot(SNR, b2, '.-', label="sub optimal thresold value", linewidth=0.5)

    plt.legend(loc=3, fontsize='x-small')
    plt.show()


plot_graph(cj)
