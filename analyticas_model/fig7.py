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

SNR = np.array([i for i in range(-5, 51, 2)])
tau = np.arange(0, 121, 5)
sij = np.random.randint(2, size=50)


def probablity(i):

    x1 = (d-r)/((4*D*i*T)**0.5)
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**0.5)
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i) for i in range(1, 7)])

Ntx = np.array([(2*lambda_0*T*(10**(i/10)))/pj[0] for i in SNR])


def Cj(Ntx, j):
    return Ntx*probablity(j+1)


cj = np.array(
    [[Cj(ntx, j) for ntx in Ntx] for j in range(0, 6)]
)

# def get_tau(Ntx):

#     c0 = Cj(Ntx, 0)
#     ans = c0 / \
#         math.log(
#             1+(c0/(sum([Cj(Ntx, i)/2 for i in range(1, L+1)]) + lambda_0*T)))
#     return ans


def p_BER(i, tau, cj):

    y = np.flipud(cj[1:])*sij[i-L-1:i-1]
    x = lambda_0*T + y.sum()
    x1 = x + cj[0]

    ans = 0.5*(sp.gammaincc(x, math.ceil(tau)) +
               1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER(cj):

    ans = [(1/(2**L))*(sum([p_BER(i, t, cj)
                            for i in range(6, len(sij)-5)])) for t in tau]

    return min(ans)


def plot_graph():

    plt.xlabel("SNR (db)")
    plt.ylabel("BER")
    plt.title("BER of actual vs theory reciver at T=30 \u0394")

    b1 = [BER(i) for i in np.transpose(cj)]
    plt.plot(SNR, b1, '.-', label="optimal thresold value", linewidth=0.5)
    plt.legend()
    plt.show()


plot_graph()
