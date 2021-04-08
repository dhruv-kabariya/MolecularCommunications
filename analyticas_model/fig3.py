import math
import numpy as np
from typing import cast
import scipy.special as sp
import random
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')
plt.grid(b=True)

lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 5  # Channel length
T2 = 50*delta_t

SNR = 30

tau = np.arange(20, 121, 5)


def probablity(i, t):

    x1 = (d-r)/((4*D*i*t)**(1/2))
    try:
        x2 = (d-r)/((4*D*(i-1)*t)**(1/2))
        return (r/d)*(sp.erfc(x1) - sp.erfc(x2))
    except:
        return (r/d)*(sp.erfc(x1))


pj = np.array([probablity(i, T) for i in range(1, 7)])
pj2 = np.array([probablity(i, T2) for i in range(1, 7)])

Ntx = (2*lambda_0*T*(10**(SNR/10)))/pj[0]
Ntx2 = (2*lambda_0*T2*(10**(SNR/10)))/pj2[0]


def Cj(Ntx, j, t):
    return Ntx*probablity(j+1, t)


cj = np.array([Cj(Ntx, j, T) for j in range(0, 6)])  # cj = Ntx*Pj
cj1 = np.array([Cj(Ntx2, j, T2) for j in range(0, 6)])

sij = np.array([random.randint(0, 1) for _ in range(50)])


def p_BER(i, tau, cj, t):

    y = np.flipud(cj[1:])*sij[i-L-1:i-1]
    x = lambda_0*t + y.sum()
    x1 = x + cj[0]

    ans = 1/2*(sp.gammaincc(x, math.ceil(tau)) +
               1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def BER(tau, cj, t):

    ans = (1/(2**L))*(sum([p_BER(i, tau, cj, t)
                           for i in range(6, len(sij)-5)]))

    return ans


def plot_graph():
    ber = [BER(t, cj, T) for t in tau]
    ber2 = [BER(t, cj1, T2) for t in tau]
    # fig, ax = plt.subplots()
    plt.yticks(ber2)

    print(ber2)
    plt.plot(tau, ber, '.-', label='Solt lenght 30', linewidth=0.5)
    plt.plot(tau, ber2, '.-', label='Solt lenght 50', linewidth=0.5)
    plt.xlabel('Thresold Value')
    plt.ylabel('BER')
    plt.title("FIG-3")
    plt.legend()
    plt.show()


plot_graph()
