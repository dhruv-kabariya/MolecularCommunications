import math
import re
from typing import cast
import scipy.special as sp

import numpy as np
import matplotlib.pyplot as plt
from scipy.special.orthogonal import legendre

plt.style.use('fivethirtyeight')
plt.grid(b=True, linewidth=0.5)


lambda_0 = 100  # ( lambda 0 )
r = 45 * 10**(-9)  # Receiver radius r
d = 500 * 10**(-9)  # Distance d
D = 4.256 * 10**(-10)  # Diffusion Coefficient D
delta_t = 9 * 10**(-6)  # Discrete Time Length
T = 30*delta_t  # Slot Length
L = 5  # Channel length

SNR = np.arange(0, 61, 1)
tau = np.arange(50, 1000, 0.5)
sij = np.random.randint(2, size=60)


sum0 = 0
sum1 = 0

global sum_cj
sum_cj = 0


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


def new_cj(Ntx, j):

    sum_cj = sum_cj + Ntx*pj[j]
    return sum_cj


cj = np.array(
    [[Cj(ntx, j) for ntx in Ntx] for j in range(1, 6)]
)


def get_tau(co, cj):

    ans = co / \
        math.log(
            1+(co/((cj.sum()/2) + lambda_0*T)))
    return ans


def new_tau(c0_n, avg1):

    t3 = c0_n/avg1
    t2 = math.log(1+t3)
    t1 = c0_n/t2

    ans1 = 0
    lg0 = avg1
    lg1 = avg1+c0_n

    for i in range(0, t1, 1):
        d1 = lg0**(i) * (math.e ** (-lg0))
        d2 = lg1**(i) * (math.e ** (-lg1))
        ans1 = ans1 + (d1 - d2)/math.factorial(i)

    return ans1


def new_tau2(c0_n, avg1, final_t):

    ans2 = 0
    lg0 = avg1
    lg1 = avg1+c0_n

    for i in range(0, final_t, 1):
        d1 = lg0**(i) * (math.e ** (-lg0))
        d2 = lg1**(i) * (math.e ** (-lg1))
        ans1 = ans1 + (d1 - d2)/math.factorial(i)

    return ans1


def p_BER(i, tau, co, cj):

    y = np.flipud(cj)*sij[i-L-1:i-1]
    x = lambda_0*T + y.sum()
    x1 = x + co

    ans = 0.5*(sp.gammaincc(x, math.ceil(tau)) +
               1 - sp.gammaincc(x1, math.ceil(tau)))

    return ans


def Q_func(lam, n):

    u1 = math.factorial(n)
    u3 = lam*n
    u2 = math.e**(-lam)*u3
    ans = u2/u1

    return ans


def BER(cj, co):

    ans = [(1/(2**L))*(sum([p_BER(i, t, co, cj)
                            for i in range(6, len(sij)-5)])) for t in tau]

    return min(ans)


def BER2(cj, co):
    ans = (1/(2**L))*(sum([p_BER(i,  get_tau(co, cj), co, cj)
                           for i in range(6, len(sij)-6)]))

    return ans


# for i in range(1, 5+1):
#     sum0 = sum0 + new_cj(Ntx[i], i)
#     sum1 = sum1 + sij[i]*new_cj(Ntx[i], i)

avg = sum0 + lambda_0*T
avg1 = avg+20
final_term = 100000

pe0 = []
pe1 = []


def new_BER():

    for i in range(1, len(SNR)):
        t = new_tau(c0[i], avg)
        pe0.append(0.5 - 0.5*(t))

        for j in range(1, len(tau)):

            lam1 = lambda_0*T + sum1
            lam2 = lam1 + c0[i]
            first = Q_func(lam1, tau[j])
            last = Q_func(lam2, tau[j])

            final = 0.5*(first+last)

            if(final < final_term):
                final_term = final
                final_tao = sij[j]

        b = new_tau2(c0[i], avg1, final_tao)
        pe1.append(0.5 - 0.5*(b))


def plot_graph(cj):

    # new_BER()
    # plt.yscale('log')
    # plt.yticks(np.arange(1, 0, 0.01))
    plt.xlabel("SNR (db)")
    plt.ylabel("BER")
    plt.title("BER of actual vs theory reciver at T=30 \u0394")

    cj = np.transpose(cj)
    b1 = [BER(cj[i], c0[i]) for i in range(len(SNR))]
    b2 = [BER2(cj[i], c0[i]) for i in range(len(SNR))]
    # print(b1)
    plt.plot(SNR, b1, '.-',
             label="optimal thresold value", linewidth=0.5)
    plt.plot(SNR, b2, '.-',
             label="sub optimal thresold value", linewidth=0.5)

    plt.legend(loc=1, fontsize='x-small')
    plt.show()


plot_graph(cj)
