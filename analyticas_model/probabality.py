import math
from typing import cast
import scipy.special as sp
import random
import matplotlib.pyplot as plt

lambda_0 = 100   #( lambda 0 )
r = 45 * 10**(-9) #Receiver radius r
d = 500 * 10**(-9)# Distance d
D = 4.256 * 10**(-10) # Diffusion Coefficient D 
delta_t = 9 *10**(-6) # Discrete Time Length
T = 30*delta_t # Slot Length 
L = 5 #Channel length
T2 = 50*delta_t

SNR = 30

tau =  [i for i in range(20,120)]

cj = [0]  #cj = Ntx*Pj
sij=  [random.randint(0,1) for i in range(50)]


def probablity(i,T):

    x1 = (d-r)/((4*D*i*T)**(1/2))
    try:
        x2 = (d-r)/((4*D*(i-1)*T)**(1/2))
        return (r/d)*(sp.erfc(x1) -  sp.erfc(x2))
    except:
        return  (r/d)*(sp.erfc(x1))

P01 = probablity(1,T) # q 
P02 = probablity(1,T2)

Ntx = (2*lambda_0*T*(10**(SNR/10)))/P01
Ntx2 = (2*lambda_0*T2*(10**(SNR/10)))/P02

def Cj(Ntx,j,T):
    return Ntx*probablity(j+1,T)


def p_BER(i,tau,Ntx,T):

    y  = [Cj(Ntx,j,T)*sij[i-j] for j in range(1,L)]
    x = lambda_0*T + sum(y)
    x1 = x + Cj(Ntx,0,T)

    ans = 1/2*(sp.gammaincc(x, math.ceil(tau)) + 1 - sp.gammaincc(x1,math.ceil(tau)))

    return ans

def BER(tau,Ntx,T):
    
    ans = (1/(2**L))*(sum([p_BER(i,tau,Ntx,T) for i in range(5, len(sij)-5)]))

    return ans



def plot_graph():
    ber = [BER(t,Ntx,T) for t in tau]
    ber2 = [ BER(t,Ntx2,T2) for t in tau ]
    plt.plot(tau, ber)
    plt.plot(tau,ber2)
    plt.show()


plot_graph()