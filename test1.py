import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as stat

# stats.norm.cdf(1.65, loc = 0, scale = 1) loc=mean and scale  =  sigma



lambda_0 = 100   #( lambda 0 )
r = 45 #Receiver radius r
d = 500 # Distance d
D = 4.256 * 10**(-10) # Diffusion Coefficient D 
delta_t = 9 # Discrete Time Length
T = 30*delta_t # Slot Length 
L = 5 #Channel length

C0 = 1 #undefine for now as i don't get it's value
#NTX is total particle


SNR = 10*math.log10(C0/(2*lambda_0*T))
NTX = 2*lambda_0*T*10
