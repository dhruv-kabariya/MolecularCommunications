import math
import numpy as np


lambda_0 = 100   #( lambda 0 )
r = 45 * 10**(-9) #Receiver radius r
d = 500 * 10**(-9)# Distance d
D = 4.256 * 10**(-10) # Diffusion Coefficient D 
delta_t = 9 *10**(-6) # Discrete Time Length
T = 30*delta_t # Slot Length 
L = 5 #Channel length


C0 = 1 #undefine for now as i don't get it's value



SNR = 10*math.log10(C0/(2*lambda_0*T))


#NTX is total particle
NTX = 2*lambda_0*T*10**(SNR/10)
