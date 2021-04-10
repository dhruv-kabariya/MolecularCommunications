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
tau = np.arange(40, 60, 0.5)
sij = np.random.randint(2, size=6)
