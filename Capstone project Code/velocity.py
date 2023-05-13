from math import cos,sin,pi
import matplotlib.pyplot as plt
import pymatgen as py
import linecache
from pymatgen.core import Structure
import numpy as np
import scipy.optimize as fitter

T=600
m = 39.948*10**-23/6.02 #particle mass(SI)
kB = 1.38*10**-23
N=200

f = open('/Users/sheena/PycharmProjects/NVT/v_absolute-200-5000.txt')
lines = f.readlines()
v = []
f.close()
for line in lines:
    v.append(float(line.split()[0]))

sample=v

n=20#number of bins


y, bin_edges = np.histogram(sample, range=(0,40),bins=n)

# bin_edges is an array of left/right edges of all bins
x = 0.5*(bin_edges[1:] + bin_edges[:-1])
# x is now an array of the bin centers

# Poisson errors are just the square root of the number of counts in each bin
error_y = np.sqrt(y)
error_y = [max(error,1) for error in error_y]


plt.errorbar(x,y,error_y, fmt='o')
#plt.show()

def model(x, A, sigma):
    '''Define your model function according to
       the Gaussian function given above.'''
    return A*(x**2)*np.exp(-(x**2)/sigma)

sig = error_y
par0 = np.array([0.002,120])
par, cov = fitter.curve_fit(model, x, y, par0, sig, absolute_sigma=True)

A=par[0]
sigma=par[1]
print(A)
print(sigma)


xfit = np.linspace(0,40,100)
y=A*(xfit**2)*np.exp(-(xfit**2)/sigma)
plt.hist(v,bin_edges)
v_t=np.linspace(0,40,100)
y_t=N*((2/pi)*(m/kB/T)**3)**0.5*v_t**2*np.exp(-m*v_t**2/2/kB/T)
plt.plot(v_t,y_t)
plt.plot(xfit,y,'r-')
plt.title('distribution of absolute velocity(T=600K)')
plt.xlabel('velocity(m/s)')
plt.ylabel('number')
plt.legend(['theory', 'simulation'], loc='upper right')


plt.show()



