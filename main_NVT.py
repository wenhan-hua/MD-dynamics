import numpy as np

N=10#Number of particles
dt=0.001#time period
total_steps=5 #total steps should be large enough to reach equilibrium

r=np.zeros((N, 3)) #create a n*3 matrix to store the positions
v=np.zeros((N, 3)) #create a n*3 matrix to store the velocities
a=np.zeros((N, 3)) #create a n*3 matrix to store the accelerations

L=20 #length of the box. The periodical boundary condition ensures that particals are restraiend wihtin the box
T=300 #Temperature for microcanonical ensemble
vmax=100 #restriction for initialzing velocities
epsilon = -0.0077*(1.602E-19) #ε in L-J potential(SI)
sigma = 4.5 #σ in L-J potential(SI)
m = 39.948E-23/6.02 #particle mass(SI)
kB = 1.38E-23 #Boltzmann constant(SI)

def initialize(): #initialize positions--uniformly place particles in a 3D box
    for n in range(N):
        for i in range(3):
            r[n][i]=np.random.uniform(0,L)

def rescale_velocities(): #rescale the velocities to maitain a constant temperature
     vSqdSum=0
     for n in range(N):
          for i in range(3):
               vSqdSum=vSqdSum+v[n][i]**2
     lamda=(3*(N-1)*T*kB/(vSqdSum*m))**0.5
     for n in range(N):
          for i in range(3):
               v[n][i]=v[n][i]*lamda

def initialize_velocities(): #randomly initialize velocities
    for n in range(N):
        for i in range(3):
            v[n][i]=np.random.uniform(-vmax,vmax)
    vCM=[0,0,0] #vCM is the velocity of the center of mass
    for n in range(N):
        for i in range(3):
            vCM[i]=vCM[i]+v[n][i]
    for i in range(3):
        vCM[i]=vCM[i]/N
    for n in range(N):
        for i in range(3):
            v[n][i]=v[n][i]-vCM[i] #velocity in the center of mass coordinate
    rescale_velocities()


def compute_accelerations():
    a=np.zeros((N, 3))
    for i in range(N-1):
        for j in range(i+1, N, 1):#add upp all two-body interactions
            rij=[0,0,0]
            rSqd=0
            for k in range(3):
                rij[k]=r[i][k]-r[j][k]
                rSqd=rSqd+rij[k]**2
            f=24*(epsilon*sigma**-1)**2*(2*(sigma**-2*rSqd)**-7)-(sigma**-2*rSqd)**-4
            #this formula is given by first order derivative of L-J potential
            for k in range(3):
                a[i][k]=a[i][k]+rij[k]*f/m
                a[j][k]=a[j][k]-rij[k]*f/m

def velocity_Verlet():
    compute_accelerations()
    for n in range(N):
        for i in range(3):
            r[n][i]=r[n][i]+v[n][i]*dt+0.5*a[n][i]*dt**2
            r[n][i]=r[n][i]%L
            v[n][i]=v[n][i]+0.5*a[n][i]*dt
    compute_accelerations()
    for n in range(N):
        for i in range(3):
            v[n][i]=v[n][i]+0.5*a[n][i]*dt

def instant_temperature(): #compute instant temperature
    #we can print the result of this function in every loop to ensure that temperature is fluctating within a small range
    sum=0
    for n in range(N):
        for i in range(3):
            sum=sum+v[n][i]**2
    return sum*m/(3*(N-1)*kB)



######-------------main loop starts here------------------######
trajectory_r=np.zeros((total_steps+1, N, 3 ))
#this matrix will recored the initial positions as well as the positions after each steps
trajectory_v=np.zeros((total_steps+1, N, 3 ))
#this matrix will recored the initial velocities as well as the velocities after each steps
initialize()
initialize_velocities()
for n in range(N):
    for k in range(3):
        trajectory_r[0][n][k] = r[n][k]
        trajectory_v[0][n][k] = v[n][k]
for i in range(total_steps):
    velocity_Verlet()
    if i%50==0:#rescale velocites every 50 steps
        rescale_velocities()
    for n in range(N):
        for k in range(3):
            trajectory_r[i+1][n][k]=r[n][k]
            trajectory_v[i+1][n][k]=v[n][k]







