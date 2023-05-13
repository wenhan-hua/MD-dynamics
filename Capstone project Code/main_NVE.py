import numpy as np
import matplotlib.pyplot as plt

N = 128   #number of simulation particles
dt = 1E-5   #simulation time period
total_steps = 500  #total steps in simulation

L = 20   #length of the box
vmax = 20   #maximum velocity of a particle in initialization
total_energy = 1E-9   #total energy for microcanonical ensemble
epsilon = -0.0077*(1.602E-19)  #ε in L-J potential(SI)
sigma = 4.5  #σ in L-J potential(SI)
m = 39.948E-27/6.02  #mass of the particle(SI)
kB = 1.38E-23  #Boltzmann constant(SI)

# initialize particle position as randomly distributed particles in a 3D box
def initializePosition(N):
    position = np.zeros((N, 3))
    for i in range(N):
        position[i] = np.random.uniform(0,L,3)

    return position

# initialize particle velocity as randomly distributed between (-vmax, vmax)
def initializeVelocity(N):
    velocity = np.zeros((N, 3))
    for i in range(N):
        velocity[i] = np.random.uniform(-vmax, vmax, 3)
    vc = velocity.sum(axis=0)/N  # central mass velocity
    for i in range(N):
        velocity[i] = velocity[i] - vc #calculate and later rescale velocity in the coordinate of central mass

    return velocity

#compute the accelarations of each particle at given position based on L-J potential
def computeAccelerations(position,N):
    a = np.zeros((N, 3))
    f = np.zeros((N,3))
    for i in range(N):
        for j in range(N):
            ri = []
            for k in range(3):
                x0 = position[i, k]
                x = position[j,k]
                ri.append((x-x0)**2)
            r = np.sqrt(sum(ri))
            for k in range(3):
                x0 = position[i, k]
                x = position[j,k]
                if i == j:
                    f[j, k] = 0
                else:
                    f[j,k] = (x0 - x) * ((sigma / r) ** 14 - 0.5 * (sigma / r) ** 8)
        for k in range(3):
            a[i] = 48 * (epsilon / sigma) ** 2 * f.sum(axis=0)
    return a

#Compute the position of particles using Verlet method
def computePosiiton(position, velocity, acceleration, dt, L):
    nextPosition = position + velocity*dt + 0.5*acceleration*(dt**2)
    position = nextPosition % L
    return position

#Compute the velocity of particles using Verlet method
def computeVelocity(position, velocity, acceleration, dt,N):
    r_nextdt = computePosiiton(position, velocity, acceleration, dt, L)
    a_nextdt = computeAccelerations(r_nextdt, N)
    return velocity + 0.5*(acceleration+a_nextdt)*dt

def instaneousTemperature(v,N):
    vsqr = []
    for i in range(N):
        for k in range(3):
            vsqr.append(v[i,k]*v[i,k])
    return sum(vsqr)*m/(3*(N-1)*kB)

def computePotentialEnergy(position,N):
    potential_energy = 0
    for i in range(N):
        for j in range(N):
            ri = []
            for k in range(3):
                x0 = position[i, k]
                x = position[j, k]
                ri.append((x - x0) ** 2)
            r = np.sqrt(sum(ri))
            if i == j :
                potential_energy = potential_energy
            else:
                potential_energy = potential_energy + 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    potential_energy = 0.5 * potential_energy
    return potential_energy

def computeKeneticEnergy(v,N):
    vSqdSum = 0
    for i in range(N):
        for k in range(3):
            vSqdSum = vSqdSum + v[i,k]*v[i,k]
    return 0.5*m*vSqdSum

#use a scaler to maintain a constant total energy in microcanonical ensemble
def rescaleVelocity(v,position,total_energy,N):
    potential_energy_total = computePotentialEnergy(position,N)
    vSqdSum = 0
    for i in range(N):
        for k in range(3):
            vSqdSum = vSqdSum + v[i,k]*v[i,k]
    lamda = np.sqrt(2*(total_energy - potential_energy_total)/ (vSqdSum * m))
    for i in range(N):
        for k in range(3):
            v[i,k] = lamda * v[i,k]
    return v


######-------------main loop starts here------------------######
### This loop will return two matrix position_output and velocity_output
### Each column records the nth step of the simulation result
### Each row records the position or velocity in three directions of the Nth particle

#0 Initialization
trajectory_r = np.zeros((total_steps+1,N,3))
trajectory_v = np.zeros((total_steps+1,N,3))
trajectory_a = np.zeros((total_steps+1,N,3))
for i in range(N):
    for k in range(3):
        trajectory_r[0][i][k] = initializePosition(N)[i][k]
        trajectory_v[0][i][k] = initializeVelocity(N)[i][k]
trajectory_v[0] = rescaleVelocity(trajectory_v[0],trajectory_r[0],total_energy,N)

print('PotentialEnergy',computePotentialEnergy(trajectory_r[0],N))
print('KeneticEnergy',computeKeneticEnergy(trajectory_v[0],N) )
#total_energy = computePotentialEnergy(trajectory_r[0],N) + computeKeneticEnergy(trajectory_v[0],N)
print('total_energy',total_energy)

#1 The simulation loop. Judge whether the system has reached the equilibrium state or not after visualization
for i in range(total_steps):
    trajectory_a[i] = computeAccelerations(trajectory_r[i], N)
    trajectory_r[i+1] = computePosiiton(trajectory_r[i], trajectory_v[i], trajectory_a[i], dt, L)
    trajectory_v[i+1] = computeVelocity(trajectory_r[i], trajectory_v[i], trajectory_a[i], dt,N)
    if i%10 == 0:     #rescale velocites every 50 steps
        trajectory_v[i+1] = rescaleVelocity(trajectory_v[i+1],trajectory_r[i+1],total_energy,N)

f = open('NVE_vx.txt','w')
v = trajectory_v[total_steps]
for i in range(N):
    f.write(f"{v[i,0]}\n")
f.close()

f = open('NVE_vy.txt','w')
v = trajectory_v[total_steps]
for i in range(N):
    f.write(f"{v[i,1]}\n")
f.close()

f = open('NVE_vz.txt','w')
v = trajectory_v[total_steps]
for i in range(N):
    f.write(f"{v[i,2]}\n")
f.close()

f = open('NVE_v.txt','w')
v = trajectory_v[total_steps]
for i in range(N):
    data = np.sqrt(v[i,0]**2 + v[i,1]**2 + v[i,2]**2)
    f.write(f"{data}\n")
f.close()
