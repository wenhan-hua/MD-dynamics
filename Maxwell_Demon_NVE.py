import numpy as np

N = 20   #number of simulation particles
dt = 0.01   #simulation time period
total_steps=5 #total steps in simulation

L = 10   #length of the box
a = 1  #length of the Maxwell demon's door
r = 0.1 #the size of a particle
total_energy = 100   #total energy for microcanonical ensemble
vmax = 0.4   #maximum velocity of a particle in initialization
epsilon = -0.0077*(1.602E-19)  #ε in L-J potential(SI)
sigma = 4.5  #σ in L-J potential(SI)
m = 39.948E-23/6.02  #mass of the particle(SI)
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
    vc = velocity.sum(axis=0)  # central mass velocity
    for i in range(N):
        velocity[i] = velocity[i] - vc #calculate and later rescale velocity in the coordinate of central mass

    return velocity

def indoor(position):
    # position is the info of a single particle, size = [3]
    if (L - r)/2 < position[0] < (L + r)/2 and (L - a)/2 < position[1] < (L + a) / 2 and (L - a)/2 < position[2] < (L + a)/2:
        return True
    return False

# Judge wether particle meet the Maxwell Demon's door
def MaxwellDemon_judgeparticle(position, velocity, N):
    # The method that Maxwell Demon judge whether a particle can pass is by compare its velocity with average one
    # If the particle is in the box in the center, it will be judged by the Maxwell Demon, if the velocity is larger than average, it will be bounded
    velocity_abs = np.zeros(N)
    
    for i in range(N):
        velocity_abs[i] = np.sqrt(velocity[i, 0] ** 2 + velocity[i, 1] ** 2 + velocity[i, 2] ** 2)
    velocity_abs_average = np.average(velocity_abs)
    for i in range(N):
        if velocity_abs[i] < velocity_abs_average and velocity[i, 0] > 0 and indoor(position[i, :]):
            velocity[i, 0] = - velocity[i, 0]
        elif velocity_abs[i] > velocity_abs_average and velocity[i, 0] < 0 and indoor(position[i, :]):
            velocity[i, 0] = - velocity[i, 0]
    return velocity

def Boundary_bounce_particle(position, velocity, N):
    for i in range(N):
        if position[i, 0] < r or L - position[i, 0] < r:
            velocity[i, 0] = - velocity[i, 0]
        if position[i, 1] < r or L - position[i, 1] < r:
            velocity[i, 1] = - velocity[i, 1]
        if position[i, 2] < r or L - position[i, 2] < r:
            velocity[i, 2] = - velocity[i, 2]
    return velocity


#compute the accelarations of each particle at given position based on L-J potential
def computeAccelerations(position,N):
    acceleration = np.zeros((N, 3))
    f = np.zeros((N,3))
    ri = []
    for i in range(N):
        for j in range(N):
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
            acceleration[i] = 48 * (epsilon / sigma) ** 2 * f.sum(axis=0)
    return acceleration

#Compute the position of particles using Verlet method
def computePosition(position, velocity, acceleration, dt, L):
    nextPosition = position + velocity*dt + 0.5*acceleration*(dt**2)
 
    return nextPosition

#Compute the velocity of particles using Verlet method
def computeVelocity(position, velocity, acceleration, dt,N):
    r_nextdt = computePosition(position, velocity, acceleration, dt, L)
    a_nextdt = computeAccelerations(r_nextdt, N)
    return MaxwellDemon_judgeparticle(r_nextdt, Boundary_bounce_particle(r_nextdt, velocity + 0.5*(acceleration+a_nextdt)*dt, N), N)

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



#use a scaler to maintain a constant total energy in microcanonical ensemble
def rescaleVelocity(v,position,N):
    potential_energy_total = computePotentialEnergy(position,N)
    vSqdSum = 0
    T = instaneousTemperature(v,N)
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
position_output = [initializePosition(N)]
velocity_output = [initializeVelocity(N)]
acceleration = computeAccelerations(position_output[0],N)

#1 The simulation loop. Judge whether the system has reached the equilibrium state or not after visualization
for i in range(total_steps):
    position_original = computePosition(position_output[i], velocity_output[i], acceleration[i], dt, L)
    position_output.append(position_original)
    velocity_original = computeVelocity(position_output[i], velocity_output[i], acceleration[i], dt,N)
    if i%50==0:     #rescale velocites every 50 steps
        velocity_output.append(rescaleVelocity(velocity_original,position_original,N))
    else:
        velocity_output.append(velocity_original)



