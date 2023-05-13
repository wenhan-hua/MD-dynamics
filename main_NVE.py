import numpy as np
import matplotlib.pyplot as plt
import quantity as qt

N = 50  #number of simulation particles
dt = 1e-10   #simulation time period
total_steps=1000 #total steps in simulation

L = 20   #length of the box
total_energy = 100   #total energy for microcanonical ensemble
vmax = 0.4   #maximum velocity of a particle in initialization
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
    vc = velocity.sum(axis=0)  # central mass velocity
    for i in range(N):
        velocity[i] = velocity[i] - vc #calculate and later rescale velocity in the coordinate of central mass

    return velocity

#compute the accelarations of each particle at given position based on L-J potential
def computeAccelerations(position,N):
    a = np.zeros((N, 3))
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
    position_original = computePosiiton(position_output[i], velocity_output[i], acceleration[i], dt, L)
    position_output.append(position_original)
    velocity_original = computeVelocity(position_output[i], velocity_output[i], acceleration[i], dt,N)
    if i%50==0:     #rescale velocites every 50 steps
        velocity_output.append(rescaleVelocity(velocity_original,position_original,N))
    else:
        velocity_output.append(velocity_original)

trajectory_r = position_output
trajectory_v = velocity_output

#visulization
def M_T(map):
    time = len(map)
    num = len(map[0])
    name_png=[]
    frames = []
    X=[]
    Y=[]
    Z=[]
    for i in range(time):
        x=[]
        y=[]
        z=[]
        for j in range(num):
            x.append(map[i][j][0])
            y.append(map[i][j][1])
            z.append(map[i][j][2])
        X.append(x)
        Y.append(y)
        Z.append(z)

    return X,Y,Z

time=len(position_output)
X,Y,Z = M_T(position_output)

for t in range(time):
    plt.figure(figsize=(16, 9))
    ax = plt.axes(projection='3d')
    # 3d contour plot
    ax.scatter3D(X[t], Y[t], Z[t])
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo['grid']['color'] = (1, 1, 1, 0)
    ax.yaxis._axinfo['grid']['color'] = (1, 1, 1, 0)
    ax.zaxis._axinfo['grid']['color'] = (1, 1, 1, 0)
    ax.set_zlim(-5, 25)
    ax.set_ylim(-5, 25)
    ax.set_xlim(-5, 25)

    # save figure with different names depend on the view
    filename = '3d/3d_picture_' + str(t) + '.png'
    plt.savefig(filename, dpi=75)

from PIL import Image
png_count = time
files = []
for t in range(time):
    seq = str(t)
    file_names = '3d/3d_picture_' + seq + '.png'
    files.append(file_names)

print(files)

# Create the frames
frames = []
for i in files:
    new_frame = Image.open(i)
    frames.append(new_frame)

    # Save into a GIF file that loops forever
frames[0].save('3d/3d_vis.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=20, loop=0)


K=qt.K_E(trajectory_v,m)
P=qt.pot(trajectory_r)
t=np.linspace(0,time,len(trajectory_r))

#print(K,'!!!!',P,'!!!!',K)

plt.close('all')

p1=plt.figure()
plt.plot(t,K)
plt.title('Kinetic energy')
plt.xlabel('time')
plt.ylabel('Kinetic energy')

p2=plt.figure()
plt.plot(t,P)
plt.title('Potential energy')
plt.xlabel('time')
plt.ylabel('Potential energy')

U=np.array(P)+np.array(K)

p3=plt.figure()
plt.plot(t,U)
plt.title('Total energy')
plt.xlabel('time')
plt.ylabel('Total energy')

plt.show()

