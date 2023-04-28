import numpy as np
import matplotlib.pyplot as plt
import quantity as qt
import os

N= 20#Number of particles
dt= 0.001#time period
total_steps=1000 #total steps should be large enough to reach equilibrium

r=np.zeros((N, 3)) #create a n*3 matrix to store the positions
v=np.zeros((N, 3)) #create a n*3 matrix to store the velocities
a=np.zeros((N, 3)) #create a n*3 matrix to store the accelerations

L= 10 #length of the box. The periodical boundary condition ensures that particals are restraiend wihtin the box
T= 300 #Temperature for microcanonical ensemble
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

time=len(trajectory_r)
X,Y,Z = M_T(trajectory_r)

for t in range(time):
    fig=plt.figure(figsize=(16, 9))
    ax = plt.axes(projection='3d')
    # 3d contour plot
    ax.scatter3D(X[t], Y[t], Z[t])
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo['grid']['color'] = (1, 1, 1, 0)
    ax.yaxis._axinfo['grid']['color'] = (1, 1, 1, 0)
    ax.zaxis._axinfo['grid']['color'] = (1, 1, 1, 0)


    # save figure with different names depend on the view
    filename = '3d/3d_picture_' + str(t) + '.png'
    plt.savefig(filename, dpi=75)
    plt.close(fig)

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
    frame = new_frame.copy()
    frames.append(frame)
    new_frame.close()

for i in files:
    os.remove(i)

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



def p(trajectory_r,L):
    #compute_accelerations()
    V=L**3
    F=m*np.array(a)
    X=[]
    for t in range(len(trajectory_r)):
        map_r=trajectory_r[t]
        f=F[t]
        x=[]
        for l in range(len(f)):
            xx=np.dot(f[l],map_r[l])
            x.append(xx)
        X.append(sum(x))
    P=(kB*N*T+(1/3)*np.array(X))/V
    return P

def c(E,T):
    kB = 1.38E-23
    T = 300
    C=(1/(kB*T**2))*(np.var(E))**2
    return C

p=p(trajectory_r,L)

p4=plt.figure()
plt.plot(t,U)
plt.title('pressure')
plt.xlabel('time')
plt.ylabel('pressure')

plt.show()

C=c(U,T)
print(C)


