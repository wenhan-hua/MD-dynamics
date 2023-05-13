import numpy as np
import matplotlib.pyplot as plt
import quantity as qt
import os


N= 300#Number of particles
dt= 1e-5 #time period
total_steps= 800 #total steps should be large enough to reach equilibrium
a_door = 10 #length of the Maxwell demon's door
radium = 0.1 #the size of a particle
ac = []

r=np.zeros((N, 3)) #create a n*3 matrix to store the positions
v=np.zeros((N, 3)) #create a n*3 matrix to store the velocities
a=np.zeros((N, 3)) #create a n*3 matrix to store the accelerations

L= 20 #length of the box. The periodical boundary condition ensures that particals are restraiend wihtin the box
T= 300 #Temperature for microcanonical ensemble
vmax=100 #restriction for initialzing velocities
epsilon = -0.0077*(1.602E-19) #ε in L-J potential(SI)
sigma = 4.5 #σ in L-J potential(SI)
m = 39.948E-27/6.02 #particle mass(SI)
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

def indoor(r):
    # position is the info of a single particle, size = [3]
    if (L - radium)/2 < r[0] < (L + radium)/2 and (L - a_door)/2 < r[1] < (L + a_door) / 2 and (L - a_door)/2 < r[2] < (L + a_door)/2:
        return True
    return False

def outdoor(r):
    # position is the info of a single particle, size = [3]
    if (L - radium)/2 < r[0] < (L + radium)/2 and (not indoor(r)):
        return True
    return False

# Judge wether particle meet the Maxwell Demon's door
def MaxwellDemon_judgeparticle(r, v, N):
    # The method that Maxwell Demon judge whether a particle can pass is by compare its velocity with average one
    # If the particle is in the box in the center, it will be judged by the Maxwell Demon, if the velocity is larger than average, it will be bounded
    v_abs = np.zeros(N)
    
    for i in range(N):
        v_abs[i] = np.sqrt(v[i, 0] ** 2 + v[i, 1] ** 2 + v[i, 2] ** 2)
    v_abs_average = np.average(v_abs)
    for i in range(N):
        if outdoor([r[i, 0], r[i, 1], r[i, 2]]):
            v[i, 0] = - v[i, 0]
        elif v_abs[i] < v_abs_average and v[i, 0] > 0 and indoor([r[i, 0], r[i, 1], r[i, 2]] ):
            v[i, 0] = - v[i, 0]
        elif v_abs[i] > v_abs_average and v[i, 0] < 0 and indoor([r[i, 0], r[i, 1], r[i, 2]]):
            v[i, 0] = - v[i, 0]

def Boundary_bounce_particle(r, v, N):
    for i in range(N):
        if r[i, 0] < radium or L - r[i, 0] < radium:
            v[i, 0] = - v[i, 0]
        if r[i, 1] < radium or L - r[i, 1] < radium:
            v[i, 1] = - v[i, 1]
        if r[i, 2] < radium or L - r[i, 2] < radium:
            v[i, 2] = - v[i, 2]


def compute_accelerations():
    a=np.zeros((N, 3))
    for i in range(N-1):
        for j in range(i+1, N, 1):#add upp all two-body interactions
            rij=[0,0,0]
            rSqd=0
            for k in range(3):
                rij[k]=r[i][k]-r[j][k]
                rSqd=rSqd+rij[k]**2
            f=24*((epsilon*sigma**-1)**2*(2*(sigma**-2*rSqd)**-7)-(sigma**-2*rSqd)**-4)
            #this formula is given by first order derivative of L-J potential
            for k in range(3):
                a[i][k]=a[i][k]+rij[k]*f/m
                a[j][k]=a[j][k]-rij[k]*f/m

def velocity_Verlet():
    compute_accelerations()
    for n in range(N):
        for i in range(3):
            r[n][i]=r[n][i]+v[n][i]*dt+0.5*a[n][i]*dt**2
            r[n][i]=r[n][i]
            v[n][i]=v[n][i]+0.5*a[n][i]*dt
    compute_accelerations()
    for n in range(N):
        for i in range(3):
            v[n][i]=v[n][i]+0.5*a[n][i]*dt
    MaxwellDemon_judgeparticle(r, v, N)
    Boundary_bounce_particle(r, v, N)


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
    c=a.tolist()
    ac.append(c)
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
    ax.plot3D([L/2, L/2, L/2, L/2, L/2],\
               [L/2 - a_door/2, L/2 + a_door/2, L/2 + a_door/2, L/2 - a_door/2, L/2 - a_door/2],\
                [L/2 + a_door/2, L/2 + a_door/2, L/2 - a_door/2, L/2 - a_door/2, L/2 + a_door/2], 'r')
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
tt=np.linspace(0,time,len(trajectory_r)-1)
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

'''print(ac)
print(len(ac))
print(len(trajectory_r),'!!!')'''

def p(trajectory_r,L):
    #compute_accelerations()
    V=L**3
    X = []
    for t in range(len(trajectory_r)-1):
        #print(t)
        b = np.array(ac[t])
        f = m * b
        map_r=trajectory_r[t]
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
plt.plot(tt,p)
plt.title('pressure')
plt.xlabel('time')
plt.ylabel('pressure')

v=trajectory_v[total_steps]
r=trajectory_r[total_steps]
v_amp = []
vx=[]
v_amp_left = []
v_left = []
v_right = []
v_amp_right = []
for t in range(len(v)):
    vx.append(v[t][0])
    v_amp.append(np.sqrt(v[t][0] ** 2 + v[t][1] ** 2 + v[t][2] ** 2))
    if r[t][0] < L / 2:
        v_amp_left.append(np.sqrt(v[t][0] ** 2 + v[t][1] ** 2 + v[t][2] ** 2))
        K_E_left = np.sum(np.array(v_amp_left) * np.array(v_amp_left) / 2 * m)
    else:
        v_amp_right.append(np.sqrt(v[t][0] ** 2 + v[t][1] ** 2 + v[t][2] ** 2))
        K_E_right = np.sum(np.array(v_amp_right) * np.array(v_amp_right) / 2 * m)

plt.figure()
plt.title('Distribution of the velocity')
plt.hist(vx)
plt.show()
plt.title('Distribution of the velocity amplitude')
plt.hist(v_amp)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(v_amp_left, bins=np.arange(0, 2300, 200), label='left box', alpha = 0.5, color= 'b')
ax.hist(v_amp_right, bins=np.arange(0, 2300, 200), label='right box', alpha = 0.5, color= 'r')
ax.set_xlim(0, 2500)
plt.legend()
plt.show()

print('The kinetic energy of left side is:{} Joule'.format(K_E_left))
print('The kinetic energy of right side is:{} Joule'.format(K_E_right))

C=c(U,T)
print(C)




