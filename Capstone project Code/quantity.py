import numpy as np

m = 39.948E-23 / 6.02  # particle mass(SI)
kB = 1.38E-23  # Boltzmann constant(SI)

def potential(r1,r2):
    epsilon = -0.0077 * (1.602E-19)  # ε in L-J potential(SI)
    sigma = 4.5  # σ in L-J potential(SI)
    r12 = np.sqrt(np.dot(r2 - r1,r2-r1))
    return 4*epsilon*((sigma/r12)**12-(sigma/r12)**6)

# The Kinetic energy
def K_E(map,m):
    time = len(map)
    K=[]
    for t in range(time):
        k_e=[]
        a = map[t]
        num = len(a)
        for n in range(num):
            P = a[n]
            PP = np.array(P)
            k_e.append((m*np.dot(PP,PP))/(2))

        K.append(sum(k_e))

    return K

# The potential energy
def pot(map):
    time = len(map)
    P=[]
    for t in range(time):
        p_e=[]
        a = map[t]
        num = len(a)
        for i in range(num):
            r1=a[i]
            for j in range(num):
                if j == i:
                    continue
                else:
                    r2 = a[j]
                    p_e.append(potential(r1,r2))

        P.append(sum(p_e))

    return P





