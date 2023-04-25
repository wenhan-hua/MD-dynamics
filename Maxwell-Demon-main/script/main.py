'''

DARIO FUMAROLA - PHYS265 FINAL PROJECT
SIMULATION OF THE MAXWELL'S DEMON
File: main.py

This program produces a visual simulation of the Maxwell's demon.
Imagine you have two chambers with gas molecules of different energy. The chambers are connected by a door, controlled by a demon. 
The demon opens and shuts the door such that fast particles can pass to one side, and slow molecules to the other side. 
This causes the two chambers to become more organized (one becomes hot and the other cold), thus the entropy decreases.
The program keeps track of temperature and kinetic energy in each of the chamber.


WORKS CITED:

Vector Addition, Collision, and Pygame generic functions adapted from Beginning Game Development with Python and Pygame (Will McGugan)

'''

    # Initial Parameters
# Import libraries necessary for the simulation
import random
import pygame
import matplotlib.pyplot as plt
import math
from drawnow import *
from blueParticle import bluePar
from redParticle import redPar


# Initialize the game
pygame.init()


# Colors
l_red = (250, 140, 140)
l_blue = (65,165,255)
l_yellow = (255,255,150)
blue = (0,0,255)
red = (255,0,0)
yellow = (255,255,0)


# Set variable for scenario
font = pygame.font.match_font('calibri')
rightEnergy= []
plt.ion()
font = pygame.font.Font(pygame.font.get_default_font(), 30)
colorBackground = (1,1,1)
(w, h) = (1000, 600)       # Width and height of window
elastic = 1                # Affect way particles rebound against each other


    # Main functions
def makeText(surf, text, x, y, color):
    """Takes positional inputs to place texts on screen."""
    surfText = font.render(str(text),True, color)
    rectText = surfText.get_rect()
    rectText.topleft = (x,y)
    surf.blit(surfText, rectText)


def measureTempR(particles):
    """Keeps track of temperature in right chamber, by counting particles in it."""
    tempR = 0                                      # Start temperature counter
    particlesR = []                                # List of particles in right chamber
    for i, particle in enumerate(particlesLyst):   # Loop through all particles
        if particle.x < 500:                       # Coordinates to be in right chamber
            tempR += particle.temp                 # Add the particle intrinsic temperature to the count
            particlesR.append(1)                   # Add particle to the list

    return round(tempR / (len(particlesR)),2)      # Return rounded result


def measureTempL(particles):
    """Keeps track of temperature in left chamber, by counting particles in it."""
    tempL = 0 
    particlesL = []                             
    for i, particle in enumerate(particlesLyst):
        if particle.x > 500:
            tempL += particle.temp
            particlesL.append(1)

    return round(tempL / (len(particlesL)),2)


def keLeft(particles):
    """Keeps track of kinetic energy in left chamber, by adding energy of particles in it."""
    keL = 0                                         # Start KE counter

    for i, particle in enumerate(particlesLyst):    # Loop through all particles      
        if particle.x < 500:                        # Coordinates to be in left chamber
            keL += particle.speed**2                # Add square of speed to kinetic energy

    return round(keL ,2)                            # Return rounded result


def keRight(particles):
    """Keeps track of kinetic energy in right chamber, by adding energy of particles in it."""
    keR = 0

    for i, particle in enumerate(particlesLyst):
        if particle.x > 500:
            keR += particle.speed**2

    return round(keR ,2)


def vectorAddition(angle1, length1, angle2, length2):
    """Standard addition of two vectors - used to determine collision between particles."""
    x = math.sin(angle1) * length1 + math.sin(angle2) * length2
    y = math.cos(angle1) * length1 + math.cos(angle2) * length2

    length = math.hypot(x, y)
    angle = 0.5 * math.pi - math.atan2(y, x)

    return (angle, length)


def collision(p1, p2):
    """Determine collision bounds between two moving particles and prevents overlap."""
    # Delta movements for inbound collisions
    dx = p1.x - p2.x
    dy = p1.y - p2.y

    dist = math.hypot(dx, dy)
    if dist < p1.size + p2.size:                  # Collision conditions
        angle = math.atan2(dy, dx) + 0.5 * math.pi
        total_mass = p1.mass + p2.mass            # Mass of the collision

        # Recalculate new angles and speeds
        (p1.angle, p1.speed) = vectorAddition(p1.angle, p1.speed * (p1.mass - p2.mass) / total_mass,
                                          angle, 2 * p2.speed * p2.mass / total_mass)
        (p2.angle, p2.speed) = vectorAddition(p2.angle, p2.speed * (p2.mass - p1.mass) / total_mass,       
                                          angle + math.pi, 2 * p1.speed * p1.mass / total_mass)
        p1.speed *= elastic                       # rebound back by elasticity
        p2.speed *= elastic

        overlap = 0.5 * (p1.size + p2.size - dist + 1)
        p1.x += math.sin(angle) * overlap
        p1.y -= math.cos(angle) * overlap
        p2.x -= math.sin(angle) * overlap
        p2.y += math.cos(angle) * overlap


# Main loop conditions
pygame.display.set_caption('Maxwell Demon Simulation')
background = pygame.image.load('background.png')
screen = pygame.display.set_mode((w, h+200))
particlesN = 50     # Number of particles
particlesLyst = [] 

# Place particles at random coordinates
for n in range(particlesN):
    size = 12
    x1 = random.randint(size, w - size)           # Red x
    y2 = random.randint(size, h - size)           # Red y
    x = random.randint(size, w-size)              # Blue x
    y = random.randint(size, h-size)              # Blue x

    # Initialize particles red
    particler = redPar(x1,y2,size)
    particler.speed = 1                              # Higher speed for cold particles
    particler.angle = random.uniform(0, math.pi*2)   # Random starting angle         
    particlesLyst.append(particler)

    # Initialize particles blue
    particleb = bluePar(x, y, size)
    particleb.speed = 0.7                            # Slower speed for cold particles
    particleb.angle = random.uniform(0, math.pi*2)   # Random starting angle         
    particlesLyst.append(particleb)


# Main loop for the simulation
playing = True
while playing:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            playing = False

    screen.fill(colorBackground)
    screen.blit(background,(0,0))


    for i, particle in enumerate(particlesLyst):
        particle.move()         # Initialize movements
        particle.rebound()       # Initialize rebounds

        for particle2 in particlesLyst[i+1:]:
            collision(particle, particle2)       # Initialize collisions
        particle.display()                       # Display particles

        # Display texts
        makeText(screen, measureTempR(particle), 300, 700, l_blue)
        makeText(screen, 'T = ', 230, 700, blue)
        makeText(screen, measureTempL(particle), 760, 700, l_red)
        makeText(screen, 'T = ', 690, 700, red)

        makeText(screen, keRight(particle), 760, 740, l_yellow)
        makeText(screen, 'KE  = ', 660, 740, yellow)

        makeText(screen, keLeft(particle), 300, 740, l_yellow)
        makeText(screen, 'KE  = ', 210, 740, yellow)
    pygame.display.flip()
# End of the script