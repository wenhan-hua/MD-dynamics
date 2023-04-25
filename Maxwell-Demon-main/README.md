# Maxwell's Demon

This repository builds a visualization of Maxwell's Demon using Pygame.

## Description

This is a thought experiment designed by James Maxwell in the 19th century that challenges the laws of thermodynamics. Imagine there are two chambers with gas molecules of different energy. The chambers are connected through a door, controlled by a demon. The demon opens and shuts the door such that fast particles can pass to one side and slow molecules to the other side. This process causes the two chambers to become more organized (one becomes hot and the other cold), thus the entropy decreases. According to the 2nd law, entropy can only increase in a closed system, violating the law. If this was true, it would be fantastic! Assume the two chambers are divided perfectly into hot and cold particles: lots of energy would be available for useful work. Once the particles get mixed, there will be no way to go back to hot and cold, so there is no energy available. If Maxwell's demon existed, this system could go back without additional work, building some sort of perpetual motion machine.

![dem](https://user-images.githubusercontent.com/86791449/130356965-b5f2f563-64ec-430c-8226-16d270b7e86b.png)

This enigma stayed unsolved for centuries and found its answer only with the coming of new information technology. If you look at the model, you could ask where does the extra work come from. It could be the demon opening and shutting the door, but we can think of a demon randomly operating the door, exerting the same amount of work. Then, the work must be done in the demon’s brain - specifically, in the information used to decide on the action of the door. The demon needs to recognize the behavior of single particles, to decide which can or cannot cross the door. From this information, the demon can create order from chaos. It follows that information is entropy, regardless of which system keeps the information - brain, computer, or particles.

This Python program produces a visualization of the two chambers randomly filled with gas and runs until the particles are divided between hot and cold ones.

![demon 1](https://user-images.githubusercontent.com/86791449/130357757-81cfce27-486f-40a6-86a0-2dee82028d50.png)


## Getting Started

### Dependencies

We will need:
* random
* math
* pygame
* matplotlib
* drawnow

You will need to pip install the last three. Also, the main file inherits from the Blue and Red classes.
Make sure you have the background PNG in the smae folder to display the two chambers.


### Executing program

Run the main script to produce the visualization. You can tweak the size, number, and speed of particles.

```
particlesN = 50     
size = 12
particleb.speed = 0.7    
particler.speed = 1                              
```

A possible improvement could be a 3-D rendering.

## Help

Make sure you have installed all dependencies and that all four files are in the same folder

## License

This project is licensed under The Unlicense - see the LICENSE.md file for details

## Version History

* 0.1
    * Initial Release
