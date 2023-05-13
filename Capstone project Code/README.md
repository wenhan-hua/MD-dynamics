In the code for our project, 
we mainly used two python libraries: numpy and matplotlib.pyplot. scipy.optimize is used for curve fitting and os is also used for file writing and deleting.

 There are four files in the code. 
 1. quantity.py is used to calculate some basic physical statistical quantities for further visualization. It is called by main_NVE.py, main_NVT.py and Maxwell_Demon_NVT.py. 
 2. The first file main_NVE.py, main_NVT.py is used for running NVE(NVT) system and corresponding visualization. 
 3. An additional file velocity.py is used for velocity distribution visualization and curve fitting. The file Maxwell_Demon_NVT.py can run the case with a Maxwell's Demon and visualize the result.