from math import sin, sinh,exp,cos,sin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def laplace_solver(psi, a, N):
    '''
    Author: Joel Adams, Date: 22.12.22
    Solving Laplace equation using the over-relaxation method

    Input:
    * psi: 2D matrix containing the initial \psi, including boundaries.
    * a: The coefficient of over-relaxation
    * N: Maximum number of iterations performed
    
    Output:
    * psi: 2D matrix of the value of \psi after (up to) N_iter iterations.
    * hist_values: (N_iter x 3) matrix that contains historical values of 3 points during % the iteration (1 in the upper half, 1 in the middle, and 1 in the lower half).

    Example use:
    >> init_psi = np.array([[0.0]*7]*7) 
    >> init_psi[6,:] = [sin(x/6)*sinh(1) for x in range(7)]
    >> init_psi[:,6] = [sin(1)*sinh(y/6) for y in range(7)] 
    >> alpha = 1.35 
    >> N = 30
    >> psi, hist = laplace_solver(init_psi, alpha, N)
    '''
    hist = [[],[],[]] # ready to be filled with the history values
    for n in range(N): # repeat N times 
        for x in range(1,6):
            for y in range(1,6): # for each coordinate
                x,y = 6-x,6-y #pan from top right corner (next to B.C. downwards and leftwards into 0's
                R = psi[x,y+1] + psi[x,y-1] + psi[x+1,y] + psi[x-1,y] - 4*psi[x,y] # perform the numerical technique
                psi[x,y] = psi[x,y] + a*R/4 # redefine that coordinates value with the new value
        hist[0].append(psi[2,1]) 
        hist[1].append(psi[3,3])
        hist[2].append(psi[4,5])# picked 3 coordinate samples to record the historical evolution
        
        if n > 2: # This is a break out test so doesnt go all the way to 30 if no need
            if abs(hist[0][-2] - hist[0][-1]) < 0.000000000001 and abs(hist[1][-2] - hist[1][-1]) < 0.000000000001 and abs(hist[2][-2] - hist[2][-1]) < 0.000000000001:
                break
    return psi, hist # after entire grid combed over 30 times, return the 'solution' and the evolution values. 


def DDD_plot():
    '''
    Create a 3D plot of psi
    
    Example use:
    >> plot = input('1: 3d \n2: contour \n3: coord evolution \n4: alpha choice\n')
    >> [DDD_plot,contour_plot, evolution_plot,alpha_plot][int(plot)-1](psi)
    '''
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(range(7), range(7))
    ax.plot_surface(X, Y, psi) # cheeky 3d plot



def contour_plot():
    '''
    Create a contour plot of psi
    
    Example use:
    >> plot = input('1: 3d \n2: contour \n3: coord evolution \n4: alpha choice\n')
    >> [DDD_plot,contour_plot, evolution_plot,alpha_plot][int(plot)-1](psi)
    '''
    plt.contourf(psi, levels = 100) 
    plt.colorbar()

def evolution_plot():
    '''
    Create a line plot of the 3  historical values

    Example use:
    >> plot = input('1: 3d \n2: contour \n3: coord evolution \n4: alpha choice\n')
    >> [DDD_plot,contour_plot, evolution_plot,alpha_plot][int(plot)-1](psi)
    '''
    plt.plot(range(len(hist[0])),hist[0]) 
    plt.plot(range(len(hist[1])),hist[1])
    plt.plot(range(len(hist[2])),hist[2])


def alpha_plot():
    '''
     Create a line plot for 4 values of alpha

    Example use:
    >> plot = input('1: 3d \n2: contour \n3: coord evolution \n4: alpha choice\n')
    >> [DDD_plot,contour_plot, evolution_plot,alpha_plot][int(plot)-1](psi)
    '''
    for items in [0.7,1.35,1.39, 1.7, 1.9]: # these are our sample alphas
        init_psi = np.array([[0.0]*7]*7)
        init_psi[:,6] = [sin(x/6)*sinh(1) for x in range(7)]
        init_psi[6,:] = [sin(1)*sinh(y/6) for y in range(7)]
        psi, hist = laplace_solver(init_psi, items, N) # repeated as above
        plt.plot(hist[2],label = str(items)) # one sample coord for each alpha
    plt.legend()  



init_psi = np.array([[0.0]*7]*7) #create a 7x7 2d array filled with floated 0's, this will be populated with our calculated laplacian
init_psi[6,:] = [sin(x/6)*sinh(1) for x in range(7)] # define the initial conditions, this is the upper edge of the square
init_psi[:,6] = [sin(1)*sinh(y/6) for y in range(7)] # this is the right edge of the square
alpha = 1.35 
N = 30
psi, hist = laplace_solver(init_psi, alpha, N) # plug these values into the laplace solver
for row in np.flip(psi):
    a = ''
    for x in np.flip(row):
        a = a+ ('%.6f ' %x)
    print(a+'\n') # this is a simple snippet that makes the solution print in an easy to read format; a la matlab


plot = input('1: 3d \n2: contour \n3: coord evolution \n4: alpha choice\n')
[DDD_plot,contour_plot, evolution_plot,alpha_plot][int(plot)-1]()
plt.show()
