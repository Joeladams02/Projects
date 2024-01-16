import numpy as np
import matplotlib.pyplot as plt
from math import exp

def solve_numerov(f,x,psi0,dpsi0):
    '''
    Author: Joel Adams, Date: 22.12.22
     Solving d^2(psi)/dx^2 = f(x) psi from x0 to x1 with boundary condition
     psi(x0) = psi0 and d(psi(x0))/dx = dpsi0

    Input:
    * f: The function to be called on the right hand side of the equation (see hint
    section), receives x and returns the value of f(x).
    * x: array of the integration, the first element is the lower bound, x0, and the last
    element is the upper bound, x1.
    * psi0: the value of \psi at x0, i.e. \psi(x0)
    
    * dpsi0: the value of \psi derivative at x0, i.e. d\psi(x0)/dx

    Output:
    * psi: array of values of \psi with each element in the array corresponding to the same
    element in x

    Example use:
    >> f = lambda x: x*x − E
    >> x = np.arange(0, x1, d)
    >> psi0 = 1
    >> dpsi0 = 0
    >> psi = solve_numerov(f, x, psi0, dpsi0)
    >> plt.plot(x, psi)
    '''
    d = x[1]-x[0] # Quick calc of the seperation
    psi = np.array([0.0]*len(x)) # Create an empty array to be populated with the solution
    psi[0] = psi0
    if psi0 == 0:
        psi[1] = d*dpsi0 + (dpsi0*f(0)*d**3)/6 #From taylor expansion for even BCs
    elif psi0 == 1:
        psi[1] = psi0 + (psi0*f(0)*d**2)/2 + (psi0*f(0)**2*d**4)/24 #From taylor expansion for odd BCs
        
    for j in range(1,len(x)-1): # Once the initial 2 values have been sorted, the Numerov program can then calculate the rest
        a = (2+5*d*d*f(x[j])/6)*psi[j]
        b = (1-(d*d/12)*f(x[j-1]))*psi[j-1]
        psi[j+1] = (a-b)/(1-d*d*f(x[j+1])/12) 
    return psi #Upon completion, return the solution


def find_oscillator_eigenvalue(E0, psi):
    '''
    Author: Joel Adams, Date: 22.12.22
    Solving for eigenvalue of d^2(psi)/dx^2 = f(x)

    Input:
    * E0: Initial guess of an eigenvalue. The program will output the nearest actual eigenvalue to this.
    * psi: The value of psi at x = x1 as a function of energy.

    Output:
    * E1: The eigenvalue closest to E0

    Example use:
    >> E0: 2.9
    >> psi = lambda E: solve_numerov(lambda x: x*x-E, x, 1, 0)[-1]
    >> 
    >> even = find_oscillator_eigenvalue(4.9,lambda E: solve_numerov(lambda x: x*x-E, x, 1, 0)[-1])
    >> odd = find_oscillator_eigenvalue(4.9,lambda E: solve_numerov(lambda x: x*x-E, x, 0, 1)[-1])
    >> print(min([odd,even]))
    '''
    for i in range(100):
        f = psi(E0)
        df = (psi(E0+0.001)-psi(E0-0.001))/0.002
        E1 = E0 - f/df
        if abs(E0 - E1) < 0.0000001:
            return E1
        else:
            E0 = E1


def hermite(arr,n):
    '''
    Author: Joel Adams, Date: 22.12.22
    Outputting the nth hermite polynomial multiplied by the exponential

    Input:
    * arr: the x values that are input to the fucntion
    * n: The nth hermite polynomial

    Output:
    * out: The array of y values
    
    Example use:
    >> arr = np.arange(0,x1,d)
    >> n = 3
    '''
    out = []
    for x in arr:
        H0 = 1
        H1 = 2*x
        for i in range(1, n):
            H2 = 2*x*H1 - 2*i*H0
            H0 = H1
            H1 = H2
        y = H1*np.exp(-x*x/2)
        out.append(y)
    if out[1] < 0:
        return [-1*y for y in out]
    else:   
        return out

d = float(input('d: ')) #Allow the user to choose the spacing between the points the function will be evaluated
x1 = float(input('x1: ')) #The range over which the function is evaluated [0,x1]5
E = float(input('E: ')) #The energy eigenvalue of the oscillator we want graphing
x = np.arange(0,x1,d) # produce the array of x values the fucntion will be evaluated on


'''
Use the specified energy eigenvalue given to calculate and plot the function for both odd and even BCs.
This is achieved by calling the numerov function i have defined above.
'''

psi_out_even = solve_numerov(lambda x: x*x-E, x, 1, 0) #even BCs - flat and non zero at x=0
psi_out_odd = solve_numerov(lambda x: x*x-E, x, 0, 1) # call the solver with odd BCs - 0 but with a gradient


'''
Given a non eigenvalue energy, find the nearest eigenvalue to it, regardless of even or odd BC.
This is achieved by implementing the Newton Raphson root finding algorithm to find the energy at which
the func is 0 at x=5. This is the defining feature of the eigenvalue as the func has to be normalisable.
'''
even = find_oscillator_eigenvalue(E,lambda E: solve_numerov(lambda x: x*x-E, x, 1, 0)[-1])
odd = find_oscillator_eigenvalue(E,lambda E: solve_numerov(lambda x: x*x-E, x, 0, 1)[-1])
print('Closest eigenvalue: %.3f' % min([odd,even])) # After calculating nearest root in both odd and even case, compare for nearest.


'''
Graph the solved function for both odd and even side by side. It is evident that only specified eigenvalues
for the correct parity give sensible results, and theyre of the form of Legendre polynomials.
'''

analytic = hermite(x,int(0.5*(E-1)))

plt.subplot(1,3,1)
plt.plot(x,psi_out_even)
plt.title('even')
plt.subplot(1,3,2)
plt.plot(x,psi_out_odd)
plt.title('odd')
plt.subplot(1,3,3)
plt.plot(x,analytic)
plt.title('analytic')
plt.show()
