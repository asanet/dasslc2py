#########################################################################
#                                                                       #
#   This is a template for the dasslc module in Python.                 #
#   Use it and modify as you wish.                                      #
#                                                                       #
#   Author: Ataide Neto                                                 #
#   email:ataide@peq.coppe.ufrj.br                                      #
#   Universidade Federal do Rio de Janeiro                              #
#   Version: 0.1-6                                                      #
#                                                                       #
#########################################################################

#t,y,yp = dasslc.solve(resfun,tspan,y0,yp0,rpar,rtol,atol,idx,ipfile,jac,display)

## Import the modules
import dasslc, time
import numpy as np
import matplotlib.pyplot as plt

SPARSE = False # Change it to True if sparse was compiled

### Defining the residual functions 

def model0(t,y,yp): #--------- Minimum of 3 input arguments
    res = np.empty(1) #------- Always allocate res as a numpy array, 
#                              even if it has len = 1.

    res[0] = yp[0] + 2*y[0] #- Declare the residual
    ires = 0 #---------------- Set ires = 0 if everything is ok

    return res, ires #-------- Beware: ires must always be returned
#                              as second output.


def model1(t,y,yp): #--------------- Just another example
    res = np.empty(2)
    res[0] = yp[0]-20*np.cos(20*t)
    res[1] = yp[1]+20*np.sin(20*t)
    return res, 0 #----------------- ires can be a literal


def model2(t,y,yp,par): #------------- Maximum of 4 input arguments
    res = np.empty(3)
         
    k1 = par[0]                     #|||||||||||||||||||||||||||||
    k2 = par[1]                     #|                           |
    CT0 = par[2]                    #| aliasing is optional,     |
    Ca = y[0]; dCa = yp[0] #--------#| but always encoraged      |
    Cb = y[1]; dCb = yp[1]          #|                           |
    Cc = y[2]; dCc = yp[2]          #|||||||||||||||||||||||||||||

    res[0] = -k1*Ca - dCa  
    res[1] = k1*Ca - k2*Cb - dCb
    res[2] = Ca + Cb + Cc - CT0
    ires = 0
    return res, ires


def model3(t,y,yp,par): #------- The parameter may be a whole class
    res = np.empty(5)
    res[0] = yp[0] - y[2]
    res[1] = yp[1] - y[3]
    res[2] = yp[2] + y[4]*y[0]
    res[3] = yp[3] + y[4]*y[1] + par.g
    if (par.dae == 3):
        res[4] = y[0]*y[0] + y[1]*y[1] - par.L*par.L
        ires = 0
    elif (par.dae == 2):
        res[4] = y[0]*y[2] + y[1]*y[3]
        ires = 0
    elif (par.dae == 1):
        res[4] = y[2]**2 + y[3]**2 - par.g*y[1] - par.L**2*y[4]
        ires = 0
    elif (par.dae == 0):
        res[4] = yp[4] + 3*y[3]*par.g/par.L**2
        ires = 0
    else:
        print("Invalid index.")
        ires = -1

    return res, ires

### Solve model0 

t0 = np.array([1, 2]) #-- Integration interval with
#                         initial and final time

y0 = np.array([1]) #----- Initial condition

# The simplest call to dasslc, with all the mandatory inputs
# and outputs. y and yp are equally spaced in all time span
t, y, yp = dasslc.solve(model0,t0,y0)

# Plot results
plt.figure(1)
plt.subplot(211)
plt.plot(t,y)
plt.ylabel('y')
plt.title('Model0 Solution')
plt.subplot(212)
plt.plot(t,yp)
plt.xlabel('time')
plt.ylabel('yp')


### Solve model1 

# The time span can also be a vector.
# In this case y and yp are returned at all values of t
t0 = np.linspace(0,1,100)

y0 = np.array([0,1]) #--- Initial condition
yp0 = np.array([1,0]) #-- Derivatives at initial condition (optional)

#-- Call with the optional yp0
t, y, yp = dasslc.solve(model1,t0,y0,yp0) 


# Plot results
plt.figure(2)
plt.plot(t,y)
plt.ylabel('y')
plt.xlabel('time')
plt.title('Model1 Solution')
plt.legend(["y1","y2"])

### Solve model2

# You can also specify only the final time. 
# In this case y and yp are equally spaced in [0 t0]
t0 = np.array([500])

y0 = np.array([1,0,0])

# If you are not passing an optional input,
# but is passing the next one, define it as None.
# Only positional arguments are supported.
yp0 = None

par = np.array([0.01,0.02,1]) #-- The optional parameter vector
atol = 1e-8 #-------------------- The absolute tolerance
rtol = 1e-6 #-------------------- The relative tolerance

# Call with optional arguments (yp0 = None)
t, y, yp = dasslc.solve(model2,t0,y0,yp0,par,rtol,atol) 

# Plot results
plt.figure(3)
plt.plot(t,y)
plt.ylabel('y')
plt.xlabel('time')
plt.title('Model2 Solution')
plt.legend(["Ca","Cb","Cc"])


### Solve model2 (steady state)

# If you specify the tspan as None,
t0 = None

# Here, y0 is the initial guess for the steady state
y0 = np.array([0,0,0])
yp0 = np.zeros(3) #--- At steady state, yp = 0

# Call with t0 = None for the steady state
t, y, yp = dasslc.solve(model2,t0,y0,yp0,par,rtol,atol) 

# Plot results
plt.figure(4)
plt.plot(0,y[0],'o')
plt.plot(1,y[1],'o')
plt.plot(2,y[2],'o')
plt.ylabel('Concentration')
plt.xlabel('Component')
plt.title('Model2 Solution (steady state)')
plt.legend(["Ca","Cb","Cc"])

#### Solve model3

#Defining the parameter class for 
class pend_par: 
    g = 9.81                    
    L = 1.0                      
    dae = 3  

t0 = np.linspace(0,50,10000)
y0 = np.array([1,0,0,0,0])
yp0 = None

# The optional parameter class initialization
par = pend_par() 
atol = 1e-10
rtol = 1e-8

# The dependent variable index vector (needed for high index DAE)
index = np.array([1,1,2,2,3]) 

t, y, yp = dasslc.solve(model3,t0,y0,yp0,par,rtol,atol,index)

# Plot results
plt.figure(5)
plt.plot(t,y)
plt.ylabel('y')
plt.xlabel('time')
plt.title('Model3 Solution')
plt.legend(["x","y","vx","vy","mu"])


### Solve model3 (with jacobian)

# The jacobian definition. See dasslc manual for more details
def jac_pend(t,y,yp,cj,par): 
    PD = np.zeros((5,5))
    PD[0][0] = cj
    PD[0][2] = -1
    PD[1][1] = cj
    PD[1][3] = -1
    PD[2][0] = y[4]
    PD[2][2] = cj
    PD[2][4] = y[0]
    PD[3][1] = y[4]
    PD[3][3] = cj
    PD[3][4] = y[1]

    if (par.dae == 3):
        PD[4][0] = 2*y[0]
        PD[4][1] = 2*y[1]
        ires = 0
    elif (par.dae == 2):
        PD[4][0] = y[2]
        PD[4][1] = y[3]
        PD[4][2] = y[0]
        PD[4][3] = y[1]
        ires = 0
    elif (par.dae == 1):
        PD[4][1] = -par.g
        PD[4][2] = 2*y[2]
        PD[4][3] = 2*y[3]
        PD[4][4] = -par.L**2
        ires = 0
    elif (par.dae == 0):
        PD[4][3] = 3*par.g/par.L**2
        PD[4][4] = cj
        ires = 0
    else:
        print("Invalid index.")
        ires = -1

    return PD, ires

# if passing the jacobian, the input file
# must be properly configured. See model3.dat file
display = False
t,y,yp = dasslc.solve(model3,t0,y0,yp0,par,rtol,atol,index,"model3.dat",jac_pend,display)

# Plot results
plt.figure(6)
plt.plot(t,y)
plt.ylabel('y')
plt.xlabel('time')
plt.title('Model3 Solution (with jacobian)')
plt.legend(["x","y","vx","vy","mu"])


if SPARSE:
    
    def model4(t,y,yp): #------------- A huge sparse system
        res = np.empty(Ns)
        #for i in range(0,Ns): #------ Avoid using for-loops 
        #   res[i] = yp[i] + y[i]      in python at all cost
        res = yp + y
        return res, 0

    Ns = 5000
    t0 = np.linspace(0,1,15)
    y0 = np.ones(Ns)
    yp0 = -y0

    #with dense algebra
    tic = time.time()
    t, y, _ = dasslc.solve(model4,t0,y0,None,None,1e-6,1e-8)
    toc = time.time() - tic

    #with sparse algebra

    tic = time.time()
    t1, y1, _ = dasslc.solve(model4,t0,y0,None,None,1e-6,1e-8,None,"model4.dat")
    toc1 = time.time() - tic

    #with sparse algebra and jacobian supplied
    def jacSparse(t,y,yp,cj): #---------------- The jacobian function
        J = np.eye(Ns,Ns) #-------------------- Just an eye matrix
        J = J + cj*J #------------------------- The true jacobian (iteration matrix)
        i = range(0,Ns) #---------------------- The rows list
        j = range(0,Ns) #---------------------- The columns list
        ires = 0
        return J, ires, i, j

    tic = time.time()
    t2, y2, _ = dasslc.solve(model4,t0,y0,None,None,1e-6,1e-8,None,"model4b.dat",jacSparse)
    toc2 = time.time() - tic

    plt.figure(7)
    plt.plot(0,toc,'o')
    plt.plot(1,toc1,'o')
    plt.plot(2,toc2,'o')
    plt.ylabel('Time (s)')
    plt.title('Model4 Dense vs. Sparse performance comparison')
    plt.legend(["Dense = %0.2f s" % toc,"Sparse = %0.2f s" % toc1,"Jac Sparse = %0.2f s" % toc2])

## Show all figures
plt.show()
