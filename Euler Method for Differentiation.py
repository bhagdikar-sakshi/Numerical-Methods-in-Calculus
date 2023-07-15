import numpy as np
import math
import matplotlib.pyplot as plt

def f(t): 
    '''This function returns the RHS of the given equation.'''
    return 0

def coefficients():
    '''This function allows to input the values of the coefficients k,c and m.'''
    m=1
    c=1
    k=1
    return k,m,c

def g(t):
    '''This function is the modification of f(t) in order to substitute for further solution.'''
    m=1
    return f(t)/m

def parameters():
    '''This function allows to set the values of the initial and final time, and time interval. '''
    delta_t=0.1
    t_0=0
    t_final=10
    return delta_t,t_0,t_final

def alpha():
    '''This function allows to set the value of alpha in order to substitute in the generalised alpha method.'''
    a=0.5
    return a

def initial_conditions():
    '''This function allows to set the initial conditions for the second order ODE.'''
    x=2
    x_dot=-1
    return x,x_dot

def eqn_matrices(a,h):
    '''This function enables the creation of the matrices required to be created in order to find the solution.'''
    k,m,c=coefficients()
    A=np.array([[0,1],[-k/m,-c/m]]) # coefficient matrix
    I=np.identity(2) # 2x2 identity matrix
    P=I+(np.multiply(((1-a)*h),A)) # matrix formed by multiplying scalars corresponding to alpha and delta_t with 
                                   # the matrix A
    Q=I-(np.multiply((a*h),A)) # matrix formed by multiplying scalars corresponding to alpha and delta_t with 
                               # the matrix A
    return A,P,Q

def euler_method():
    '''This function is the actual implementation of the genralised alpha Euler method to find the solution 
    of the given ODE.'''
    
    h,t_i,t_f=parameters()
    a=alpha()
    x_0,v_0=initial_conditions()
    
    n=math.floor((t_f-t_i)/h)
    t=np.linspace(t_i,t_f,n)

    A,P,Q=eqn_matrices(a,h)

    R=np.linalg.inv(Q)

    X=[] # array to store the values of x
    for k in range(n):
        X.append(0)

    V=[] # array to store the values of v
    for k in range(n):
        V.append(0)

    X[0]=x_0
    V[0]=v_0

    for i in range(n-1): # actual matrix manipulation to find successive values of x and v
        C=(np.multiply(a*h,np.array([[0],[g(t[i+1])]])))+(np.multiply((1-a)*h,np.array([[0],[g(t[i])]]))) 
        # Here C is the matrix of the non-coupled terms
        [[X[i+1]],[V[i+1]]]=np.matmul(R,(np.matmul(P,np.array([[X[i]],[V[i]]]))+C))
        
    return X,V,t

def x(t): 
    '''This function is the analytical solution obtained for the given differential equation.'''
    return np.exp(-0.5*t)*(((3)**(0.5))*np.cos((((3)**(0.5))/2)*t)+np.sin((((3)**(0.5))/2)*t))

def graphs(t,X):
    '''This function plots the graphs for the calculations performed and the analytical solution.'''
    p=x(t)
    plt.plot(t,X)
    plt.plot(t,p)
    plt.legend(["Obtained Solution","Analytical Solution"])
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.title("Graph of x vs t for a=0.5, h=0.1")
    plt.show()   
    
    return 0

X_soln,V_soln,t_span=euler_method()
graphs(t_span,X_soln)