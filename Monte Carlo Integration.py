import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, dblquad,tplquad
import time


def f1(x):
    return np.exp((-1)*np.sin(x))

def f2(x,y):
    return (2*(x**2))+(y**2)

def f3(x,y,z):
    return (y*np.sin(x))+(z*np.cos(x))

def monte_carlo_integration_1D(a,b,n):
    x=np.zeros(n)
    for i in range(len(x)):
        x[i]=random.uniform(a,b)
 
    sum=0.0
    for i in range(n):
        sum=sum+f1(x[i])

    I=((b-a)/float(n))*sum
        
    return I

def monte_carlo_integration_2D(ax,bx,ay,by,n):
    x=np.zeros(n)
    for i in range(len(x)):
        x[i]=random.uniform(ax,bx)
    y=np.zeros(n)
    for i in range(len(y)):
        y[i]=random.uniform(ay,by)
    
    sum=0.0
    for i in range(n):
        sum=sum+f2(x[i],y[i])

    I=(((bx-ax)*(by-ay))/float(n))*sum
        
    return I

def monte_carlo_integration_3D(ax,bx,ay,by,az,bz,n):
    x=np.zeros(n)
    for i in range(len(x)):
        x[i]=random.uniform(ax,bx)
    y=np.zeros(n)
    for i in range(len(y)):
        y[i]=random.uniform(ay,by)
    z=np.zeros(n)
    for i in range(len(z)):
        z[i]=random.uniform(az,bz)
        
    sum=0.0
    for i in range(n):
        sum=sum+f3(x[i],y[i],z[i])
        
    I=(((bx-ax)*(by-ay)*(bz-az))/float(n))*sum
    
    return I

def area(a,b,n):
    for i in range(n):
        x=np.zeros(n)
        for i in range(len(x)):
            x[i]=random.uniform(a,b)
    
            sum=0.0
        for j in range(n):
            sum=sum+f1(x[j])

        I=((b-a)/float(n))*sum
        
        s=[]
        for k in range(n):
            s.append(I)
        
    return s

x=np.arange(-100,100,0.01)
y=f1(x)
plt.title("Graph of the function f(x)")
plt.plot(x,y)
plt.show()

a=area(1,4,1000)
plt.hist(a,density=True,bins=10)
plt.title("Areas obtained on integration")
plt.show()

    


a=time.time() 
I1=monte_carlo_integration_1D(a=1,b=4,n=1000)
b=time.time()
print("One- dimensional Monte Carlo Integration= ",I1)
print("Time= ",b-a)

s=time.time()
I_1=quad(f1,1,4) [0]
t=time.time()
print("Single Integration= ",I_1)
print("Time= ",t-s)


c=time.time()
I2=monte_carlo_integration_2D(ax=0,bx=2,ay=0,by=1,n=1000)
d=time.time()
print("Two-dimensional Monte Carlo Integration= ",I2)
print("Time= ",d-c)

p=time.time()
I_2=dblquad(f2,0,2,0,1)[0]
q=time.time()
print("Double Integration= ",I_2)
print("Time= ",q-p)

e=time.time()
I3=monte_carlo_integration_3D(ax=0,bx=np.pi,ay=0,by=1,az=-1,bz=1,n=1000)
f=time.time()
print("Three-dimensional Monte Carlo Integration= ",I3)
print("Time= ",f-e)
        
u=time.time()
I_3=tplquad(f3,0,np.pi,0,1,-1,1)[0]
v=time.time()
print("Triple Integration= ",I_3)
print("Time= ",v-u)


fig=plt.figure()
ax=plt.axes(projection='3d')
x=np.arange(-1*np.pi,np.pi,0.1)
y=np.arange(-1*np.pi,np.pi,0.1)
X,Y=np.meshgrid(x,y)
ax.plot_surface(X,Y,f2(X,Y))
plt.show() 