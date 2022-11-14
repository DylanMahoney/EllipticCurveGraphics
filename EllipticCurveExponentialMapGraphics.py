#This code computes and visualizes a homomorphism (in particular, a type of exponential map) between the complex numbers (as an additive group) and an elliptic curve over the complex numbers
import math
import numpy as np
import matplotlib.pyplot as plt
import random
from decimal import *

modes = ['first_component,second_component,poles'] #Different kinds of graph it can make

a=-1 #The coefficient a in the Weierstrass form of the elliptic curve

b = 0 #The coefficient a in the Weierstrass form of the elliptic curve
n = 12 #Log base 2 of how many terms to go to in the sequence when calculating the value of the homomorphism

Px = 1 #The x-component of the point which we're treating as the identity of the elliptic curve

Py = 0 #The y-component of the point which we're treating as the identity of the elliptic curve
how_far_out = 40 #The graph will be on the square in the complex plane centered at the origin with side length double this variable (with the sides of the square parallel to the imaginary and real axes)
mode='first_component'

#Checking that the point P is on the elliptic curve
print("%s should be equal to %s" % (Py**2,Px**3+a*Px+b))

#Check that the elliptic curve is non-singular
print(str(4*a**3 + 27*b**2))

def lengthofvector(point):
  return float(abs(point[0])**2 + abs(point[1])**2)
def addtwopoints(point1,point2):
  #Let's make sure we're dealing with the right data type
  x1 = complex(point1[0])
  y1 = complex(point1[1])
  
  x2 = complex(point2[0])
  y2 = complex(point2[1])

  #Comes from of elliptic curve operation
  if (x1,y1) == (0,math.inf):
    return point2
  if (x2,y2) == (0,math.inf):
    return point1
  if (x1,y1) == (x2,y2):
    return doublepoint(point1)
  if x1 == x2:
    return (0,math.inf)
  m = (y2-y1)/(x2-x1)
  newx = m**2 - x1 - x2
  return (newx,m*(x1-newx)-y1)

def doublepoint(point):
  #Comes from definition of elliptic curve operation
  x = point[0]
  y = point[1]

  if y==0 or y==math.inf or y==0+0j:
    return (0,math.inf)
  m = (3*x**2+a)/(2*y)
  newx = m**2-2*x
  #Modify elliptic curve operation so that P is the identity instead of infinity
  return addtwopoints(((newx, m*(x-newx)-y)),(Px,-Py))

def multiplypointbytwotothen(point,n):
  #Recursive function based on doubling a point
  if n==0:
    return point
  point = doublepoint(point)
  return multiplypointbytwotothen(point,n-1)

def f(z):
  #Given a complex number z, give the number in C^2 corresponding to a displacement of z in the tangent space of E at P
  return (Px+2*Py*z,Py+(3*Px**2+a)*z)

def phi(point):
  #Subtracts P from "point" (w.r.t the elliptic curve operation)
  x = point[0]
  y = point[1]
  if (x,y) == (0,math.inf):
    return (Px,-Py)
  if (complex(x),complex(y))==(complex(Px),complex(Py)):
    return (0,math.inf)
  if (complex(x),complex(y))==(complex(Px),complex(-Py)):
    return doublepoint((x,y))
  if complex(x)==complex(Px):
    return (0,math.inf)
  return (((y+Py)**2-(x+Px)*(x-Px)**2)/((x-Px)**2),(-(y+Py)**3+x*(y+Py)*(x-Px)**2+2*Px*(y+Py)*(x-Px)**2+Py*(x-Px)**3)/((x-Px)**3))

def homomorphism(z,i):
  #Divide the complex number by 2^i, move a distance of z from P in the tangent plane of E at P, multiplying the resulting point by itself 2^i times with respect to the elliptic curve operation (even though it's not quite on the elliptic curve the error is small) then changing coordinates so that the point P goes to infinity to match the convention for the Weierstrass p-function
  return phi(multiplypointbytwotothen(f(z/2**i),i))

def checkhomomorphism(z,w):
  #Just the definition of a homomorphism
  print("%s should be equal to %s" % (homomorphism(z+w,n),addtwopoints(homomorphism(z,n),homomorphism(w,n))))

def randomlycheckhomomorphism(min,max):
  #Choose some pairs of random complex numbers in a specified domain, check whether the homomorphism is actually a homomorphism
  #Susceptible to floating point arithmetic errors and difficulties with how many terms in the homomorphism sequence are really necessary for a good approximation
  randoms = [0,0,0,0]
  for i in range(4):
    randoms[i] = random.uniform(min,max)
  z = randoms[0] + randoms[1]*1j
  w = randoms[2] + randoms[3]*1j
  print(z,w)
  checkhomomorphism(z,w)

def randomlycheckderivative():
  #Choose some pairs of random complex numbers in a specified domain, check whether the derivative of the first component of the homomorphism is approximately double the second component (as we would expect)
  #Susceptible to floating point arithmetic errors and difficulties with how many terms in the homomorphism sequence are really necessary for a good approximation
  randoms = [0,0,0,0]
  randoms[0] = random.uniform(-2,2)
  randoms[1] = random.uniform(-2,2)
  randoms[2] = random.uniform(-0.005,0.005)
  randoms[3] = random.uniform(-0.005,0.005)
  startingPoint = randoms[0] + 1j*randoms[1]
  endingPoint = startingPoint + randoms[2] + 1j*randoms[3]
  z = homomorphism(startingPoint,n)
  w = homomorphism(endingPoint,n)
  approximateDerivative = (w[0]-z[0])/(randoms[2]+1j*randoms[3])
  print("%s should do a good job of approximating %s" % (2*z[1],approximateDerivative))

#Do the checks
print("Homomorphism checks:")
for i in range(10):
  randomlycheckhomomorphism(-1,1)
print("Derivative checks:")
for i in range(4):
  randomlycheckderivative()
#See what the homomorphism gives for some values I chose
testvalues = [-0.01,0.01,-0.5,0.5,-1,1,0+1j,0-1j,1.99,2]
for z in testvalues:
  print("z = %s" % z)
  print(homomorphism(z,n))

#Create a graph of the homomorphism
grid_length = 400
#Generate grid points and compute relevant number at each grid point
xs,xe,rx, ys,ye,ry = -how_far_out,how_far_out,grid_length, -how_far_out,how_far_out,grid_length
x,y = np.ogrid[xs:xe:1j*rx, ys:ye:1j*ry]
z = x - 1j*y
for i in range(grid_length):
  for j in range(grid_length):
    input = z[i][j]
    if mode=='poles':
      #We impose a cutoff so that the graph is more readable
      z[i][j] = min(lengthofvector(homomorphism(input,n)),1e6)
    elif mode=='first_component':
      z[i][j] = complex(homomorphism(input,n)[0])
    elif mode=='second_component':
      z[i][j] = complex(homomorphism(input,n)[1])

z = z.T #So that the real and imaginary axes are oriented correctly
if mode != 'poles':
  #Assign an rgba value based on the complex number
  r = (np.cos(np.angle(z))+1)/2
  g = (np.cos(np.angle(z) + 2*np.pi/3)+1)/2
  bb = (np.cos(np.angle(z) + 4*np.pi/3)+1)/2
  aa = np.arctan(np.abs(z).astype('float'))

#Make the graph
fig, ax = plt.subplots()
if mode=='poles':
  ax = plt.imshow(z.T.astype('float'),cmap='jet')
  plt.title('abs(exp_1(z))^2 + abs(exp_2(z))^2')
else:
  ax = plt.imshow(np.dstack((r,g,bb,aa)).astype('float'))
  if mode == 'first_component':
    plt.title('arg(exp_1(z))')
  if mode == 'second_component':
    plt.title('arg(exp_2(z))')

plt.xlabel('Re(z)')
plt.ylabel('Im(z)')

x = np.linspace(0, grid_length, 5)
plt.xticks(x,how_far_out*(x-(grid_length/2))/(grid_length/2))
plt.yticks(np.flip(x),-how_far_out*(np.flip(x)-(grid_length/2))/(grid_length/2))

#Save the graph
print("a=%sb=%sr=%sn=%s%s_new" % (a,b,how_far_out,n,mode))
plt.savefig(fname = "a=%sb=%sr=%sn=%s%s_new" % (a,b,how_far_out,n,mode))
