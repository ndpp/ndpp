#!/usr/bin/env sage
from sage.all import *

# Set the variables.
Nl = 11
(fc0, fc1, fc2, fc3, fc4, fc5,fc6,fc7,fc8,fc9,fc10,A) = var('fc0','fc1','fc2','fc3','fc4','fc5','fc6','fc7','fc8','fc9','fc10','A')
(fc11, fc12, fc13, fc14, fc15,fc16,fc17,fc18,fc19,fc20) = var('fc11','fc12','fc13','fc14','fc15','fc16','fc17','fc18','fc19','fc20')
#(fc0, fc1, fc2, fc3, fc4, fc5,fc6,A) = var('fc0','fc1','fc2','fc3','fc4','fc5','fc6','A')
umin = var('umin')
assume(A > 1)
assume(umin >= -1)
assume(umin <= 1)
assume(umin - 1 < 0)

fc = [fc0, fc1, fc2, fc3, fc4, fc5,fc6,fc7,fc8,fc9,fc10,fc11,fc12,fc13,fc14,fc15,fc16,fc17,fc18,fc19,fc20]
#fc = [fc0, fc1, fc2]
Nc = len(fc)
u = var('u')
w = (u**2 - 1 + u * sqrt(u**2+A**2-1))/A

# Initialize the Legendres, store them in Pw, Pu
Pw = [None for lcm in xrange(Nc)]
Pu = [None for l in xrange(Nl)]
for lcm in xrange(Nc):
  Pw[lcm] = gen_legendre_P(lcm,0,w)
for l in xrange(Nl):
  Pu[l] = gen_legendre_P(l,0,u)

# Set the change of variable term (dw/dmu)
#dwdu = (2.0*u + sqrt(u**2+A**2-1.0) + (u**2 / sqrt(u**2+A**2-1.0))) / A
dwdu = diff(w,u)

# Calculate Fc(w)
print('Getting CM Functions')
Fc = 0.0
for lcm in xrange(Nc):
  Fc = Fc + (0.5 * (2.0 * lcm + 1.0) * fc[lcm] * Pw[lcm])
  print('Completed Moment ' + str(lcm))


# Initialize the solution space
fl = [None for l in xrange(Nl)]   # the lab system moments
Fl = 0.0         # the lab system function

# Get the lab system function
Fl = Fc * dwdu
print('Getting Lab Moments')
# Get the lab system moments
for l in xrange(Nl):
  fl[l] = SR(integral(Fl*Pu[l],u,umin,1.0)).full_simplify().expand()
  print('Completed Moment ' + str(l))

print('Getting Jacobian')
# Calculate the jacobian matrix and print the results.    
outFile = "mu_min.txt"
f = open(outFile, 'w')
import sympy
jacobian = [[None for lcm in xrange(Nc)] for l in xrange(Nl)]
jacobian_sympy = [[None for lcm in xrange(Nc)] for l in xrange(Nl)]
for l in xrange(Nl):
  for lcm in xrange(Nc):
    jacobian[l][lcm] = fl[l].coeff(fc[lcm]).full_simplify()
    jacobian_sympy[l][lcm]=jacobian[l][lcm]._sympy_()
    f.write(str(jacobian[l][lcm]))
    f.write('\n')
  print(jacobian[l])
  f.write('\n')
f.close()
