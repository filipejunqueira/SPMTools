import numpy as np


f = np.array([5,5,10,15,20,25,30,40,50,70,90])
z = np.array([0,1,2,3,4,5,6,6.1,4,9,19])


dnum = np.ediff1d(f)
dden = np.ediff1d(z)
print(dnum - dnum[1:])
print(dden)
print(dnum/dden)



