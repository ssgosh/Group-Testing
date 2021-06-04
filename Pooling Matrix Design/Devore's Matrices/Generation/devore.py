# This code has been reproduced from https://github.com/pandey-tushar/Compressed-Sensing

import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)

def poly(n,prime,r,p):
    ind = int(n/p)
    rem = n%p
    ans = rem
    i = 1
    while(ind != 0 and i <= r):
        ans+=(ind%p)*(prime**i)
        i+=1
        ind = int(ind/p)    
    return ans%p

def devore(p,r):
    if r>=p:
        print('r should be less than p')
        return None
    else:
        matrix = np.zeros((p*p, p**(r+1)),dtype='int')    
        for i in range(p**(r+1)):
            for j in range(p):
                val = poly(i,j,r,p)
                row_num = val + j*p
                matrix[row_num][i] = 1
        return matrix

mat = devore(7,2)
print(mat)
