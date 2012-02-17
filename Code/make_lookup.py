#!/usr/bin/env python

import numpy as num
from itertools import permutations

K = 100*num.ones((3,3,3,3))

count = 0
for i in range(0,3):
    for j in range(0,3):
        for k in range(0,3):
            for l in range(0,3):
                update = 0
                for p in permutations([i,j,k,l]):
                    if(K[p] > 21):
                        K[p] = 6+count
                        update += 1
                if(update > 0):
                    count += 1

print 'unsigned int D[][3] = {{0,1,2},{1,3,4},{2,4,5}};'

print 'unsigned int K[][3][3][3] = { \\'
for i in range(0,3):
    print '{ \\'
    for j in range(0,3):
        print '{ \\'
        for k in range(0,3):
            print '{', int(K[i,j,k,0]), ',', int(K[i,j,k,1]), ',', int(K[i,j,k,2]), '}, \\'
        print '}, \\'
    print '}, \\'
print '};'

