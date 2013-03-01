#!/usr/bin/env python
from itertools import permutations

var = {0:"n.gx", 1:"n.gy", 2:"n.gz"}

def block1():
    print """void make_C_block1(unsigned int r,
		      const DiffusionEncodingDirection &n,
		      vnl_matrix<double> &C)
{"""

    s = dict()
    l = []
    for i in range(0,3):
	for j in range(0,3):
	    index = None
	    for p in permutations([i,j]):
		if p in s:
		    index = p
		    break
	    if index == None:
		s[(i,j)] = 1
		l.append((i,j))
	    else:
		s[index] +=1

    c = 0
    codestr = "    C[%(row)s][%(col)i] = -%(factor)s*%(t1)s*%(t2)s;"
    for e in l:
	print codestr % {'row':'r',
			 'col':c,
			 'factor':s[e],
			 't1':var[e[0]],
			 't2':var[e[1]]}
	c += 1

    codestr = "    C[%(row)s][%(col)i] = 0;"
    while c < 21:
	print codestr % {'row':'r',
			 'col':c}
	c += 1

    print "}"

def block2():
    print """void make_C_block2(unsigned int r,
		      const DiffusionEncodingDirection &n,
		      vnl_matrix<double> &C)
{"""

    c = 0
    codestr = "    C[%(row)s][%(col)i] = 0;"
    while c < 6:
	print codestr % {'row':'r',
			 'col':c}
	c += 1

    s = dict()
    L = []
    for i in range(0,3):
	for j in range(0,3):
	    for k in range(0,3):
		for l in range(0,3):
		    index = None
		    for p in permutations([i,j,k,l]):
			if p in s:
			    index = p
			    break
		    if index == None:
			s[(i,j,k,l)] = 1
			L.append((i,j,k,l))
		    else:
			s[index] +=1

    codestr = "    C[%(row)s][%(col)i] = -%(factor)s*%(t1)s*%(t2)s*%(t3)s*%(t4)s;"
    for e in L:
	print codestr % {'row':'r',
			 'col':c,
			 'factor':s[e],
			 't1':var[e[0]],
			 't2':var[e[1]],
			 't3':var[e[2]],
			 't4':var[e[3]]}
	c += 1

    print "}"

def block3():

    print """void make_C_block3(unsigned int r,
				const DiffusionEncodingDirection &n,
				const double bmax,
				vnl_matrix<double> &C)
{"""

    s = dict()
    l = []
    for i in range(0,3):
	for j in range(0,3):
	    index = None
	    for p in permutations([i,j]):
		if p in s:
		    index = p
		    break
	    if index == None:
		s[(i,j)] = 1
		l.append((i,j))
	    else:
		s[index] +=1

    c = 0
    codestr = "    C[%(row)s][%(col)i] = (-3/bmax)*%(factor)s*%(t1)s*%(t2)s;"
    for e in l:
	print codestr % {'row':'r',
			 'col':c,
			 'factor':s[e],
			 't1':var[e[0]],
			 't2':var[e[1]]}
	c += 1

    s = dict()
    L = []
    for i in range(0,3):
	for j in range(0,3):
	    for k in range(0,3):
		for l in range(0,3):
		    index = None
		    for p in permutations([i,j,k,l]):
			if p in s:
			    index = p
			    break
		    if index == None:
			s[(i,j,k,l)] = 1
			L.append((i,j,k,l))
		    else:
			s[index] +=1

    codestr = "    C[%(row)s][%(col)i] = %(factor)s*%(t1)s*%(t2)s*%(t3)s*%(t4)s;"
    for e in L:
	print codestr % {'row':'r',
			 'col':c,
			 'factor':s[e],
			 't1':var[e[0]],
			 't2':var[e[1]],
			 't3':var[e[2]],
			 't4':var[e[3]]}
	c += 1

    print "}"


block1()
block2()
block3()
