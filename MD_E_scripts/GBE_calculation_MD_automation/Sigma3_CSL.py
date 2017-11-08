from __future__ import division
import numpy as np
import pandas as pd
import os
import GBpy
from shutil import copy, move
from numpy.random import randint
from GBpy import lattice
from GBpy import bp_basis
from GBpy import find_csl_dsc
import inspect
import pickle
from itertools import product
# import matplotlib.pyplot as plt
# import GBpy
# from GBpy import lattice
# from GBpy import bp_basis
# from GBpy import find_csl_dsc
# import inspect
# import pickle
# from itertools import product

#------Sigma-3 CSL Unit Cell Calculator---------------------
def Sigma3CSL(left,right):
	fcc_rotation_location =  os.path.join(os.path.dirname(inspect.getfile(GBpy)), 'pkl_files', 'cF_Id_csl_common_rotations.pkl')
	with open(fcc_rotation_location) as infile:
		fcc_rotations = pickle.load(infile)
	# misorientation for sigma 3 boundary
	sig3 = fcc_rotations['3']['N'][0] / fcc_rotations['3']['D'][0]
	# create lattice instance and multiply by lattice parameter of nickel
	Ni = lattice.Lattice()
	a_Ni = 3.52;
	a_Al = 4.050; #4.032, LEA (Liu, Ercolessi, Adams) potential on NIST
	a_ref = 1.;
	Ni.l_g_go *= a_Al
	basis = [[ 0. ,  0.5,  0.5],
	[ 0.5,  0. ,  0.5],
	[ 0.5,  0.5,  0. ]]
	# get primitive CSL vectors and convert them to orthogonal basis
	csl_p, __ = find_csl_dsc.find_csl_dsc(basis, sig3)
	print csl_p
	csl_o = np.dot(Ni.l_g_go, csl_p).T
	# take boundary plane from orientation matrix, use it to get
	# 2D CSL in boundary plane, and convert to orthogonal basis
	boundary_plane = left[0]
	gb =  bp_basis.gb_2d_csl(boundary_plane, sig3, basis, 'normal_go', 'g1')
	#print "gb =",gb
	csl_2d = np.dot(Ni.l_g_go, gb[0]).T
	#print "csl_2d =",csl_2d

	# normalize orientation matrix to make it a rotation matrix,
	# and rotate both CSL and 2D CSL into crystal frame

	rotation_matrix = [list(np.array(i) / np.linalg.norm(i)) for i in left]
	#print rotation_matrix

	csl_2d = [np.dot(rotation_matrix, i) for i in csl_2d]
	csl = [np.dot(rotation_matrix, i) for i in csl_o]
	# iterate through the vertices of the CSL unit cell
	# and use the highest and lowest x coordinates
	# to calculate the period in the x direction

	xlo, xhi = 0, 0

	for i in product([0, 1], [0, 1], [0, 1]):
	    csl_vertex = np.sum(np.array(csl).T * i, axis=1)
	    xlo = min(xlo, csl_vertex[0])
	    xhi = max(xhi, csl_vertex[0])
	    
	x_period = xhi - xlo

	# iterate through the vertices of the 2D CSL unit cell
	# and use the highest and lowest y and z coordinates
	# to calculate the period in the y and z directions

	ylo, yhi, zlo, zhi = 0, 0, 0, 0

	for i in product([0, 1], [0, 1]):
	    csl_2d_vertex = np.sum(np.array(csl_2d).T * i, axis=1)
	    ylo = min(ylo, csl_2d_vertex[1])
	    yhi = max(yhi, csl_2d_vertex[1])
	    zlo = min(zlo, csl_2d_vertex[2])
	    zhi = max(zhi, csl_2d_vertex[2])
	
	y_period = yhi - ylo
	z_period = zhi - zlo

	return x_period, y_period, z_period

# grain1_rot_mat = [[11, 8, 5],
#                   [1, -2, 1],
#                   [6, -2, -10]]
# grain2_rot_mat = [[11, 8, -5],
#                   [-1, 2, 1],
#                   [6, -2, 10]]
# grain1_rot_mat = [[-1, 10, 11],
#                   [7, 4, 3],
#                   [1, -1, -1]]
# grain2_rot_mat = [[29, -1, -34],
#                   [19, 7, 16],
#                   [1, -5, 1]]
grain1_rot_mat = [[4, 1, 1],
                  [0, 1, -1],
                  [-2, 4, 4]]
grain2_rot_mat = [[1, 1, 0],
                  [1, -1, 0],
                  [0, 0, -6]]
g1_rot_mat = [[5, 5, 2],
                  [1, -1, 0],
                  [2, 2, -10]]
g2_rot_mat = [[2, 1, 1],
                  [0, -1, 1],
                  [6, -6, 6]]
x1, y1, z1 = Sigma3CSL(grain1_rot_mat,grain2_rot_mat)
print x1, y1, z1