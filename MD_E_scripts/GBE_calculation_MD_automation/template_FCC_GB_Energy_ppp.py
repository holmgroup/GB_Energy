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

#------ CSL Unit Cell Calculator---------------------
def CSL(Sigma_str,left,right):
	# INPUT: 
	# Sigma_str: Sigma value as string, e.g. '3' for Sigma-3
	# left, right: Unnormalized orientation matrices for grains as 3x3 matrix, with rows as orthogonal crystallographic directions along simulation box [h k l]
	# Express left,right as lists, as below:
	# left =[[11, 8, 5],
	# 	  [1, -2, 1],
	# 	  [6, -2, -10]]
	# right = [[11, 8, -5],
	# 	  [-1, 2, 1],
	# 	  [6, -2, 10]]

	# OUTPUT: CSL unit cell lengths along x,y,z, orthogonal directions of simulation box 

	fcc_rotation_location =  os.path.join(os.path.dirname(inspect.getfile(GBpy)), 'pkl_files', 'cF_Id_csl_common_rotations.pkl')
	with open(fcc_rotation_location) as infile:
		fcc_rotations = pickle.load(infile)
	# misorientation for sigma 3 boundary
	sig3 = fcc_rotations[Sigma_str]['N'][0] / fcc_rotations[Sigma_str]['D'][0]
	# create lattice instance and multiply by lattice parameter of nickel
	Ni = lattice.Lattice()
	a_Ni = 3.52;
	a_Al = 4.032; #LEA (Liu, Ercolessi, Adams) potential on NIST
	a_ref = 1.;
	Ni.l_g_go *= a_ref
	basis = [[ 0. ,  0.5,  0.5],
	[ 0.5,  0. ,  0.5],
	[ 0.5,  0.5,  0. ]]
	# get primitive CSL vectors and convert them to orthogonal basis
	csl_p, __ = find_csl_dsc.find_csl_dsc(Ni.l_g_go, sig3)
	csl_o = np.dot(Ni.l_g_go, csl_p).T
	# take boundary plane from orientation matrix, use it to get
	# 2D CSL in boundary plane, and convert to orthogonal basis
	boundary_plane = left[0]
	gb =  bp_basis.gb_2d_csl(boundary_plane, sig3, Ni.l_g_go, 'normal_go', 'g1')
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

#------READING IN POTENTIAL FILES AND MATERIALS DATA ------------

# IMPORT data about interatomic potential. See text file for proper format.
# Add multiple rows to text file to consider multiple potentials. 

metaldata = pd.read_table('Params_FCC_Cu.txt',skiprows=3,delim_whitespace = True)
msort = metaldata

filler_Tm = 1600 #Filler Tm for those potentials where precise melting point is unknown (can be measured)
msort['T_m'].fillna(filler_Tm,inplace=True)

#------CSL Unit Cell Calculator---------------------

#Define left/right grain orientations, unnormalized
#rows = [h k l] along orthogonal simulation box reference frame

g1 = [[1, 1, 1],
	  [1, -1, 0],
	  [1, 1, -2]]
g2 = [[1, 1, 1],
	  [-1, 1, 0],
	  [-1, -1, 2]]


# Sigma value as string
# Calculate CSL unit cell, assuming boundary plane normal along x
Sigma_str = '3';
x1, y1, z1 = CSL(Sigma_str,g1,g2)
print x1, y1, z1

BP_str = '111' #119 in Olmsted Survey
file_string = ''.join(['Cu_E_out_Sigma_',Sigma_str,'_',BP_str])

# Define interatomic potential calls as used in lammps 
potential_strings = []
for i in range(0,len(msort)):
    curr_elt = msort['Element'][i]
    if curr_elt == 'Ni2':
        curr_pot = 'pair_coeff * * {}.eam {}'.format(curr_elt, 1)
    else: 
        curr_pot = 'pair_coeff * * {0}.eam {0}'.format(curr_elt)
    potential_strings.append(curr_pot)


#create directory, if one doesn't already exist
if not os.path.isdir(file_string):
   os.mkdir(file_string)
n_nodes = 1
submit_script_lines = ['#!/bin/bash',
                        '#SBATCH --job-name={}'.format(file_string),
                        '#SBATCH --output={}.stdout'.format(file_string),
                        '#SBATCH --ntasks={}'.format(n_nodes*32),
                        '#SBATCH --partition=batch',
                        '#SBATCH --exclusive',
                        '#SBATCH --mail-type=END',
                        '#SBATCH --mail-user=ichesser@andrew.cmu.edu',
                        # openmpi and LAMMPS require gcc
                        'module load gcc',
                        # the 2014-06-28 version of LAMMPS requres openmpi instead of mpich2
                        'module load openmpi',
                        # the 2014-06-28 version of LAMMPS is the only one with the ECO force 
                        'module load lammps/2014.06.28']
                        # line to execute scripts will be addeded iteratively

#-----------GENERATE LAMMPS SCRIPTS --------------
# Loop over potentials (materials)
for i_mat in range(0,len(msort)): #loop for material type
	elt = msort['Element'][i_mat] #string, e.g. 'Al'
	a = msort['a'][i_mat] #lattice parameter (angstroms)
	Tm = msort['T_m'][i_mat] #Melting T (K)
	num_x_periods = 10 #30,need to have enough repeats in this direction to get good mobility statistics
	xlen = x1
	num_y_periods = 5 #10,should be 10 to match Jon's thesis
	ylen = y1  #you can play around with correction factors here for CSL calculations
	num_z_periods = 5 #5, should be 5 to match Jon's thesis
	zlen = z1  

	# print temp
	lat_param = a
	atom_delete_cutoff = 0.5 * a #cutoff for deleting atoms that are too close
	name_string = ''.join([elt,'_'])
	xstal_string = ''.join([elt,'_',Sigma_str,'_',BP_str]) #to get xstallography specific files
	x_boundary = 0 #origin at half of box size

	count = 0

  # print count
  #define curr_name
	curr_name = ''.join([name_string,file_string])
	curr_restart = ''.join(['restart_',curr_name,'.equil'])
	min_name = ''.join([name_string,'0K_',str(count)])
	lammps_script_lines1 = [# create new logfile and set global parameters
              
              'log log-{}.lammps'.format(curr_name),
              'units metal',
              'atom_style atomic',
              'dimension 3',
              'boundary p p p',
              'variable a equal {}'.format(lat_param),
              'variable xcsl equal {}'.format(x1),
              'variable ycsl equal {}'.format(y1),
              'variable zcsl equal {}'.format(z1),
              'variable nx equal {}'.format(num_x_periods),
              'variable ny equal {}'.format(num_y_periods),
              'variable nz equal {}'.format(num_z_periods),
              'variable xdim equal ${a}*${nx}*${xcsl}',
              'variable ydim equal ${a}*${ny}*${ycsl}',
              'variable zdim equal ${a}*${nz}*${zcsl}',
              'variable xlo equal -${xdim}/2',
              'variable xhi equal ${xdim}/2',
              'variable ylo equal -${ydim}/2',
              'variable yhi equal ${ydim}/2',
              'variable zlo equal -${zdim}/2',
              'variable zhi equal ${zdim}/2',

              # define and create simulation box

              'region total block ${xlo} ${xhi} ${ylo} ${yhi} ${zlo} ${zhi}',
              'create_box 1 total',

              # set interatomic potential
  
              'pair_style eam/alloy']
	lammps_script_potential = potential_strings[i_mat]
	lammps_script_lines2 = [
	              # define regions for different grains and create atoms
	              'region grain1 block ${xlo} -0.001 ${ylo} ${yhi} ${zlo} ${zhi}',
	              'lattice fcc {} orient x {} {} {} orient y {} {} {} orient z {} {} {}'.format(lat_param,*np.ravel(g1)),
	              'create_atoms 1 region grain1',
	              'region grain2 intersect 2 total grain1 side out',
	              'lattice fcc {} orient x {} {} {} orient y {} {} {} orient z {} {} {}'.format(lat_param,*np.ravel(g2)),
	              'create_atoms 1 region grain2',
	              
	              #Delete overlapping atoms

	              'delete_atoms overlap {} all all'.format(round(1.6*lat_param/3.52,2)),  #value is somewhat arbitrary

	              #First minimization
	              'compute eng all pe/atom',
	              'compute eatoms all reduce sum c_eng',
	              'compute MinAtomEnergy all reduce min c_eng',
	              'compute csym all centro/atom fcc',

	              'thermo_style    custom step lx ly lz press pxx pyy pzz c_eatoms c_MinAtomEnergy',
	              # 'thermo_modify lost ignore flush yes',
	              'thermo                  1000',   #outputs data every 100 steps
	              'min_style cg',
	              'minimize 1.0e-8 1.0e-10 10000 10000',

	              #Second Minimization: minimizing pressure on simulation box
	              'fix 1 all box/relax x 0.1 y 0.1 z 0.1', #can only fix pressures in periodic dimensions
	              'min_style cg',
	              'minimize 1e-8 1e-10 10000 10000',
                  'variable Eo equal -3.540*atoms',
				  'variable Ef equal "c_eatoms"',
				  'variable Cf equal 1.60217657e-16',
				  'variable A equal (${ydim}*${zdim}*1e-20)',
				  'variable GBE equal ((${Ef}-${Eo})*${Cf})/(2*${A})',
				  'print "GB energy is ${GBE} mJ/m^2"',
	              'dump dump1 all custom 1 {}.dump id type x y z c_eng c_csym'.format(min_name),
	              'run 1',
	              'unfix 1',
	              'undump dump1'
              	]
    # write lammps script to file
	curr_runline = 'mpirun lammps -i {}.in'.format(curr_name)
	submit_script_lines.append(curr_runline)
	with open('{}.in'.format(curr_name), 'w') as outfile:
		outfile.write('\n'.join(lammps_script_lines1))
		outfile.write('\n')
		outfile.write((lammps_script_potential))
		outfile.write('\n')
		outfile.write('\n'.join(lammps_script_lines2))
		move('{}.in'.format(curr_name), file_string)

	copy('EAM/{}.eam'.format(elt), file_string) #get filepath right here!

# write lines to generate submission script on Hippolyta
with open('{}.submit'.format(file_string), 'w') as outfile:
    outfile.writelines('\n'.join(submit_script_lines))
move('{}.submit'.format(file_string), file_string)





