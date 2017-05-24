#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

# this file contains all the functions to construct the essentials vectors and matrices common with all 
# different algorithms
# language: python 3.4
# 30 march 2017
# Hugo Martin

############################### imports ########################################################################
import numpy as np
import ctypes as ct





############################### build_D_U_B_T ##################################################################

# this function build the vector vec_D with distances in the normal directions, and the vector Uof regular 
# and smooth forces on a time step. it uses the function vector_D and vector_F writted in the matrices.c file
def build_D_U_B_T(q_n, number_discs, number_contacts, discs_radius, side_length, number_walls,\
                             gravity_constant, slope, time_step, mass_matrix, u_n):

	size = 3*number_discs
	vec_D = np.zeros(number_contacts)
	vec_D_temp = np.zeros(number_contacts)
	vec_F = np.zeros(size)
	mat_B = np.zeros(2*number_discs*number_contacts)
	mat_T = np.zeros(3*number_discs*number_contacts)

	clib = ct.cdll.LoadLibrary('./matrices.so') # definition of c library
	vec_type = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')

	clib.vector_D.argtypes = (vec_type, vec_type, ct.c_uint, ct.c_uint, ct.c_double, ct.c_double, ct.c_uint)
	clib.vector_F.argtypes = (vec_type, ct.c_uint, ct.c_double, ct.c_double)
	clib.matrix_B_T_transp.argtypes = (vec_type, vec_type, vec_type, ct.c_uint, ct.c_uint, ct.c_double, ct.c_double, ct.c_uint)

	clib.vector_D(vec_D, q_n, number_discs, number_contacts, discs_radius, side_length, number_walls)
	clib.vector_F(vec_F, number_discs, gravity_constant, slope)
	clib.matrix_B_T_transp(mat_B, mat_T, q_n, number_discs, number_contacts, discs_radius, side_length, number_walls)

	vec_D = np.concatenate((vec_D,np.zeros(number_contacts)),axis=0)
	
	mat_B = -np.reshape(mat_B,[2*number_discs, number_contacts]) # we are getting the mat B transp

	mat_B = (np.concatenate((np.concatenate((mat_B, np.zeros([2*number_discs,number_contacts])),axis=1),\
                              np.zeros([number_discs,2*number_contacts])),axis=0)).T	

	mat_T = np.reshape(mat_T,[number_contacts,size]) # we are getting the mat T transp

	mat_T = (np.concatenate(((np.zeros([3*number_discs,number_contacts]),mat_T.T)),axis=1)).T

	vec_U = np.sum([u_n,np.dot(time_step,vec_F)],axis=0)

	q_n_temp = np.sum([q_n, np.dot(time_step, vec_U) ],axis=0)

	clib.vector_D(vec_D_temp, q_n_temp, number_discs, number_contacts, discs_radius, side_length, number_walls)

	#low_values_indices = np.absolute(vec_D) < 1e-15  # Where values are low
	#vec_D[low_values_indices] = 0.
	#low_values_indices = np.absolute(vec_U) < 1e-15  # Where values are low
	#vec_U[low_values_indices] = 0.
	#low_values_indices = np.absolute(mat_B) < 1e-15  # Where values are low
	#mat_B[low_values_indices] = 0.
	#low_values_indices = np.absolute(mat_T) < 1e-15  # Where values are low
	#mat_T[low_values_indices] = 0.
	#low_values_indices = np.absolute(vec_D_temp) < 1e-15  # Where values are low
	#vec_D_temp[low_values_indices] = 0.


	return vec_D, vec_U, mat_B, mat_T, vec_D_temp, vec_F
























