#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

####################### This file contains the function to test a solution from mosek ##########################
# language: python 3.4
# 10 april 2017
# Hugo Martin

################################################################################################################
############################### Imports ########################################################################

import numpy as np
import ctypes as ct


############################### test_post_mosek ################################################################
def test_post_mosek(lambda_n, mat_A, vec_U, time_step, mass_matrix_i, q_n, number_discs, number_contacts, discs_radius, side_length, number_walls):

	test = int(1)
	vecteur_force = -np.dot(mat_A.T,lambda_n)
	u_n = np.sum([vec_U, np.dot(time_step, np.dot(mass_matrix_i,vecteur_force))], axis=0)
	q_n_temp = np.sum([q_n , np.dot(time_step,u_n)],axis=0)
	vec_D_temp = np.zeros(number_contacts)

	clib = ct.cdll.LoadLibrary('./matrices.so') # definition of c library
	vec_type = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')

	clib.vector_D.argtypes = (vec_type, vec_type, ct.c_uint, ct.c_uint, ct.c_double, ct.c_double, ct.c_uint)
	clib.vector_D(vec_D_temp, q_n_temp, number_discs, number_contacts, discs_radius, side_length, number_walls)
	
	for dist in vec_D_temp: # on parcours le vecteur distance temporaire et si l'une des composantes est négative alors on dit que la solution est mauvaise
		if dist < -1.e-8:
			print( vec_D_temp )
			wait = input("PRESS ENTER TO CONTINUE.")
			test = int(0)

	return test 