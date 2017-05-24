#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

####################### This file contains all model functions #################################################
# language: python 3.4
# 06 march 2017
# Hugo Martin

############################### imports ########################################################################
import numpy as np
import scipy.optimize as sc
import mosek as mo
import itertools as it
from conversion_for_mosek import conv_mosek
from mosek_function import use_mosek

####################### model_of_february_2017 #################################################################
# the key point of this model is the set of binary vectors z to take into account the orthogonality between the 
# two vectors of normal forces and noninterpenetration constraint

# When we are in this function we are sure there exists a constraint which is wrong

def model_of_february_2017(number_discs, number_contacts, vec_D, vec_U, mat_B, mat_T, mass_matrix_i,time_step, mu_friction, epsilon, print_info, vec_D_temp):

	one = np.ones(1)[0]
	zero = np.zeros(1)[0] # to get double precision (but it's not important actually)
	set_Z = np.asarray(list(it.product([one,zero],repeat=number_contacts))) # we create all vectors z
	sol_set_value = [] # set of value of the quadratic function for solution
	SOL_SET = [] # set of solution vectors
	
	# usefull definition without a depedency in mat_Z
	mat_A = (mat_B+mat_T) # definition of mat_A
	M_A = np.dot(mass_matrix_i,mat_A.T)	
	T_u2 = np.sum([np.dot(time_step,np.dot(mat_B,vec_U)),-vec_D],axis=0)
	K = np.dot(time_step*time_step,np.dot(mat_B,M_A))	
	temp_T = np.dot(mat_T,np.dot(mass_matrix_i,mat_A.T)) 
	T_U = np.dot(mat_T,vec_U)

	# operators for linear constraints without depedency in mat_Z

	LIN = np.asarray(K[:number_contacts,:]) # interpenetration LIN
	BIGD = np.asarray(T_u2[:number_contacts]) # interpenetration BIGD 
	LIN = np.concatenate((LIN, np.concatenate((np.dot(mu_friction,\
		np.eye(number_contacts)), -np.eye(number_contacts)), axis=1) ),axis=0) # Coulomb1 LIN
	BIGD = np.concatenate((BIGD,np.zeros(number_contacts)),axis=0) # Coulomb1 BIGD
	LIN = np.concatenate((LIN, np.concatenate((np.dot(mu_friction,\
		np.eye(number_contacts)), np.eye(number_contacts)), axis=1) ),axis=0) # Coulomb2 LIN
	BIGD = np.concatenate((BIGD,np.zeros(number_contacts)),axis=0) # Coulomb2 BIGD

	for vec_Z in set_Z: # we take a binary vector z in set_Z

		vec_Z_test = int(1) # par defaut je veux bien travailler
		for elem in range(number_contacts):
			if((vec_D_temp[elem]<=1e-15)and(vec_Z[elem])==one): # mais si une contrainte et violee et qu en plus on a desactive le contact dans vec_Z
					vec_Z_test = int(0) # alors je sais qu il ne faut pas travailler

		if(vec_Z_test==1): # si il faut que je travaille
		
			mat_ZT = np.diag(np.concatenate((np.zeros(number_contacts), np.sum([np.ones(number_contacts),-vec_Z],axis=0)),axis=0))
			temp_Z = np.concatenate((np.sum([np.ones(number_contacts),-vec_Z],axis=0),np.zeros(number_contacts)),axis=0)
			temp_B = np.dot(mat_B.T,(np.diag(np.concatenate((np.sum([np.ones(number_contacts),-vec_Z],axis=0),np.zeros(number_contacts)),axis=0))).T)

			number_of_one = int(np.dot(vec_Z,np.ones(number_contacts))) # I take the number of ones in vec_Z

			# operators for objective function (with a depedency in mat_Z)   { min 1/2 x.Qx + c.x + f }
			Q = np.dot(time_step*time_step,\
						np.dot(epsilon,  np.dot(np.dot(temp_B.T,mass_matrix_i),temp_B))  \
						+ np.dot(np.dot(temp_T.T,mat_ZT.T),temp_T)   \
						 ) # with a depedency in mat_Z
	
			c = -np.dot(time_step, np.dot(np.dot(temp_T.T,mat_ZT),T_U)) # with a depedency in mat_Z
			f = np.dot(0.5,np.dot(np.dot(mat_ZT,T_U),T_U)) # with a depedency in mat_Z may be not necessary
	
			# operators for linear constraints with depedency in mat_Z
			LIN_final = LIN
			BIGD_final = BIGD
			minus_one_constraint = int(0) # it s a artafact variable to take into account of the possibility to have one constraint remove
			if number_of_one<number_contacts:
				if number_of_one==int(0): # all z components are 0
					minus_one_constraint = int(1)
					LIN_final = np.concatenate((LIN_final, np.matrix(-np.dot( K.T, temp_Z))),axis=0) # ortho 2 for w LIN
					BIGD_final = np.concatenate((BIGD_final,-np.dot(temp_Z, T_u2)*np.ones(1)),axis=0) # ortho 2 for w BIGD
				else: # there are both zeros and ones
					LIN_final = np.concatenate((LIN_final, np.matrix(np.concatenate((-vec_Z,np.zeros(number_contacts)),axis=0))),axis=0) # ortho 1 for lambda LIN		
					BIGD_final = np.concatenate((BIGD_final,np.zeros(1)),axis=0) # ortho 1 for lambda BIGD
					LIN_final = np.concatenate((LIN_final, np.matrix(-np.dot( K.T, temp_Z))),axis=0) # ortho 2 for w LIN
					BIGD_final = np.concatenate((BIGD_final,-np.dot(temp_Z, T_u2)*np.ones(1)),axis=0) # ortho 2 for w BIGD
			else: # all z components are one 
				minus_one_constraint = int(1)
				LIN_final = np.concatenate((LIN_final, np.matrix(np.concatenate((-vec_Z,np.zeros(number_contacts)),axis=0))),axis=0) # ortho 1 for lambda LIN		
				BIGD_final = np.concatenate((BIGD_final,np.zeros(1)),axis=0) # ortho 1 for lambda BIGD 

			# conversion of elements 
			sizeQ = np.size(Q,axis=0)
			ligneLIN = np.size(LIN_final,axis=0)
			colonne_LIN = np.size(LIN_final,axis=1)
	
			asub, aval, qsubi, qsubj, qval = conv_mosek(Q, LIN_final, sizeQ, ligneLIN, colonne_LIN)

			lambda_n_temp = use_mosek(asub, aval, qsubi, qsubj, qval, c, BIGD_final, number_contacts, -minus_one_constraint, print_info)

			if lambda_n_temp!='no_solution':
				sol_set_value.append(np.dot(0.5,np.dot(lambda_n_temp,np.dot(Q,lambda_n_temp)))+np.dot(c,lambda_n_temp)+f)
				SOL_SET.append(lambda_n_temp)
		
		################### END OF LOOP	

	sol_set_value = np.asarray(sol_set_value, dtype=np.float64)
	min_value = np.min(sol_set_value)

	if min_value>=0:
		argmin = np.argmin(sol_set_value)
		lambda_n = np.asarray(SOL_SET[argmin])
	else:
		print("NO SOLUTION FOUND !")
		return 'no_solution'

	vecteur_force = -np.dot(mat_A.T,lambda_n)
	u_n = np.sum([vec_U, np.dot(time_step, np.dot(mass_matrix_i,vecteur_force))], axis=0)

	return u_n, vecteur_force
