#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

####################### This file contains the model function of march #########################################
# language: python 3.4
# 06 april 2017
# Hugo Martin

################################################################################################################
############################### Imports ########################################################################

import numpy as np
import itertools as it
from conversion_for_mosek import conv_mosek
from mosek_function import use_mosek
from test_post_mosek import test_post_mosek


################################ model_of_march_2017 ############################################################
# 

def model_of_march_2017(number_discs, number_contacts, vec_D, vec_U, mat_B, mat_T, mass_matrix_i,time_step,\
				 mu_friction, epsilon, print_info, vec_D_temp, q_n, discs_radius, side_length, number_walls,\
				analytic_solution, file_number, choix_pour_ana_8):

	number_of_solution_test_ana_8 = -1

	set_Z = np.asarray(list(it.product([1.,0.],repeat=number_contacts))) # we create all vectors z
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

	
	LIN = np.concatenate((np.dot(mu_friction,\
		np.eye(number_contacts)), -np.eye(number_contacts)), axis=1) # Coulomb1 LIN
	BIGD = np.zeros(number_contacts) # Coulomb1 BIGD
	LIN = np.concatenate((LIN, np.concatenate((np.dot(mu_friction,\
		np.eye(number_contacts)), np.eye(number_contacts)), axis=1) ),axis=0) # Coulomb2 LIN
	BIGD = np.concatenate((BIGD,np.zeros(number_contacts)),axis=0) # Coulomb2 BIGD

	for vec_Z in set_Z: # we take a binary vector z in set_Z
		
		# if vec_D_temp[elem] <= 1e-15 I am sure that I need to activate the contact so vec_Z = 1.
		# if vec_Z_i = 1 then the contact number i is activated so we want : # mat_ZT_i=1 in objective function,
		# temp_Z_for_K_i = -1 for K WITH A MINUS,
		# temp_Z_for_lambda_i = -0 for lambda WITH A MINUS


		number_of_one = int(np.dot(vec_Z,np.ones(number_contacts))) # I take the number of ones in vec_Z

		# test pour savoir si vec_Z est a premiere vue bien construit
		vec_Z_test = int(1) # par defaut je veux bien travailler
		for elem in range(number_contacts):
			if((vec_D_temp[elem] <= 1e-15)and(vec_Z[elem]) == 0.): # mais si une contrainte et violee et qu en plus on a desactive le contact dans vec_Z
					vec_Z_test = int(0) # alors je sais d avance qu il ne faut pas travailler


		vec_Z_test = 1 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(vec_Z_test==1): # si il faut que je travaille (c'est a dire qu a premiere vue, vec_Z est satisfaisant au regard de vec_D_temp)
 
			# some usefull operators
			mat_ZT = np.diag(np.concatenate((np.zeros(number_contacts), vec_Z),axis=0)) # A matrix which is used in the objective function
			temp_Z_for_K = -np.concatenate((vec_Z, np.zeros(number_contacts)),axis=0) # A vector used in linear constraint for Ortho. constraint again K WITH A MINUS
			temp_Z_for_lambda = -np.concatenate((np.sum([np.ones(number_contacts),-vec_Z],axis=0)\
						,np.zeros(number_contacts)), axis=0) # A vector used in linear constraint for Ortho. constraint again lambda WITH A MINUS
			temp_B = np.dot(mat_B.T,(np.diag(np.concatenate((vec_Z, np.zeros(number_contacts)),axis=0))).T) # used in objective function

			


			# operators for objective function (with a depedency in mat_Z)   { min 1/2 x.Qx + c.x + f }
			Q = np.dot(time_step*time_step,\
						np.dot(epsilon,  np.dot(np.dot(temp_B.T,mass_matrix_i),temp_B))  \
						+\
					 np.dot(np.dot(temp_T.T,mat_ZT.T),temp_T)   \
						 ) # with a depedency in mat_Z
	
			c = -np.dot(time_step, np.dot(np.dot(temp_T.T,mat_ZT),T_U)) # with a depedency in mat_Z
			f = np.dot(0.5,np.dot(np.dot(mat_ZT,T_U),T_U)) # with a depedency in mat_Z may be not necessary
	



			# operators for linear constraints with depedency in mat_Z
			LIN_final = LIN
			BIGD_final = BIGD

			# here we have the non interpenetration constraint only for activated contacts
			# je me balade dans vec_Z, lorsque je trouve z_i = 1 alors j'affiche la ligne i dans K 
			# (je vais donc avoir un nombre de contraintes egales a number_of_one pour cette partie de LIN)
			for i in range(number_contacts):
				if(vec_Z[i]==1.):
					LIN_final = np.concatenate((LIN_final, np.matrix(K[i,:])), axis=0) # je rajoute la ligne i de K interpenetration LIN
					BIGD_final = np.concatenate((BIGD_final, np.dot(T_u2[i],np.ones(1))), axis=0) # interpenetration BIGD interpenetration BIGD 

			minus_one_constraint = int(0) # it s a artafact variable to take into account of the possibility to have one constraint remove
			if number_of_one<number_contacts:
				if number_of_one==int(0): # all z components are 0 so all contacts are non activated and I need just the lambda constraint
					minus_one_constraint = int(1)
					LIN_final = np.concatenate((LIN_final, np.matrix(temp_Z_for_lambda)), axis=0) # ortho 1 for lambda LIN		
					BIGD_final = np.concatenate((BIGD_final, np.zeros(1)),axis=0) # ortho 1 for lambda BIGD
				else: # there are both zeros and ones
					LIN_final = np.concatenate((LIN_final, np.matrix(temp_Z_for_lambda)), axis=0) # ortho 1 for lambda LIN		
					BIGD_final = np.concatenate((BIGD_final, np.zeros(1)),axis=0) # ortho 1 for lambda BIGD
					LIN_final = np.concatenate((LIN_final, np.matrix(np.dot(K.T, temp_Z_for_K))),axis=0) # ortho 2 for w LIN
					BIGD_final = np.concatenate((BIGD_final, np.dot(temp_Z_for_K, T_u2)*np.ones(1)),axis=0) # ortho 2 for w BIGD
			else: # all z components are one so all contacts are activated I need just the K constraint
				minus_one_constraint = int(1)
				LIN_final = np.concatenate((LIN_final, np.matrix(np.dot(K.T, temp_Z_for_K))),axis=0) # ortho 2 for w LIN
				BIGD_final = np.concatenate((BIGD_final, np.dot(temp_Z_for_K, T_u2)*np.ones(1)),axis=0) # ortho 2 for w BIGD
	

			# conversion of elements 
			sizeQ = np.size(Q,axis=0)
			ligneLIN = np.size(LIN_final,axis=0)
			colonne_LIN = np.size(LIN_final,axis=1)
	
			asub, aval, qsubi, qsubj, qval = conv_mosek(Q, LIN_final, sizeQ, ligneLIN, colonne_LIN)

			
			# resolution of the problem
			lambda_n_temp = use_mosek(asub, aval, qsubi, qsubj, qval, c, BIGD_final, number_contacts, -minus_one_constraint, print_info)

			# In this algorithme, all problems defined by vec_Z have a solution
			if lambda_n_temp!='no_solution':

				# test post solution
				test = test_post_mosek(lambda_n_temp, mat_A, vec_U, time_step, mass_matrix_i, q_n,\
							 number_discs, number_contacts, discs_radius, side_length, number_walls)

				if(test==1): # !!!!!!!!!!!!!!
					sol_set_value.append(np.dot(0.5,np.dot(lambda_n_temp,np.dot(Q,lambda_n_temp)))+np.dot(c,lambda_n_temp)+f)
					SOL_SET.append(lambda_n_temp)
					if(analytic_solution==7)or(analytic_solution==8):
						number_of_solution_test_ana_8 = number_of_solution_test_ana_8 + 1
						np.savetxt("./temporary_files/nonuniquess_" + str(file_number+number_of_solution_test_ana_8) + ".dat", np.concatenate((np.matrix(sol_set_value[-1]),np.matrix(vec_Z)),axis=1) , fmt='%.15e')
				print("\n\n\n La valeur du test : ",test,"\n\n Le vecteur z est : ",vec_Z,"\n\n Le vecteur sol_set_value est : ",sol_set_value,"\n\n")# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				wait = input("PRESS ENTER TO CONTINUE.")

		
		
################### END OF LOOP ################################################################################	

	sol_set_value = np.asarray(sol_set_value, dtype=np.float64)
	min_value = np.min(sol_set_value)
	argmin = np.argmin(sol_set_value)

	print('\n\n\n\n',sol_set_value,'\n\n\n\n')
	print("La valeur de argmin est : ",argmin,"\n\n")
	lambda_n = np.asarray(SOL_SET[argmin])
	if(analytic_solution==8):
		lambda_n = np.asarray(SOL_SET[choix_pour_ana_8])

	vecteur_force = -np.dot(mat_A.T,lambda_n)
	u_n = np.sum([vec_U, np.dot(time_step, np.dot(mass_matrix_i,vecteur_force))], axis=0)


	print('\n\n La valeur de la norme de la vitesse est : ',np.sqrt(np.dot(u_n,u_n)),'\n\n')

	return u_n, vecteur_force
