#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

# language: python 3.4
# 31 march 2017
# Hugo Martin


# this file contains all initialization functions


############################### imports ########################################################################
import math
import numpy as np


############################### ana_sol_param_ini ##############################################################

# initialisation of simulation parameters if we considere an analitic solution. the different situations are 
# described in the DC_execution.py file with definition of the variable analytic_solution.
def ana_sol_param_ini(analytic_solution,one_time_step,time_step):
	discs_radius = 1.
	discs_mass = 3.
	gravity_constant = 10.
	side_length = 10.
	xmin = -2*discs_radius-1
	ymin = -3*discs_radius
	xmax = side_length-2*discs_radius+1
	mu_friction = 0.5
	slope = 0.
	number_walls = 1
	number_discs = 2

	if analytic_solution == 0: # free fall
		final_time = 3.
		number_discs = 1
		slope = (math.pi)/8.
		ymax = 7*discs_radius
	elif (analytic_solution == 1 or analytic_solution == 2): # rolling disc with slip and wihtout
		final_time = 2.
		number_discs = 1
		ymax = 2*discs_radius
		if analytic_solution == 1:
			slope = (math.pi)/3 # with slip
		else:
			slope = (math.pi)/8 # without slip
	elif (analytic_solution == 3 or analytic_solution == 4): # static and dynamic 3xdiscs
		number_discs = 3
		ymax = 4*discs_radius
		if analytic_solution == 3:
			final_time = 2.0
			mu_friction = 0.8
		else:
			final_time = 2.0
			mu_friction = 0.2
	elif (analytic_solution == 5 or analytic_solution == 6): # nonorthogonal and othogonal binar collision
		final_time = 4.
		number_discs = int(2)
		if analytic_solution == 6:
			gravity_constant = 0.
		ymax = 4*discs_radius

	elif (analytic_solution == 7):
		final_time = time_step
		number_discs = int(3)
		gravity_constant = 0.
		ymax = 4*discs_radius

	elif(analytic_solution==8):
		final_time = time_step
		number_discs = int(3)
		gravity_constant = 0.
		ymax = 4*discs_radius
	
	if one_time_step==1:
		final_time = time_step

	return discs_radius, number_discs, discs_mass, gravity_constant, final_time, slope, mu_friction,\
            number_walls, side_length, xmin, ymin, xmax, ymax



############################### variables_initialization #######################################################
# initialisation of some variables which are common with all simulation 

def variables_initialization(number_discs,number_walls,discs_mass,discs_radius):
	current_time = 0.
	number_contacts = int(number_discs*(number_discs-1)/2+number_walls*number_discs)  
	J = 0.5*discs_mass*discs_radius*discs_radius

	temp = np.ones(3*number_discs)
	temp[:2*number_discs] = discs_mass*temp[:2*number_discs]
	temp[2*number_discs:3*number_discs] = J*temp[2*number_discs:3*number_discs]

	mass_matrix = np.diag(temp)
	mass_matrix_i = np.diag(1./temp)

	file_number = 0

	return current_time, number_contacts, mass_matrix, mass_matrix_i, file_number



############################### position_and_velocities_initialization #########################################
# all is in the title

def position_and_velocities_initialization(discs_radius, number_discs, number_discs_side, analytic_solution,\
                                                   side_length):
	q_n = np.zeros(3*number_discs)
	u_n = np.zeros(3*number_discs)

	if analytic_solution == 0: # free fall
		q_n[0]= 0.
		q_n[1] = 6*discs_radius
          
	elif (analytic_solution == 1 or analytic_solution == 2): # rolling disc with slip and wihtout
		q_n[0] = 0.
		q_n[1] = -discs_radius
          
	elif (analytic_solution == 3 or analytic_solution == 4): # static and dynamic 3xdiscs
		q_n[0] = 0.5*side_length/3.
		q_n[1] = -discs_radius
		q_n[2] = q_n[0]+2*discs_radius
		q_n[3] = q_n[1]
		q_n[4] = q_n[0]+discs_radius
		q_n[5] = q_n[1]+discs_radius*np.sqrt(3)
                    
	elif (analytic_solution == 5 or analytic_solution == 6): # nonorthogonal and othogonal binar collision
		q_n[0] = 0.5*side_length/3.
		if(analytic_solution==5):
			q_n[1] = 3.-discs_radius/2.
		else:
			q_n[1] = 0.
		q_n[2] = q_n[0]+5*discs_radius
		q_n[3] = q_n[1]+discs_radius

		u_n[0] = 4.

		if analytic_solution == 6:
			q_n[3] = q_n[1]
			u_n[0] = 2.

	elif (analytic_solution == 7):
		q_n[0] = 0.5*side_length/3.
		q_n[1] = 1.6*discs_radius
		q_n[2] = q_n[0] #+2.3*discs_radius
		q_n[3] = q_n[1]+2.2*discs_radius
		q_n[4] = q_n[2]
		q_n[5] = q_n[1]-2.2*discs_radius

		u_n[3] = -10.
		#u_n[2*number_discs+2] = -2.

	elif (analytic_solution == 8):
		refx = 0.5*side_length/3.
		refy = 2.*discs_radius
		q_n[0] = refx+0.4*discs_radius
		q_n[1] = refy
		q_n[2] = refx+1.6*discs_radius
		q_n[3] = refy
		q_n[4] = refx+discs_radius
		q_n[5] = refy+discs_radius*np.sqrt(3)
		

	else: # discs are put in a box configuration
		pas = 2*discs_radius+discs_radius/6. 
		varx = 0. 
		varz = -pas # temporary variables to create the box

		for i in range(1,number_discs+1):
		
			if (i-1)%number_discs_side == 0:
				varz = varz + pas
				varx = 0. # pour pouvoir incrementer en z

			q_n[2*i-2] = varx+0.03*((i-1)%2)
			q_n[2*i-1] = varz+0.05*((i-1)%2) # initialisation avec une petite perturbation
			varx = varx + pas # incrementation en x

	return q_n, u_n



































