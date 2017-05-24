#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

################################################################################################################
############################### Imports ########################################################################

from initialization import *
from matrices_build import *
from display import current_time_diplay
from model_of_february_2017 import *
from model_of_march_2017 import *
from make_data_files import *
import math
import sys

################################################################################################################
############################### execution script for contact dynamics simulation ###############################
################################################################################################################
# language: python 3.4
# 16 May 2017
# Hugo Martin







############################### general variables for simulation type ##########################################

analytic_solution = 8 # analytic solution type {-1:not analytic solution, 0:free fall, 1:rolling disc with slip,
                      # 2:rolling disc without slip, 3:3xdiscs static, 4:3xdiscs dynamic, 5:nonorthogonal binar 
                      # collision, 6:orthogonal binar collision, 7:nonuniqueness of solution, 8:nonuniqueness with 3 beads}
choice_for_ana_8 = 3 # 0, 1, 2 or 3 the diffrent solutions you can take

one_time_step = 0 # if we want simulate over just one time step {1:yes, 0:no}
print_info = 0 # if we want to print informations in the resolution step
epsilon = 1.e-10 # the epsilon present in most of models
algo = 2 # select the algorithme : 1 for a loop on the all infeasible problems, 2 for a loop on all feasible problems
data_files = 1 # 1 if create data files, 0 otherwise
display_forces = 1 # 1 if display forces in png files, 0 otherwise

# Parameters to compare with modygs code.
last_file_number = 20000000000 # This variable is exactly the number of the last iteration in the execution.
iteration_step = 1 # The arithmetic reason for the data files.



############################### default parameters initialisation ##############################################

########### default physical parameters intialisation
# default physical parameters if there is no analytic solution selected

discs_radius = 1. # discs radius
number_discs = 4 # number of discs
discs_mass = 3. # discs mass
gravity_constant = 10. # gravity constant
time_step = 1.e-1 # time step
final_time = 3 # final time of simulation
slope = 0. # slope between the bottom and the horizontal direction
mu_friction = 0.5 # friction coefficient in Coulomb s law

########### default domain parameters initialisation
# default domain parameters if there is no analytic solution selected

number_walls = int(1) # 1=< number_walls =< 3 number of walls in the order (bottom, left, right)
number_discs_side = int(3) # number of discs by side if the domain is a box (number_walls = 3)
side_length = 7. # length of each box side

########### domain display parameters
xmin = -2*discs_radius-1.
ymin = -3*discs_radius
xmax = side_length-2*discs_radius+1.
ymax = side_length-2*discs_radius








































################################################################################################################
############################### Code ###########################################################################
################################################################################################################


############################### program launch #################################################################
print('')
print('')
print("                      *** contact dynamics program launched ***")
print('')

############################### parameters initialisation for analytic solutions ###############################
if (analytic_solution>8 or analytic_solution<-1):
	sys.exit("> WARNING : ANALYTIC SOLUTION DOESN'T EXIST ! PROGRAM STOPPED !")
if analytic_solution!=-1:
	discs_radius, number_discs, discs_mass, gravity_constant, final_time, slope, mu_friction,\
	number_walls, side_length, xmin, ymin, xmax, ymax\
                = ana_sol_param_ini(analytic_solution,one_time_step,time_step) # initialisation function
	print(">  analytic solution is ON")
else:
	print(">  analytic solution is OFF")


############################### initialization of essential variables ##########################################
current_time, number_contacts, mass_matrix, mass_matrix_i, file_number\
                = variables_initialization(number_discs,number_walls,discs_mass,discs_radius)
print(">  variables and constants are initialized")

############################### initialization of positions and velocities #####################################
q_n, u_n = position_and_velocities_initialization(discs_radius, number_discs, number_discs_side,\
                                    analytic_solution, side_length)
print(">  positions and velocities are initialized")

############################## first matrices construction #####################################################
vec_D, vec_U, mat_B, mat_T, vec_D_temp, vec_F = build_D_U_B_T(q_n, number_discs, number_contacts, discs_radius, side_length,\
               number_walls, gravity_constant, slope, time_step, mass_matrix, u_n)

############################### display of initial time and some informations ##################################
vecteur_force_i = np.zeros(3*number_discs)
vecteur_force = np.zeros(3*number_discs)
gravite = [vec_F[0], vec_F[1]]
current_time_diplay(vecteur_force_i, gravite, q_n, number_discs,side_length,discs_radius,number_walls,xmin,\
                              xmax,ymin,ymax,file_number,display_forces)

if(data_files==1):
	make_data_files(file_number, q_n) # making data files

print("\n>\tSome parameters :\n\t- Configuration is :",analytic_solution,"\n\t- One time step is :",one_time_step)
rep = 'y'
if((last_file_number*time_step < final_time)and(one_time_step==0)): # some checks on the time parameters
	message = "\n>\tWarning ! One_time_step is 0 and the final time calculated by ' Last_iteration*time_step ' : "+str(last_file_number*time_step)+" is lower than the ' final_time ' variable : "+str(final_time)+". Be sure of what you are doing. Do you want to continue ? (y/n)\n\t\t> "
	rep = input(message)
elif((last_file_number*time_step > final_time)and(one_time_step==0)):
	message = "\n>\tWarning ! One_time_step is 0 and the final time calculated by ' Last_iteration*time_step ' : "+str(last_file_number*time_step)+" is higher than the ' final_time ' variable : "+str(final_time)+". Be sure of what you are doing. Do you want to continue ? (y/n)\n\t\t> "
	rep = input(message)

if(rep=='n'):
	sys.exit(">\tExecution stopped.")

############################### update time and file_number ####################################################
current_time = current_time + time_step
file_number = file_number + 1
print('')

############################### time loop ######################################################################
# START of time loop
print(">  start of the time loop")
print('')
print(">  current time is : 0.000000")
while ((current_time <= final_time)and(file_number <= last_file_number)):

######## matrices construction
	vec_D, vec_U, mat_B, mat_T, vec_D_temp, vec_F = build_D_U_B_T(q_n, number_discs, number_contacts,\
					 discs_radius, side_length, number_walls, gravity_constant, slope,\
					 time_step, mass_matrix, u_n) # calculus of vector D, vector forces, matrix B and T

######## resolution of optimization problem

	work_to_do = int(0)
	for elem in range(number_contacts):
		if(vec_D_temp[elem]<=1e-15):
			work_to_do = int(1) # est-ce que je dois travailler ?

	if(work_to_do==1): # oui je le dois

		if(algo==1):
			u_n, vecteur_force = model_of_february_2017(number_discs, number_contacts, vec_D, vec_U, mat_B, mat_T, mass_matrix_i, time_step,\
					mu_friction, epsilon, print_info, vec_D_temp) # optimization problem resolution
		elif(algo==2):
			u_n, vecteur_force = model_of_march_2017(number_discs, number_contacts, vec_D, vec_U, mat_B, mat_T, mass_matrix_i, time_step,\
				mu_friction, epsilon, print_info, vec_D_temp, q_n, discs_radius, side_length,\
						 number_walls, analytic_solution, file_number, choice_for_ana_8) # optimization problem resolution
	else: # non pas besoin
		u_n = vec_U

######## updating of q_n
	q_n = np.sum([q_n , np.dot(time_step,u_n)],axis=0) # update of q_n

######## display
	if((data_files==1)and( (file_number % iteration_step) == 0 )):
		current_time_diplay(vecteur_force, vec_F, np.asarray(q_n),number_discs,side_length,discs_radius,number_walls,xmin,\
                             xmax,ymin,ymax,file_number,display_forces) # display of current time		
		make_data_files(file_number, q_n) # making data files


######## other things
	vecteur_force = vecteur_force_i
	print(">  current time is : %f" % current_time)
	current_time = current_time + time_step # update of current_time
	file_number = file_number + 1 # update of file_number
     

# END of time loop



############################### execution end ##################################################################
print('')
print("                        *** program end without any error ***")
print('')