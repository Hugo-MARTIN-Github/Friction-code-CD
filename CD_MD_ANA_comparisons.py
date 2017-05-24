#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

####################### load data files of both python and modygs codes and compare them with analytic solution
# language: python 3.4
# 10 april 2017
# Hugo Martin

################################################################################################################
############################### Imports ########################################################################

import matplotlib.pyplot as plt
import numpy as np
import math 
import matplotlib.patches as mpatches
from initialization import *

################################################################################################################
# A LITTLE SCRIPT IS BELOW THE FUNCTIONS.







# number_last_file is the number of the last file to compare with
# the arithmetic reason of the increasing of files

def compare(number_last_file, reason, discs_radius, side_length, number_walls, xmin, xmax, ymin, ymax, just_ana, display):

	for file_number in np.arange(0 ,number_last_file+reason ,reason , dtype=int):

		string_file_number = str(file_number)

		if(just_ana==0):
			# extracting python code data
			python_file_name = "./temporary_files/python_data_" + string_file_number + ".dat"
			python_positions = np.loadtxt(python_file_name) 
			number_discs = int(np.size(python_positions)/3)

			# extracting modygs data
			modygs_file_name = "./temporary_files/modygs_data_" + string_file_number + ".dat" # names definitions
			modygs_positions_temp = np.matrix(np.loadtxt(modygs_file_name)) # data extractions
			modygs_positions = np.concatenate((\
					np.asarray(modygs_positions_temp[:,1:3]),\
					np.asarray(modygs_positions_temp[:,6:7])), axis=1) # getting informations you need
			modygs_positions = np.reshape(modygs_positions,(1,3*number_discs))
		
		# creating analytic data
		real_time = file_number*time_step
		q_analytic = analytic(analytic_solution, real_time, slope, mu_friction, gravity_constant, q_0, u_0, discs_radius)
		np.savetxt("./temporary_files/analytic_data_" + str(file_number) + ".dat", q_analytic, fmt='%.15e')

		#print(' python\n',python_positions,'\n modygs\n',modygs_positions,'\n\n',modygs_positions[0,2])
		#print(1/0.)

		
		# At this stage, we've got our data well arranged. Each row is a disc.
		# We can launch the display.
		if(display==1):
			display_comparison(q_analytic, python_positions, modygs_positions, discs_radius, side_length, number_discs,\
				 number_walls, xmin, xmax, ymin, ymax, file_number)




############################### diplay_comparison ##############################################################
def display_comparison(analytic_positions, python_positions, modygs_positions, discs_radius, side_length, number_discs,\
			 number_walls, xmin, xmax, ymin, ymax, file_number):

	fig, ax = plt.subplots(1,1)
	R = discs_radius
	L = side_length

	# looping on data
	for i in range(number_discs):
		
		# making a blue Python code disc
		qi_py = python_positions[2*i]
		qip_py = python_positions[2*i+1]
		theta_py = python_positions[2*number_discs+i]
		
		x_0_py = discs_radius*math.cos(theta_py)
		y_0_py = discs_radius*math.sin(theta_py)
		
		plt.gca().add_patch(mpatches.Circle((qi_py, qip_py), discs_radius, facecolor="white", edgecolor="b", linewidth=1))
		plt.plot([qi_py,qi_py+x_0_py],[qip_py,qip_py+y_0_py],'b')

		# making a red Modygs code disc
		qi_mo = modygs_positions[0,2*i]
		qip_mo = modygs_positions[0,2*i+1]
		theta_mo = modygs_positions[0,2*number_discs+i]
		
		x_0_mo = discs_radius*math.cos(theta_mo)
		y_0_mo = discs_radius*math.sin(theta_mo)
		
		plt.gca().add_patch(mpatches.Circle((qi_mo, qip_mo), discs_radius, facecolor="white", edgecolor="red", linewidth=1))
		plt.plot([qi_mo,qi_mo+x_0_mo],[qip_mo,qip_mo+y_0_mo],'red')

		# making a green analytic code disc
		qi_ana = analytic_positions[2*i]
		qip_ana = analytic_positions[2*i+1]
		theta_ana = analytic_positions[2*number_discs+i]
		
		x_0_ana = discs_radius*math.cos(theta_ana)
		y_0_ana = discs_radius*math.sin(theta_ana)
		
		plt.gca().add_patch(mpatches.Circle((qi_ana, qip_ana), discs_radius, facecolor="white", edgecolor="g", linewidth=1))
		plt.plot([qi_ana,qi_ana+x_0_ana],[qip_ana,qip_ana+y_0_ana],'g')

	plt.axis('equal') # to get nice circles

	# Walls
	ax.plot([-2*R, L-2*R],[-2*R, -2*R],'lightgray', label='Walls') # bottom
	if number_walls>=3:
		ax.plot([L-2*R, L-2*R],[-2*R, L-2*R],'lightgray'); # plot right wall
	if number_walls>=2:
		ax.plot([-2*R, -2*R],[-2*R, L-2*R],'lightgray'); # plot left wall

	# legende
	lightgray_patch = mpatches.Patch(color='lightgray', label='Walls')
	blue_patch = mpatches.Patch(color='b', label='Python code', linewidth=1)
	red_patch = mpatches.Patch(color='red', label='Modygs')
	g_patch = mpatches.Patch(color='g', label='Analytic', linewidth=1)
	plt.legend(handles=[blue_patch, lightgray_patch, red_patch, g_patch])
   
	# making and saving
	ax.set(xlabel='x',ylabel='y',xlim=(xmin,xmax),ylim=(ymin,ymax))
	plt.savefig("./outputs/comparison_" + str(file_number) + ".png")
	
	plt.close(fig)




############################### analytic_solution ##############################################################

def analytic(analytic_solution, real_time, slope, mu_friction, gravity, q_0, u_0, discs_radius):

	if(analytic_solution == 1 or analytic_solution == 2):
		
		if (math.tan(slope) <= 3.*mu_friction):

			x = 1/3*gravity*math.sin(slope)*real_time*real_time+q_0[1]+u_0[1]*real_time+discs_radius
			z = q_0[2]+u_0[2]*real_time-discs_radius
			theta = -1/(3*discs_radius)*gravity*math.sin(slope)*real_time*real_time+q_0[3]+u_0[3]*real_time
			return [x,z,theta]
		else:

			x = 1/2*gravity*(math.sin(slope)-mu_friction*math.cos(slope))*real_time*real_time+q_0[1]+u_0[1]*real_time+discs_radius
			z = q_0[2]+u_0[2]*real_time-discs_radius
			theta = -mu_friction/discs_radius*gravity*math.cos(slope)*real_time*real_time; 
			return [x,z,theta]

		

	return [100.,100.,100]
		















# script's below






























############################## comparison script ###############################################################
# General parameters
analytic_solution = 1 # See execution_dc.py file to know what analytic_solution means.
one_time_step = 0 
time_step = 1.e-2 # take it as in execution_dc.py file
just_ana = 1 # just creates ana files
display = 0 # creates display png files

# Parameters to compare with modygs code.
last_file_number = 200
iteration_step = 10




# other things
number_discs_side = 3


print("Starting the comparison:\n")

discs_radius, number_discs, discs_mass, gravity_constant, final_time, slope, mu_friction,\
	number_walls, side_length, xmin, ymin, xmax, ymax\
                = ana_sol_param_ini(analytic_solution,one_time_step,time_step) # setting variable parameters

q_0, u_0 = position_and_velocities_initialization(discs_radius, number_discs, number_discs_side, analytic_solution,\
                                                   side_length)

compare(last_file_number,iteration_step, discs_radius, side_length, number_walls, xmin, xmax, ymin, ymax, just_ana, display) # compare

print("End of comparison.")