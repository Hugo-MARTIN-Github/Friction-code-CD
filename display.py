#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

# language: python 3.4
# 31 march 2017
# Hugo Martin

# this file contains the display function to create outputs

############################### imports ########################################################################
import matplotlib.pyplot as plt
import numpy as np
import math 
import matplotlib.patches as mpatches


############################### current_time_diplay ############################################################
def current_time_diplay(vecteur_force, gravite, q_n,number_discs,side_length,discs_radius,\
					number_walls,xmin,xmax,ymin,ymax,file_number, display_forces):

	fig, ax = plt.subplots(1,1)
	R = discs_radius
	L = side_length
	affiche_poids = 0
	affiche_force_contact = 0

	for i in range(number_discs):
		qi = q_n[2*i]
		qip = q_n[2*i+1]
		theta = q_n[2*number_discs+i]

		if display_forces==1:
			fi = vecteur_force[2*i]/10.
			fip = vecteur_force[2*i+1]/10.
			fipp = vecteur_force[2*number_discs+i]
			
		else:
			fi = 0.
			fip = 0.
			fipp = 0.

		gi = gravite[0]/10.
		gip = gravite[1]/10.
		
		normg = math.sqrt(np.sum([gi*gi,gip*gip],axis=0))
		norm = math.sqrt(np.sum([fi*fi,fip*fip],axis=0))
		
		x_0 = discs_radius*math.cos(theta)
		y_0 = discs_radius*math.sin(theta)
		
		plt.gca().add_patch(mpatches.Circle((qi, qip), discs_radius, facecolor="white", edgecolor="b", linewidth=1))
		plt.plot([qi,qi+x_0],[qip,qip+y_0],'b')

		bool1 = normg>=1.e-14
		bool2 = norm>=1.e-14

		if(bool1):
			plt.arrow(qi, qip, gi, gip , width=0.1, color='sandybrown')
			affiche_poids = 1			
	
		if(bool2):
			plt.arrow(qi, qip, fi, fip ,width=0.1,color='mediumseagreen')
			plt.arrow(qi, qip-discs_radius, fipp, 0. ,width=0.05, color='r')
			affiche_force_contact = 1


	plt.axis('equal')

	ax.plot([-2*R, L-2*R],[-2*R, -2*R],'lightgray', label='Walls') # bottom

	if number_walls>=3:
		ax.plot([L-2*R, L-2*R],[-2*R, L-2*R],'lightgray'); # plot right wall
     
	if number_walls>=2:
		ax.plot([-2*R, -2*R],[-2*R, L-2*R],'lightgray'); # plot left wall



	# legende
	lightgray_patch = mpatches.Patch(color='lightgray', label='Walls')
	blue_patch = mpatches.Patch(color='b', label='Discs', linewidth=1)
	red_patch = mpatches.Patch(color='mediumseagreen', label='Normal forces / 10.')
	mediumseagreen_patch = mpatches.Patch(color='r', label='Friction forces')
	sandybrown_patch = mpatches.Patch(color='sandybrown', label='Weight / 10.')


	if(affiche_poids or affiche_force_contact):
		if((affiche_force_contact) and (affiche_poids)):
			plt.legend(handles=[red_patch, blue_patch, mediumseagreen_patch, lightgray_patch, sandybrown_patch])
		elif(affiche_force_contact):
			plt.legend(handles=[red_patch, blue_patch, mediumseagreen_patch, lightgray_patch])
		else:
			plt.legend(handles=[blue_patch, lightgray_patch, sandybrown_patch])
	else:
		plt.legend(handles=[blue_patch, lightgray_patch])
   
	ax.set(xlabel='x',ylabel='y',xlim=(xmin,xmax),ylim=(ymin,ymax))
	plt.savefig("./outputs/CD_plot_" + str(file_number) + ".png")
	plt.close(fig)