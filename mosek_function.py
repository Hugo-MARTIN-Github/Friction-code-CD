#! /home/hugo/anaconda3/envs/py34/bin/ipython

# the function to sole the optimization problem
# language: python 3.4
# 06 march 2017
# Hugo Martin

import sys
import os
import numpy as np

import mosek

# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:
inf = 0.0

# Define a stream printer to grab output from MOSEK
def streamprinter(text):
	sys.stdout.write(text)
	sys.stdout.flush()

# We might write everything directly as a script, but it looks nicer
# to create a function.
def use_mosek(asub, aval, qsubi, qsubj, qval, c, BIGD, number_contacts, number_of_one, print_info):
	# Open MOSEK and create an environment and task
	# Make a MOSEK environment
	with mosek.Env () as env:
		# Attach a printer to the environment
		env.set_Stream (mosek.streamtype.log, streamprinter)
		# Create a task
		with env.Task() as task:

			#task.putdouparam( param, param value) # value to control the resolution precision
			#print ('default value for parameter =',param)			

			task.set_Stream (mosek.streamtype.log, streamprinter) # attach a stream printer to the task
			numvar = 2*number_contacts # number of variables
			numcon = np.size(BIGD) # number of constraints

			# Set up and input bounds and linear coefficients
			bkc   = [ mosek.boundkey.lo ]*numcon # types of linear constraints
			blc   = [ 0.0 ]*numcon # lower values for linear constraints
			buc   = [ inf ]*numcon # upper values for linear constraints

			bkx   = [ mosek.boundkey.lo ]*number_contacts + [ mosek.boundkey.fr ]*number_contacts # types of bound constraints 
			blx   = [ 0.0 ]*number_contacts + [ inf ]*number_contacts # lower values for bound constraints
			bux   = [ inf ]*numvar # upper values for bound constraints

			# Appends
			# Append 'numcon' empty constraints.
			# The constraints will initially have no bounds.
			task.appendcons(numcon)

			# Append 'numvar' variables.
			# The variables will initially be fixed at zero (x=0).
			task.appendvars(numvar)


			# Inputs
			# Input linear objects of constraints and objective 
			for j in range(numvar):
				# Set the linear term c_j in the objective.
				task.putcj(j,c[j])

				# Set the bounds on variable j
				# blx[j] <= x_j <= bux[j]
				task.putbound(mosek.accmode.var,j,bkx[j],blx[j],bux[j])

				# Input column j of A
				task.putacol( j, # Variable (column) index.
					asub[j], # Row index of non-zeros in column j.
					aval[j]) # Non-zero Values of column j.

			# Input bound constraints
			for i in range(numcon):
				# Set the constant vector in linear constraints
				blc[i] = BIGD[i] 
				task.putbound(mosek.accmode.con,i,bkc[i],blc[i],buc[i])

			# Input the quadratic term in objective function
			task.putqobj(qsubi,qsubj,qval)

			# Input the objective sense (minimize/maximize)
			task.putobjsense(mosek.objsense.minimize)


			############## Solve and print
			# Optimize
			task.optimize()
			# Print a summary containing information
			# about the solution for debugging purposes
			if print_info==1:
				task.solutionsummary(mosek.streamtype.msg)
			# prosta = task.getprosta(mosek.soltype.itr) dont'see the interest
			solsta = task.getsolsta(mosek.soltype.itr)

			# Output a solution
			xx = [0.]*numvar
			task.getxx(mosek.soltype.itr,
				xx)

			if solsta == mosek.solsta.optimal or solsta == mosek.solsta.near_optimal:
				if print_info==1:
					print("Optimal solution: %s" % xx)
				return xx

			if print_info==1:
				if solsta == mosek.solsta.dual_infeas_cer:
					print("Primal or dual infeasibility.\n")
				elif solsta == mosek.solsta.prim_infeas_cer:
					print("Primal or dual infeasibility.\n")
				elif solsta == mosek.solsta.near_dual_infeas_cer:
					print("Primal or dual infeasibility.\n")
				elif  solsta == mosek.solsta.near_prim_infeas_cer:
					print("Primal or dual infeasibility.\n")
				elif mosek.solsta.unknown:
					print("Unknown solution status")
				else:
					print("Other solution status")


	return 'no_solution'
















