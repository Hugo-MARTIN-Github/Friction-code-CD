#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

################################################################################################################
############################### Imports ########################################################################

import numpy as np

################################################################################################################
# language: python 3.4
# 15 May 2017
# Hugo Martin

############################### make_data_files ###############################################################
# All is in the function title

def make_data_files(file_number, position_vector):

	np.savetxt("./temporary_files/python_data_" + str(file_number) + ".dat", position_vector, fmt='%.15e')
