#! /home/hugo/anaconda3/envs/py34/bin/ipython
# -*- coding: cp1252 -*-

####################### This file contains function to build objects for mosek #################################
# language: python 3.4
# 30 march 2017
# Hugo Martin

############################### imports ########################################################################
import numpy as np
import ctypes as ct



############################### conv_mosek #####################################################################
# this function build objects for mosek from optimization problem operators 
def conv_mosek(Q, LIN, sizeQ, ligneLIN, colonne_LIN):


######################## definition of c library and types
	cclib = ct.cdll.LoadLibrary('./matrices.so') # definition of c library

	vec_type = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
	cclib.build_asub_and_aval.argtypes = (vec_type, ct.c_uint, ct.c_uint)
	cclib.build_qsubi_qsubj_and_qval.argtypes = (vec_type, ct.c_uint)
	

######################## writing of aval and asub in files
	sizze = colonne_LIN*ligneLIN	
	LIN = np.reshape(LIN.T,[1,sizze])
	LINN = np.zeros(sizze)
	for i in range(sizze):
		LINN[i] = LIN[0,i] #convert 2darray LIN in tuple 1d

	cclib.build_asub_and_aval(LINN, ligneLIN, colonne_LIN) # I write aval and asub in files


######################## writing qsubi qsubj and qval in files
	size = sizeQ*sizeQ
	Q = np.reshape(Q,[1,size])
	QQ = np.zeros(size)
	for i in range(size):
		QQ[i] = Q[0,i] #convert 2darray LIN in tuple 1d
	cclib.build_qsubi_qsubj_and_qval(QQ, sizeQ)



######################## reading of asub and aval files
	asub = []
	aval = []
	asub_file = open("./temporary_files/asub_file.dat", 'r')
	aval_file = open("./temporary_files/aval_file.dat",'r')
	for line in asub_file:
		temp = []
		for numberString in line.split():
			temp.append(int(numberString))
		
		asub.append(temp)

	for line in aval_file:
		temp = []
		for numberString in line.split():
			temp.append(float(numberString.replace(',','.')))
		
		aval.append(temp)


######################## reading of qsubi, qsubj and qval files
	qsubi_file = open("./temporary_files/qsubi_file.dat", 'r')
	qsubj_file = open("./temporary_files/qsubj_file.dat","r")	
	qval_file = open("./temporary_files/qval_file.dat",'r')

	qsubi =	np.fromfile(qsubi_file, dtype=int, sep=' ').tolist()
	qsubj = np.fromfile(qsubj_file, dtype=int, sep=' ').tolist()

	qval = []
	for line in qval_file:
		for numberString in line.split():
			qval.append(float(numberString.replace(',','.')))


######################### closes
	asub_file.close()
	aval_file.close()
	qsubi_file.close()
	qsubj_file.close()
	qval_file.close()

	return asub, aval, qsubi, qsubj, qval