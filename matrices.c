// this file is for the c code for the functions to build matrices and vectors
// language: C
// 30 march 2017
// Hugo Martin

//////////////////////////////////// includes //////////////////////////////////////////////////////////////////
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


//////////////////////////////////// vector_D //////////////////////////////////////////////////////////////////
// this function builds the relative distances vector D in the normal directions
// WARNING ! The vector Dist MUST BE ALREADY DEFINED WITH NULL COMPONENTS EVERYWHERE
int vector_D(double *Dist, double const *q, unsigned const number_discs, unsigned const number_contacts, double const discs_radius, double const side_length, unsigned const number_walls){
	
        unsigned i, j, temp, N = number_discs;
	double L = side_length, R = discs_radius;  

        int k=-1;
	for(i=0;i<N;++i){
		for(j=i+1;j<N;++j){
                        temp = k+j-i;
			Dist[temp] = sqrt(((q[2*i]-q[2*j])*(q[2*i]-q[2*j]))+((q[2*i+1]-q[2*j+1])*(q[2*i+1]-q[2*j+1])))-2*R;
			if((Dist[temp]<=1.0e-15)&&(Dist[temp]>=-1.0e-15)){
				Dist[temp]= 0.;
			};
		};
		
		if(i<N){
			temp = k+N-i;
			Dist[temp] = sqrt((q[2*i+1]+2*R)*(q[2*i+1]+2*R))-R;//distance par rapport au sol
			if((Dist[temp]<=1.0e-15)&&(Dist[temp]>=-1.0e-15)){
				Dist[temp]= 0.;
			};
			if(number_walls>=3){
				Dist[temp+2] = sqrt((q[2*i]-(L-2*R))*(q[2*i]-(L-2*R)))-R;//distance par rapport au mur de droite
				if((Dist[temp+2]<=1.0e-15)&&(Dist[temp+2]>=-1.0e-15)){
					Dist[temp+2]= 0.;
				};
			};
			if(number_walls>=2){
				Dist[temp+1] = sqrt((q[2*i]+2*R)*(q[2*i]+2*R))-R;//distance par rapport au mur de gauche
				if((Dist[temp+1]<=1.0e-15)&&(Dist[temp+1]>=-1.0e-15)){
					Dist[temp+1]= 0.;
				};

			}
			k = temp+number_walls-1;			
		}; 
	};

	return 0;
};


//////////////////////////////////// vector_F //////////////////////////////////////////////////////////////////
// this function builds the relative distances vector D in the normal directions
// WARNING ! The vector Forces MUST BE ALREADY DEFINED WITH NULL COMPONENTS EVERYWHERE
int vector_F(double *Forces, unsigned const number_discs, double const gravity, double const slope){
	
	unsigned i;

	for(i=0;i<number_discs;++i){		
		Forces[2*i] = gravity*sin(slope);
		Forces[2*i+1] = -gravity*cos(slope);
	};
	
	return 0;
};


//////////////////////////////////// matrix_B //////////////////////////////////////////////////////////////////
// this function builds the matrix constraints B^t
// WARNING ! WE CALCULATE THE TRANSPOSE OF MATRICE B ACTUALLY
// WARNING ! THE MATRIX (actually it's a vector here) mat_B_transpose MUST BE ALREADY DEFINED WITH NULL COMPONENTS EVERYWHERE
int matrix_B_T_transp(double *mat_B_transpose, double *matrix_T, double const *q, unsigned const number_discs, unsigned const number_contacts, double const discs_radius, double const side_length, unsigned const number_walls){
	
	unsigned i, j, temp, N = number_discs, size = number_contacts, size2 = 3*N;
	double norm, L = side_length, R = discs_radius, E_ij1, E_ij2, qi, qip, qj, qjp, dR = 2.*R;  
	int k=-1;

	for(i=0;i<N;++i){
		for(j=i+1;j<N;++j){
                        temp = k+j-i;

			qi = q[2*i];
			qip = q[2*i+1];
			qj = q[2*j];
			qjp = q[2*j+1];
			
			norm = sqrt(((qi-qj)*(qi-qj))+((qip-qjp)*(qip-qjp)));

			mat_B_transpose[2*i*size+temp] = (qi-qj)/norm;
			mat_B_transpose[(2*i+1)*size+temp] = (qip-qjp)/norm;

                        E_ij1 = -mat_B_transpose[2*i*size+temp];
                        E_ij2 = -mat_B_transpose[(2*i+1)*size+temp];

			mat_B_transpose[2*j*size+temp] = E_ij1;
			mat_B_transpose[(2*j+1)*size+temp] = E_ij2;

                        matrix_T[temp*size2+2*i]  = E_ij2;
                        matrix_T[temp*size2+2*i+1] = -E_ij1;                               

                        matrix_T[temp*size2+2*j] = -E_ij2;
                        matrix_T[temp*size2+2*j+1] = E_ij1; 

                        matrix_T[temp*size2+2*N+i] = -R;
                        matrix_T[temp*size2+2*N+j] = -R;

		};
		
		if(i<N){
			temp = k+N-i;

			mat_B_transpose[(2*i+1)*size+temp] = (qip+dR)/sqrt((qip+dR)*(qip+dR)); // derivee par rapport au sol
                        
			E_ij1 = 0;
                        E_ij2 = -mat_B_transpose[(2*i+1)*size+temp];

                        matrix_T[temp*size2+2*i]  = E_ij2;
                        matrix_T[temp*size2+2*i+1] = -E_ij1;
                        matrix_T[temp*size2+2*N+i] = -R;
			
                        if(number_walls>=2){
                	        mat_B_transpose[2*i*size+temp+1] = (qi+dR)/sqrt((qi+dR)*(qi+dR)); // derivee par rapport au mur de gauche
                                
				E_ij1 = -mat_B_transpose[2*i*size+temp+1];
                                E_ij2 = 0;

                  	        matrix_T[(temp+1)*size2+2*i]  = E_ij2;
                        	matrix_T[(temp+1)*size2+2*i+1] = -E_ij1;
                      		matrix_T[(temp+1)*size2+2*N+i] = -R;

                        };
                        if(number_walls>=3){
                	        mat_B_transpose[2*i*size+temp+2] = (qi-(L-dR))/sqrt((qi-(L-dR))*(qi-(L-dR))); // derivee par rapport au mur de droite
                        	
				E_ij1 = -mat_B_transpose[2*i*size+temp+2];
                                E_ij2 = 0;

                  	        matrix_T[(temp+2)*size2+2*i]  = E_ij2;
                        	matrix_T[(temp+2)*size2+2*i+1] = -E_ij1;
                      		matrix_T[(temp+2)*size2+2*N+i] = -R;

			};

                        k = temp+number_walls-1;			
		}; 
	};

	return 0;
};



//////////////////////////////////// build_asub_and_aval ///////////////////////////////////////////////////////
// all in the title, conversion for one of mosek arguments for linear operator
// it's writing aval and asub arguments for mosek in text files
int build_asub_and_aval(const double * LIN, unsigned const number_lines, unsigned const number_colonnes){


	unsigned i, j, lines = number_colonnes, columns = number_lines; // we exchange lines and columns because we will work on the transpose of LIN
	char* asub_name = "./temporary_files/asub_file.dat"; 
	char* aval_name = "./temporary_files/aval_file.dat";
	double temp;

	FILE* asub_file = fopen(asub_name,"w"); // to stock asub array
	FILE* aval_file = fopen(aval_name,"w"); // to stock aval array
	if((asub_file==NULL)||(aval_file==NULL)){
		printf("\t PROBLEM IN THE ASUB OR AVAL CONSTRUCTION !\n");
		return 1; 
	};


	for(i=0;i<lines;++i){ // I am looking for all the lines of the transpose of LIN operator
		for(j=0;j<columns;j++){
			temp = LIN[i*columns+j];
			//if((temp>=1.0e-15)||(temp<=-1.0e-15)){
			if(temp!=0){			
				fprintf(asub_file, "%d ", j); // write the column number in asub_file
				fprintf(aval_file, "%.16e ", temp); // write the value in aval_file
			};
		};
		fprintf(asub_file,"\n"); 
		fprintf(aval_file,"\n"); // we end the line so change line in both files
	};
	
	
	fclose(asub_file);
	fclose(aval_file);
	
	return 0;
};



//////////////////////////////////// build_qsubi_qsubj_and_qval ////////////////////////////////////////////////
// all in the title, conversion for one of mosek arguments for quadratic operator
// it's writing qsubi, qsubj and qval arguments for mosek in text files
int build_qsubi_qsubj_and_qval(const double * Q, unsigned const sizeQ){


	unsigned i, j, size = sizeQ;
	char* qsubi_name = "./temporary_files/qsubi_file.dat"; 
	char* qsubj_name = "./temporary_files/qsubj_file.dat";
	char* qval_name = "./temporary_files/qval_file.dat";
	double temp;


	FILE* qsubi_file = fopen(qsubi_name,"w"); // to stock qsubi array
	FILE* qsubj_file = fopen(qsubj_name,"w"); // to stock qsubj array
	FILE* qval_file = fopen(qval_name,"w"); // to stock qval array


	if((qsubi_file==NULL)||(qval_file==NULL)||(qsubj_file==NULL)){
		printf("\t PROBLEM IN THE QSUBI, QVAL OR QSUBJ CONSTRUCTION !\n");
		return 1; 
	};


	for(i=0;i<size;++i){ // I run on matrix Q (which is not a matrix)
		for(j=i;j<size;j++){
			temp = Q[i*size+j];
			//if((temp>=1.0e-15)||(temp<=-1.0e-15)){
			if(temp!=0.){			
				fprintf(qsubi_file, "%d ", j); // write the column number in asub_file
				fprintf(qsubj_file, "%d ", i); // write the column number in asub_file
				fprintf(qval_file, "%.16e ", temp); // write the value in aval_file
			};
		};
	};
	
	fclose(qsubi_file);
	fclose(qval_file);
	fclose(qsubj_file);
	
	return 0;
};





















































