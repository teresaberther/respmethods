#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES // for C
#include <math.h>

typedef struct Array Array;

struct Array {
	int len;
	double *data;
};

Array read_data(void){

	FILE *fp = fopen("testarray.csv", "r");
    /*if (!fp) {
        perror("open");
        return NULL;
    }*/	

	// first pass, count the lines	
	int n = 0;
	double tmp;
	while (fscanf(fp, "%lf", &tmp)==1){
		n++;
	}

	// allocate memory
	double *data = (double*)malloc(n* sizeof(double));
	/*if (data==NULL){
		return NULL;
	}*/

	// second pass
	rewind(fp);
	for (int j=0; j<n; j++){
		fscanf(fp, "%lf", &data[j]);		
	}	

	// create structure
	Array my_array;
	my_array.data = data;
	my_array.len = n;

	return my_array;	
}

int main(void) {

 	// double *arr = read_data();
	Array arr = read_data();

	/*
	printf("%d\n", arr.len);
	printf("%lf\n", arr.data[6]);
	*/

	// find peaks & troughs
	int *null_der = (int*)malloc(arr.len* sizeof(int));
	int acc_ = 0;

	// first pass: assign 1 and -1 to peaks and troughs, respectively
	for (int j=1; j<arr.len-1; j++){

		int peakhere = ((arr.data[j]>arr.data[j+1]) & (arr.data[j]>arr.data[j-1]));
		int troughhere = -1*((arr.data[j]<arr.data[j+1]) & (arr.data[j]<arr.data[j-1]));
		null_der[j] = peakhere + troughhere;

		acc_ = acc_ + peakhere - troughhere;

	}

	// second pass: store the indexes of events (peaks and troughs)
	acc_ = 0;
	int *idx_der = (int*)malloc(arr.len* sizeof(int));

	for (int j=0; j<arr.len; j++){

		if (null_der[j]){

			idx_der[acc_] = j;
			acc_ ++;
				
		}
		
	}

	// last pass: interpolate phase
	double *phase_vect = (double*)malloc(arr.len*sizeof(double));

	double this_entry = NAN;
	int acc_event = 0;
	int npoints = 0;
	double acc_phase = NAN;

	phase_vect[0] = NAN;
	for (int j=1; j<arr.len; j++){

		if (null_der[j]){

			this_entry = M_PI*(1-null_der[j])/(2);
			
			if (acc_event<acc_){

				npoints = idx_der[acc_event+1]-idx_der[acc_event];
				acc_phase = (M_PI/(double)npoints); // *(-null_der[j]); // still needs to get the right p ratio here, now testing purpose only				
				printf("%d\n", null_der[j]);
				printf("%lf\n", acc_phase);
				
				acc_event ++;
				
			}

			else {

				acc_phase = NAN;
			}
			
		}

		if (null_der[j-1]==-1){

			this_entry = -M_PI;
			
		}

		else {

			this_entry = this_entry+acc_phase;			

		}

		phase_vect[j] = this_entry;		

	}


	// write (just for debug purpose)
	FILE *out = fopen("out.dat", "w");
	for (int j=0; j<arr.len; j++){
		fprintf(out, "%g\n", phase_vect[j]);
	}
	fclose(out);
	 		
    return 0;
    
}


/*



*/
