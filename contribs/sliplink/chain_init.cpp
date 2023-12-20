# include <math.h>
#include "slink.h"

//extern double b;
//extern long int *in_ran1;
//extern int n_ch;
//extern double ***XYZ;
//extern double box[3];
//extern int n_m;

void chain_init(Bead *bead, double b,long int *in_ran1, int n_ch, int n_m, double box[3] )
{
int i,k,nc,m,id;
double w1,w2,lw1,du1,u1;
double grnd1,grnd2,grnd3;
double sigma2;
double ***XYZinit;

sigma2=(double) b/sqrt(3.0);
XYZinit=(double ***)malloc(3*sizeof(double **));

for (k=0;k<3;k++) {
	XYZinit[k]=(double **)malloc(n_m*sizeof(double *));
	for (m=0;m<n_m;m++) {
	XYZinit[k][m]=(double *)malloc(n_ch*sizeof(double ));
		    }
		  }

    grnd1=ran1(in_ran1);
    grnd2=ran1(in_ran1);
    grnd3=ran1(in_ran1);

	for (nc=0;nc<n_ch;nc++) {
    		grnd1=ran1(in_ran1);
    		grnd2=ran1(in_ran1);
    		grnd3=ran1(in_ran1);
    		XYZinit[0][0][nc]=grnd1*box[0]; 
    		XYZinit[1][0][nc]=grnd2*box[1]; 
    		XYZinit[2][0][nc]=grnd3*box[2]; 

		for (m=1;m<n_m;m++) {	
    			w1=ran1(in_ran1);
    			w2=ran1(in_ran1);
    			lw1=log(w1);
    			du1=sqrt(-2.0*lw1);
    			u1=du1*cos(2.0* M_PI *w2);

			XYZinit[0][m][nc]=XYZinit[0][m-1][nc]+u1*sigma2;

    			w1=ran1(in_ran1);
    			w2=ran1(in_ran1);
    			lw1=log(w1);
    			du1=sqrt(-2.0*lw1);
    			u1=du1*cos(2.0*M_PI*w2);

			XYZinit[1][m][nc]=XYZinit[1][m-1][nc]+u1*sigma2;

    			w1=ran1(in_ran1);
    			w2=ran1(in_ran1);
    			lw1=log(w1);
    			du1=sqrt(-2.0*lw1);
    			u1=du1*cos(2.0*M_PI*w2);

			XYZinit[2][m][nc]=XYZinit[2][m-1][nc]+u1*sigma2;
		}
	}
		id=-1;
		for (nc=0;nc<n_ch;nc++) {
			for (m=0;m<n_m;m++) {
			id++;
			bead[id].index = id;
			bead[id].id_inf = id-1; 
                        if( m == 0 ) { bead[id].id_inf = -1; }
			bead[id].id_sup = id+1;
                        if( m == (n_m-1) ) { bead[id].id_sup = -1; }
			bead[id].id_ch=nc;
			bead[id].XYZ[0]=XYZinit[0][m][nc];
			bead[id].XYZ[1]=XYZinit[1][m][nc];
			bead[id].XYZ[2]=XYZinit[2][m][nc];
			if ((m==(n_m-1))) {
			bead[id].terminal = -1;
			//bead[id].id_sup = -1;
			} else {
			bead[id].terminal=0; 
			}
			bead[id].index_in_chain = m;
			}	
}
	for (i=0;i<n_m*n_ch;i++) {
//	fprintf(stdout,"Bead id=%d has chain index:%d\n",bead[i].index,bead[i].index_in_chain);
	}
	for (k=0;k<3;k++) { 
		for (m=0;m<n_m;m++) {
		free(XYZinit[k][m]);
		}	
		free(XYZinit[k]);
		}
}

