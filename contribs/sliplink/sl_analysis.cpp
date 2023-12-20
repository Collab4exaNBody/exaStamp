# include <math.h>
# include <malloc.h>
# include "slink.h"

void sl_analysis(Bead *bead, double ***XYZ0, double *MSD, double **COM0, double *P2, double *Rbmean, double ***Rb0, double *endtoend, double *Ree2, double *DRee, double **Ree0, int n_ch, int n_m, int istep) {

int l,m,j,jj,dN;
int i,k,nc;
int tmp;
double tmpx,tmpy,tmpz;
int *sort_index;
double dree_tmp,msd_tmp,drb_tmp,com_tmp;
double ***Rb,**Ree, **COM; 
double ***XYZan,**XYZsort;
double *diff;
FILE *reefile,*beadcoor;

//beadcoor=fopen("beadcoor.dat","a+");
//reefile=fopen("ree.dat","a");
//fprintf(reefile,"istep=%d\n",istep);
diff=(double *)malloc(3*sizeof(double ));
COM=(double **)malloc(3*sizeof(double *));
//COM0=(double **)malloc(3*sizeof(double *));
        for (k=0;k<3;k++) {
                COM[k]=(double *)malloc(n_ch*sizeof(double ));
 //               COM0[k]=(double *)malloc(n_ch*sizeof(double ));
                }

sort_index=(int *)malloc((n_m*n_ch)*sizeof(int ));
XYZsort=(double **)malloc((n_m*n_ch)*sizeof(double *));
for (i=0;i<n_m*n_ch;i++) XYZsort[i]=(double *)malloc(3*sizeof(double ));

XYZan=(double ***)malloc(3*sizeof(double **));
for (k=0;k<3;k++) {
	XYZan[k]=(double **)malloc(n_m*sizeof(double *));
	for (j=0;j<n_m;j++) {
		XYZan[k][j]=(double *)malloc(n_ch*sizeof(double ));	
	}
}
Rb=(double ***)malloc(3*sizeof(double **));
//Rb0=(double ***)malloc(3*sizeof(double **));
        for (k=0;k<3;k++) {
                Rb[k]=(double **)malloc((n_m-1)*sizeof(double *));
 //               Rb0[k]=(double **)malloc((n_m-1)*sizeof(double *));
                for (j=0;j<(n_m-1);j++) {
                        Rb[k][j]=(double *)malloc((n_ch)*sizeof(double ));
  //                      Rb0[k][j]=(double *)malloc((n_ch)*sizeof(double ));
                }
        }
//Ree0=(double **)malloc(3*sizeof(double *));
Ree=(double **)malloc(3*sizeof(double *));
        for (i=0;i<3;i++) {
 //       Ree0[i]=(double *)malloc(n_ch*sizeof(double ));
        Ree[i]=(double *)malloc(n_ch*sizeof(double ));
    	}

for (k=0;k<3;k++) {
	for (nc=0;nc<n_ch;nc++) {
	Ree[k][nc]=0.0;
}
}
/*DRee[nc]=0.0;
P2[nc]=0.0;
MSD[nc]=0.0;
for (k=0;k<3;k++) {
COM[k][nc]=0.0;
}  */
/* Calculate end-to-end vector */
/* Center of Mass */
i=-1;
	//fprintf(stdout,"Starting reading....\n");
	for (i=0;i<n_m*n_ch;i++) {
	XYZsort[i][0]=bead[i].XYZ[0];
	XYZsort[i][1]=bead[i].XYZ[1];
	XYZsort[i][2]=bead[i].XYZ[2];
	//fprintf(stdout,"%lf %lf %lf %lf %lf %lf\n",XYZsort[i][0],XYZsort[i][1],XYZsort[i][2], bead[i].XYZ[0],bead[i].XYZ[1],bead[i].XYZ[2]); 
	sort_index[i]=bead[i].index;
				} 
	//fprintf(stdout,"Reading ok.\n");
//	fprintf(stdout,"Starting reading....\n");
	for (i=0;i<n_m*n_ch;i++) {
		for (j=0;j<n_m*n_ch;j++) {
		if (sort_index[j] > sort_index[i]) {
		tmp = sort_index[i];
		tmpx = XYZsort[i][0];
		tmpy = XYZsort[i][1];
		tmpz = XYZsort[i][2];
		XYZsort[i][0] = XYZsort[j][0];
		XYZsort[i][1] = XYZsort[j][1];
		XYZsort[i][2] = XYZsort[j][2];
		sort_index[i]=sort_index[j];
		XYZsort[j][0] = tmpx;
		XYZsort[j][1] = tmpy;
		XYZsort[j][2] = tmpz;
		sort_index[j]=tmp;
					   }
			}
		}	
	/*for (i=0;i<n_m*n_ch;i++) 
	fprintf(stdout,"%d %lf %lf %lf\n",sort_index[i],XYZsort[i][0],XYZsort[i][1],XYZsort[i][2]); */ 

//	fprintf(stdout,"Sorting ok.\n");
i=-1;
 	for (nc=0;nc<n_ch;nc++) {
 		for (l=0;l<n_m;l++) {
			i++;
			XYZan[0][l][nc]=XYZsort[i][0];
			XYZan[1][l][nc]=XYZsort[i][1];
			XYZan[2][l][nc]=XYZsort[i][2];
		//	if (istep%10==0) 
		//	fprintf(beadcoor,"%lf  %lf  %lf\n",XYZan[0][l][nc],XYZan[1][l][nc],XYZan[2][l][nc]);
			}
		}	

for (nc=0;nc<n_ch;nc++) {
//nc=bead[id].id_ch;
DRee[nc]=0.0;
P2[nc]=0.0;
MSD[nc]=0.0;
}
for (k=0;k<3;k++) {
for (nc=0;nc<n_ch;nc++) {
COM[k][nc]=0.0;
}
}
for (nc=0;nc<n_ch;nc++) {
for (j=0;j<n_m;j++) {
COM[0][nc]=COM[0][nc]+XYZan[0][j][nc];
COM[1][nc]=COM[1][nc]+XYZan[1][j][nc];
COM[2][nc]=COM[2][nc]+XYZan[2][j][nc];
}
}
for (nc=0;nc<n_ch;nc++) {
COM[0][nc]=COM[0][nc]/(double) (n_m);
COM[1][nc]=COM[1][nc]/(double) (n_m);
COM[2][nc]=COM[2][nc]/(double) (n_m);
//fprintf(stdout,"COMx=%lf,COMy=%lf,COMz=%lf\n",COM[0][nc],COM[1][nc],COM[2][nc]);
}
/* End to end vector */

for (nc=0;nc<n_ch;nc++) {
//for (k=0;k<3;k++) 
        Ree[0][nc]=XYZan[0][n_m-1][nc]-XYZan[0][0][nc];
                Ree[1][nc]=XYZan[1][n_m-1][nc]-XYZan[1][0][nc];
                        Ree[2][nc]=XYZan[2][n_m-1][nc]-XYZan[2][0][nc];
			// fprintf(reefile,"Ree=%lf %lf %lf\n",Ree[0][nc],Ree[1][nc],Ree[2][nc]);
	}

	for (nc=0;nc<n_ch;nc++) {
                      Ree2[nc]=Ree[0][nc]*Ree[0][nc]+Ree[1][nc]*Ree[1][nc]+Ree[2][nc]*Ree[2][nc];
			
			//if (istep%10==0)
			// fprintf(reefile,"Ree2=%lf\n",Ree2[nc]);	
		}
/*	for (nc=0;nc<n_ch;nc++)
		fprintf(reefile,"Ree2=%lf\n",Ree2[nc]); */
 /* Store initial values for Ree, XYZ and COM */
if (istep==0) {
for (nc=0;nc<n_ch;nc++) {
    for (l=0;l<n_m;l++) {
    XYZ0[0][l][nc]=XYZan[0][l][nc];
    XYZ0[1][l][nc]=XYZan[1][l][nc];
    XYZ0[2][l][nc]=XYZan[2][l][nc];
}
}
}
if (istep==0) {
for (nc=0;nc<n_ch;nc++) {
    Ree0[0][nc]=Ree[0][nc];
    Ree0[1][nc]=Ree[1][nc];
    Ree0[2][nc]=Ree[2][nc];
    COM0[0][nc]=COM[0][nc];
    COM0[1][nc]=COM[1][nc];
    COM0[2][nc]=COM[2][nc];
}
}
/* End-to-end vector */
for (nc=0;nc<n_ch;nc++) { 
   endtoend[nc]=(XYZan[0][n_m-1][nc]-XYZan[0][0][nc])*(XYZan[0][n_m-1][nc]-XYZan[0][0][nc])+(XYZan[1][n_m-1][nc]-XYZan[1][0][nc])*(XYZan[1][n_m-1][nc]-XYZan[1][0][nc])+(XYZan[2][n_m-1][nc]-XYZan[2][0][nc])*(XYZan[2][n_m-1][nc]-XYZan[2][0][nc]);
//if (istep%10==0)
//fprintf(reefile,"Endtoend of chain %d is:%lf\n",nc,endtoend[nc]);
}
/* End to end vector autocorrelation */
	
for (nc=0;nc<n_ch;nc++) { 
        dree_tmp=(Ree0[0][nc]*Ree[0][nc]+Ree0[1][nc]*Ree[1][nc]+Ree0[2][nc]*Ree[2][nc])/(Ree0[0][nc]*Ree0[0][nc]+Ree0[1][nc]*Ree0[1][nc]+Ree0[2][nc]*Ree0[2][nc]);
        DRee[nc]=dree_tmp;
}
/* Calculate the bond vector auto-correlation function */

for (nc=0;nc<n_ch;nc++) { 
Rbmean[nc]=0.0;
                for (m=0;m<n_m-1;m++) {
                        Rb[0][m][nc]=XYZan[0][m+1][nc]-XYZan[0][m][nc];
                        Rb[1][m][nc]=XYZan[1][m+1][nc]-XYZan[1][m][nc];
                        Rb[2][m][nc]=XYZan[2][m+1][nc]-XYZan[2][m][nc];
        if (istep==0) {
                Rb0[0][m][nc]=Rb[0][m][nc];
                Rb0[1][m][nc]=Rb[1][m][nc];
                Rb0[2][m][nc]=Rb[2][m][nc];
                }
                drb_tmp=(Rb0[0][m][nc]*Rb[0][m][nc]+Rb0[1][m][nc]*Rb[1][m][nc]+Rb0[2][m][nc]*Rb[2][m][nc])/(Rb0[0][m][nc]*Rb0[0][m][nc]+Rb0[1][m][nc]*Rb0[1][m][nc]+Rb0[2][m][nc]*Rb0[2][m][nc]);
                P2[nc]=P2[nc]+drb_tmp;
                Rbmean[nc]=Rbmean[nc]+(Rb[0][m][nc]*Rb[0][m][nc]+Rb[1][m][nc]*Rb[1][m][nc]+Rb[2][m][nc]*Rb[2][m][nc]);
                }
}

/* Calculate chain endtoend distance 
 * */
/* Calculate MSD for COM */
for (nc=0;nc<n_ch;nc++) {
        diff[0]=COM[0][nc]-COM0[0][nc];
        diff[1]=COM[1][nc]-COM0[1][nc];
        diff[2]=COM[2][nc]-COM0[2][nc];
          msd_tmp=diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
        MSD[nc]=MSD[nc]+msd_tmp;

}
for (i=0;i<n_m*n_ch;i++) { 
	free(XYZsort[i]);
	}
	free(XYZsort);
for (k=0;k<3;k++) {
	for (m=0;m<n_m;m++) {
	free(XYZan[k][m]);
	}
	free(XYZan[k]);
	}
	free(XYZan);
for (k=0;k<3;k++) 
free(Rb[k]);
free(Rb);    
free(diff);
for (k=0;k<3;k++) {
free(Ree[k]);
free(COM[k]);
}
free(Ree);
free(COM);
free(sort_index);
//fclose(beadcoor);
}                                 
