# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include <time.h>
# include <assert.h>
# include <sys/time.h>

# include <map>
# include <random>
# include <iostream>

# include "slink.h"

std::map< int , int > g_bead_index;

#define anint(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))
#define aint(x) ((x)>0 ? floor(x) : ceil((x)))
#define nint(x) ((x)>0 ? (int)((x)+0.5) : (int)((x)-0.5)

# define RMAX  32767


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (1.2E-07)
#define RNMX (1.0-EPS)

double ran1(long *idum)
    {
        int j;
        long k;
        double temp;
        static long iy=0;
        static long iv[NTAB];
        if ( (*idum<=0) || (iy==0) ) {
            if (-(*idum)<1) {
              *idum=1;
            } else {
              *idum=-(*idum);
            }
            for (j=NTAB+7;j>=0;j--) {
              k=*idum/IQ;
              *idum=IA*(*idum-k*IQ)-IR*k;
              if (*idum<0) *idum+=IM;
              if (j<NTAB) iv[j]=*idum;
            }
            iy=iv[0];
        }
        k=*idum/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum<0) *idum+=IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j]=*idum;
        if ( (temp=AM*iy ) > RNMX ) {
            return RNMX;
        } else {
            return temp;
        }
    }

#undef IA
#undef IM 
#undef AM 
#undef IQ 
#undef IR 
#undef NTAB 
#undef NDIV 
#undef EPS 
#undef RMNX 

double my_normal_distribution(long *in_ran1)
{
  double w1=ran1(in_ran1);
  double w2=ran1(in_ran1);
  double lw1=log(w1);
  double du1 = sqrt(-2.0*lw1);
  double u1 = du1*cos(2.0*M_PI*w2);
  return u1;
}

void test_normal_distribution()
{
  long in_ran1 = -3333;
  
  int N = 100000;
  double mu = 0.0;
  double sigma = 0.0;
  for(int i=0;i<N;i++)
  {
    double X = my_normal_distribution( &in_ran1 );
    mu += X;
    sigma += fabs(X);
  }
  mu /= N;
  sigma /= N;
  std::cout <<"internal method : mu="<<mu<<", sigma="<<sigma<<std::endl;
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0.0,1.0};

  mu = 0.0;
  sigma = 0.0;
  for(int i=0;i<N;i++)
  {
    double X = d(gen);
    mu += X;
    sigma += fabs(X);
  }
  mu /= N;
  sigma /= N;
  std::cout <<"c++ std normal distribution : mu="<<mu<<", sigma="<<sigma<<std::endl;
}

double random_number(void) {
  // generateur uniforme librairie C
     double temp ;
     temp = (double) rand() / RAND_MAX ;
     return temp;
    }
  
double gaussienne (void) {

  double t, r1, r2, res1, res2 ;

  t  = random_number();
  r1 = random_number();

  if ( r1 == 0 ) r1 = 1e-9 ;

  r2 =  random_number();

  if (t < .5){
    res1 = sqrt(-2.*log(r1))*cos(2.*M_PI*r2);
    return (res1);
  } else {
    res2 = sqrt(-2.*log(r1))*sin(2.*M_PI*r2);
    return (res2);
  }
}


int main()
{
test_normal_distribution();

Bead *bead;
SlipLink *sl;
long int *in_ran1;
int *ran_list;
int i,zzz,n_ch,nc;
int n_e,n_m;
int *destroy;
int maxstep,istep;
int output;
double tmax;
double vol;
double dt;
double ksi,ksis;
double b,kBT;
double Ns,Rgir2;
double box[3];
double pi;
double h,h1,tauzero;
double dt1;
double kf,ks;
double rho;
double *endtoend;
double *Rbmean,*P2,*MSD;
double *Ree2,*DRee;
double ***XYZ0,**COM0;
double **Ree0,***Rb0;
double valtot,msd_tot,DRee_tot;
double Ree2_tot,Rbmean_tot,P2_tot;
int *num_SL;
//FILE *slcoorfile,*coorfile,*oldcoor,*fslfile,*rousefile,*nslfile;
int npairs,id;
int j,k,l;
int m,n;
int Nperm,randi,randj;
char fName[34];
int zf;
FILE *initfile,*inputfile,*end2endfile,*reefile,*p2file,*msdfile,*numslfile,*beadidfile;

/*Declaration of function pointers */
in_ran1=(long *)malloc(sizeof(long));
/* Initialize seed*/ 
*in_ran1=-3333;

/* Define model parameters for the creation of slip-links. */
maxstep=1000000;
kBT=1.0;
b=1.0;
Ns=0.5;
n_e=4;
ksi=1.0;
rho=1.0;
zf=12;
h1=0.01;
output=1000;
n_ch=100;
n_m=32;
zzz=n_ch*n_m/n_e;
//zzz=2;
printf("Status ok\n");
/* sprintf(fName,"init_conf_%d.dat",zf);
	initfile=fopen(fName,"r");
	if (initfile==NULL)
	{
		printf("Can't open %s\n",fName);
		exit(1);
	}
	fprintf(stdout,"Reading file %s\n",fName); */ 
/* Initialization */
num_SL=(int *)malloc(n_ch*sizeof(int ));
destroy=(int *)malloc(zzz*sizeof(int ));
//age_SL=(int *)malloc(zzz*sizeof(int ));
endtoend=(double *)malloc(n_ch*sizeof(double ));
Ree2=(double *)malloc(n_ch*sizeof(double ));
DRee=(double *)malloc(n_ch*sizeof(double ));
P2=(double *)malloc(n_ch*sizeof(double ));
Rbmean=(double *)malloc(n_ch*sizeof(double ));
MSD=(double *)malloc(n_ch*sizeof(double ));

tauzero=(ksi*b*b)/(3.0*M_PI*M_PI*kBT);
h=h1*tauzero;
tmax=maxstep*h1;
npairs=zzz/2;
istep=0;

//fscanf(initfile,"%lf %lf %lf\n",&box[0],&box[1],&box[2]);
/*Open output files used to print data */
//initfile=fopen("init_conf.dat","w");
end2endfile=fopen("end2end.dat","w");
reefile=fopen("ree.dat","w");
p2file=fopen("p2.dat","w");
msdfile=fopen("msdcom.dat","w");
beadidfile=fopen("beadid.dat","w");  
  //*Allocate memory */
//COM=(double **)malloc(3*sizeof(double *));
COM0=(double **)malloc(3*sizeof(double *));
	for (k=0;k<3;k++) {
//		COM[k]=(double *)malloc(n_ch*sizeof(double ));
		COM0[k]=(double *)malloc(n_ch*sizeof(double ));
		}
Rb0=(double ***)malloc(3*sizeof(double **));
//Rb=(double ***)malloc(3*sizeof(double **));
	for (k=0;k<3;k++) {
		Rb0[k]=(double **)malloc((n_m-1)*sizeof(double *));
//		Rb[k]=(double **)malloc((n_m-1)*sizeof(double *));
		for (j=0;j<(n_m-1);j++) {
			Rb0[k][j]=(double *)malloc((n_ch)*sizeof(double ));
//			Rb[k][j]=(double *)malloc((n_ch)*sizeof(double ));
		}
	}		
    Ree0=(double **)malloc(3*sizeof(double *));
    //Ree=(double **)malloc(3*sizeof(double *));
    for (i=0;i<3;i++) {
//	Ree[i]=(double *)malloc(n_ch*sizeof(double ));	
	Ree0[i]=(double *)malloc(n_ch*sizeof(double ));	
    }
    XYZ0=(double ***)malloc(3*sizeof(double **));
    for(k=0;k<3;k++) {
    	XYZ0[k]=(double **)malloc(n_m*sizeof(double *));
    	for (j=0;j<n_m;j++) {
		XYZ0[k][j]=(double *)malloc(n_ch*sizeof(double ));
        }
    }
bead=(struct Bead*)malloc((n_ch*n_m)*sizeof(struct Bead));
sl=(struct SlipLink*)malloc(zzz*sizeof(struct SlipLink));

/* Create the initial configuration of the polymer melt */
vol=(double) (n_ch*n_m)*pow(b,3)/rho;
box[0]=exp(1.0/3.0*log(vol));
box[1]=box[0];
box[2]=box[0];
initfile=fopen("init_conf.dat","w"); 
//
//inputfile=fopen("input.dat","w");
	chain_init(bead, b,in_ran1,n_ch,n_m, box);
	id=-1;
  	for (m=0;m<n_ch;m++) {
		for (l=0;l<n_m;l++) {
		id++;
		fprintf(initfile,"%lf %lf %lf\n",bead[id].XYZ[0],bead[id].XYZ[1],bead[id].XYZ[2]); 
		}
	}  
	for (id=0;id<(n_ch*n_m);id++) {
		fscanf(initfile,"%lf %lf %lf\n",&bead[id].XYZ[0],&bead[id].XYZ[1],&bead[id].XYZ[2]);
	//	fprintf(inputfile,"%lf %lf %lf\n",bead[id].XYZ[0],bead[id].XYZ[1],bead[id].XYZ[2]);
	}  
fclose(initfile);


//fclose(inputfile);
/* Initialize bead map */
// for (i=0;i<n_m*n_ch;i++) g_bead_index[i]=i;

/* Permutate the elements of the map */

Nperm=1200;
for (i=0;i<Nperm;i++) {
randi=rand()%(n_m*n_ch);
randj=rand()%(n_m*n_ch);
std:: swap( bead[randi] , bead[randj] );
}


for (i=0;i<n_m*n_ch;i++) { g_bead_index[i] = -1;  }

for (i=0;i<n_m*n_ch;i++)
{
	assert( bead[i].index>=0 && bead[i].index < (n_m*n_ch) );
	g_bead_index[ bead[i].index] = i;
	fprintf(beadidfile,"i=%d, bead[i].index=%d, g_bead_index[bead[i].index]=%d \n",i,bead[i].index, g_bead_index[bead[i].index]);
}
fclose(beadidfile);

for (i=0;i<n_m*n_ch;i++) { assert( g_bead_index[i] != -1);  }

sl_init(sl,bead, b,Ns,zzz,n_e,n_ch, n_m, in_ran1, num_SL);


//fclose(inputfile);
/* Calculate the radius of gyration */
Rgir2=(n_m)*(b*b)/6.0;

for (k=0;k<3;k++) {
	for (m=0;m<(n_m-1);m++) {
		for (n=0;n<n_ch;n++) {
		Rb0[k][m][n]=0.0;	
		}
	}
}
/* Loop over time */
for (istep=0;istep<maxstep;istep++) {
/* Initialize and reset accumulators */
for (k=0;k<3;k++) {
	for (m=0;m<(n_m-1);m++) {
		for (n=0;n<n_ch;n++) {
		}
	}
}
for(m=0;m<n_ch;m++) Ree2[m]=0.0;
valtot=0.0;
msd_tot=0.0;
DRee_tot=0.0;
Ree2_tot=0.0;
Rbmean_tot=0.0;
P2_tot=0.0;
for (nc=0;nc<n_ch;nc++) endtoend[nc]=0.0;
	for (m=0;m<n_m*n_ch;m++) {
		for (k=0;k<3;k++) {
		bead[m].F_SL[k]=0.0;
		bead[m].F_SP[k]=0.0;
		}
	}
for (j=0;j<zzz;j++) destroy[j]=0;

/* Loop over polymer monomers */
/* Store position from previous timestep */
        for (id=0;id<n_m*n_ch;id++) {
            	bead[id].XYZold[0]=bead[id].XYZ[0]; 
            	bead[id].XYZold[1]=bead[id].XYZ[1]; 
                bead[id].XYZold[2]=bead[id].XYZ[2];
            }
/* Calculation of Rouse forces */
		force_rouse(bead,b,kBT,n_ch, n_m, istep);
/* Calculation of slip-link forces */
		force_sliplink(bead,sl,b,kBT,n_ch,n_m,Ns,zzz, istep);
/* Perform static and dynamic analysis for the polymer chains */
		sl_analysis(bead, XYZ0,MSD,COM0,P2, Rbmean, Rb0, endtoend, Ree2, DRee, Ree0, n_ch, n_m, istep);

 /* End to end vector */

/*	if ((istep%output)==0) {
	endtoend[nc]=(XYZ[0][n_m-1][nc]-XYZ[0][0][nc])*(XYZ[0][n_m-1][nc]-XYZ[0][0][nc])+(XYZ[1][n_m-1][nc]-XYZ[1][0][nc])*(XYZ[1][n_m-1][nc]-XYZ[1][0][nc])+(XYZ[2][n_m-1][nc]-XYZ[2][0][nc])*(XYZ[2][n_m-1][nc]-XYZ[2][0][nc]);
} */
/* Beads position evolution */
        move_beads(bead, sl, h,ksi,kBT,in_ran1,n_ch,n_m, istep);

/* Slip-link position evolution */
        move_sl(bead, sl,b,h,kBT,ksi,Ns,n_ch,n_m,in_ran1,zzz,destroy);

/* end of loop over polymer beads */
/* Slip-link renewal algorithm:constraint release */
        sl_destroy(destroy, zzz);
        sl_renew(b,bead, sl, Ns,zzz,n_ch, istep,n_m, in_ran1,destroy);
/* Update of the sliplink list and sliplink pointers */
        update_sl(bead, sl, zzz, n_m, n_ch, num_SL, istep);
/* Write data to output files */

for (j=0;j<n_ch;j++) {
valtot+=endtoend[j];
DRee_tot+=DRee[j];
Ree2_tot+=Ree2[j];
Rbmean_tot+=Rbmean[j];
P2_tot+=P2[j];
msd_tot+=MSD[j];
}
	if (istep%output==0) {
	fprintf(end2endfile,"%d %lf\n",istep,valtot/(double) (n_ch));
	fprintf(reefile,"%lf %lf %lf\n",istep*h1,DRee_tot/(double) n_ch,Ree2_tot/(double) (n_ch));
	fprintf(p2file,"%lf %lf %lf\n",istep*h1,P2_tot/((double) (n_ch))/((double) (n_m-1)),Rbmean_tot/((double) (n_ch))/((double) (n_m-1)));
	fprintf(msdfile,"%lf %lf\n",istep*h1,msd_tot/((double) (n_ch)));
	}
}  
fclose(reefile);
fclose(p2file);
fclose(end2endfile);
fclose(msdfile);
/* Print the number of SLs as well as the chain and the monomer they belong to at the end of the run */
/* Free memory */
	free(Ree0);
for (k=0;k<3;k++) {
	for (m=0;m<(n_m-1);m++) {
	//	free(Rb[k][m]);
		free(Rb0[k][m]);
	}
//	free(Rb[k]);
	free(Rb0[k]);
	}
//	free(Rb);
	free(Rb0); 
//printf("Status 3 ok\n");
	free(num_SL);
	free(in_ran1);
	free(endtoend);
	free(destroy);
	free(DRee);
	free(Ree2);
	free(P2);
	free(Rbmean);
	free(MSD);
//	free(diff);
return 0;
} 


