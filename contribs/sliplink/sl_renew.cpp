# include <math.h>
# include "slink.h"
# include <map>

#define anint(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))
#define aint(x) ((x)>0 ? floor(x) : ceil((x)))
#define nint(x) ((x)>0 ? (int)((x)+0.5) : (int)((x)-0.5)

void sl_renew(double b,Bead *bead, SlipLink *sl, double Ns, int zzz, int n_ch, int istep, int n_m, long int *in_ran1,int *destroy ) {

int ic,nnc,i,j,l,jj,ii,k,txjj,txjjj;
int id,idjj;
int guess;
double w1,w2,lw1,du1,u1;
double rand_ncnew,rand_nncnew;
double rand_idch_sl;
double frac;
double sigma3;
int nc_new,nnc_new;
int bead_index, sup_index, id_sup;
FILE *renew;
sigma3=(sqrt(Ns)*b)/sqrt(3.0);
//renew=fopen("renew.dat","a+");
//fprintf(renew,"istep=%d\n",istep);
	
        for (j=0;j<zzz;j++) {
        /* Destruction-recreation of the paired sliplink: binary correspondence between sliplinks j& jj */
	jj=j+1;
        if (j%2==1) jj=j-1;
        if (destroy[j]==1) {
//	fprintf(renew,"SL %d destroyed!\n",j);
        /* Choose a new monomer among the n_ch chains */
	frac = drand48();
	do {
        id=rand()%(n_ch*n_m);
	} while ( bead[g_bead_index[id]].terminal==-1) ;
	sl[j].id_left_bead=id;
        /* Slip-link j recreated  */
	bead_index=g_bead_index[id];
	id_sup = bead[bead_index].id_sup;
	sup_index = g_bead_index[id_sup];
	sl[j].truncxj=bead[bead_index].index_in_chain;
        sl[j].xj=(double) (sl[j].truncxj)+frac;
	frac=sl[j].xj - (double) (sl[j].truncxj);	
   		for (k=0;k<3;k++) {
        	sl[j].XYZsj[k]=bead[bead_index].XYZ[k]+(frac)*(bead[sup_index].XYZ[k]-bead[bead_index].XYZ[k]);
	//	  fprintf(renew,"new sj coordinates for j=%d are:%lf %lf %lf\n",j,sl[j].XYZsj[0],sl[j].XYZsj[1],sl[j].XYZsj[2]);
		}
          	for (k=0;k<3;k++) {
                        w1=ran1(in_ran1);
                        w2=ran1(in_ran1);
                        lw1=log(w1);
                        du1=sqrt(-2.0*lw1);
                        u1=du1*cos(2.0*M_PI*w2);
                        sl[j].XYZaj[k]=sl[j].XYZsj[k]+u1*sigma3;
        	}
        rand_idch_sl=rand()%(n_ch);

        nnc_new=rand_idch_sl;
        
  /* Slip-link re-created on chain ends of chain nnc_new or at the end of the previous chain */

        frac=drand48();
	w1=drand48();	
        if (w1-0.5>0.000001) {
	idjj=(nnc_new)*(n_m);
	sl[jj].id_left_bead=idjj;
	bead_index=g_bead_index[idjj];
	id_sup = bead[bead_index].id_sup;
	sup_index=g_bead_index[id_sup];
	sl[jj].truncxj=bead[bead_index].index_in_chain;
        sl[jj].xj=sl[jj].truncxj+frac;
        } else {
		if (nnc_new==0)
		idjj=n_m-2;
		else  idjj=(nnc_new)*(n_m)-2; 
	sl[jj].id_left_bead=idjj;
	bead_index=g_bead_index[idjj];
	id_sup = bead[bead_index].id_sup;
	sup_index=g_bead_index[id_sup];
	sl[jj].truncxj=bead[bead_index].index_in_chain;
        sl[jj].xj=sl[jj].truncxj+frac;
        }

   for (k=0;k<3;k++) {
        sl[jj].XYZsj[k]=bead[bead_index].XYZ[k]+(frac)*(bead[sup_index].XYZ[k]-bead[bead_index].XYZ[k]);
	//	  fprintf(renew,"new sj coordinates for jj=%d are:%lf %lf %lf\n",jj,sl[jj].XYZsj[0],sl[jj].XYZsj[1],sl[jj].XYZsj[2]);
}
            for (k=0;k<3;k++) {
                        w1=ran1(in_ran1);
                        w2=ran1(in_ran1);
                        lw1=log(w1);
                        du1=sqrt(-2.0*lw1);
                        u1=du1*cos(2.0*M_PI*w2);

                        sl[jj].XYZaj[k]=sl[jj].XYZsj[k]+u1*sigma3;

            }
        } /*endif destroy[j] */
        }
//fclose(renew);
}

