# include <math.h>
# include "slink.h"
# include <map>

void force_sliplink(Bead *bead, SlipLink *sl, double b, double kBT,int n_ch, int n_m, double Ns, int zzz, int istep)
{
int id,i,j,k,l;
double cte8,m0,m1,m2;
FILE *fsl,*idfile;
int bead_index,id_sup,sup_index;
int bead_size=n_m*n_ch;
//fsl=fopen("fsl.dat","a+");
cte8=3.0*kBT/(Ns*b*b);
//fprintf(fsl,"istep=%d\n",istep);
for (j=0;j<zzz;j++) {
	
	id=sl[j].id_left_bead;
	bead_index=g_bead_index[id];
	id_sup=bead[bead_index].id_sup;
	sup_index = g_bead_index[id_sup];
        bead[bead_index].F_SL[0]=bead[bead_index].F_SL[0]+cte8*(1.0-(sl[j].xj-(double)(sl[j].truncxj)))*(sl[j].XYZaj[0]-sl[j].XYZsj[0]);
        bead[bead_index].F_SL[1]=bead[bead_index].F_SL[1]+cte8*(1.0-(sl[j].xj-(double)(sl[j].truncxj)))*(sl[j].XYZaj[1]-sl[j].XYZsj[1]);
        bead[bead_index].F_SL[2]=bead[bead_index].F_SL[2]+cte8*(1.0-(sl[j].xj-(double)(sl[j].truncxj)))*(sl[j].XYZaj[2]-sl[j].XYZsj[2]);
        m0=bead[sup_index].F_SL[0]+cte8*(sl[j].xj-(double) (sl[j].truncxj))*(sl[j].XYZaj[0]-sl[j].XYZsj[0]);
        m1=bead[sup_index].F_SL[1]+cte8*(sl[j].xj-(double) (sl[j].truncxj))*(sl[j].XYZaj[1]-sl[j].XYZsj[1]);
        m2=bead[sup_index].F_SL[2]+cte8*(sl[j].xj-(double) (sl[j].truncxj))*(sl[j].XYZaj[2]-sl[j].XYZsj[2]);
        bead[sup_index].F_SL[0]=m0;
        bead[sup_index].F_SL[1]=m1;
        bead[sup_index].F_SL[2]=m2;
	//if (istep%100==0) 
        //fprintf(fsl," j=%d, bead id= %d  %lf  %lf  %lf \n",j,id,bead[bead_index].F_SL[0],bead[bead_index].F_SL[1],bead[bead_index].F_SL[2]);

//	fprintf(fsl," %lf  %lf %lf \n",sl[j].XYZaj[0]-sl[j].XYZsj[0],sl[j].XYZaj[1]-sl[j].XYZsj[1], sl[j].XYZaj[2]-sl[j].XYZsj[2]);

}
//fclose(fsl);
}
