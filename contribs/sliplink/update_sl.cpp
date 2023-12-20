# include <math.h>
#include "slink.h"
# include <map>

#define anint(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))
#define aint(x) ((x)>0 ? floor(x) : ceil((x)))
#define nint(x) ((x)>0 ? (int)((x)+0.5) : (int)((x)-0.5)


void update_sl(Bead *bead, SlipLink *sl, int zzz, int n_m, int n_ch, int *num_SL, int istep)
{
int k,i,j,txjj;
int ic,id;
int bead_index,sup_index;
int id_sup;
FILE *update_sl;
//update_sl=fopen("update.dat","w");
for (ic=0;ic<n_ch;ic++) num_SL[ic]=0;
for (j=0;j<zzz;j++) {
ic=sl[j].chainj;
num_SL[ic]=num_SL[ic]+1;
}


/*Update sliplink coordinates */
        for (j=0;j<zzz;j++) {
        ic=sl[j].chainj;
        txjj=(int) (aint(sl[j].xj));
        sl[j].truncxj=txjj;
	id=sl[j].id_left_bead;
        bead_index = g_bead_index[id];
	id_sup= bead[bead_index].id_sup;
	sup_index = g_bead_index[id_sup];
                for (k=0;k<3;k++) {
        sl[j].XYZsj[k]=bead[bead_index].XYZ[k]+(sl[j].xj-(double) (sl[j].truncxj))*(bead[sup_index].XYZ[k]-bead[bead_index].XYZ[k]);
        }
}
//fclose(update_sl);
}
