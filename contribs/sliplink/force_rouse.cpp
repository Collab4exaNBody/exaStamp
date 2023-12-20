# include <math.h>
# include "slink.h"
# include <map>

void force_rouse(Bead *bead, double b, double kBT,int n_ch, int n_m, int istep)
{
int j,k,i;
double cte7;
FILE *frouse;
int bead_index,id_sup,sup_index;
cte7=3.0*kBT/(b*b);
//frouse=fopen("frouse.dat","w");
//fprintf(stdout,"istep=%d\n",istep);
        int bead_size = n_m*n_ch;

        for (i=0;i<bead_size;i++) {
	bead_index=i; //g_bead_index[i]; 
        id_sup = bead[bead_index].id_sup;
        sup_index = g_bead_index[ id_sup ];
//	fprintf(stdout,"i=%d,bead_index=%d,id_sup=%d,sup_index=%d\n",i,bead_index,id_sup,sup_index);
	if (bead[bead_index].terminal==0) {
	//bead_index=g_bead_index[i]; 
        //id_sup = bead[bead_index].id_sup;
        //sup_index = g_bead_index[ id_sup ];
        for (k=0;k<3;k++) {
        bead[bead_index].F_SP[k] = bead[bead_index].F_SP[k] + ( bead[sup_index].XYZ[k] - bead[bead_index].XYZ[k] )*cte7;
	//fprintf(stdout,"%lf\n",bead[i].F_SP[k]);
        bead[sup_index].F_SP[k] = bead[sup_index].F_SP[k] - ( bead[sup_index].XYZ[k] - bead[bead_index].XYZ[k] )*cte7;
}
}
}
//	fprintf(frouse,"%lf %lf %lf\n",bead[11].F_SP[0], bead[11].F_SP[1], bead[11].F_SP[2]);
//fclose(frouse);
}
