# include <math.h>
# include "slink.h"
# include <assert.h>
# include <map>
# include <cmath>

void move_sl(Bead *bead, SlipLink *sl, double b, double h, double kBT, double ksi, double Ns, int n_ch, int n_m, long int *in_ran1, int zzz, int *destroy)
{
int i,j,l,txjj,tid;
double cte4,sigma4,ksis;
double w1,w2,lw1,u1,du1,xjtmp;
int bead_index,id_sup,sup_index;
FILE *slmove;

//slmove=fopen("slmove.dat","w");
ksis=0.1*ksi;
cte4=3.0*kBT/(Ns*b*b*ksis);
sigma4=sqrt((2.0*h*kBT)/ksis);
/* Langevin equation for slip-links */
for (j=0;j<zzz;j++) {
        int oldtruncxj = static_cast<int>( std::floor(sl[j].xj) );
        assert( oldtruncxj == sl[j].truncxj );

        tid=sl[j].id_left_bead;
	bead_index=g_bead_index[tid];
	id_sup = bead[bead_index].id_sup;	
	sup_index = g_bead_index[id_sup];

	assert( sl[j].truncxj == bead[bead_index].index_in_chain );

        w1=ran1(in_ran1);
        w2=ran1(in_ran1);
        lw1=log(w1);
        du1=sqrt(-2.0*lw1);
        u1=du1*cos(2.0*M_PI*w2);
//fprintf(slmove,"For SL %d tid=%d, sup_index=%d\n",j,tid, sup_index);
//fprintf(slmove,"For SL %d ajmsjx=%lf, ajmsjy=%lf, ajmsjz=%lf\n",j,sl[j].XYZaj[0]-sl[j].XYZsj[0], sl[j].XYZaj[1]-sl[j].XYZsj[1], sl[j].XYZaj[2]-sl[j].XYZsj[2]);
        xjtmp=sl[j].xj+cte4*h*((bead[sup_index].XYZold[0]-bead[bead_index].XYZold[0])*(sl[j].XYZaj[0]-sl[j].XYZsj[0])
        +                      (bead[sup_index].XYZold[1]-bead[bead_index].XYZold[1])*(sl[j].XYZaj[1]-sl[j].XYZsj[1])

        +                      (bead[sup_index].XYZold[2]-bead[bead_index].XYZold[2])*(sl[j].XYZaj[2]-sl[j].XYZsj[2]))+sigma4*u1;
        sl[j].xj=xjtmp;

	int newtruncxj = static_cast<int>( std::floor(sl[j].xj) );
	if( newtruncxj < oldtruncxj )
	{
		sl[j].id_left_bead = bead[bead_index].id_inf;
	}
	if( newtruncxj > oldtruncxj )
	{
		sl[j].id_left_bead = bead[bead_index].id_sup;
	}
}
/* Destroy=1 if a slip-link leaves the chain */
 for (j=0;j<zzz;j++) {
        if ((sl[j].xj<0.0)||(sl[j].xj>(double) (n_m-1))) 
	{
	destroy[j]=1;
	//fprintf(slmove,"For j=%d, destroy[j]=%d\n",j,destroy[j]);
	}
        }
//fclose(slmove);
}

