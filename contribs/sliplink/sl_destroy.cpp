# include <math.h>
# include "slink.h"

 void sl_destroy(int *destroy, int zzz)
{
int j,q,jj;

        for (j=0;j<zzz;j++) {
        /* Destruction-recreation of the paired SL : Binary correspondance between SLs j and jj */
        q=j+1;

        if ((destroy[q]==1)&&(destroy[q-1]==1)) {
//	fprintf(stdout,"SL %d destroyed!\n",q);
        destroy[q-1]=0;
	}
        }
}

