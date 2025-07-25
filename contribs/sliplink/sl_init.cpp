/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

# include <math.h>
# include <map>
#include <stdio.h>
#include "slink.h"

#define anint(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))
#define aint(x) ((x)>0 ? floor(x) : ceil((x)))
#define nint(x) ((x)>0 ? (int)((x)+0.5) : (int)((x)-0.5)


void sl_init(SlipLink *sl, Bead *bead, double b, double Ns, int zzz, int n_e, int n_ch, int n_m,long int *in_ran1, int *num_SL )
{
int i,ij,jf,j,jjj,k,l,ii,nc,nnc,q,m;
int id,idjj;
double sigma3;
double frac;
double ranxj,rand_idch_sl;
int nbead;
int bead_index,id_sup,sup_index;
double w1,w2,lw1,du1,u1;
FILE *slcoorfile,*slidfile;

/* Initialisation of SLIPLINKS  */
sigma3=(sqrt(Ns)*b)/sqrt(3.0);
//slcoorfile=fopen("slcoor.dat","w");
for (m=0;m<n_ch;m++) num_SL[m]=0;

/* Initial construction of SLs */

for (j=0;j<zzz/2;j++) {
        q=2*j+1;

/* Chain chosen randomly */

   do
   {
   	id = rand()%(n_m*n_ch);
   }
   while( bead[ g_bead_index[id] ].terminal==-1 );

   bead_index = g_bead_index[id];
   id_sup = bead[bead_index].id_sup;
   sup_index = g_bead_index[id_sup];

   frac=drand48();
   sl[q].id_left_bead = id;
   nc = bead[bead_index].id_ch;
   sl[q].chainj = nc;
   num_SL[nc]++;

   sl[q].truncxj = bead[bead_index].index_in_chain;
   sl[q].xj = (double) (sl[q].truncxj)+frac;
   
	frac=sl[q].xj-(double) (sl[q].truncxj);
   	for (k=0;k<3;k++) {
        sl[q].XYZsj[k] = bead[bead_index].XYZ[k] + frac * ( bead[sup_index].XYZ[k] - bead[bead_index].XYZ[k] );
	
        }
   	for (k=0;k<3;k++) {
                        w1=ran1(in_ran1);
                        w2=ran1(in_ran1);
                        lw1=log(w1);
                        du1=sqrt(-2.0*lw1);
                        u1=du1*cos(2.0*M_PI*w2);
                        sl[q].XYZaj[k]=sl[q].XYZsj[k]+u1*sigma3;
        }

   do
   {
   	idjj = rand()%(n_m*n_ch);
   }
   while( bead[ g_bead_index[idjj]].terminal==-1 );

    bead_index = g_bead_index[idjj];
    id_sup = bead[bead_index].id_sup;
    sup_index = g_bead_index[id_sup];
   sl[q-1].id_left_bead=idjj;
   nnc=bead[bead_index].id_ch;
   sl[q-1].chainj=nnc;
   num_SL[nnc]++;
   sl[q-1].truncxj=bead[bead_index].index_in_chain;
   sl[q-1].xj=(double) (sl[q-1].truncxj)+frac;
	frac=sl[q-1].xj-(double) (sl[q-1].truncxj);
   	for (k=0;k<3;k++) {
        sl[q-1].XYZsj[k] = bead[bead_index].XYZ[k]+ frac*(bead[sup_index].XYZ[k]-bead[bead_index].XYZ[k]);

        }
   	for (k=0;k<3;k++) {
                        w1=ran1(in_ran1);
                        w2=ran1(in_ran1);
                        lw1=log(w1);
                        du1=sqrt(-2.0*lw1);
                        u1=du1*cos(2.0*M_PI*w2);

                        sl[q-1].XYZaj[k]=sl[q-1].XYZsj[k]+u1*sigma3;
        }

}
        for (nc=0;nc<n_ch;nc++) num_SL[nc]=0;

        for (j=0;j<zzz;j++) {
        nc=sl[j].chainj;
        num_SL[nc]=num_SL[nc]+1;
        }
        for (j=0;j<zzz;j++) {
        nc=sl[j].chainj;
        }
}

