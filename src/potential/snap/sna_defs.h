#pragma once

#include <cmath>

#define ULIST_IDIM      idxu_max
#define DULIST_IDIM     3
#define SUM_INT_SQR(i)      ( ( (i) * ((i)+1) * (2*(i)+1) ) / 6 )
#define IDXU_BLOCK(i)       SUM_INT_SQR(i)

#ifndef SNA_NO_ARRAY_ACCESSOR_MACROS

#define SNA_ARRAY(a,i)    a [i]

#define ULIST_R_JJ(i)   SNA_ARRAY( ulist_r_ij_jj , (i) )
#define ULIST_I_JJ(i)   SNA_ARRAY( ulist_i_ij_jj , (i) )
#define ULISTTOT_R(i)   SNA_ARRAY(ulisttot_r,i)
#define ULISTTOT_I(i)   SNA_ARRAY(ulisttot_i,i)
#define DULIST_R(j,i)   SNA_ARRAY(dulist_r, (j) * DULIST_IDIM + (i) )
#define DULIST_I(j,i)   SNA_ARRAY(dulist_i, (j) * DULIST_IDIM + (i) )
#define ZLIST_R(i)      SNA_ARRAY(zlist_r,i)
#define ZLIST_I(i)      SNA_ARRAY(zlist_i,i)
#define BLIST(i)        SNA_ARRAY(blist,i)
#define YLIST_R(i)      SNA_ARRAY(ylist_r,i)
#define YLIST_I(i)      SNA_ARRAY(ylist_i,i)

#define BZERO(i)            bzero[i]
#define IDXZ(i)             idxz[i]
#define IDXZ_ALT(i)         idxz_alt[i]
#define IDXB(i)             idxb[i]
#define ROOTPQARRAY(p,q)    rootpqarray[((p)-1)*twojmax+(q)-1] /// = sqrt(static_cast<double>(p)/(q)) // rootpqarray[i][j] // dims = twojmax x twojmax [1;twojmax] => indices [0;twojmax-1] => side = twojmax
#define CGLIST(i)           cglist[i]
#define IDXCG_BLOCK(k,j,i)  idxcg_block[((k)*(twojmax+1)+(j))*(twojmax+1)+(i)]  // idxcg_block[k][j][i] // dims = (twojmax+1) x (twojmax+1) x (twojmax+1)
#define IDXZ_BLOCK(k,j,i)   idxz_block[((k)*(twojmax+1)+(j))*(twojmax+1)+(i)] // idxz_block[k][j][i]  // dims = (twojmax+1) x (twojmax+1) x (twojmax+1)
#define IDXB_BLOCK(k,j,i)   idxb_block[((k)*(twojmax+1)+(j))*(twojmax+1)+(i)] // idxb_block[k][j][i]  // dims = (twojmax+1) x (twojmax+1) x (twojmax+1)

#endif

