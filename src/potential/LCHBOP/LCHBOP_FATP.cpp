#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/algorithm.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>

#include "lchbop_FATP.h"
//#include "grid_particle_nbtbl.h"
#include "lchbop_utils.h"
#include "lchbop_common.h"



using std::sqrt;
//using std::min;

void subdsig(Vec3d, Vec3d*);
inline Mat3d tensorsum(Mat3d&,double (&t)[3][3]);

//------------------------------------------------------------------------//


void FATP::subFcnj_CH(int ispcc,int ispch,int ifracnb,
int ipmrnb,double rij,Vec3d sigij,double vsmra,vector<uint8_t>& ispc,
LNBStr* nbsij,double& fcnij,double* pccij,FATPLoc& Ftbls,NBStr& nbs)
//void FATP::subFcnj_CH(uint8_t ispcc,uint8_t ispch,int ifracnb,
//int ipmrnb,double rij,Vec3d sigij,double vsmra,LNBStr* nbsij,
//double& fcnij,double* pccij,FATPLoc& Ftbls,NBStr& nbs)
{
  int nnij3=0,iwk[3],nnij[2][4][2];  // initialisation of nnij3 needs te be checked
  int nijfl[2],nchij[2];
  Vec3d dfrc;                     // size needs to be adapted ?!!!!! //

//  FILE *VIRFRC2;    // for testing only
//  VIRFRC2=fopen("TMPDAT/virfrc2.dat","a");    // for testing only
//  size_t iptest=-1;

  double snki1,dsnki1,snki2,dsnki2;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
//CC Comments:
//CC	 + for ij=0, ip side is calculated	 
//CC	 + for ij=1, jp side is calculated		 
//CC	 + nijfull = number full neighbors of i[j] not equal to j[i] 
//CC	 + nijfrac = number SR-fractional neighbors of i[j] not equal to j[i]	
//CC	 + lstcnb	= list of neighbor indices of conneted neighbors 
//CC	 + nij		= number connected neighbors k to atom i
//CC	 + nji		= number connected neighbors l to atom j
//CC	 + array nnij[2][nij][2]
//CC		 - nnij[ij][nij][0] = first index of configurations for atom ijp[ip]=ip/jp with NIJ connected atoms k/l not equal to jp/ip
//CC		 - nnij[ij][nij][1] = last  index of configurations for atom ijp[ip]=ip/jp with NIJ connected atoms k/l not equal to jp/ip
//CC	 + nbr		= number of connection branches
//CC		 - E.g., if nijfull=2 and nijfrac=3, there are 8 possible primary branches of connectivity configurations: 
//CC										[11000],[11100],[11010],[11001],[11110],[11101],[11011],[11111]
//CC 		         however, all connectivity configurations with 3 full neighbours can be taken together, 
//CC                     as for nij>=3 Fconj does not depend on Ncj (fcj0=fcj1), so that there are only 2 relevant connectivity configurations, 
//CC                     namely: [11000] and [11100]+[11010]+[11001]+[11110]+[11101]+[11011]+[11111]
//CC Nota bene !!!
//CC In Stamp the components of snstor are defined as:
//CC 0) snstor[0] = Sn
//CC 1) snstor[1] = dSn/dr
//CC 2) snstor[2] = Sn-1
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//

//C Initialisations
  fcnij=0.;
  pccij[0]=pccij[1]=0.;
  
//C Compute configuration-tree & PCC for both ijp[0] and ijp[1]

  for(int ij=0;ij<2;ij++)
  {
    int inbj=nbsij[ij].inbj;
    size_t ip=nbsij[ij].nbtbl[0];
//    size_t jp=nbsij[ij].nbtbl[nbsij[ij].inbj];    // for testing only
    int nijfull=nbsij[ij].nnbfl-2+ifracnb;
    ijp[ij]=ip;
    nijfl[ij]=nijfull;
    
    if(nijfull<3) // ip is unsaturated: connectivity tree required
    {
//C precompute prwij ((1.-snik)-product for frac. neigh's; snstor2=(snik-1.))
      double prwijmin=1.;
      for(int inb=nbsij[ij].nnbfr-1;inb>=nbsij[ij].nnbfl;inb--)
      {
        if(inb==inbj)continue;                   // inbj has to be initialized !!!!! //
      	prwijmin*=-nbsij[ij].snstor[inb][0];
        assert(inb<nbsij[ij].nnbsrm+1);
      	Ftbls.prwij[inb]=prwijmin;
      }
      Ftbls.prwij[nbsij[ij].nnbfr]=1.;   // what is going on here ????? //
//      Ftbls.prwij[nbsij[ij].nnbfr+1]=1.;   // what is going on here ????? //
      
// nnij stores min/max indices of tree portions containing branches with 
// fixed numbers of connected neigh. (0 for min, 1 for max)
      for(int k=0;k<4;k++){for(int l=0;l<2;l++)nnij[ij][k][l]=0;}
      nnij[ij][nijfull][0]=1;   // index change for nnij done  //
      nnij[ij][nijfull][1]=1;
      int nij=nijfull;
      int nijfrac=nbsij[ij].nnbfr-nbsij[ij].nnbfl-ifracnb+ipmrnb;  // are ifracnb and ipmrnb defined ?!!!!!!! //
//      if(ip==4&&nbsij[ij].nbtbl[inbj]==86)
//      {
//        printf("ij,ip,jp,nijfull,nijfrac,ifracnb :%6d%6d%6d%6d%6d%6d\n",ij,ip,nbsij[ij].nbtbl[inbj],nijfull,nijfrac,ifracnb);
//        std::cin.get();
//      }  

      for(int i=1;i<=nijfrac;i++) 
      {
      	int nijp1=nij+1;
      	if(nijp1==3) 
        {
      	  nnij3=nnij[ij][2][1]+1;
      	  nnij[ij][3][0]=nnij3;
      	  nnij[ij][3][1]=nnij3;
      	  wijstor[ij][nnij3]=0.;   // index change for wijstor done  //
      	  break;
      	}
      	nnij[ij][nijp1][0]=nnij[ij][nij][1]+1;
      	nnij[ij][nijp1][1]=nnij[ij][nij][1]+imat[i][nijfrac];  // check imat !!! //
      	nij=nijp1;
      }

// calculate Nelijbr[1] = # of avail. e- after removing all full neigh. contribs
      int icnb=0;
      int nbr=1;     // first branch; nbr=number of connectivity branches
      nijbr[nbr]=nijfull; // number of connected neigh. in branch 1
      wijbr[nbr]=1.;      // weight of branch 1
      
// list of connected neigh. for the first branch
//      if((ip==8&&jp==342)||(jp==8&&ip==342))
//      {
//        printf("OK1; ij,ip,jp,nijfull :%6d%6d%6d%6d\n",ij,ip,jp,nijfull);
//      }
      for(int inb=1;inb<nbsij[ij].nnbfl;inb++)  // Is this correct ? !!!!!!!! //
      {
        if(inb==inbj)continue;                   // inbj has to be initialized !!!!! //
//        size_t kp=nbsij[ij].nbtbl[inb];
//      	int ispck=nbsij[ij].ispc[inb];
      	icnb=icnb+1;
        assert(icnb<3);
      	lstcnbbr[nbr][icnb]=inb;     // index change for lstcnnbr done //
//        if(ij==1&&nbr==6&&inb>=nbsij[ij].nnbsrm)
//        {
//          printf("OK1; ij,nbr,icnb,lstcnbbr :%6d%6d%6d%6d\n",ij,nbr,icnb,lstcnbbr[nbr][icnb]);
//        }  
      }

      if(nijfull==0)  // Case Nij=0 for first branch
      {
      	Nelijbr[nbr]=4.;         // 4 e- by default;  index change for dNelijbrdik and dNelijbrdki done  //			 
//        dNelijbrdik[nbr][0]=0.  // CAREFUL: for now we store all force contribs in the ispcc component // ????? //
//        dNelijbrdik[nbr][1]=0.  // CAREFUL: for now we store all force contribs in the ispcc component // ????? //
      	dNelijbrdki[nbr][0]=0.;	 // since there is no more specie-dependence in derivatives
      	dNelijbrdki[nbr][1]=0.;	 // since there is no more specie-dependence in derivatives
      } 
      else if(nijfull==1) // Case Nij=1 for first branch
      {
        int inb=lstcnbbr[nbr][1];
        int ispck=nbsij[ij].ispc[inb];
//        uint8_t ispck=nbsij[ij].ispc[inb];
#ifdef F1NEW
        if(ispck==ispch)
        {
          Nelijbr[nbr]=3.0;
//          dNelijbrdik[nbr][0]=0.;    // CAREFUL: for now we store all force contribs in the ispcc component
//          dNelijbrdik[nbr][1]=0.;    // CAREFUL: for now we store all force contribs in the ispcc component
          dNelijbrdki[nbr][0]=0.;    // since there is no more specie-dependence in derivatives
          dNelijbrdki[nbr][1]=0.;    // since there is no more specie-dependence in derivatives
//          printf("1) No specy dependence \n");
        }
        else
        {  
          size_t kp=nbsij[ij].nbtbl[inb];
          double Nki=nbs.crd[kp]-nbsij[ij].snstor[inb][1];
          double NCki=nbs.crdspc[kp][ispcc]-nbsij[ij].snstor[inb][1];
          double dnelij,dnelijdnki[2]; 
          F1new(ispcc,ispch,Nki,NCki,dnelij,dnelijdnki);
//          if((ip==8&&jp==342)||(jp==8&&ip==342))
//          printf("1) ispck,Nki,NCki,dnelij,dnelijdnki : %6d%16.10lf%16.10lf%16.10lf%16.10lf%16.10lf\n",
//          ispck,Nki,NCki,dnelij,dnelijdnki[0],dnelijdnki[1]);
          Nelijbr[nbr]=2.0+dnelij;
          dNelijbrdki[nbr][0]=dnelijdnki[0];    // since there is no more specie-dependence in derivatives
          dNelijbrdki[nbr][1]=dnelijdnki[1];    // since there is no more specie-dependence in derivatives
        }
#else
        if(ispck==ispch)   // ispch needs to defined here !!!!!! //
        {
          snki=1.;   // H neigh. always contribute for 1 e- (no force)
          dsnki=0.;  //
        } 
        else 
        {
      	  size_t kp=nbsij[ij].nbtbl[inb];
          double Nki=nbs.crd[kp]-nbsij[ij].snstor[inb][1];	// reduced coordination of k wrt i
          SqdSqpup(Nki,1.,3.,0.5,p2sm1,snki,dsnki);
        }
        Nelijbr[nbr]=2.+snki;      // C neigh. contribution depends on reduced coordination = sigma-bond (1) + partial pi-bond (1-snik)
//        printf("1) ij,nbr,Nelijbr : %6d%6d%16.10lf\n",ij,nbr,Nelijbr[nbr]);
//        dNelijbrdik[nbr][0]=0.;	   // CAREFUL: for now we store all force contribs in the ispcc component
//        dNelijbrdik[nbr][1]=0.;    // CAREFUL: for now we store all force contribs in the ispcc component
        dNelijbrdki[nbr][0]=dsnki; // since there is no more specie-dependence in derivatives
        dNelijbrdki[nbr][1]=0.;    // since there is no more specie-dependence in derivatives
#endif
      } 
      else if(nijfull==2)  //CASE Nij=2 for first branch
      {
        int inb=lstcnbbr[nbr][1];  // contrib. of 1st connected neigh.
      	int ispck=nbsij[ij].ispc[inb];
//        uint8_t ispck=nbsij[ij].ispc[inb];
        if(ispck==ispch)
        {
          snki1=1.;   // H neigh. always contribute for 1 e- (no force)
          dsnki1=0.;
        } 
        else 
        {
      	  size_t kp=nbsij[ij].nbtbl[inb];
          double Nki1=nbs.crd[kp]-nbsij[ij].snstor[inb][1];	// reduced coordination of k wrt i
          SqdSqpup(Nki1,2.,3.,1.,p2sm2,snki1,dsnki1);
//          SqdSqpup(Nki1,2.,3.,1.,-3.,snki1,dsnki1);
        }
        
        inb=lstcnbbr[nbr][2];   // contrib. of 2nd connected neigh.
      	ispck=nbsij[ij].ispc[inb];
        if(ispck==ispch)
        {
          snki2=1.;             // H neigh. always contribute for 1 e- (no force)
          dsnki2=0.;            // 
        } 
        else 
        {
      	  size_t kp=nbsij[ij].nbtbl[inb];
          double Nki2=nbs.crd[kp]-nbsij[ij].snstor[inb][1];	// reduced coordination of k wrt i
          SqdSqpup(Nki2,2.,3.,1.,p2sm2,snki2,dsnki2);
//          SqdSqpup(Nki2,2.,3.,1.,-3.,snki2,dsnki2);
        }
        
        Nelijbr[nbr]=(4.+snki1+snki2)/3.;  // C neigh. contributions depend on reduced coordinations
        dNelijbrdki[nbr][0]=oneby3*dsnki1; //	since there is no more specie-dependence in derivatives
        dNelijbrdki[nbr][1]=oneby3*dsnki2; //since there is no more specie-dependence in derivatives
//        if((ip==8&&jp==342)||(jp==8&&ip==342))
//        {
//          printf("1) ij,nbr,Nelijbr : %6d%6d%16.10lf%16.10lf%16.10lf\n",ij,nbr,Nelijbr[nbr],dNelijbrdki[nbr][0],dNelijbrdki[nbr][1]);
//        }
      }

// run over frac. neigh to build further branches (connected/disconnected tree)
      int inbjm=nbsij[ij].nnbfr;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)   // loop over fractional neighbours
      {
        if(inb==inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inb];
      	int ispck=nbsij[ij].ispc[inb];
//        uint8_t ispck=nbsij[ij].ispc[inb];
      	double snik=nbsij[ij].snstor[inb][1];
        int nbrold=nbr;          // store current branch number (nbr may be incremented within the loop)
//        uint16_t nbrold=nbr;          // store current branch number (nbr may be incremented within the loop)
        for(int ibr=1;ibr<=nbrold;ibr++) // loop over all existing branches
//        for(uint16_t ibr=1;ibr<=nbrold;ibr++) // loop over all existing branches
        {
//          if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//          {
//            printf("OK1; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//          }  
//          if(nbr>6)printf("OKOKOK2a; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
          double wij=wijbr[ibr];      // store "parent" branch weight
          wijbr[ibr]=wij*(1.-snik);   // current fractional neighbour disconnected for parent branch: decreasing branch weight
          int nij=nijbr[ibr];     // # of connected neighbours of parent branch
          int nijp1=nij+1;        // # of connected neighbours of potential new branch

          if(nijp1==1)  // CASE Nij==1
          {
//            if(nbr>6)printf("OKOK1a; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
            nbr++;                    // new branch number
            assert(nbr<nbsij[ij].nbrm);
            nijbr[nbr]=nijp1;         // store the # of connected neigh. of the new branch
            lstcnbbr[nbr][1]=inb;     // 1st connected neigh. of the new list
            wijbr[nbr]=wij*snik;      // weight of new branch (parent weight * cutoff_inb)
#ifdef F1NEW
            if(ispck==ispch)
            {
              Nelijbr[nbr]=3.0;
              dNelijbrdki[nbr][0]=0.;  // since there is no more specie-dependence in derivatives
              dNelijbrdki[nbr][1]=0.;  // since there is no more specie-dependence in derivatives
//              if((ip==8&&jp==342)||(jp==8&&ip==342))
//              printf("2) No specy dependence \n");
            }
            else
            {  
              double Nki=nbs.crd[kp]-nbsij[ij].snstor[inb][1];
              double NCki=nbs.crdspc[kp][ispcc]-nbsij[ij].snstor[inb][1];
              double dnelij,dnelijdnki[2]; 
              F1new(ispcc,ispch,Nki,NCki,dnelij,dnelijdnki);
//              if((ip==8&&jp==342)||(jp==8&&ip==342))
//              printf("2) ispck,Nki,NCki,dnelij,dnelijdnki : %6d%16.10lf%16.10lf%16.10lf%16.10lf%16.10lf\n",
//              ispck,Nki,NCki,dnelij,dnelijdnki[0],dnelijdnki[1]);
              Nelijbr[nbr]=2.0+dnelij;
              dNelijbrdki[nbr][0]=dnelijdnki[0];    // since there is no more specie-dependence in derivatives
              dNelijbrdki[nbr][1]=dnelijdnki[1];    // since there is no more specie-dependence in derivatives
            }
//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK2; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  
#else
            if(ispck==ispch)
            {
              snki=1.;         // H neigh. always contribute for 1 e- (no force)
              dsnki=0.;
            } 
            else
            {
              double Nki=nbs.crd[kp]-snik;   // reduced coordination of k wrt i
              SqdSqpup(Nki,1.,3.,0.5,p2sm1,snki,dsnki);
            }

            Nelijbr[nbr]=2.+snki;       // C neigh. contribution depends on reduced coordination
//            printf("3) ij,nbr,Nelijbr : %6d%6d%16.10lf\n",ij,nbr,Nelijbr[nbr]);
//            dNelijbrdik[nbr][0]=0.;     // CAREFUL: for now we store all force contribs in the ispcc component
//            dNelijbrdik[nbr][1]=0.;     // CAREFUL: for now we store all force contribs in the ispcc component
            dNelijbrdki[nbr][0]=dsnki;  // since there is no more specie-dependence in derivatives
            dNelijbrdki[nbr][1]=0.;     // since there is no more specie-dependence in derivatives
//            if(nbr>6)printf("OKOK1b; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
#endif
          }
          else if (nijp1==2)  // CASE Nij==2
          {
//            if(nbr>6)printf("OKOK2a; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
            nbr++;                          // new branch number
            nijbr[nbr]=nijp1;               // store the # of connected neigh. of the new branch
            int inb0=lstcnbbr[ibr][1];  // new list of connected neigh.= parent list + inb 
            lstcnbbr[nbr][1]=inb0;          //
            lstcnbbr[nbr][2]=inb;           //
//          if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//          {
//            printf("OK3; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//          }  
//            printf("OK2a; ij,nbr,icnb,lstcnbbr :%6d%6d%6d%6d\n",ij,nbr,1,lstcnbbr[nbr][1]);
//            printf("OK2b; ij,nbr,icnb,lstcnbbr :%6d%6d%6d%6d\n",ij,nbr,2,lstcnbbr[nbr][2]);
            wijbr[nbr]=wij*snik;            // weight of new branch (parent weight * cutoff_inb)
                
// contrib. of 1st connected neigh.
      	    int ispck=nbsij[ij].ispc[inb0];
//            uint8_t ispck=nbsij[ij].ispc[inb0];
//            printf("inb0,ispck : %6d%6zu\n",inb0,ispck);
            if(ispck==ispch)
            {
              snki1=1.;              // H neigh. always contribute for 1 e- (no force)
              dsnki1=0.;
            } 
            else 
            {
      	      size_t kp0=nbsij[ij].nbtbl[inb0];                    // recover the ID of the 1st connected neigh.
//              printf("inb0,kp0 : %6d%6zu\n",inb0,kp0);
              double Nki1=nbs.crd[kp0]-nbsij[ij].snstor[inb0][1];  // reduced coordination of k1 wrt i
              SqdSqpup(Nki1,2.,3.,1.,p2sm2,snki1,dsnki1);
//              SqdSqpup(Nki1,2.,3.,1.,-3.,snki1,dsnki1);
            }
                  
// contrib. of 2nd connected neigh.
      	    ispck=nbsij[ij].ispc[inb];
 //           printf("inb,ispck,ispch : %6d%6d%6d\n",inb,ispck,ispch);
            if(ispck==ispch) 
            {
              snki2=1.;                 // H neigh. always contribute for 1 e- (no force)
              dsnki2=0.;
            } 
            else 
            {
//              size_t kp=nbsij[ij].nbtbl[inb];  // recover the ID of the 1st connected neigh.
//              printf("inb,kp : %6d%6zu\n",inb,kp);
              double Nki2=nbs.crd[kp]-snik;    // reduced coordination of k2 wrt; is kp initialized ?!!!!!! //
              SqdSqpup(Nki2,2.,3.,1.,p2sm2,snki2,dsnki2);
//              SqdSqpup(Nki2,2.,3.,1.,-3.,snki2,dsnki2);
            }

            Nelijbr[nbr]=(4.+snki1+snki2)/3.;  // C neigh. contributions depend on reduced coordinations
//            printf("4) ij,nbr,Nelijbr : %6d%6d%16.10lf\n",ij,nbr,Nelijbr[nbr]);
//            printf("4) snki1,snki2    : %16.10lf%16.10lf\n",snki1,snki2);
            dNelijbrdki[nbr][0]=oneby3*dsnki1; // since there is no more specie-dependence in derivatives
            dNelijbrdki[nbr][1]=oneby3*dsnki2; //          since there is no more specie-dependence in derivatives
//            if(nbr>6)printf("OKOK2b; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK4; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  
//            if((ip==8&&jp==342)||(jp==8&&ip==342))
//            {
//              printf("2) ij,nbr,Nelijbr : %6d%6d%16.10lf%16.10lf%16.10lf\n",ij,nbr,Nelijbr[nbr],dNelijbrdki[nbr][0],dNelijbrdki[nbr][1]);
//            }
          } 
          else   // CASE Nij>=3 (add to nnij3 branch, no influence of Nconjij) 
          {
//            if(nbr>6)printf("OKOK3a; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
            double wijk=wij;	        // save prefactor for the force and pcc
            wij*=snik;	                // weight of the new branch
            wijstor[ij][nnij3]+=wij;    // global weight of branches with nij>=3 (no dependance on Nconjij, so can be accumulated)
//            if(nbr>6)printf("OKOK3aa; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK5a; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  
            
// contrib. of inb to the forces & virial
            double pf=wijk*nbsij[ij].snstor[inb][2];
            Vec3d dfrc=pf*nbsij[ij].sigstor[inb];
            Ftbls.frc2[ij][inb]=Ftbls.frc2[ij][inb]+dfrc;
#ifdef VIRIELLOCAL
      	    double rik=nbsij[ij].rstor[inb];
            Ftbls.virielfrc2[ij][inb]=Ftbls.virielfrc2[ij][inb]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
      	    double rik=nbsij[ij].rstor[inb];
            Ftbls.virielfrc2[ij][0]=Ftbls.virielfrc2[ij][0]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
                  
// contrib. of the previous inb neighbours (already considered for connection/disconnection)
            for(int icnb=nijfull+1;icnb<=2;icnb++)  // distinguish between connected and disconnected frac. neigh. using iscnb table
            {
              int inbp=lstcnbbr[ibr][icnb];
              assert(inbp!=inbj);
              assert(inbp<nbsij[ij].nnbsrm);
              Ftbls.iscnb[inbp]=1;
            }

            assert(inb<nbsij[ij].nnbsrm);
            for(int inbp=nbsij[ij].nnbfl;inbp<inb;inbp++)
            { 
              if(inbp==inbj){continue;}
//              int ispckp=nbsij[ij].ispc[inbp];
//              uint8_t ispckp=nbsij[ij].ispc[inbp];
              double snikfl=nbsij[ij].snstor[inbp][Ftbls.iscnb[inbp]]; // (Sn-1) stored in snstorij[2][inbp][ij]
              double wijkloc=wij/snikfl;                               // (Sn-1) instead of (1-Sn) to account for the 
              double pf=wijkloc*nbsij[ij].snstor[inbp][2];  //  "- sign" coming from the derivation of (1-Sn)
              Vec3d dfrc=pf*nbsij[ij].sigstor[inbp];
              Ftbls.frc2[ij][inbp]=Ftbls.frc2[ij][inbp]+dfrc;
#ifdef VIRIELLOCAL
              double rikp=nbsij[ij].rstor[inbp];
              Ftbls.virielfrc2[ij][inbp]=Ftbls.virielfrc2[ij][inbp]+tensor(rikp*nbsij[ij].sigstor[inbp],dfrc);
#else 
#ifdef VIRIEL
              double rikp=nbsij[ij].rstor[inbp];
              Ftbls.virielfrc2[ij][0]=Ftbls.virielfrc2[ij][0]+tensor(rikp*nbsij[ij].sigstor[inbp],dfrc);
#endif
#endif
              Ftbls.iscnb[inbp]=0;
            }

//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK5c; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  
                  
 //PCC correction
            double wijkpcc=wijk;
            if(inb+1==inbj&&inbj<inbjm){wijkpcc*=Ftbls.prwij[inb+2];}
            else{wijkpcc*=Ftbls.prwij[inb+1];}
            double wijpcc=wijkpcc*snik;
//            if(ij==1&&ip==19)
//            {
//              printf("ij,inb,inbj,nbsij.nnbfr : %6d%6d%6d%6d\n",ij,inb,inbj,nbsij[ij].nnbfr);
//              printf("wijk,wijkpcc,wijpcc : %16.10lf%16.10lf%16.10lf\n",wijk,wijkpcc,wijpcc);
//            }
            for(int k=0;k<2;k++)nchij[k]=0;
            for(int icnb=1;icnb<=2;icnb++) 
            {
              int inbp=lstcnbbr[ibr][icnb];
      	      int ispckp=nbsij[ij].ispc[inbp];
//              uint8_t ispckp=nbsij[ij].ispc[inbp];
              nchij[ispckp]++;
            }
            nchij[ispck]++;
            double pccijcnf=pcc[nchij[ispcc]][nchij[ispch]];
            pccij[ij]+=wijpcc*pccijcnf;

//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK5d; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  

// PCC forces & virial
            double pf0=pccijcnf*vsmra;
            pf=wijkpcc*pf0*nbsij[ij].snstor[inb][2];
            dfrc=pf*nbsij[ij].sigstor[inb];
//            if(kp==iptest||ip==iptest)
//            {
//              printf("1 kp,ip,dfrc+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",kp,ip,dfrc.x,dfrc.y,dfrc.z);
//            }  
            nbs.force[kp]=nbs.force[kp]+dfrc;
            nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
            nbs.viriel[kp]=nbs.viriel[kp]-tensor(0.5*rik*nbsij[ij].sigstor[inb],dfrc);
            nbs.viriel[ip]=nbs.viriel[ip]-tensor(0.5*rik*nbsij[ij].sigstor[inb],dfrc);
#else 
#ifdef VIRIEL
            nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
            
            for(int icnb=nijfull+1;icnb<=2;icnb++) // distinguish between connected and disconnected frac. neigh. using iscnb table
            { 
              int inbp=lstcnbbr[ibr][icnb];
              Ftbls.iscnb[inbp]=1;
            }

//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK5e; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  

            for(int inbp=nbsij[ij].nnbfl;inbp<nbsij[ij].nnbfr;inbp++) // loop over previous inb neigh (already considered for connect/disconnect)
            {
              if(inbp==inbj||inbp==inb)continue;				// skip inb
      	      size_t kpalt=nbsij[ij].nbtbl[inbp];
      	      double snikfl=nbsij[ij].snstor[inbp][Ftbls.iscnb[inbp]];
              Ftbls.iscnb[inbp]=0;
              double wijkpcc=wijpcc/snikfl;
              double pf=wijkpcc*pf0*nbsij[ij].snstor[inbp][2];
              Vec3d dfrc=pf*nbsij[ij].sigstor[inbp];
//              if(kpalt==iptest||ip==iptest)
//              {
//                printf("2 kpalt,ip,dfrc+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",kpalt,ip,dfrc.x,dfrc.y,dfrc.z);
//              }  
              nbs.force[kpalt]=nbs.force[kpalt]+dfrc;
              nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
              double rikp=nbsij[ij].rstor[inbp];
              nbs.viriel[kpalt]=nbs.viriel[kpalt]-tensor(0.5*rikp*nbsij[ij].sigstor[inbp],dfrc);
              nbs.viriel[ip]=nbs.viriel[ip]-tensor(0.5*rikp*nbsij[ij].sigstor[inbp],dfrc);
#else 
#ifdef VIRIEL
              double rikp=nbsij[ij].rstor[inbp];
              nbs.viriel[0]=nbs.viriel[0]-tensor(rikp*nbsij[ij].sigstor[inbp],dfrc);
#endif
#endif
            }
//            if(nbr>6)printf("OKOK3b; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
//            if(ij==1&&ip==188&&inbj==3&&nbr>=6)
//            {
//              printf("OK5; ibr,nbr,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ibr,nbr,inb,lstcnbbr[6][1]);
//            }  
          }
//          if(nbr>6)printf("OKOKOK2b; ij,nbr,lstcnbbr[6][2] : %6d%6d%6d\n",ij,nbr,lstcnbbr[6][2]);
        }
//        if(ij==1&&ip==188&&inbj==3)
//        {
//          printf("nbsij[ij].nnbsrm :%6d\n",nbsij[ij].nnbsrm);
//          printf("ip,inbj,inb,lstcnbbr[6][1] :%6d%6d%6d%6d\n",ip,inbj,inb,lstcnbbr[6][1]);
//        }  
      }
//connected/disconnected tree built

//      if(ij==1&&ip==188&&inbj==3)
//      {
//        printf("ip,inbj,lstcnbbr[6][1] :%6d%6d%6d\n",ip,inbj,lstcnbbr[6][1]);
//      }  
//      for(int k=1;k<=nij;k++)
//      {
//        printf("ij,ind,k,ibr :%6d%6d%6d%6d\n",ij,ind,k,ibr);
//        lstcnb[ij][ind][k]=lstcnbbr[ibr][k];    // change of indices done for lstcnb //
//      }  
      

//until there wij, Mij, and Nelij data are stored in tables indexed along the connection/disconnection process.
//  this is not very convenient since we want to treat the cases with same # of connected neigh. jointly.
//  so we transfer these tables to global tables, arranging the different branches with respect to their number of connected neigh.
      iwk[0]=nnij[ij][0][0];	   // iwk[nij] stores the next destination index for each nij-branch (with nij connected neigh.)
      iwk[1]=nnij[ij][1][0];	   // see previous comment for nnij
      iwk[2]=nnij[ij][2][0];
      for(int ibr=1;ibr<=nbr;ibr++)  // run over all the branches (including ibr=0 for 2 full neigh. case)
//      for(uint16_t ibr=1;ibr<=nbr;ibr++)  // run over all the branches (including ibr=0 for 2 full neigh. case)
      {
        int nij=nijbr[ibr];            // get # of connected neigh. for current branch
        int ind=iwk[nij];	           // update destination index
        wijstor[ij][ind]=wijbr[ibr];       // transfer branch weight
        Nelijstor[ij][ind]=Nelijbr[ibr];   // transfer Nelij and derivatives;  index change for Nelijstor done  //
//        printf("ij,ind,ibr,Nelij :%6d%6d%6d%16.10lf\n",ij,ind,ibr,Nelijbr[ibr]);

        for(int k=1;k<=nij;k++)
        {
//          printf("ij,ip,jp :%6d%6d%6d\n",ij,ip,nbsij[ij].nbtbl[nbsij[ij].inbj]);
//          printf("nnbfl,nnbfr  :%6d%6d\n",nbsij[ij].nnbfl-1,nbsij[ij].nnbfr-nbsij[ij].nnbfl);
          if(lstcnbbr[ibr][k]>=nbsij[ij].nnbsrm)
          {
            printf("ij,ind,k : %6d%6d%6d\n",ij,ind,k);
            printf("ibr,lstcnbbr[ibr][k] : %6d%6d\n",ibr,lstcnbbr[ibr][k]);
          }  
          lstcnb[ij][ind][k]=lstcnbbr[ibr][k];    // change of indices done for lstcnb //
//          printf("ibr,ij,ind,k,lstcnb :%6d%6d%6d%6d%6d\n",ibr,ij,ind,k,lstcnb[ij][ind][k]);
//          printf("ibr,k,lstcnbbr      :%6d%6d%6d\n",ibr,k,lstcnbbr[ibr][k]);
        }  
        for(int k=0;k<2;k++)
        {
//          dNelijdikstor[ij][ind][k]=dNelijbrdik[ibr][k];            //  transfer force contribs.
          dNelijdkistor[ij][ind][k]=dNelijbrdki[ibr][k];
        }
        iwk[nij]++;             //  update destination index for nij-branches
      }

//PCC correction
      for(int ibr=1;ibr<=nbr;ibr++)  //run over all the branches
//      for(uint16_t ibr=1;ibr<=nbr;ibr++)  //run over all the branches
      {
        for(int k=0;k<2;k++)nchij[k]=0;
        int nij=nijbr[ibr];
//compute pcc
        for(int icnb=1;icnb<=nij;icnb++)
        {
          int inb=lstcnbbr[ibr][icnb];
          int ispck=nbsij[ij].ispc[inb];
//          uint8_t ispck=nbsij[ij].ispc[inb];
          nchij[ispck]++;
        }
        double wij=wijbr[ibr];
        double pccijcnf=pcc[nchij[ispcc]][nchij[ispch]];    // is the order of pcc indices okay ? !!!!! //
        pccij[ij]+=wij*pccijcnf;
          
//table of connection flags
        for(int icnb=nijfull+1;icnb<=nij;icnb++)
        {
          int inb=lstcnbbr[ibr][icnb];
          Ftbls.iscnb[inb]=1;
        }
          
//compute local force contribs.
        double pf0=pccijcnf*vsmra;
        for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==inbj){continue;}                     // Is this correct ? !!!!!!! //
          double wijk=wij/nbsij[ij].snstor[inb][Ftbls.iscnb[inb]];       // value stored in snstorij[2][ip][ij]
          double pf=wijk*pf0*nbsij[ij].snstor[inb][2];
          Vec3d dfrc=pf*nbsij[ij].sigstor[inb];
          Ftbls.frc1[0][inb]=Ftbls.frc1[0][inb]+dfrc;
#ifdef VIRIELLOCAL
          double rik=nbsij[ij].rstor[inb];
          Ftbls.virielfrc1[0][inb]=Ftbls.virielfrc1[0][inb]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#else 
#ifdef VIRIEL
          double rik=nbsij[ij].rstor[inb];
          Ftbls.virielfrc1[0][0]=Ftbls.virielfrc1[0][0]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
          Ftbls.iscnb[inb]=0;
        }
      }

//PCC forces (COULD BE MERGED TO THE PREVIOUS PART FOR SIMPLICITY, THOUGH LESS EFFICIENT)
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inb];
//        if(kp==iptest||ip==iptest)
//        {
//          printf("3 kp,ip,dfrc1+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",
//          kp,ip,Ftbls.frc1[0][inb].x,Ftbls.frc1[0][inb].y,Ftbls.frc1[0][inb].z);
//        }  
        nbs.force[kp]=nbs.force[kp]+Ftbls.frc1[0][inb];
        nbs.force[ip]=nbs.force[ip]-Ftbls.frc1[0][inb];
        Ftbls.frc1[0][inb]={0.,0.,0.};
#ifdef VIRIELLOCAL
        nbs.viriel[kp]=nbs.viriel[kp]-0.5*Ftbls.virielfrc1[0][inb];
        nbs.viriel[ip]=nbs.viriel[ip]-0.5*Ftbls.virielfrc1[0][inb];
        Ftbls.virielfrc1[0][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
      }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
      nbs.viriel[0]=nbs.viriel[0]-Ftbls.virielfrc1[0][0];
      Ftbls.virielfrc1[0][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif
    } 
    else // IJP[IJ] IS SATURATED
    {
//only one saturated branch, no contribution of any potential frac. neigh.
      for(int l=0;l<2;l++){for(int k=0;k<3;k++)nnij[ij][k][l]=0;}
      nnij[ij][3][0]=1;
      nnij[ij][3][1]=1;
      wijstor[ij][1]=1.;
        
//PCC correction
      if(nijfull==3)
      {
        double wij=1.;
        for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==inbj){continue;}
          wij*=-nbsij[ij].snstor[inb][0];
        }
        
        for(int k=0;k<2;k++)nchij[k]=0;
//compute pcc
        for(int inb=1;inb<nbsij[ij].nnbfl;inb++)
        {
          if(inb==inbj){continue;}
      	  int ispck=nbsij[ij].ispc[inb];
          nchij[ispck]++;
        }
        double pccijcnf=pcc[nchij[ispcc]][nchij[ispch]];
        pccij[ij]+=wij*pccijcnf;
          
//compute pcc forces
        double pf0=pccijcnf*vsmra;
        for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==inbj){continue;}
          size_t kp=nbsij[ij].nbtbl[inb];
          double wijk=wij/nbsij[ij].snstor[inb][0];
          double pf=wijk*pf0*nbsij[ij].snstor[inb][2];
          Vec3d dfrc=pf*nbsij[ij].sigstor[inb];
//          if(kp==iptest||ip==iptest)
//          {
//            printf("4 kp,ip,dfrc+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",kp,ip,dfrc.x,dfrc.y,dfrc.z);
//          }  
          nbs.force[kp]=nbs.force[kp]+dfrc;
          nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
          double rik=nbsij[ij].rstor[inb];
          nbs.viriel[kp]=nbs.viriel[kp]-tensor(0.5*rik*nbsij[ij].sigstor[inb],dfrc);
          nbs.viriel[ip]=nbs.viriel[ip]-tensor(0.5*rik*nbsij[ij].sigstor[inb],dfrc);
#else 
#ifdef VIRIEL
          double rik=nbsij[ij].rstor[inb];
          nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
        }
      }
    }
  }

//Compute FCONJij
  fcnij=0.;
  int nijfull=std::min(nijfl[0],3);    // if nijfull and/or njifull >3, use index 3 in 'fconj' table.
  int njifull=std::min(nijfl[1],3);    // if nijfull and/or njifull >3, use index 3 in 'fconj' table.
//  double wijsum=0.0;    // for testing only
//  size_t ip=nbsij[0].nbtbl[0];    // for testing only
//  size_t jp=nbsij[0].nbtbl[nbsij[0].inbj];    // for testing only
  for(int nij=nijfull;nij<=3;nij++)               // Double-loop over nij (# of connected neigh of ijp[0])
  {
    if(nnij[0][nij][1]==0)continue;           //  (skip if at. ijp[0] has no branch with nij connected neigh.)
    for(int inij=nnij[0][nij][0];inij<=nnij[0][nij][1];inij++) //  use nnij[ij][nij][0;1] to run over the appropriate portion of the tables
//    for(uint16_t inij=nnij[0][nij][0];inij<=nnij[0][nij][1];inij++) //  use nnij[ij][nij][0;1] to run over the appropriate portion of the tables
    {
      double wij=wijstor[0][inij];
      double Nelij=Nelijstor[0][inij];
//      wijsum+=wij;
      
//local force contribs for ij=0
      frcsFcnj(0,inij,nij,nijfull,wij,dNelijdkistor[0][inij][0],
      dNelijdkistor[0][inij][1],ispc,nbsij,Ftbls,nbs);  // to be checked ? !!!!!! //
      
      for(int nji=njifull;nji<=3;nji++)       // Double-loop over nji (# of connected neigh of ijp[1])
      {
      	if(nnij[1][nji][1]==0)continue;    //  (skip if at. ijp[1] has no branch with nji connected neigh.)
        for(int inji=nnij[1][nji][0];inji<=nnij[1][nji][1];inji++) //  use nnij[ij][nji][0;1] to run over the appropriate portion of the tables
//        for(uint16_t inji=nnij[1][nji][0];inji<=nnij[1][nji][1];inji++) //  use nnij[ij][nji][0;1] to run over the appropriate portion of the tables
        {
          double wji=wijstor[1][inji];
          double wijwji=wij*wji;
          double fcnijcnf=0.0,dFcnjdNelij=0.,dFcnjdNelji=0.;
          
          if(nij==0&&nji==0)             // CASES INDPDT OF Ncj (fncj[0]==fcnj[1])
          {
            fcnijcnf=fcnj[nij][nji][0];
            fcnij+=wijwji*fcnijcnf;
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342)
printf("     fcnj 00 : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,0.,fcnj[nij][nji][0],fcnj[nij][nji][1],wijwji*fcnijcnf,fcnij);
#endif
          }
          else if(nij==3||nji==3)       // CASES INDPDT OF Ncj (fncj[0]==fcnj[1])
          {
            fcnijcnf=fcnj[nij][nji][0];
            fcnij+=wijwji*fcnijcnf;
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342)
printf("     fcnj 33 : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,0.,fcnj[nij][nji][0],fcnj[nij][nji][1],wijwji*fcnijcnf,fcnij);
#endif
          } 
          else if (nij==0) // OTHER CASES WITH NIJ==0
          {
            double Nelji=Nelijstor[1][inji];
            double Nel=Nelij+Nelji;
            double dNcjdNel=1./dxmel2[nij][nji];
            double Ncj=(Nel-xmelmin2[nij][nji])*dNcjdNel;
            double prefNcj=(cf1+(1.-cf1)*Ncj)*Ncj;
            double dprefNcj=cf1+2.*(1.-cf1)*Ncj;
            double fcnj1m0=fcnj[nij][nji][1]-fcnj[nij][nji][0];
            fcnijcnf=fcnj[nij][nji][0]+prefNcj*fcnj1m0;
            fcnij+=wijwji*fcnijcnf;
            dFcnjdNelji=fcnj1m0*dprefNcj*dNcjdNel;
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342)
printf("     fcnj 0x : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,Ncj,fcnj[nij][nji][0],fcnj[nij][nji][1],wijwji*fcnijcnf,fcnij);
//printf("     fcnj 0x : %5zu%5zu%5zu%5zu%18.12lf%18.12lf%18.12lf\n",nij,nji,inij,inji,Nelij,Nelji,Nel);
#endif
          } 
          else if(nji==0)  // OTHER CASES WITH NIJ==0
          {
            double Nelji=4.;
            double Nel=Nelij+Nelji;
            double dNcjdNel=1./dxmel2[nij][nji];
            double Ncj=(Nel-xmelmin2[nij][nji])*dNcjdNel;
            double prefNcj=(cf1+(1.-cf1)*Ncj)*Ncj;
            double dprefNcj=cf1+2.*(1.-cf1)*Ncj;
            double fcnj1m0=fcnj[nij][nji][1]-fcnj[nij][nji][0];
            fcnijcnf=fcnj[nij][nji][0]+prefNcj*fcnj1m0;
            fcnij+=wijwji*fcnijcnf;
            dFcnjdNelij=fcnj1m0*dprefNcj*dNcjdNel;
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342)
printf("     fcnj x0 : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,Ncj,fcnj[nij][nji][0],fcnj[nij][nji][1],wijwji*fcnijcnf,fcnij);
#endif
          }
          else  // DEFAULT
          {
            double Nelji=Nelijstor[1][inji];
            double Nel=Nelij+Nelji;
            double dNcjdNel=1./dxmel2[nij][nji];
            double Ncj=(Nel-xmelmin2[nij][nji])*dNcjdNel;
            double prefNcj=(cf1+(1.-cf1)*Ncj)*Ncj;
            double dprefNcj=cf1+2.*(1.-cf1)*Ncj;
            double fcnj1m0=fcnj[nij][nji][1]-fcnj[nij][nji][0];
            fcnijcnf=fcnj[nij][nji][0]+prefNcj*fcnj1m0;
            dFcnjdNelij=fcnj1m0*dprefNcj*dNcjdNel;
            dFcnjdNelji=dFcnjdNelij;
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342){
printf("     fcnj xx : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,Ncj,fcnj[nij][nji][0],fcnj[nij][nji][1],wijwji*fcnijcnf,fcnij+wijwji*fcnijcnf);
printf("     fcnj xx : %6d%6d%18.12lf\n",nij,nji,Nel);}
#endif

// antibonding
//            if(nij!=nji)
//            {
              double Az=0.,dAzdz1=0.,dAzdz2=0.;
              AdAz(nij,nji,Nelij,Nelji,Az,dAzdz1,dAzdz2);   // to be checked !!!!!!!!double  //
              fcnijcnf+=Az;
              dFcnjdNelij+=dAzdz1;
              dFcnjdNelji+=dAzdz2;
#ifdef VERBOSE
printf("     fcnj az : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,Az,wijwji*Az,fcnij+wijwji*fcnijcnf);
#endif
//            }
          
//torsion
            if(nij==2&&nji==2)
            {
              for(int k=0;k<2;k++)inbcnb[0][k]=lstcnb[0][inij][k+1];
              for(int k=0;k<2;k++)inbcnb[1][k]=lstcnb[1][inji][k+1];
              
//              double Tyz=0.,dTyzdy=0.,dTyzdz=0.;
//              torsion1(ijp[0],ijp[1],Ncj,wijwji,vsmra,rij,sigij,Tyz,dTyzdy,dTyzdz,
//              nbsij,nbs.force);   // to be checked !!!!!! //
              double Tyz=0.,dTyzdysq=0.,dTyzdz=0.;
              torsion2(ijp,Ncj,wijwji,vsmra,rij,sigij,Tyz,dTyzdysq,dTyzdz,
              nbsij,nbs.force,nbs.viriel);
//              std::cin.get();
              fcnijcnf+=Tyz;
              dFcnjdNelij+=dTyzdz*dNcjdNel;  // The force due to y=COSijkl is computed inside the torsion routine
              dFcnjdNelji+=dTyzdz*dNcjdNel;  // -> only the force due to z=Ncj must included in the fcnij
#ifdef VERBOSE_FCNJ
if(ip==8&&jp==342)
printf("     fcnj ty : %5d%5d%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf%18.12lf\n",
ijp[0],ijp[1],wij,wji,wijwji,Tyz,wijwji*Tyz,fcnij+wijwji*fcnijcnf);
#endif
            }
            fcnij+=wijwji*fcnijcnf;
          } 
        
//local force contribs for ij=1
//          printf("inji,nji,njifull : %6d%6d%6d\n",inji,nji,njifull);
          frcsFcnj(1,inji,nji,njifull,wji,dNelijdkistor[1][inji][0],
          dNelijdkistor[1][inji][1],ispc,nbsij,Ftbls,nbs);   // to be checked !!!!!!!! //
        
//global forces, both sides
          frcsFcnjFull(ijp,nij,nji,wij,wji,wijwji,fcnijcnf,vsmra,
          dFcnjdNelij,dFcnjdNelji,nbsij,Ftbls,nbs.force,nbs.viriel);   // to be checked !!!!!!!! //

//Reset local force contributions frc1 & frc3 for ij=1
          for(int inb=nbsij[1].nnbfl;inb<nbsij[1].nnbfr;inb++)
          {
            Ftbls.frc1[1][inb]={0.,0.,0.};
          }
#ifdef VIRIELLOCAL
          for(int inb=nbsij[1].nnbfl;inb<nbsij[1].nnbfr;inb++)
          {
            Ftbls.virielfrc1[1][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
          }
#else 
#ifdef VIRIEL         
          Ftbls.virielfrc1[1][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif            
#endif            

          for(int inb=1;inb<=nfrc3[1];inb++)
          {
            Ftbls.frc3[1][inb]={0.,0.,0.};
          }	
#ifdef VIRIELLOCAL
          for(int inb=1;inb<=nfrc3[1];inb++)
          {
            Ftbls.virielfrc3[1][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
          }	
#else 
#ifdef VIRIEL         
          Ftbls.virielfrc3[1][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif            
#endif            
        }
      }

//Reset local force contributions frc1 & frc3 for ij=0
      for(int inb=nbsij[0].nnbfl;inb<nbsij[0].nnbfr;inb++)
      {
      	Ftbls.frc1[0][inb]={0.,0.,0.};
      }	
#ifdef VIRIELLOCAL
      for(int inb=nbsij[0].nnbfl;inb<nbsij[0].nnbfr;inb++)
      {
        Ftbls.virielfrc1[0][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
      }	
#else 
#ifdef VIRIEL
        Ftbls.virielfrc1[0][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif

      for(int inb=1;inb<=nfrc3[0];inb++)
      {
      	Ftbls.frc3[0][inb]={0.,0.,0.};
      }	
#ifdef VIRIELLOCAL
      for(int inb=1;inb<=nfrc3[0];inb++)
      {
        Ftbls.virielfrc3[0][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
      }	
#else 
#ifdef VIRIEL
      Ftbls.virielfrc3[0][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif
    }
  }

//  if(fabs(wijsum-1.0)>1e-10)
//  {
//    printf("wijsum :%16.10%lf\n",wijsum);
//    std::cin.get();
//  }

//Reset local force and viriel contributions to zero
  for(int ij=0;ij<2;ij++)
  {
    for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
    {
      Ftbls.frc2[ij][inb]={0.,0.,0.};
    }
  }
#ifdef VIRIELLOCAL
  for(int ij=0;ij<2;ij++)
  {
    for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
    {
      Ftbls.virielfrc2[ij][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    }
  }
#else 
#ifdef VIRIEL
  for(int ij=0;ij<2;ij++)
  {
    Ftbls.virielfrc2[ij][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  }
#endif            
#endif            

//printf("   fcnjcc:  %24.12lf\n",fcnij);
//printf("   pccij:   %24.12lf\n",pccij[0]);
//printf("   pccji:   %24.12lf\n",pccij[1]);

#ifdef INFO
/* Fin de la fonction */
FinCPUFonction(NumFonc_SUBFCNJ_CHlchbop,"SUBFCNJ_CHlchbop",tdebutFonction);
#endif

/*  if( VIRFRC2 != nullptr )
  {
    fclose(VIRFRC2);
  }
*/

}

//========================================================================//

void FATP::F1new(int ispcc,int ispch,double Nki,double NCki,
double& dnelij,double* dnelijdnki)
//void FATP::F1new(uint8_t ispcc,uint8_t ispch,double Nki,double NCki,
//double& dnelij,double* dnelijdnki)
{
  if(fabs(Nki-NCki)<1e-12)
  {
    if(Nki>1.0)
    {
      if(Nki<3.0)
      {
        double Nkim3=Nki-3.0;
        double qv1=1.0-p2sm1*Nkim3;
        double qv2=0.50*(Nki-1.0);
        double qv=qv1*qv2;
        double dqvdnc=-p2sm1*qv2+0.50*qv1;
        double dsqv;
        SqdSqup(qv,0.0,1.0,1.0,dnelij,dsqv);
        dnelijdnki[ispcc]=dqvdnc*dsqv;
        dnelijdnki[ispch]=0.0;
      }  
      else
      {
        dnelij=1.0;
        dnelijdnki[0]=0.0;
        dnelijdnki[1]=0.0;
      }
    }  
    else
    {
      dnelij=0.0;
      dnelijdnki[0]=0.0;
      dnelijdnki[1]=0.0;
    }
  }  
  else if(NCki<1e-12)
  {
    double qv=Nki-2.0;
    double dsqv;
    SqdSqup(qv,0.0,1.0,1.0,dnelij,dsqv);
    dnelijdnki[ispcc]=0.0;
    dnelijdnki[ispch]=dsqv;
  }  
  else
  {
    if(Nki<3.0)
    {
      double fcki=NCki/Nki;
      double qv2num=Nki-2.0+fcki;
      if(qv2num>0.0)
      {
        double Nkim3=Nki-3.0;
        double onebyNkisq=1.0/(Nki*Nki);
        double qv1=1.0-p2sm1*fcki*Nkim3;
        double qv2den=1.0+fcki;
        double qv2=qv2num/qv2den;
        double onebyqv2densq=1.0/(qv2den*qv2den);
        double qv=qv1*qv2;
        double dfckidnc=(Nki-NCki)*onebyNkisq;
        double dfckidnh=-NCki*onebyNkisq;
        double dqv1dnc=-p2sm1*(dfckidnc*Nkim3+fcki);
        double dqv1dnh=-p2sm1*(dfckidnh*Nkim3+fcki);
        double dqv2dnc=((1.0+dfckidnc)*qv2den-qv2num*dfckidnc)*onebyqv2densq;
        double dqv2dnh=((1.0+dfckidnh)*qv2den-qv2num*dfckidnh)*onebyqv2densq;
        double dqvdnc=dqv1dnc*qv2+qv1*dqv2dnc;
        double dqvdnh=dqv1dnh*qv2+qv1*dqv2dnh;
        double dsqv;
        SqdSqup(qv,0.0,1.0,1.0,dnelij,dsqv);
        dnelijdnki[ispcc]=dqvdnc*dsqv;
        dnelijdnki[ispch]=dqvdnh*dsqv;
      }
      else
      {
        dnelij=0.0;
        dnelijdnki[0]=0.0;
        dnelijdnki[1]=0.0;
      } 
    }
    else
    {
      dnelij=1.0;
      dnelijdnki[0]=0.0;
      dnelijdnki[1]=0.0;
    }
  }
//  printf("Nki,NCki          :%14.8lf%14.8lf\n",Nki,NCki);
//  printf("dnelij,dnelijdnki :%14.8lf%14.8lf\n",dnelij,dnelijdnki);
} 

//========================================================================

#ifdef F1NEW
void FATP::frcsFcnj(int ij,int inij,int nij,int nijfull,double wij,
double dNelijdki1,double dNelijdki2,vector<uint8_t>& ispc,LNBStr* nbsij,
FATPLoc& Ftbls,NBStr& nbs)
#else
void FATP::frcsFcnj(int ij,int inij,int nij,int nijfull,double wij,
double dNelijdki1c,double dNelijdki2c,vector<uint8_t> ispc,LNBStr* nbsij,
FATPLoc& Ftbls,NBStr& nbs)
#endif
{
// special case nij>=3: force contributions already computed (tables frc2 & virielfrc2)
  if(nij>=3)return;
  
// table of connection flags
  for(int icnb=nijfull+1;icnb<=nij;icnb++)
  {
    int inb=lstcnb[ij][inij][icnb];
//    printf("OK1; ij,inij,icnb,inb : %6d%6d%6d%6d\n",ij,inij,icnb,inb);
    Ftbls.iscnb[inb]=1;
  }
  
// Forces due to weight factor wij (only frac. neigh. contribute)
  for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
  {
    if(inb==nbsij[ij].inbj){continue;}
    double wijk=wij/nbsij[ij].snstor[inb][Ftbls.iscnb[inb]];
    double pf=wijk*nbsij[ij].snstor[inb][2];
    Vec3d dfrc=pf*nbsij[ij].sigstor[inb];
    Ftbls.frc1[ij][inb]=Ftbls.frc1[ij][inb]+dfrc;
#ifdef VIRIELLOCAL
    double rik=nbsij[ij].rstor[inb];
    Ftbls.virielfrc1[ij][inb]=Ftbls.virielfrc1[ij][inb]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#else 
#ifdef VIRIEL
    double rik=nbsij[ij].rstor[inb];
    Ftbls.virielfrc1[ij][0]=Ftbls.virielfrc1[ij][0]+tensor(rik*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
    Ftbls.iscnb[inb]=0;
  }	

//special case nij==0: all force contributions have been computed
  if(nij==0)return;
  
// Forces due to Ncj 
// Nelik & Nelki are treated separately because they depend on distinct variables
  int ifrc3=0; 	// initial number of force contributions to Nelik
//  uint16_t ifrc3=0; 	// initial number of force contributions to Nelik
  nfrc3[ij]=0; 	//
  
// these local tables are convenient to loop over the connected neighbours
#ifdef F1NEW
  double dNelijdkiabs[3],dNelijdkiloc[2][3];
  if(nij==1)
  {
    dNelijdkiloc[0][1]=dNelijdki1;
    dNelijdkiloc[1][1]=dNelijdki2;
    dNelijdkiabs[1]=fabs(dNelijdki1)+fabs(dNelijdki2);
  }  
  else
  {
    dNelijdkiloc[0][1]=dNelijdkiloc[1][1]=dNelijdki1;
    dNelijdkiloc[0][2]=dNelijdkiloc[1][2]=dNelijdki2;
    dNelijdkiabs[1]=fabs(dNelijdki1); // is fabs different from abs ?! //
    dNelijdkiabs[2]=fabs(dNelijdki2);
  }  
#else
  double dNelijdkiloc[3]={0.,dNelijdki1c,dNelijdki2c};
#endif

// loop over all connected neigh. (1 or 2)
  for(int icnb=1;icnb<=nij;icnb++) 
  {
    int inb=lstcnb[ij][inij][icnb];
    if(inb>=nbsij[ij].nnbsrm)
    {
      printf("BUG; nnbsrm : %6d\n",nbsij[ij].nnbsrm);
      printf("ij,inij,icnb,inb : %6d%6d%6d%6d\n",ij,inij,icnb,inb);
      abort();
    }  
    size_t kp=nbsij[ij].nbtbl[inb];
    int nnbfrkp=nbs.nnbfrac[kp]-nbs.nnbfull[kp];
//    if(inb>=nbsij[ij].nnbfl){nnbfrkp=nnbfrkp--;}
    if(nbsij[ij].snstor[inb][1]<1.){nnbfrkp--;}
    
// NELKI dependance (cutoff functions within Nki)
#ifdef F1NEW
    if(nnbfrkp>0&&(dNelijdkiabs[icnb]>1e-12))
#else
    if(nnbfrkp>0&&(fabs(dNelijdkiloc[icnb])>1e-12||fabs(dNelijdkiloc[icnb])>1e-12)) // is fabs different from abs ?! //
#endif
    {
//      int ispck=nbsij[ij].ispc[inb];
//      uint8_t ispck=nbsij[ij].ispc[inb];
      int ifrck=ifrc3+1;           // force index: ifrck is the index for atom kp 
//      uint16_t ifrck=ifrc3+1;           // force index: ifrck is the index for atom kp 
      ifrc3=ifrck;			// then ifrc3 is incremented for all -non-ip- frac. neigh. of kp
      assert(ifrc3!=0);
      assert(ifrc3<nbsij[ij].nnbsrm*nbsij[ij].nnbsrm);
      Ftbls.lstfrc3[ij][ifrck]=kp;		// index change for lstfrc3 done //
//      if(kp==4)
//      {
//        printf("ij,ip,jp,kp,ifrck : %6d%6d%6d%6d%6d\n",ij,nbsij[ij].nbtbl[0],nbsij[1-ij].nbtbl[0],kp,ifrck);
//      }
      for(size_t inbk=nbs.nnbfull[kp];inbk<nbs.nnbfrac[kp];inbk++)
      {
        size_t mp=nbs.nbtbl[inbk];
        if(mp==nbsij[ij].nbtbl[0])continue;
//        if(mp==ijp[ij])continue;
        ifrc3++;
        assert(ifrc3!=0);
        assert(ifrc3<nbsij[ij].nnbsrm*nbsij[ij].nnbsrm);
        Ftbls.lstfrc3[ij][ifrc3]=mp;
//        if(mp==4)
//        {
//          printf("ij,ip,jp,mp,ifrc3 : %6d%6d%6d%6d%6d\n",ij,nbsij[ij].nbtbl[0],nbsij[1-ij].nbtbl[0],mp,ifrc3);
//        }
#ifdef F1NEW
        int ispcm=ispc[mp];
//        int ispcm=nbs.ispc[mp];
//        uint8_t ispcm=nbs.ispc[mp];
        if(fabs(dNelijdkiloc[ispcm][icnb])<1e-12)continue;
        double pf=dNelijdkiloc[ispcm][icnb]*nbs.snstor[inbk][2];
#else
        if(fabs(dNelijdkiloc[icnb])<1e-12)continue;  //NOT ANYMORE ! now everything is stored in the ispcc component - the tables may be simplified
        double pf=dNelijdkiloc[icnb]*nbs.snstor[inbk][2];   // NOT ANYMORE ! now everything is stored in the ispcc component - the tables may be simplified
#endif
        Vec3d dfrc=pf*nbs.sigstor[inbk];
        assert(ifrc3!=0);
        assert(ifrck!=0);
        assert(ifrc3<nbsij[ij].nnbsrm*nbsij[ij].nnbsrm);
        assert(ifrck<nbsij[ij].nnbsrm*nbsij[ij].nnbsrm);
        Ftbls.frc3[ij][ifrc3]=Ftbls.frc3[ij][ifrc3]+dfrc;
        Ftbls.frc3[ij][ifrck]=Ftbls.frc3[ij][ifrck]-dfrc;
#ifdef VIRIELLOCAL
        double rkm=nbs.rstor[inbk];
        Ftbls.virielfrc3[ij][ifrc3]=Ftbls.virielfrc3[ij][ifrc3]-tensor(0.5*rkm*nbs.sigstor[inbk],dfrc);
        Ftbls.virielfrc3[ij][ifrck]=Ftbls.virielfrc3[ij][ifrck]-tensor(0.5*rkm*nbs.sigstor[inbk],dfrc);
#else 
#ifdef VIRIEL
        double rkm=nbs.rstor[inbk];
//        printf("kp,mp,icnb,ispcm    :%6d%6d%6d%6d\n",kp,mp,icnb,ispcm);
//        printf("rkm,dNelijdkiloc,pf :%16.10lf%16.10lf%16.10lf\n",rkm,dNelijdkiloc[ispcm][icnb],pf);
//        printf("dfrc                :%16.10lf%16.10lf%16.10lf\n",dfrc.x,dfrc.y,dfrc.z);
//        printf("sigstor             :%16.10lf%16.10lf%16.10lf\n",nbs.sigstor[inbk].x,nbs.sigstor[inbk].y,nbs.sigstor[inbk].z);
//        std::cin.get();
        Ftbls.virielfrc3[ij][0]=Ftbls.virielfrc3[ij][0]-tensor(rkm*nbs.sigstor[inbk],dfrc);
#endif
#endif
      }
    }
  }
  nfrc3[ij]=ifrc3; // total number of Nelik force contributions for ijp[ij]
}

//------------------------------------------------------------------------//
void FATP::frcsFcnjFull(size_t ijp[2],int nij,int nji,double wij,double wji,
double wijwji,double fcnijcnf,double vsmra,double dFcnj1, double dFcnj2,
LNBStr* nbsij,FATPLoc& Ftbls,vector<Vec3d>& force,vector<Mat3d>& viriel)
{
//  int ij;
//  int k,l;
//  int ip,kp;
//  int ifrc3;
//  int inb;
//  int nijloc[2];
//  
//  double pf,pf0,pf1,hpf;
//  double dfrc[3];
//  double wijloc[2],dFcnjloc[2];

//  size_t iptest=-1; 
  
  double pf0=fcnijcnf*vsmra;
  
// local tables used to loop over ij-sides
  int nijloc[2]={nij,nji};
  double wijloc[2]={wij,wji};
  double dFcnjloc[2]={dFcnj1,dFcnj2};

//loop over ij sides
  for(int ij=0;ij<2;ij++)
  {
    size_t ip=ijp[ij];
    if(nijloc[ij]<3)
    {
//contribs from wij/wji
      double pf=wijloc[1-ij]*pf0;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==nbsij[ij].inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inb];
        Vec3d dfrc=pf*Ftbls.frc1[ij][inb];
//        if(kp==iptest||ip==iptest)
//        {
//          printf("5 kp,ip,dfrc+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",kp,ip,dfrc.x,dfrc.y,dfrc.z);
//        }  
        force[kp]=force[kp]+dfrc;
        force[ip]=force[ip]-dfrc;
#ifdef VIRIELLOCAL
        double hpf=0.5*pf;
        viriel[kp]=viriel[kp]-hpf*Ftbls.virielfrc1[ij][inb];
        viriel[ip]=viriel[ip]-hpf*Ftbls.virielfrc1[ij][inb];
#endif
      }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
      viriel[0]=viriel[0]-pf*Ftbls.virielfrc1[ij][0];
#endif
#endif
    
//contribs from Nelij/Nelji
      if(fabs(dFcnjloc[ij])>1e-12)
      {
        double pf1=wijwji*dFcnjloc[ij]*vsmra;
        for(int ifrc3=1;ifrc3<=nfrc3[ij];ifrc3++)
//        for(uint16_t ifrc3=1;ifrc3<=nfrc3[ij];ifrc3++)
        {
          size_t kp=Ftbls.lstfrc3[ij][ifrc3];
//          if(kp==iptest)
//          {
//            printf("6 kp,ifrc3,nfrc3 : %6d%6d%6d\n",kp,ifrc3,nfrc3[ij]);
//            printf("6 kp,frc3+ : %6d%18.12lf%18.12lf%18.12lf\n",
//            kp,Ftbls.frc3[ij][ifrc3].x,Ftbls.frc3[ij][ifrc3].y,Ftbls.frc3[ij][ifrc3].z);
//          }  
          force[kp]=force[kp]+pf1*Ftbls.frc3[ij][ifrc3];
#ifdef VIRIELLOCAL
//          Mat3d virloc=pf1*Ftbls.virielfrc3[ij][ifrc3];
//          printf("kp,pf1 :%6d%16.10lf\n",kp,pf1);
//          printf("%16.10lf%16.10lf%16.10lf\n",virloc.m11,virloc.m21,virloc.m31);
//          printf("%16.10lf%16.10lf%16.10lf\n",virloc.m12,virloc.m22,virloc.m32);
//          printf("%16.10lf%16.10lf%16.10lf\n",virloc.m13,virloc.m23,virloc.m33);
//          std::cin.get();
          viriel[kp]=viriel[kp]+pf1*Ftbls.virielfrc3[ij][ifrc3];
#endif
        }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
        viriel[0]=viriel[0]+pf1*Ftbls.virielfrc3[ij][0];
#endif
#endif
      }
    }
    else
    {	
//only contribs from wij/wji for nij/nji=3
      double pf=wijloc[1-ij]*pf0;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==nbsij[ij].inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inb];
      	Vec3d dfrc=pf*Ftbls.frc2[ij][inb];
//        if(kp==iptest||ip==iptest)
//        {
//          printf("7 kp,ip,dfrc+- : %6d%6d%18.12lf%18.12lf%18.12lf\n",kp,ip,dfrc.x,dfrc.y,dfrc.z);
//        }  
      	force[kp]=force[kp]+dfrc;
      	force[ip]=force[ip]-dfrc;
#ifdef VIRIELLOCAL
        double hpf=0.5*pf;
        viriel[kp]=viriel[kp]-hpf*Ftbls.virielfrc2[ij][inb];
        viriel[ip]=viriel[ip]-hpf*Ftbls.virielfrc2[ij][inb];
#endif
      }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
//      if(ij==0&&ijp[0]==4&&ijp[1]==86)
//      {
//        printf("ij,ip,jp,wij,wji,pf : %6d%6d%6d%16.10lf%16.10lf%16.10lf\n",ij,ijp[0],ijp[1],wijloc[ij],wijloc[1-ij],pf);
//        std::cin.get();
//      }
      viriel[0]=viriel[0]-pf*Ftbls.virielfrc2[ij][0];
#endif
#endif
    }
  }
}

//------------------------------------------------------------------------//
//
//void felmin(double var1,double var2,double& varmin,double& dvarmindvar1,double& dvarmindvar2)
//{
//  double del,dvar,varloc;
//  
//  del=1./6.;
//  dvar=var1-var2;
//  varloc=dvar/del;
//  
//  if(varloc>0.&&varloc<1.)
//  {
//    varmin=var1-varloc*(2.-varloc)*dvar;
//    dvarmindvar2=varloc*(4.-3.*varloc);
//    dvarmindvar1=1.-dvarmindvar2;
//  } 
//  else if(dvar<=0.)
//  {
//    varmin=var1;
//    dvarmindvar1=1.;
//    dvarmindvar2=0.;
//  }
//  else
//  {
//    varmin=var2;
//    dvarmindvar1=0.;
//    dvarmindvar2=1.;
//  } 
//}
//
//------------------------------------------------------------------------//
//
//void felmax(double var1,double var2,double& varmin,double& dvarmindvar1,double& dvarmindvar2)
//{
//  double del,dvar,varloc;
//  del=1./6.;
//  dvar=var2-var1;
//  varloc=dvar/del;
//  
//  if(varloc>0. && varloc<1.)
//  {
//    varmin=var1+varloc*(2.-varloc)*dvar;
//    dvarmindvar2=varloc*(4.-3.*varloc);
//    dvarmindvar1=1.-dvarmindvar2;
//  } 
//  else if(dvar<=0.)
//  {
//    varmin=var1;
//    dvarmindvar1=1.;
//    dvarmindvar2=0.;
//  } 
//  else 
//  {
//    varmin=var2;
//    dvarmindvar1=0.;
//    dvarmindvar2=1.;
//  } 
//}
//
//------------------------------------------------------------------------//

void FATP::subPCH(int ij,int ispcc,int ispch,int ifracnb,int ipmrnb,
double vsmra,double& pchij,LNBStr* nbsij,FATPLoc& Ftbls,NBStr& nbs)
//void FATP::subPCH(int ij,uint8_t ispcc,uint8_t ispch,int ifracnb,
//int ipmrnb,double vsmra,double& pchij,LNBStr* nbsij,FATPLoc& Ftbls,NBStr& nbs)
{
  int nchij[2];
//  int nijfl[2],nchij[2];
  Vec3d dfrc;

//#ifdef INFO
///* Initialisation */
//double tdebutFonction;
//InitialisationCPUFonction(NumFonc_SUBPCHlchbop,"SUBPCHlchbop",&tdebutFonction);
//#endif
//
////{/*TESTFORCES*/pchij=0.;printf("PCH desactive\n");return;}/*CANCELS OUT PCH*/
//
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
////CC Comments:
////CC   + for ij=0, ip side is calculated 
////CC   + for ij=1, jp side is calculated 
////CC   + nijfull = number full neighbors of i[j] not equal to j[i]
////CC   + nijfrac = number fractional neighbors of i[j] not equal to j[i]
////CC   + ncnb		= number of connected neighbors
////CC   + lstcnb	= lst of neighbor indices of conneted neighbors
////CC   + nij		 = number connected neighbors k to atom i
////CC   + nji		 = number connected neighbors l to atom j
////CC   + nbr		 = number of connection possibilities/branches
////CC  	- E.g., if nijfull=2 and nijfrac=3, there are 8 possible primary 
////CC      branches of connectivity configurations: [11000],[11100],[11010],
////CC      [11001],[11110],[11101],[11011],[11111]
////CC   	  However, all connectivity configurations with 3 full neighbours 
////CC      can be taken together, as for nij equal or larger then 3 the 
////CC      conjugation number does not depend on coordinations of 
////CC   	  the connected neighbours, so that there are only 2 relevant 
////        connetivity configurations, namely: [11000] and 
////        [11100]+[11010]+[11001]+[11110]+[11101]+[11011]+[11111]
////CC   + nchij[ispci] = number of type-ispci-neighbours of ip other than jp
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
//
////for(j=0;j<2;j++){for(k=0;k<nnbm;k++){for(l=0;l<3;l++){for(m=0;m<3;m++)Ftbls.virielfrc1[m][l][k][j]=0.;}}}
  
  size_t ip=nbsij[ij].nbtbl[0];
  int inbj=nbsij[ij].inbj;
//  size_t inbj=nbsij[ij].inbj;
  int nijfull=nbsij[ij].nnbfl-2+ifracnb;
//  int nijfrac=nbsij[ij].nnbfr-nbsij[ij].nnbfl-ifracnb+ipmrnb;  // are ifracnb and imrnb defined ?!!!!!!! //
//  nijfl[ij]=nijfull;
  pchij=0.;

  if(nijfull<4)
  {
//precompute prwij ( (1.-snik)-product for neigh. above inb - careful, snstor2=(snik-1.) )
    double prwijmin=1.;
    for(int inb=nbsij[ij].nnbfr-1;inb>=nbsij[ij].nnbfl;inb--)
    {
      if(inb==inbj)continue;                   // inbj has to be initialized !!!!! //
      prwijmin*=-nbsij[ij].snstor[inb][0];
      Ftbls.prwij[inb]=prwijmin;
    }
//    Ftbls.prwij[nbsij[ij].nnbfr+1]=1.;  // what is going on here ?????
    Ftbls.prwij[nbsij[ij].nnbfr]=1.;  // wouldn't this be better ?????
    
    if(nijfull==3)
    {
      double wij=prwijmin;
      for(int k=0;k<2;k++)nchij[k]=0;
      for(int inb=1;inb<nbsij[ij].nnbfl;inb++)
      {
        if(inb==inbj)continue;                   // inbj has to be initialized !!!!! //
      	int ispcj=nbsij[ij].ispc[inb];
//        uint8_t ispcj=nbsij[ij].ispc[inb];
      	nchij[ispcj]++;
      }
      double pchijcnf=pch[nchij[ispcc]][nchij[ispch]];
      pchij+=wij*pchijcnf;
      double pf0=pchijcnf*vsmra;
      for(int inbp=nbsij[ij].nnbfl;inbp<nbsij[ij].nnbfr;inbp++)
      {
        if(inbp==nbsij[ij].inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inbp];
      	double wijk=wij/nbsij[ij].snstor[inbp][0];
      	double pf=wijk*pf0*nbsij[ij].snstor[inbp][2];
      	dfrc=pf*nbsij[ij].sigstor[inbp];
      	nbs.force[kp]=nbs.force[kp]+dfrc;
      	nbs.force[ip]=nbs.force[ip]-dfrc;

#ifdef VIRIELLOCAL
       double hrik=0.5*nbsij[ij].rstor[inbp];
       nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);	
       nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);
#else  
#ifdef VIRIEL
       double rik=nbsij[ij].rstor[inbp];
       nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inbp],dfrc);	
#endif
#endif
      }
    } 
    else
    {
      int inbjm=nbsij[ij].nnbfr;
      int ncnb=0;
      for(int inb=1;inb<nbsij[ij].nnbfl;inb++)
      {
        if(inb==inbj)continue;
        lstcnbbr[1][++ncnb]=inb;
      }
      
      int nbr=1;
//      uint16_t nbr=1;
      nijbr[1]=nijfull;
      wijbr[1]=1.;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
        double snik=nbsij[ij].snstor[inb][1];
        int nbrold=nbr;
//        uint16_t nbrold=nbr;
        for(int ibr=1;ibr<=nbrold;ibr++)
//        for(uint16_t ibr=1;ibr<=nbrold;ibr++)
        {
          double wij=wijbr[ibr];
          wijbr[ibr]=wij*(1.-snik);
          int nij=nijbr[ibr];
          int nijp1=nij+1;
          if(nijp1<3)
          {
            nbr++;
            nijbr[nbr]=nijp1;
            for(int k=1;k<=nij;k++)lstcnbbr[nbr][k]=lstcnbbr[ibr][k];
            lstcnbbr[nbr][nijp1]=inb;
            wijbr[nbr]=wij*snik;
          }
          else
          {
            double wijk=wij;
            if(inb+1==inbj&&inbj<inbjm){wijk*=Ftbls.prwij[inb+2];}
            else{wijk*=Ftbls.prwij[inb+1];}
//            double wijk=wij*Ftbls.prwij[inb+1];
            wij=wijk*snik;
            size_t kp=nbsij[ij].nbtbl[inb];
            for(int k=0;k<2;k++)nchij[k]=0;
            for(int icnb=1;icnb<=2;icnb++)
            {
              int inbp=lstcnbbr[ibr][icnb];
      	      int ispcj=nbsij[ij].ispc[inbp];
//              uint8_t ispcj=nbsij[ij].ispc[inbp];
              nchij[ispcj]++;
            }
      	    int ispcj=nbsij[ij].ispc[inb];
//            uint8_t ispcj=nbsij[ij].ispc[inb];
            nchij[ispcj]++;
            double pchijcnf=pch[nchij[ispcc]][nchij[ispch]];
            pchij+=wij*pchijcnf;
            double pf0=pchijcnf*vsmra;
            double pf=wijk*pf0*nbsij[ij].snstor[inb][2];
            dfrc=pf*nbsij[ij].sigstor[inb];
            nbs.force[kp]=nbs.force[kp]+dfrc;
            nbs.force[ip]=nbs.force[ip]-dfrc;
// sign of accumulation in virial term MUST match the sign of accumulation in "ip"  force => loosy explanation, dig more

#ifdef VIRIELLOCAL
            double hrik=0.5*nbsij[ij].rstor[inb];
            nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);	
            nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
            double rik=nbsij[ij].rstor[inb];
            nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inb],dfrc);	
#endif
#endif

            for(int icnb=nijfull+1;icnb<=2;icnb++)
            {
              int inbp=lstcnbbr[ibr][icnb];
              Ftbls.iscnb[inbp]=1;
            }
            for(int inbp=nbsij[ij].nnbfl;inbp<inb;inbp++)
            {
              if(inbp==nbsij[ij].inbj){continue;}
              size_t kp=nbsij[ij].nbtbl[inbp];
              wijk=wij/nbsij[ij].snstor[inbp][Ftbls.iscnb[inbp]];
              pf=wijk*pf0*nbsij[ij].snstor[inbp][2];
              dfrc=pf*nbsij[ij].sigstor[inbp];
              nbs.force[kp]=nbs.force[kp]+dfrc;
              nbs.force[ip]=nbs.force[ip]-dfrc;

#ifdef VIRIELLOCAL
              double hrik=0.5*nbsij[ij].rstor[inbp];
              nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);	
              nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);
#else
#ifdef VIRIEL
              double rik=nbsij[ij].rstor[inbp];
              nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inbp],dfrc);	
#endif
#endif
              Ftbls.iscnb[inbp]=0;
            }
            for(int inbp=inb+1;inbp<nbsij[ij].nnbfr;inbp++)
            {
              if(inbp==nbsij[ij].inbj){continue;}
      	      size_t kp=nbsij[ij].nbtbl[inbp];
              wijk=wij/nbsij[ij].snstor[inbp][0];
              pf=wijk*pf0*nbsij[ij].snstor[inbp][2];
              dfrc=pf*nbsij[ij].sigstor[inbp];
              nbs.force[kp]=nbs.force[kp]+dfrc;
              nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
              double hrik=0.5*nbsij[ij].rstor[inbp];
              nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);	
              nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inbp],dfrc);
#else
#ifdef VIRIEL
              double rik=nbsij[ij].rstor[inbp];
              nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inbp],dfrc);	
#endif
#endif
            }
          }
        }
      }

      for(int ibr=1;ibr<=nbr;ibr++)
//      for(uint16_t ibr=1;ibr<=nbr;ibr++)
      {
        int nij=nijbr[ibr];
        double wij=wijbr[ibr];
        for(int k=0;k<2;k++)nchij[k]=0;
        for(int icnb=1;icnb<=nij;icnb++)
        {
      	  int inb=lstcnbbr[ibr][icnb];
      	  int ispcj=nbsij[ij].ispc[inb];
//          uint8_t ispcj=nbsij[ij].ispc[inb];
      	  nchij[ispcj]++;
      	}
      	double pchijcnf=pch[nchij[ispcc]][nchij[ispch]];
      	pchij+=wij*pchijcnf;
      	double pf0=pchijcnf*vsmra;
      	for(int icnb=nijfull+1;icnb<=nij;icnb++)
        {
      	  int inb=lstcnbbr[ibr][icnb];
      	  Ftbls.iscnb[inb]=1;
      	}
        for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==nbsij[ij].inbj){continue;}
          double wijk=wij/nbsij[ij].snstor[inb][Ftbls.iscnb[inb]];
          double pf=wijk*pf0*nbsij[ij].snstor[inb][2];
          dfrc=pf*nbsij[ij].sigstor[inb];               // TAKE CARE
          Ftbls.frc1[0][inb]=Ftbls.frc1[0][inb]+dfrc;   // TAKE CARE
//          Ftbls.frc1[0][inb]=Ftbls.frc1[0][inb]+pf*nbsij[ij].sigstor[inb];   // TAKE CARE
#ifdef VIRIELLOCAL
          double rinb=nbsij[ij].rstor[inb];
          Ftbls.virielfrc1[0][inb]=Ftbls.virielfrc1[0][inb]+tensor(rinb*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
          double rinb=nbsij[ij].rstor[inb];
          Ftbls.virielfrc1[0][0]=Ftbls.virielfrc1[0][0]+tensor(rinb*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
          Ftbls.iscnb[inb]=0;
        }
      }
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==nbsij[ij].inbj){continue;}
      	size_t kp=nbsij[ij].nbtbl[inb];
        nbs.force[kp]=nbs.force[kp]+Ftbls.frc1[0][inb];
        nbs.force[ip]=nbs.force[ip]-Ftbls.frc1[0][inb];
        Ftbls.frc1[0][inb]={0.,0.,0.};
// sign of accumulation in virial term MUST match the sign of accumulation in "ip"  force => loosy explanation, dig more
#ifdef VIRIELLOCAL
        nbs.viriel[kp]=nbs.viriel[kp]-0.5*Ftbls.virielfrc1[0][inb];
        nbs.viriel[ip]=nbs.viriel[ip]-0.5*Ftbls.virielfrc1[0][inb];
        Ftbls.virielfrc1[0][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
        Ftbls.iscnb[inb]=0;
      }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
        nbs.viriel[0]=nbs.viriel[0]-Ftbls.virielfrc1[0][0];
        Ftbls.virielfrc1[0][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif
    }
  }

//printf("   pch:     %24.12lf\n",pchij);
//
#ifdef INFO
/* Fin de la fonction */
FinCPUFonction(NumFonc_SUBPCHlchbop,"SUBPCHlchbop",tdebutFonction);
#endif

} //end function subPCH

//-----------------------------------------------------------------------//

void FATP::subPHC(int ij,int ifracnb,double vsmra,double& phcij,
LNBStr* nbsij,NBStr& nbs)
{
  Vec3d dfrc;
  double phcijloc=0.;

//  int k,l,m;
//  int ip,jp,kp;
//  
//  double phcijloc;
//  double pf,pf0;
//  double dfrc[3];
//
//#ifdef INFO
///* Initialisation */
//double tdebutFonction;
//InitialisationCPUFonction(NumFonc_SUBPHClchbop,"SUBPHClchbop",&tdebutFonction);
//#endif
//
////{/*TESTFORCES*/phcij=0.;printf("PHC desactive\n");return;}/*CANCELS OUT PHC*/
//
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
////CC Comments:
////CC  + for ij=0, ip side is calculated	 
////CC  + for ij=1, jp side is calculated		 
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
//
////printf("phc: %d %d\n",nnbij[ij],nnbsrij[ij]);

  if(nbsij[ij].nnbfr>1)
  {
//    if(nbsij[ij].nnbfl>2)
//    {
//      double phcijloc=0.;
////      double pf=0.;
//    }
//    else
    if(nbsij[ij].nnbfl-2+ifracnb<1)
    {
      size_t ip=nbsij[ij].nbtbl[0];
      int inbj=nbsij[ij].inbj;
//      size_t inbj=nbsij[ij].inbj;
      phcijloc=phc0;
      for(int inb=1;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
        phcijloc*=-nbsij[ij].snstor[inb][0];       // Take Care !!!!!!! //
//        size_t kp=nbsij[ij].nbtbl[inb];
//        if(kp!=jp){phcijloc*=(1.-nbsij[ij].snstor[inb][1]);}
      }
    
      double pf0=phcijloc*vsmra;
      double pf;
      for(int inb=1;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
        size_t kp=nbsij[ij].nbtbl[inb];
        {
          if(nbsij[ij].snstor[inb][0]<-1e-12)
          {pf=pf0/nbsij[ij].snstor[inb][0]*nbsij[ij].snstor[inb][2];}else{pf=0.;}
//          {pf=pf0/(nbsij[ij].snstor[inb][0]-1.)*nbsij[ij].snstor[inb][1];}else{pf=0.;}
          dfrc=pf*nbsij[ij].sigstor[inb];
          nbs.force[kp]=nbs.force[kp]+dfrc;
          nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
          double hrik=0.5*nbsij[ij].rstor[inb];
          nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);	
          nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
          double rik=nbsij[ij].rstor[inb];
          nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inb],dfrc);	
#endif
#endif
        }
      }
    }
//    else
//    {
//      phcijloc=0.;
////      double pf=0.;
//    }
  }
  else
  {
    phcijloc=phc0;
//    double pf=0.;
  }

  phcij=phcijloc;

//  if(nbsij[ij].nbtbl[0]==92&&nbsij[1-ij].nbtbl[0]==42)
//  {
//    printf("ip,jp,phcij : %6zu%6zu%16.10lf\n",nbsij[ij].nbtbl[0],nbsij[1-ij].nbtbl[0],phcij);
//    printf("nnbfl,nnbfr,ifracnb : %6zu%6zu%6d\n",nbsij[ij].nnbfl,nbsij[ij].nnbfr,ifracnb);
//  }
//  std::cin.get();

//#ifdef INFO
///* Fin de la fonction */
//FinCPUFonction(NumFonc_SUBPHClchbop,"SUBPHClchbop",tdebutFonction);
//#endif
} // end function subPHC

//-----------------------------------------------------------------------//

void FATP::subPHH(int ifracnb,double vsmra,double& phhij,
LNBStr* nbsij,NBStr& nbs)
{
  Vec3d dfrc;

//  int ij;
//  int k,l,m;
//  int ip,jp,kp;
//  int nijfull;
//  int flag;
//  
//  double phhijloc;
//  double pf,pf0;
//  double dfrc[3];
//
//#ifdef INFO
///* Initialisation */
//double tdebutFonction;
//InitialisationCPUFonction(NumFonc_SUBPHHlchbop,"SUBPHHlchbop",&tdebutFonction);
//#endif
//
////{/*TESTFORCES*/phhij=0.;printf("PHH desactive\n");return;}/*CANCELS OUT PHH*/
//
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
////CC Promotion function PHH for hydrogen
////CC  - Only non-zero and positive if atom i, which is an H atom, 
////CC    has no full neighbors among atoms not equal to j  
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//
////Initialisation
  phhij=0.;
  
  for(int ij=0;ij<2;ij++)
  {
    double phhijloc=0.5*pHH0;
    if(nbsij[ij].nnbfr>2)
    {
      int nijfull=nbsij[ij].nnbfl-2+ifracnb;
//      if(nijfull>0)
//      {
//        uint8_t flag=1;
////        double phhijloc=0.;
//      }
//      else
      if(nijfull==0)
      {
//        flag=2;
        size_t ip=nbsij[ij].nbtbl[0];
//        size_t jp=nbsij[1-ij].nbtbl[0];
        int inbj=nbsij[ij].inbj;
        for(int inb=1;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==inbj){continue;}
//          size_t kp=nbsij[ij].nbtbl[inb];
          assert(nbsij[ij].nbtbl[inb]!=nbsij[1-ij].nbtbl[0]);
          phhijloc*=-nbsij[ij].snstor[inb][0];
//          phhijloc*=(1.-nbsij[ij].snstor[inb][0]);
//          if(kp!=jp){phhijloc*=(1.-nbsij[ij].snstor[inb][0]);}
        }
        
        double pf0=phhijloc*vsmra;
        double pf;
        for(int inb=1;inb<nbsij[ij].nnbfr;inb++)
        {
          if(inb==inbj){continue;}
          size_t kp=nbsij[ij].nbtbl[inb];
          assert(kp!=nbsij[1-ij].nbtbl[0]);
          if(nbsij[ij].snstor[inb][0]<-1e-12)
          {pf=pf0/nbsij[ij].snstor[inb][0]*nbsij[ij].snstor[inb][2];}else{pf=0.;}
          dfrc=pf*nbsij[ij].sigstor[inb];
          nbs.force[kp]=nbs.force[kp]+dfrc;
          nbs.force[ip]=nbs.force[ip]-dfrc;
#ifdef VIRIELLOCAL
          double hrik=0.5*nbsij[ij].rstor[inb];
          nbs.viriel[kp]=nbs.viriel[kp]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);	
          nbs.viriel[ip]=nbs.viriel[ip]-tensor(hrik*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
          double rik=nbsij[ij].rstor[inb];
          nbs.viriel[0]=nbs.viriel[0]-tensor(rik*nbsij[ij].sigstor[inb],dfrc);	
#endif
#endif
        }
      }
      else{phhijloc=0.;}
    }
    phhij+=phhijloc;
  }

////printf("   phh:     %24.12lf\n",phhij);
//
//#ifdef INFO
///* Fin de la fonction */
//FinCPUFonction(NumFonc_SUBPHHlchbop,"SUBPHHlchbop",tdebutFonction);
//#endif

}

//-----------------------------------------------------------------------//

void FATP::subFcnj_HH(int ifracnb,double vsmra,double& fcnij,
LNBStr* nbsij,FATPLoc& Ftbls,NBStr& nbs)
{
////  int flag;
//
//#ifdef INFO
///* Initialisation */
//double tdebutFonction;
//InitialisationCPUFonction(NumFonc_SUBFCNJ_HHlchbop,"SUBFCNJ_HHlchbop",&tdebutFonction);
//#endif
//
////{/*TESTFORCES*/fcnij=0.;return;}
//    
//// HH-energy promotion only if ip & jp have no other full neighbour than 
//// each other (matches increased stability of isolated H2 moelcule)

  if(nbsij[0].nnbfl-2+ifracnb>0||nbsij[1].nnbfl-2+ifracnb>0)
  {
    fcnij=0.;
  }
  else
  {
    fcnij=FcnjH2;
    
    for(int ij=0;ij<2;ij++)
    {
      int inbj=nbsij[ij].inbj;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
        double onemsn=-nbsij[ij].snstor[inb][0];
//        double onemsn=(1.-nbsij[ij].snstor[inb][0]);
        fcnij*=onemsn;
        double pf=nbsij[ij].snstor[inb][2]/onemsn;
        Vec3d dfrc=pf*nbsij[ij].sigstor[inb];   // TAKE CARE
        Ftbls.frc1[ij][inb]=Ftbls.frc1[ij][inb]-dfrc;   // TAKE CARE
//        Ftbls.frc1[ij][inb]=Ftbls.frc1[ij][inb]-pf*nbsij[ij].sigstor[inb];   // TAKE CARE
#ifdef VIRIELLOCAL
        double rinb=nbsij[ij].rstor[inb];
        Ftbls.virielfrc1[ij][inb]=Ftbls.virielfrc1[ij][inb]-tensor(rinb*nbsij[ij].sigstor[inb],dfrc);
#else
#ifdef VIRIEL
        double rinb=nbsij[ij].rstor[inb];
        Ftbls.virielfrc1[ij][0]=Ftbls.virielfrc1[ij][0]-tensor(rinb*nbsij[ij].sigstor[inb],dfrc);
#endif
#endif
      }	
    }	
    
    double pf=fcnij*vsmra;
    for(int ij=0;ij<2;ij++)
    {
      int inbj=nbsij[ij].inbj;
      for(int inb=nbsij[ij].nnbfl;inb<nbsij[ij].nnbfr;inb++)
      {
        if(inb==inbj){continue;}
        size_t ip=nbsij[ij].nbtbl[0];
        size_t lp=nbsij[ij].nbtbl[inb];
        Vec3d dfrc=pf*Ftbls.frc1[ij][inb];
        nbs.force[ip]=nbs.force[ip]-dfrc;
        nbs.force[lp]=nbs.force[lp]+dfrc;
        Ftbls.frc1[ij][inb]={0.,0.,0.};
#ifdef VIRIELLOCAL
        double hpf=0.5*pf;
        nbs.viriel[ip]=nbs.viriel[ip]-hpf*Ftbls.virielfrc1[ij][inb];
        nbs.viriel[lp]=nbs.viriel[lp]-hpf*Ftbls.virielfrc1[ij][inb];
        Ftbls.virielfrc1[ij][inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
      }	
#ifdef VIRIEL
#ifndef VIRIELLOCAL
        nbs.viriel[0]=nbs.viriel[0]-pf*Ftbls.virielfrc1[ij][0];
        Ftbls.virielfrc1[ij][0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif
    }	
  }

////printf("   fcnjhh:  %24.12lf\n",fcnij);
//
//#ifdef INFO
///* Fin de la fonction */
//FinCPUFonction(NumFonc_SUBFCNJ_HHlchbop,"SUBFCNJ_HHlchbop",tdebutFonction);
//#endif

}

//-----------------------------------------------------------------------//

void FATP::AdAz(int nij,int nji,double xNelij,double xNelji,
double& Az,double& dAzdz1,double& dAzdz2)
{
  double a2,deltael,onebyoneptendeltabs,dAzddelta;
  
#ifdef INFO
/* Initialisation */
double tdebutFonction;
InitialisationCPUFonction(NumFonc_AZDAZ_Clchbop,"AZDAZ_Clchbop",&tdebutFonction);
#endif
  
//{/*TESTFORCES*/*Az=0.;*dAzdz1=0.;*dAzdz2=0.;printf("AZ desactive\n");return;}/*CANCELS OUT AZ CONTRIBS*/
  
  a2=aab2[nij][nji];
  deltael=xNelij-xNelji;
  onebyoneptendeltabs=1./(1.+aab1*fabs(deltael));

//  printf("nij,nji,a2,deltael  :%6d,%6d%16.10lf%16.10lf\n",nij,nji,a2,deltael);
  
  Az=a2*deltael*deltael*onebyoneptendeltabs;
  dAzddelta=a2*deltael*onebyoneptendeltabs*(1.+onebyoneptendeltabs);
  dAzdz1=dAzddelta;
  dAzdz2=-dAzddelta;

/* Fin de la fonction */
#ifdef INFO
FinCPUFonction(NumFonc_AZDAZ_Clchbop,"AZDAZ_Clchbop",tdebutFonction);
#endif

}

//-----------------------------------------------------------------------//

void FATP::torsion1(size_t ip,size_t jp,double Ncj,double wijwji,double vsmra,
double rij,Vec3d sigij,double& Tyz,double& dTyzdy,double& dTyzdz,
LNBStr* nbsij,vector<Vec3d>& force)
{
//  int ij;
//  int k,l,kp,lp;
//  int inb,icnb,ixyz;
  
//  double sqrt3=1.732050807569;
//  double rik1,rik2;
//  double wminabs,pf,pf0,pf00,pf01,pf02,oneppf,onempf,dpf,hpf,
//  double sigijwmin,dsigijwmin,sigijdwmin,wmindwmin,wminsq;
//  double vec1vec2;
//  double sigik1[3],sigik2[3];
//  double wmin[3],wpls[3],vec1[3],vec2[3],vec3[3],vec4[3];
//  double y,ysq,tijktjil,tijksqtjilsq,titj,ybytijksq,ybytjilsq;
//  double dy[3],tijksq[2],tijk[3][2],dtijk[3][3][3][2];
//  double Tyzloc,dTyzdyloc,dTyzdzloc;
//  double virieltijk[3][3][3][3][2],dviriel[3][3];

//  double sqrt3=1.732050807569;
//  double rik1,rik2;
//  double wminabs,pf,pf0,pf00,pf01,pf02,oneppf,onempf,dpf,hpf,
//  double sigijwmin,dsigijwmin,sigijdwmin,wmindwmin,wminsq,vec1vec2;
//  double y,ysq,tijktjil,tijksqtjilsq,titj,ybytijksq,ybytjilsq;
//  double Tyzloc,dTyzdyloc,dTyzdzloc;
//  double dy[3],tijksq[2];
//  Vec3d sigik1,sigik2;
//  Vec3d wmin,wpls,vec1,vec2,vec3,vec4;
//  Vec3d tijk[2],dsigijk[3],dtijk[2][3][3];

  double dy[3],tijksq[2];
  Vec3d dfrc,tijk[2],dsigijk[3],dtijk[2][3][3];

#ifdef INFO
/* Initialisation */
double tdebutFonction;
InitialisationCPUFonction(NumFonc_TORSIONlchbop,"TORSIONlchbop",&tdebutFonction);
#endif

//!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
//!CC Comments:
//!CC  + inbcnb[2][2] = neighbours indexes of connected neighbours of
//!CC                   ip and jp other than jp and ip
//!CC  + sigij        = unit Vec3d (rj-ri)/rij
//!CC                   (when entering this routine sigij=(ri-rj)/rij !!!!)
//!CC  + sigik1       = unit Vec3d (rk1-ri)/rik1, with k1 first neighbour
//!CC                   other than i
//!CC  + sigik2       = unit Vec3d (rk2-ri)/rik2, with k2 first neighbour
//!CC                   other than i
//!CC  + dsigij[3]    = Vec3d array containing derivatives sigij w.r.t. rj
//!CC  + wmin         = Vec3d containing sigik1-sigik2
//!CC  + wpls         = Vec3d containing sigik1+sigik2
//!CC  + tijk[2]      = Vec3d array containing torion vectors tijk and tjil
//!CC  + dtijksq[2]   = array containing squared lengths of torsion vectors
//!CC  + dtijk[2][3][3] =
//!CC       Vec3d array containing derivatives with respect to x,y,z (first index)
//!CC       with respect to the position rj, rik1 and rik2 (third index)
//!CC       of the three components of the tijk (second index)
//!CC       for both torsion vectors (first index)
//!CC  + vec1,vec2,vec3,vec4 and dy are Vec3d workvectors
//!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!

  for(int ij=0;ij<2;ij++)
  {
    subdsig(sigij,dsigijk); // dsigijk is a 3x3 matrix but not of type Mat3d for the time being !!!!! //
//    printf("dsigijk  :%16.10lf%16.10lf%16.10lf\n",dsigijk[0].x,dsigijk[0].y,dsigijk[0].z);
//    printf("dsigijk  :%16.10lf%16.10lf%16.10lf\n",dsigijk[1].x,dsigijk[1].y,dsigijk[1].z);
//    printf("dsigijk  :%16.10lf%16.10lf%16.10lf\n",dsigijk[2].x,dsigijk[2].y,dsigijk[2].z);
    int inb=inbcnb[ij][0];
    double rik1=nbsij[ij].rstor[inb];
    Vec3d sigik1=nbsij[ij].sigstor[inb];
    inb=inbcnb[ij][1];
    double rik2=nbsij[ij].rstor[inb];
    Vec3d sigik2=nbsij[ij].sigstor[inb];
    sigij=-sigij;
    Vec3d wmin=sigik1-sigik2;
    Vec3d wpls=sigik1+sigik2;

//    printf("ip,jp   :%6zu%6zu\n",ip,jp);
//    printf("sigij   :%16.10lf%16.10lf%16.10lf\n",sigij.x,sigij.y,sigij.z);
//    printf("sigik1  :%16.10lf%16.10lf%16.10lf\n",sigik1.x,sigik1.y,sigik1.z);
//    printf("sigik2  :%16.10lf%16.10lf%16.10lf\n",sigik2.x,sigik2.y,sigik2.z);
//    printf("wmin    :%16.10lf%16.10lf%16.10lf\n",wmin.x,wmin.y,wmin.z);
//    printf("wpls    :%16.10lf%16.10lf%16.10lf\n",wpls.x,wpls.y,wpls.z);
  
    double sigijwmin=dot(sigij,wmin);
    double wminsq=dot(wmin,wmin);
    double wminabs=sqrt(wminsq);
    double pf00=sqrt3/wminabs;
    double pf=pf00*sigijwmin;

    Vec3d vec1=cross(sigij,wmin);
    Vec3d vec2=cross(sigij,wpls);
    tijk[ij]=vec1+pf*vec2;
    tijksq[ij]=dot(tijk[ij],tijk[ij]);
//    tijksq[ij]=tijk[0][ij]*tijk[0][ij]+tijk[1][ij]*tijk[1][ij]+tijk[2][ij]*tijk[2][ij];

//    printf("pf    :%16.10lf\n",pf);
//    printf("vec1  :%16.10lf%16.10lf%16.10lf\n",vec1.x,vec1.y,vec1.z);
//    printf("vec2  :%16.10lf%16.10lf%16.10lf\n",vec2.x,vec2.y,vec2.z);
//    printf("tijk  :%16.10lf%16.10lf%16.10lf\n",tijk[ij].x,tijk[ij].y,tijk[ij].z);
  
//      pf0=pf00;						//dtijk/drij
//      for (ixyz=0;ixyz<3;ixyz++) {
//	for (k=0;k<3;k++) vec3[k]=dsigijk[k][ixyz]/rij;
//	VECPROD(vec3,wmin,vec1);
//	VECPROD(vec3,wpls,vec4);
//	SCLPROD(vec3,wmin,&dsigijwmin);
//	dpf=pf0*dsigijwmin;
//	for (k=0;k<3;k++) dtijk[ixyz][0][k][ij]=vec1[k]+pf*vec4[k]+dpf*vec2[k];
//	for (l=0;l<3;l++) { for (k=0;k<3;k++) virieltijk[ixyz][l][0][k][ij]=dtijk[ixyz][0][k][ij]*sigij[l]*rij;}
//	for (l=0;l<3;l++) { for (k=0;k<3;k++) virieltijk[ixyz][l][0][k][ij]=dtijk[ixyz][0][k][ij]*sigij[l]*rij;}
//      printf("ixyz,dtijk0 :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,dtijk[ixyz][0][ij],dtijk[ixyz][0][1][ij],dtijk[ixyz][0][2][ij]);
//	}

//dtijk/drij
    double pf0=pf00;
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=dsigijk[ixyz]/rij;  // can we make a operator returning a column/row of a matrix ????? //
      vec1=cross(vec3,wmin);
      Vec3d vec4=cross(vec3,wpls);
      double dsigijwmin=dot(vec3,wmin);
      double dpf=pf0*dsigijwmin;
      dtijk[ij][0][ixyz]=vec1+pf*vec4+dpf*vec2;   // dtijk is a tensor of rank 4 with dimensions 2x3x3x3,
                                                  // defined as an 3D array of 2x3x3 Vec3d objects,
                                                  // to be declared as: Vec3d dtijk[2][3][3] //
//    for(l=0;l<3;l++){for(k=0;k<3;k++)virieltijk[ixyz][l][0][k][ij]=dtijk[ixyz][0][k][ij]*sigij[l]*rij;}
//      printf("ixyz,pf,dpf   :%6d%16.10lf%16.10lf\n",ixyz,pf,dpf);
//      printf("ixyz,vec3     :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,vec3.x,vec3.y,vec3.z);
//      printf("ixyz,vec1     :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,vec1.x,vec1.y,vec1.z);
//      printf("ixyz,pf*vec4  :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,pf*vec4.x,pf*vec4.y,pf*vec4.z);
//      printf("ixyz,dpf*vec2 :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,dpf*vec2.x,dpf*vec2.y,dpf*vec2.z);
//      printf("ixyz,dtijk0   :%6d%16.10lf%16.10lf%16.10lf\n",ixyz,dtijk[ij][0][ixyz].x,dtijk[ij][0][ixyz].y,dtijk[ij][0][ixyz].z);
    }
   
// dtijk/drik1
    subdsig(sigik1,dsigijk);
    double oneppf=1.+pf;
    double pf01=pf00;
    double pf02=pf/wminsq;
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=dsigijk[ixyz]/rik1;   // can we make a operator returning a column/row of a matrix ????? //
      vec1=cross(sigij,vec3);
      double sigijdwmin=dot(sigij,vec3);
      double wmindwmin=dot(wmin,vec3);
      double dpf=pf01*sigijdwmin-pf02*wmindwmin;
      dtijk[ij][1][ixyz]=oneppf*vec1+dpf*vec2;
//    for(l=0;l<3;l++){for(k=0;k<3;k++)virieltijk[ixyz][l][1][k][ij]=dtijk[ixyz][1][k][ij]*sigik1[l]*rik1;}
    }

// dtijk/drik2
    subdsig(sigik2,dsigijk);
    double onempf=1.-pf;
    pf01=pf00;
    pf02=pf/wminsq;
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=-dsigijk[ixyz]/rik2;  // can we make a operator returning a column/row of a matrix ????? //
      vec1=cross(sigij,vec3);
      double sigijdwmin=dot(sigij,vec3);
      double wmindwmin=dot(wmin,vec3);
      double dpf=pf01*sigijdwmin-pf02*wmindwmin;
      dtijk[ij][2][ixyz]=onempf*vec1+dpf*vec2;
//    for (l=0;l<3;l++){for(k=0;k<3;k++)virieltijk[ixyz][l][2][k][ij]=dtijk[ixyz][2][k][ij]*sigik2[l]*rik2;}
    }
  }
//  printf("tijk0 :%16.10lf%16.10lf%16.10lf\n",tijk[0].x,tijk[0].y,tijk[0].z);
//  printf("tijk1 :%16.10lf%16.10lf%16.10lf\n",tijk[1].x,tijk[1].y,tijk[1].z);

  double tijktjil=dot(tijk[0],tijk[1]);
  double tijksqtjilsq=tijksq[0]*tijksq[1];
//  double ysq=tijktjil*tijktjil/tijksqtjilsq;
  double y=tijktjil/sqrt(tijksqtjilsq);		//COSIJKL
//  printf("tijk :%16.10lf%16.10lf%16.10lf\n",tijk[0].x,tijk[0].y,tijk[0].z);
//  printf("tjil :%16.10lf%16.10lf%16.10lf\n",tijk[1].x,tijk[1].y,tijk[1].z);
//  printf("y    :%16.10lf\n",y);
//  printf("z    :%16.10lf\n",Ncj);
  TdTyz1(y,Ncj,Tyz,dTyzdy,dTyzdz);
//  printf("Tyz    :%16.10lf\n",Tyz);
//  printf("dTyzdy :%16.10lf\n",dTyzdy);
//  printf("dTyzdz :%16.10lf\n",dTyzdz);

  double pf=wijwji*dTyzdy*vsmra;
//  hpf=-0.5*pf;
  double titj=sqrt(tijksqtjilsq);

//  size_t ip=nbsij[0].nbtbl[0];
//  size_t jp=nbsij[1].nbtbl[0];
  double ybytijksq=y/tijksq[0];                       // Force due to vector Tijk
  Vec3d vec1=tijk[1]/titj-ybytijksq*tijk[0];                // dCOSIJKL/dTijk
//  printf("pf1,pf0 :%16.10lf%16.10lf\n",1./titj,-ybytijksq);
  for(int ixyz=0;ixyz<3;ixyz++)
  {
    Vec3d vec2=dtijk[0][0][ixyz];
//    printf("ixyz,pf,ybytijksq :%6d%16.10lf%16.10lf\n",ixyz,pf,ybytijksq);
//    printf("vec1 :%16.10lf%16.10lf%16.10lf\n",vec1.x,vec1.y,vec1.z);
//    printf("vec2 :%16.10lf%16.10lf%16.10lf\n",vec2.x,vec2.y,vec2.z);
    double vec1vec2=dot(vec1,vec2);
    dy[ixyz]=pf*vec1vec2;     // how to solve this
//    for(l=0;l<3;l++)
//    { 
//      for(k=0;k<3;k++)vec2[k]=virieltijk[ixyz][l][0][k][0];
//      sclprod(vec1,vec2,vec1vec2);
//      dviriel[l][ixyz]=hpf*vec1vec2;
//    }
  }
  dfrc={dy[0],dy[1],dy[2]};
//  printf("jp,ip,dfrc_1 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",jp,ip,dfrc.x,dfrc.y,dfrc.z);
  force[jp]=force[jp]+dfrc;    // i-j bond, atom j    // This should normally work (?), I geuss //
  force[ip]=force[ip]-dfrc;    // i-j bond, atom i
//  force[jp]=force[jp]+dy;                           // i-j bond, atom j
//  force[ip]=force[ip]-dy;                           // i-j bond, atom i
//  for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][jp]+=dviriel[l][k];}  // i-j bond, atom j
//  for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][ip]+=dviriel[l][k];}  // i-j bond, atom i

  for(int icnb=0;icnb<2;icnb++)
  {
    int inb=inbcnb[0][icnb];
    size_t kp=nbsij[0].nbtbl[inb];
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec2=dtijk[0][icnb+1][ixyz];
      double vec1vec2=dot(vec1,vec2);
      dy[ixyz]=pf*vec1vec2;
//      for(l=0;l<3;l++)
//      { 
//      	for (k=0;k<3;k++) vec2[k]=virieltijk[ixyz][l][icnb+1][k][0];
//      	sclprod(vec1,vec2,vec1vec2);
//      	dviriel[l][ixyz]=hpf*vec1vec2;
//      }
    }
    dfrc={dy[0],dy[1],dy[2]};
//    printf("kp,ip,dfrc_2 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",kp,ip,dfrc.x,dfrc.y,dfrc.z);
    force[kp]=force[kp]+dfrc;      // i-j bond, atom j    // This should normally work (?), I geuss //
    force[ip]=force[ip]-dfrc;      // i-j bond, atom i
//    for(k=0;k<3;k++)force[k][kp]+=dy[k];                               // i-k bond, atom k
//    for(k=0;k<3;k++)force[k][ip]-=dy[k];                               // i-k bond, atom i
//    for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][kp]+=dviriel[l][k];} // i-k bond, atom k
//    for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][ip]+=dviriel[l][k];} // i-k bond, atom i
  }

  double ybytjilsq=y/tijksq[1];                // Force due to vector Tjil
  vec1=tijk[0]/titj-ybytjilsq*tijk[1];         // dCOSIJKL/dTjil
  for(int ixyz=0;ixyz<3;ixyz++)
  {
    Vec3d vec2=dtijk[1][0][ixyz];
    double vec1vec2=dot(vec1,vec2);
    dy[ixyz]=pf*vec1vec2;
//    for(l=0;l<3;l++)
//    { 
//      for (k=0;k<3;k++) vec2[k]=virieltijk[ixyz][l][0][k][1];
//      sclprod(vec1,vec2,vec1vec2);
//      dviriel[l][ixyz]=hpf*vec1vec2;
//    }
  }
  dfrc={dy[0],dy[1],dy[2]};
//  printf("ip,jp,dfrc_3 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",ip,jp,dfrc.x,dfrc.y,dfrc.z);
  force[ip]=force[ip]+dfrc;        // i-j bond, atom j    // This should normally work (?), I geuss //
  force[jp]=force[jp]-dfrc;        // i-j bond, atom i
//  for(k=0;k<3;k++)force[k][ip]+=dy[k];                               // j-i bond, atom i
//  for(k=0;k<3;k++)force[k][jp]-=dy[k];                               // j-i bond, atom j
//  for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][ip]+=dviriel[l][k];} // j-i bond, atom i
//  for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][jp]+=dviriel[l][k];} // j-i bond, atom j

  for(int icnb=0;icnb<2;icnb++)
  {
    int inb=inbcnb[1][icnb];
    size_t lp=nbsij[1].nbtbl[inb];
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec2=dtijk[1][icnb+1][ixyz];
      double vec1vec2=dot(vec1,vec2);
      dy[ixyz]=pf*vec1vec2;
//      for(l=0;l<3;l++)
//      { 
//      	for(k=0;k<3;k++)vec2[k]=virieltijk[ixyz][l][icnb+1][k][1];
//      	sclprod(vec1,vec2,vec1vec2);
//      	dviriel[l][ixyz]=hpf*vec1vec2;
//      }
    }
    dfrc={dy[0],dy[1],dy[2]};
//    printf("lp,jp,dfrc_4 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",lp,jp,dfrc.x,dfrc.y,dfrc.z);
    force[lp]=force[lp]+dfrc;        // i-j bond, atom j    // This should normally work (?), I geuss //
    force[jp]=force[jp]-dfrc;        // i-j bond, atom i
//    for(k=0;k<3;k++)force[k][lp]+=dy[k];                               // j-l bond, atom l
//    for(k=0;k<3;k++)force[k][jp]-=dy[k];                               // j-l bond, atom j
//    for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][lp]+=dviriel[l][k];} // j-l bond, atom l
//    for(l=0;l<3;l++){for(k=0;k<3;k++)viriel[l][k][jp]+=dviriel[l][k];} // j-l bond, atom j
  }

#ifdef INFO
/* Fin de la fonction */
FinCPUFonction(NumFonc_TORSIONlchbop,"TORSIONlchbop",tdebutFonction);
#endif

}  // end torsion1

//----------------------------------------------------------------------//

void FATP::torsion2(size_t ijp[2],double Ncj,double wijwji,double vsmra,
double rij,Vec3d sigij,double& Tyz,double& dTyzdysq,double& dTyzdz,
LNBStr* nbsij,vector<Vec3d>& force,vector<Mat3d>& viriel)
{
  double sqrt3=1.732050807569;
  double rik1[2],rik2[2],wminsq[2],wminabs[2],dysqdtijk[2];
  double stijk[2],dstijkdw[2];
//  double stijk[2],stijk12,dstijkdw[2];
  double invtijkabs[2],dy[3];
  double wminlow,wminhigh,invwminhighlow;
  Vec3d sigik1[2],sigik2[2],wmin[2],tijk[2],dsigijk[3],dtijk[2][3][3];

//  int ij;
//  int k,l,kp,lp;
//  int inb,icnb,ixyz;
//  double pf,pf0,pf00,pf01,pf02,oneppf,onempf,dpf,hpf;
//  double sigijwmin,dsigijwmin,sigijdwmin,wmindwmin;
//  double vec1vec2;
//  double y0,y,ysq,dysqds1s2,dysqdtijk[2];
//  double dy[3],tijktjil,y0bytijk;
#ifdef VIRIEL
  Vec3d virieltijk[2][3][3][3];
  double dviriel[3][3];
#endif

//#ifdef INFO
///* Initialisation */
//double tdebutFonction;
//InitialisationCPUFonction(NumFonc_TORSIONlchbop,"TORSIONlchbop",&tdebutFonction);
//#endif

//!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
//!CC Comments:
//!CC   + inbcnb[2][2] = neighbours indexes of connected neighbours of
//!CC                   ip and jp other than jp and ip
//!CC   + sigij[3]    = unit vector (rj-ri)/rij
//!CC                   (when entering this routine sigij=(ri-rj)/rij !!!!)
//!CC   + sigik1[3]   = unitvector (rk1-ri)/rik1, with k1 first neighbour
//!CC                   other than i
//!CC   + sigik2[3]   = unitvector (rk2-ri)/rik2, with k2 first neighbour
//!CC                   other than i
//!CC   + dsigij[3][3] = derivatives of of sigij w.r.t. rj
//!CC   + wmin[3]     = sigik1-sigik2
//!CC   + wpls[3]     = sigik1+sigik2
//!CC   + dsigij[3][3] = derivatives of of sigij w.r.t. rj
//!CC   + tijk[3][2]   = torsion vectors tijk and tjil
//!CC   + dtijksq[2]  = squared length of the torsion vectors
//!CC   + dtijk[3][3][3][2] =
//!CC        derivatives with respect to x,y,z (first index)
//!CC        with respect to the position rj, rik1 and rik2 (second index)
//!CC        of the three components of the tijk (third index)
//!CC        for both torsion vectors (fourth index)
//!CC   + vec1,vec2,vec3,vec4 and dy are workvectors
//!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!

//C Boundaries for cutoff function on wminabs
//  wminlow=0.005;wminhigh=0.020;
  wminlow=0.05;wminhigh=0.20;
  invwminhighlow=1./(wminhigh-wminlow);

//  printf("ip,jp,sigij :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",ijp[0],ijp[1],sigij.x,sigij.y,sigij.z);
//  std::cin.get();

//pre-calculation of the wmin and their cutoff functions
//  leave routine if one of the contribs is 0 (no torsion).
  for(int ij=0;ij<2;ij++) 
  {
    int inb=inbcnb[ij][0];
    rik1[ij]=nbsij[ij].rstor[inb];
    sigik1[ij]=nbsij[ij].sigstor[inb];
    inb=inbcnb[ij][1];
    rik2[ij]=nbsij[ij].rstor[inb];
    sigik2[ij]=nbsij[ij].sigstor[inb];
    wmin[ij]=sigik1[ij]-sigik2[ij];
    wminsq[ij]=dot(wmin[ij],wmin[ij]);
    wminabs[ij]=sqrt(wminsq[ij]);
    
  //  add a correction term to avoid singuarity when i, k1 and k2 are aligned
    if(wminabs[ij]<=wminlow) //  wmin[ij] too small, torsion canceled
    {
      Tyz=0.;
      dTyzdysq=0.;
      dTyzdz=0.;
//      printf("OK1 in torsion2; wminabs[ij] = %16.10lf\n",wminabs[ij]);
      return;
    } 
    else if(wminabs[ij]<wminhigh) // correction function activated
    {
//      printf("OK2 in torsion2 \n");
      SqdSqup(wminabs[ij],wminlow,wminhigh,invwminhighlow,stijk[ij],dstijkdw[ij]);
    }
    else   // wmin[ij] large enough, correction function not needed
    {
      stijk[ij]=1.;
      dstijkdw[ij]=0.;
    }
  }
//  stijk12=stijk[0]*stijk[1];

// compute vectors tijk and derivatives dtijk
  for(int ij=0;ij<2;ij++)
  {
// derivative tensor of i-j vector should be computed before i-j is modified
    subdsig(sigij,dsigijk); // dsigijk is a 3x3 matrix but not of type Mat3d for the time being !!!!! //
    
// compute tijk and its inverse module
// needs wpls and various prefactors (wmin has been computed and stored previously)
    sigij=-sigij;
    Vec3d wpls=sigik1[ij]+sigik2[ij];
    double sigijwmin=dot(sigij,wmin[ij]);
    double pf00=sqrt3/wminabs[ij];
    double pf=pf00*sigijwmin;
    Vec3d vec1=cross(sigij,wmin[ij]);
    Vec3d vec2=cross(sigij,wpls);
    tijk[ij]=vec1+pf*vec2;
    invtijkabs[ij]=1./sqrt(dot(tijk[ij],tijk[ij]));

//    printf("ij,pf,invtijkabs[ij] :%6d%16.10lf%16.10lf\n",ij,pf,invtijkabs[ij]);
//    printf("sigik1 :%16.10lf%16.10lf%16.10lf\n",sigik1[ij].x,sigik1[ij].y,sigik1[ij].z);
//    printf("sigik2 :%16.10lf%16.10lf%16.10lf\n",sigik2[ij].x,sigik2[ij].y,sigik2[ij].z);
//    printf("wmin   :%16.10lf%16.10lf%16.10lf\n",wmin[ij].x,wmin[ij].y,wmin[ij].z);
//    printf("wpls   :%16.10lf%16.10lf%16.10lf\n",wpls.x,wpls.y,wpls.z);
//    printf("tijk   :%16.10lf%16.10lf%16.10lf\n",tijk[ij].x,tijk[ij].y,tijk[ij].z);
//    std::cin.get();
    
// compute dtijk/drij
// dtijk is a tensor of rank 4 with dimensions 2x3x3x3,
// supposed to be defined as an 3D array of 2x3x3 Vec3d objects,
// to be declared as: Vec3d dtijk[2][3][3] //
    double pf0=pf00;
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=dsigijk[ixyz]/rij;   // can we make a operator returning a column (or row) of a matrix ????? //
      vec1=cross(vec3,wmin[ij]);
      Vec3d vec4=cross(vec3,wpls);
      double dsigijwmin=dot(vec3,wmin[ij]);
      double dpf=pf0*dsigijwmin;
      dtijk[ij][0][ixyz]=vec1+pf*vec4+dpf*vec2;   
#ifdef VIRIEL
      virieltijk[ij][0][ixyz][0]=rij*sigij.x*dtijk[ij][0][ixyz];
      virieltijk[ij][0][ixyz][1]=rij*sigij.y*dtijk[ij][0][ixyz];
      virieltijk[ij][0][ixyz][2]=rij*sigij.z*dtijk[ij][0][ixyz];
#endif
    }
    
//  compute dtijk/drik1
    subdsig(sigik1[ij],dsigijk);
    double oneppf=1.+pf;
    double pf01=pf00;
    double pf02=pf/wminsq[ij];
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=dsigijk[ixyz]/rik1[ij];  // can we make a operator returning a column (or row) of a matrix ????? //
      vec1=cross(sigij,vec3);
      double sigijdwmin=dot(sigij,vec3);
      double wmindwmin=dot(wmin[ij],vec3);
      double dpf=pf01*sigijdwmin-pf02*wmindwmin;
      dtijk[ij][1][ixyz]=oneppf*vec1+dpf*vec2;
#ifdef VIRIEL
      virieltijk[ij][1][ixyz][0]=rik1[ij]*sigik1[ij].x*dtijk[ij][1][ixyz];
      virieltijk[ij][1][ixyz][1]=rik1[ij]*sigik1[ij].y*dtijk[ij][1][ixyz];
      virieltijk[ij][1][ixyz][2]=rik1[ij]*sigik1[ij].z*dtijk[ij][1][ixyz];
#endif
    }
  
// dtijk/drik2
    subdsig(sigik2[ij],dsigijk);
    double onempf=1.-pf;
    pf01=pf00;
    pf02=pf/wminsq[ij];
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec3=-dsigijk[ixyz]/rik2[ij];   // can we make a operator returning a column (or row) of a matrix ????? //
      vec1=cross(sigij,vec3);
      double sigijdwmin=dot(sigij,vec3);
      double wmindwmin=dot(wmin[ij],vec3);
      double dpf=pf01*sigijdwmin-pf02*wmindwmin;
      dtijk[ij][2][ixyz]=onempf*vec1+dpf*vec2;
#ifdef VIRIEL
      virieltijk[ij][2][ixyz][0]=rik2[ij]*sigik2[ij].x*dtijk[ij][2][ixyz];
      virieltijk[ij][2][ixyz][1]=rik2[ij]*sigik2[ij].y*dtijk[ij][2][ixyz];
      virieltijk[ij][2][ixyz][2]=rik2[ij]*sigik2[ij].z*dtijk[ij][2][ixyz];
#endif
    }
  }

//####################################################################
//COSIJKLsq

  double tijktjil=dot(tijk[0],tijk[1]);
  double y0=tijktjil*invtijkabs[0]*invtijkabs[1];
  double dysqds1s2=0.;
  double ysq=y0*y0;
  if(stijk[0]*stijk[1]==1.) // no contributions from cutoff functions
  {
    dysqdtijk[0]=2.*y0*invtijkabs[0];
    dysqdtijk[1]=2.*y0*invtijkabs[1];
  }
  else //  contributions from cutoff functions, need to update dtijk
  {
    double stijk12=stijk[0]*stijk[1];
    ysq=1.+stijk12*(ysq-1.);
    dysqds1s2=(y0*y0-1.);
    dysqdtijk[0]=stijk12*2.*y0*invtijkabs[0];
    dysqdtijk[1]=stijk12*2.*y0*invtijkabs[1];
  }
//####################################################################

//  printf("tijk: %lf %lf %lf\n",tijk[0][0],tijk[1][0],tijk[2][0]);
//  printf("tjil: %lf %lf %lf\n",tijk[0][1],tijk[1][1],tijk[2][1]);
//  printf("ysq: %lf\n",ysq);
//  printf("z: %lf\n",Ncj);

  TdTyz2(ysq,Ncj,Tyz,dTyzdysq,dTyzdz);
//  printf("ysq,Tyz,dTyzdysq,dTyzdz : %16.10lf%16.10lf%16.10lf%16.10lf\n",ysq,Tyz,dTyzdysq,dTyzdz);
//  std::cin.get();

// compute atomic forces and virial
  for(int ij=0;ij<2;ij++)
  {
    double pf=wijwji*vsmra*dTyzdysq*dysqdtijk[ij];
#ifdef VIRIELLOCAL
    double hpf=-0.5*pf;
#else 
#ifdef VIRIEL
    double hpf=-pf;
#endif      
#endif      
// dCOSIJKL/dTijk
    double y0bytijk=y0*invtijkabs[ij];
    Vec3d vec1=invtijkabs[1-ij]*tijk[1-ij]-y0bytijk*tijk[ij];
//    printf("pf1,pf0 :%16.10lf%16.10lf\n",invtijkabs[1-ij],-y0bytijk);
    
// derivatives wrt. i-j bond
    for(int ixyz=0;ixyz<3;ixyz++)
    {
      Vec3d vec2=dtijk[ij][0][ixyz];
//      printf("ixyz,pf,y0bytijk :%6d%16.10lf%16.10lf\n",ixyz,pf,y0bytijk);
//      printf("vec1 :%16.10lf%16.10lf%16.10lf\n",vec1.x,vec1.y,vec1.z);
//      printf("vec2 :%16.10lf%16.10lf%16.10lf\n",vec2.x,vec2.y,vec2.z);
//      for(k=0;k<3;k++)vec2[k]=dtijk[ixyz][0][k][ij];
      double vec1vec2=dot(vec1,vec2);
      dy[ixyz]=pf*vec1vec2;     // how to solve this
#ifdef VIRIEL
      for(int l=0;l<3;l++)
      { 
        vec2=virieltijk[ij][0][ixyz][l];
        double vec1vec2=dot(vec1,vec2);
        dviriel[l][ixyz]=hpf*vec1vec2;
      }
#endif
    }
    Vec3d dfrc={dy[0],dy[1],dy[2]};
//    printf("jp,ip,dfrc_1 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",ijp[1-ij],ijp[ij],dfrc.x,dfrc.y,dfrc.z);
    force[ijp[ij]]=force[ijp[ij]]-dfrc;       // i-j bond, atom j
    force[ijp[1-ij]]=force[ijp[1-ij]]+dfrc;   // i-j bond, atom i
#ifdef VIRIELLOCAL
    viriel[ijp[ij]]=tensorsum(viriel[ijp[ij]],dviriel);
    viriel[ijp[1-ij]]=tensorsum(viriel[ijp[1-ij]],dviriel);
#else 
#ifdef VIRIEL
    viriel[0]=tensorsum(viriel[0],dviriel);
#endif
#endif

// derivatives wrt. i-k1/k2 bonds
    for(int icnb=0;icnb<2;icnb++)
    {
      int inb=inbcnb[ij][icnb];
      size_t kp=nbsij[ij].nbtbl[inb];
      for(int ixyz=0;ixyz<3;ixyz++)
      {
        Vec3d vec2=dtijk[ij][icnb+1][ixyz];
        double vec1vec2=dot(vec1,vec2);
        dy[ixyz]=pf*vec1vec2;
#ifdef VIRIEL
        for(int l=0;l<3;l++)
        { 
          vec2=virieltijk[ij][icnb+1][ixyz][l];
          double vec1vec2=dot(vec1,vec2);
          dviriel[l][ixyz]=hpf*vec1vec2;
        }
#endif
      }
      Vec3d dfrc={dy[0],dy[1],dy[2]};
//      printf("kp,ip,dfrc_2 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",kp,ijp[ij],dfrc.x,dfrc.y,dfrc.z);
      force[kp]=force[kp]+dfrc;            // i-j bond, atom j
      force[ijp[ij]]=force[ijp[ij]]-dfrc;  // i-j bond, atom i
#ifdef VIRIELLOCAL
      viriel[kp]=tensorsum(viriel[kp],dviriel);
      viriel[ijp[ij]]=tensorsum(viriel[ijp[ij]],dviriel);
#else 
#ifdef VIRIEL
      viriel[0]=tensorsum(viriel[0],dviriel);
#endif
#endif
    }

// stijk contribution
    if(dstijkdw[ij]!=0.)
    {
// global prefactor
      double pf0=wijwji*vsmra*dTyzdysq*dysqds1s2*stijk[1-ij]*dstijkdw[ij]/wminabs[ij];
    
// i-k1 contribution
      subdsig(sigik1[ij],dsigijk);
      double pf01=pf0/rik1[ij];
      for(int ixyz=0;ixyz<3;ixyz++)
      {
        Vec3d vec3=dsigijk[ixyz];	// + sign BECAUSE dwmin=dsigik1-dsigik2
        double wmindwmin=dot(wmin[ij],vec3);
        dy[ixyz]=pf01*wmindwmin;
      }
      int inb=inbcnb[ij][0];
      size_t kp=nbsij[ij].nbtbl[inb];
      Vec3d dfrc={dy[0],dy[1],dy[2]};
//      printf("kp,ip,dfrc_3 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",kp,ijp[ij],dfrc.x,dfrc.y,dfrc.z);
      force[kp]=force[kp]+dfrc;            // i-j bond, atom j
      force[ijp[ij]]=force[ijp[ij]]-dfrc;  // i-j bond, atom i
#ifdef VIRIELLOCAL
      viriel[kp]=viriel[kp]+tensor(sigik1[ij],rik1[ij]*dfrc);		//  i-k1 bond, atom k1
      viriel[ijp[ij]]=viriel[ijp[ij]]+tensor(sigik1[ij],rik1[ij]*dfrc);      //  i-k1 bond, atom i
#else 
#ifdef VIRIEL
      viriel[0]=viriel[0]+tensor(sigik1[ij],2.0*rik1[ij]*dfrc);        //  i-k1 bond, atom i
#endif
#endif

// i-k2 contribution
      subdsig(sigik2[ij],dsigijk);
      pf01=pf0/rik2[ij];
      for(int ixyz=0;ixyz<3;ixyz++)
      {
        Vec3d vec3=-dsigijk[ixyz];	// - sign BECAUSE dwmin=dsigik1-dsigik2
        double wmindwmin=dot(wmin[ij],vec3);
        dy[ixyz]=pf01*wmindwmin;
      }
      inb=inbcnb[1][ij];
      kp=nbsij[ij].nbtbl[inb];
      dfrc={dy[0],dy[1],dy[2]};
//      printf("kp,ip,dfrc_4 :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",kp,ijp[ij],dfrc.x,dfrc.y,dfrc.z);
      force[kp]=force[kp]+dfrc;            // i-j bond, atom j
      force[ijp[ij]]=force[ijp[ij]]-dfrc;  // i-j bond, atom i
#ifdef VIRIELLOCAL
      viriel[kp]=viriel[kp]+tensor(sigik2[ij],rik2[ij]*dfrc);		//  i-k2 bond, atom k2
      viriel[ijp[ij]]=viriel[ijp[ij]]+tensor(sigik2[ij],rik2[ij]*dfrc);      //  i-k1 bond, atom i
#else 
#ifdef VIRIEL
      viriel[0]=viriel[0]+tensor(sigik2[ij],2.0*rik2[ij]*dfrc);        //  i-k2 bond, atom i
#endif
#endif
    }
  }

#ifdef INFO
/* Fin de la fonction */
FinCPUFonction(NumFonc_TORSIONlchbop,"TORSIONlchbop",tdebutFonction);
#endif

}  // end torsion2

//==========================================================================
//------------------------------------------------------------------------//
//
//void SUBDSIG(double *sig)
//{
//int i1,i2;
//
//for (i1=0;i1<3;i1++) {
//	for (i2=i1;i2<3;i2++) {
//		dsigijk[i2][i1]=-sig[i2]*sig[i1];
//		dsigijk[i1][i2]=dsigijk[i2][i1];
//	}
//}
//for (i1=0;i1<3;i1++) dsigijk[i1][i1]=dsigijk[i1][i1]+1.;
//
//}
//
//------------------------------------------------------------------------//

inline void subdsig(Vec3d sig, Vec3d* dsigijk)
{
  dsigijk[0].x=1.0-sig.x*sig.x;
  dsigijk[1].y=1.0-sig.y*sig.y;
  dsigijk[2].z=1.0-sig.z*sig.z;

  dsigijk[0].y=-sig.x*sig.y;
  dsigijk[0].z=-sig.x*sig.z;
  dsigijk[1].z=-sig.y*sig.z;

  dsigijk[1].x=dsigijk[0].y;
  dsigijk[2].x=dsigijk[0].z;
  dsigijk[2].y=dsigijk[1].z;
  
//  for (i1=0;i1<3;i1++)
//  for (i2=i1;i2<3;i2++)
//  {
//    dsigijk[i2][i1]=-sig[i2]*sig[i1];
//    dsigijk[i1][i2]=dsigijk[i2][i1];
//  }
//  for(i1=0;i1<3;i1++)dsigijk[i1][i1]=dsigijk[i1][i1]+1.;
}

//-----------------------------------------------------------------------//

void FATP::TdTyz1(double y, double z,double& Tyz,double& dTyzdy,double& dTyzdz)
{
//{/*TESTFORCES*/Tyz=0.;dTyzdy=0.;dTyzdz=0.;printf("TYZ desactive\n");return;}/*CANCELS OUT Tyz CONTRIBS*/

  double ysq=y*y;
  double ytr=ysq*y;
  double yqu=ytr*y;
  
  double omysq=1.-ysq;
  double domysq=-2.*y;

  if(z<1.e-12)
  {
    double xsi1=At10+At11*z;
    double xsi2=At20+At21*z;
    double tau01=xsi1+xsi2*yqu;
    double dtau01dy=4.*xsi2*ytr;
    double dtau01dz=At11+At21*yqu;
    double omysq=1.-ysq;
    double domysq=-2.*y;
    Tyz=tau01*omysq;
    dTyzdy=dtau01dy*omysq+tau01*domysq;
    dTyzdz=dtau01dz*omysq;
//    Tyz=tau0;
//    dTyzdy=dtau0dy;
//    dTyzdz=dtau0dz;
  }
  else if
  (z>1.-1.e-12)
  {
    double tau11=Bt1+Bt2*ysq+Bt3*yqu;
    double tau12=1.+Bt4*ysq;
    double dtau11dy=2.*Bt2*y+4.*Bt3*ytr;
    double dtau12dy=2.*Bt4*y;
    Tyz=tau11/tau12*omysq;
    dTyzdy=(dtau11dy*omysq+tau11*domysq-Tyz*dtau12dy)/tau12;
    dTyzdz=0.;
  }
  else
  {
    double xsi1=At10+At11*z;
    double xsi2=At20+At21*z;
    double tau01=xsi1+xsi2*yqu;
    double dtau01dy=4.*xsi2*ytr;
    double dtau01dz=At11+At21*yqu;
    double omysq=1.-ysq;
    double domysq=-2.*y;
    double tau0=tau01*omysq;
    double dtau0dy=dtau01dy*omysq+tau01*domysq;
    double dtau0dz=dtau01dz*omysq;
    
    double tau11=Bt1+Bt2*ysq+Bt3*yqu;
    double tau12=1.+Bt4*ysq;
    double dtau11dy=2.*Bt2*y+4.*Bt3*ytr;
    double dtau12dy=2.*Bt4*y;
    double tau1=tau11/tau12*omysq;
    double dtau1dy=(dtau11dy*omysq+tau11*domysq-tau1*dtau12dy)/tau12;

    double zsq=z*z;
    double omztr=1.-zsq*z;
    Tyz=omztr*tau0+zsq*tau1;
    dTyzdy=omztr*dtau0dy+zsq*dtau1dy;
    dTyzdz=-3.*zsq*tau0+omztr*dtau0dz+2.*z*tau1;
  }

//printf("\t\t\ttyz: (%18.12lf;%18.12lf)   %18.12lf %18.12lf  %18.12lf %18.12lf  %18.12lf\n",ysq,z,omztr*(1.+zpzsq*tau2),tau0,ztr,tau1,Tyz);

}

//==========================================================================

void FATP::TdTyz2(double ysq, double z,double& Tyz,double& dTyzdysq,double& dTyzdz)
{
//{/*TESTFORCES*/*Tyz=0.;*dTyzdy=0.;*dTyzdz=0.;printf("TYZ desactive\n");return;}/*CANCELS OUT Tyz CONTRIBS*/
  double yqu=ysq*ysq;
  double omysq=1.-ysq;
  
  if(z<1.e-12)
  {
    double xsi1=At10+At11*z;
    double xsi2=At20+At21*z;
    double tau01=xsi1+xsi2*yqu;
    double dtau01dysq=2.*xsi2*ysq;
    double dtau01dz=At11+At21*yqu;
    Tyz=tau01*omysq;
    dTyzdysq=dtau01dysq*omysq-tau01;
    dTyzdz=dtau01dz*omysq;
  }
  else if(z>1.-1.e-12)
  {
    double tau11=Bt1+Bt2*ysq+Bt3*yqu;
    double dtau11dysq=Bt2+2.*Bt3*ysq;
    double tau12=1.+Bt4*ysq;
    double dtau12dysq=Bt4;
    Tyz=tau11/tau12*omysq;
    dTyzdysq=(dtau11dysq*omysq-tau11-Tyz*dtau12dysq)/tau12;
    dTyzdz=0.;
  }
  else
  {
    double zsq=z*z;
    double omztr=1.-zsq*z;
    double xsi1=At10+At11*z;
    double xsi2=At20+At21*z;
    double tau01=xsi1+xsi2*yqu;
    double dtau01dysq=2.*xsi2*ysq;
    double dtau01dz=At11+At21*yqu;
    double tau0=tau01*omysq;
    double dtau0dysq=dtau01dysq*omysq-tau01;
    double dtau0dz=dtau01dz*omysq;
  
    double tau11=Bt1+Bt2*ysq+Bt3*yqu;
    double dtau11dysq=Bt2+2.*Bt3*ysq;
    double tau12=1.+Bt4*ysq;
    double dtau12dysq=Bt4;
    double tau1=tau11/tau12*omysq;
    double dtau1dysq=(dtau11dysq*omysq-tau11-tau1*dtau12dysq)/tau12;
  
    Tyz=omztr*tau0+zsq*tau1;
    dTyzdysq=omztr*dtau0dysq+zsq*dtau1dysq;
    dTyzdz=-3.*zsq*tau0+omztr*dtau0dz+2.*z*tau1;
  }
//printf("\t\t\ttyz: (%18.12lf;%18.12lf)     %18.12lf %18.12lf       %18.12lf %18.12lf     %18.12lf\n",ysq,z,omztr*(1.+zpzsq*tau2),tau0,ztr,tau1,*Tyz);
}

//-----------------------------------------------------------------------//

inline Mat3d tensorsum(Mat3d& t1,double (&t2)[3][3])
{
  return Mat3d{t1.m11+t2[0][0],t1.m12+t2[0][1],t1.m13+t2[0][2],
               t1.m21+t2[1][0],t1.m22+t2[1][1],t1.m23+t2[1][2],
               t1.m31+t2[2][0],t1.m32+t2[2][1],t1.m33+t2[2][2]};
}

//-----------------------------------------------------------------------//

