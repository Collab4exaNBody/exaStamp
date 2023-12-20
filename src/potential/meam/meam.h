/// @file 
/// @brief Definition of MEAM  potential
/// @author Raphael Pratt (c) CEA 2016

#pragma once

#include <cmath>
#include <string>
#include <cassert>

#include <exanb/core/basic_types.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro

#include "meam_parameters.h"

#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_math.h>

/// @brief Shorcut to swap to variables
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}


using exanb::Vec3d;


/// @brief MEAM potential where force is computed without 
///
/// This is a MEAM potential developed to calculate the properties of tin near melting curve
struct MeamPotential
{
  MeamPotential() = default;
  MeamPotential(const MeamPotential&) = default;
  MeamPotential(MeamPotential &&) = default;
  MeamPotential& operator = (const MeamPotential&) = default;
  MeamPotential& operator = (MeamPotential &&) = default;

  /// @brief Constructor
  /// @param [in] rcut Cutoff radius
  /// @param [in] parameters Parameters
  MeamPotential(double rcut, const MeamParameters& parameters)
    : p(parameters) 
    {
    using onika::cuda::max;
    // ne sert pas car le polynome de lissage assure la nullité a rcut
    //this->initializeCutoffVariables(); 
        
    double tmp = max(rcut,parameters.rc*parameters.Cmax/sqrt(4.*(parameters.Cmax-1.)));
    // tmp = rcut le vrai
    p.rd = tmp-p.rp;
    p.Rcut = tmp;

    this->ComputeMatrix();              // Compute the matrix fill by coeff of the smoothing polynomial
    
    assert( p.rc == p.Rcut );
  }

  /// @brief Destructor (nothing to do)
  ~MeamPotential() = default;

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  ONIKA_HOST_DEVICE_FUNC inline void fEmbed(double rho, double& f, double& df) const
  {
    rho/=p.Z;                         // to compute Frho
    double coeff=p.A*p.E0;		        //coeff  = A.Ecoh
    double coeff1=coeff*log(rho);	//coeff1 = A.Ecoh.ln(rho)

    //f = A . Ecoh . roh . ln(rho)
    f=coeff1*rho;

    //df= A . Ecoh . ( 1 + ln(rho) )
    df=coeff+coeff1;
  }

  /// @brief Inverse a square matrix (same function as in STAMP)
  /// @param [in,out] a Square matrix
  /// @param [in] n Size of the matrix
  template<int n>
  void InversionMatriceCarree(double a[n][n])
  {
//    static constexpr size_t n = 4;
  	int indxc[n], indxr[n], ipiv[n];
  	int i,icol,irow,j,k,l,ll;	
  	double big,dum,pivinv,temp;
  	//double	ar[n][n];
  	
      icol=-1;
      irow=-1;

//  	indxc = new int[n];
//  	indxr = new int[n];
//  	ipiv = new int[n];

/*
  	ar = new double*[n];
  	for(i=0;i<n;i++){
  		ar[i] = new double[n];
  	}	
*/
  	
  	for(i=0;i<n;i++)
  		indxc[i] = indxr[i] = ipiv[i] = 0;
  			
  	for (j=0;j<n;j++) ipiv[j]=0;
  	for (i=0;i<n;i++) {
  		big=0.0;
  		for (j=0;j<n;j++)
  			if (ipiv[j] != 1)
  				for (k=0;k<n;k++) {
  					if (ipiv[k] == 0) {
  						if (fabs(a[j][k]) >= big) {
  							big=fabs(a[j][k]);
  							irow=j;
  							icol=k;
  						}
  					}
  				}
  		++(ipiv[icol]);
  		if (irow != icol) {
  			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
  		}
  		
  		indxr[i]=irow;
  		indxc[i]=icol;
  		if (a[icol][icol] == 0.0) {
              std::abort(); //exit(0);
  		}
  		
  		pivinv=1.0/a[icol][icol];
  		
  		a[icol][icol]=1.0;
  		for (l=0;l<n;l++) a[icol][l] *= pivinv;
  		for (ll=0;ll<n;ll++)
  			if (ll != icol) {
  				dum=a[ll][icol];
  				a[ll][icol]=0.0;
  				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
  			}
  	}
  	for (l=n-1;l>=0;l--) {
  		if (indxr[l] != indxc[l])
  			for (k=0;k<n;k++)
  				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  	}
  }

  /// @brief Calculate the coefficients of the lifting polynome of Emu and rho_i
  ///
	///
  inline void ComputeMatrix() 
  {
    //assert( cutoffRadius == p.Rcut );

    p.rc = p.Rcut; //cutoffRadius; // doit etre vérifié dans ExaStamp/MEAM
    double rd = p.rc-p.rp;

    //matrix allocation
    //p.m=new double*[5];
    //for(unsigned int j=0;j<5;j++) p.m[j] =new double[4];

    double mat[4][4]; //=new double*[4];
    //for(unsigned int j=0;j<4;j++) mat[j] =new double[4];

    double p_m[5][4];
    for(unsigned int i=0;i<5;i++){
		  for(unsigned int j=0;j<4;j++){
			  p_m[i][j] = 0;		
		  }
	  }

    double u[4];
    double nbeta[4];    

    mat[0][0] = 1.;
    mat[0][1] = p.rc;	
    mat[0][2] = p.rc*p.rc;		
    mat[0][3] = p.rc*p.rc*p.rc;			

    mat[1][0] = 0.;
    mat[1][1] = 1;	
    mat[1][2] = 2.*p.rc;		
    mat[1][3] = 3.*p.rc*p.rc;			
		
    mat[2][0] = 1.;
    mat[2][1] = rd;	
    mat[2][2] = rd*rd;		
    mat[2][3] = rd*rd*rd;
	
    mat[3][0] = 0.;
    mat[3][1] = 1;	
    mat[3][2] = 2.*rd;		
    mat[3][3] = 3.*rd*rd;

    if(p.rp != 0.) InversionMatriceCarree(mat);

    u[0]=0.;u[1]=0.;

    //fill negative beta vector
    nbeta[0]=-p.beta0;
    nbeta[1]=-p.beta1;
    nbeta[2]=-p.beta2;
    nbeta[3]=-p.beta3;

    double tmp0;
	  double	as,as2,as3,texp,dasdrij;
		
    // cut function        
    // functions f(r) depend on the distance, they are permuted by a smoothing polynomial P(r) degree 3 if r>rd
    // such as P(rdMEAM) = f(rdMEAM), P'(rdMEAM)=f'(rdMEAM), P(rCut)=0 and P'(rCut)=0 

    //for polynome Emu
	  as = p.alpha * (rd*p.r0 -1.);
	  dasdrij = p.alpha*p.r0;
	  as2 = as * as;
	  as3 = as2 * as;
	  
	  texp = exp(-as);
	  
	  u[2] = -p.E0 * texp * (1. + as + (p.delta / p.r0) * as3 / rd );
	  u[3] = -u[2] * dasdrij - p.E0 * texp * ( dasdrij + 3. * (p.delta / p.r0) * as2 / rd * dasdrij - (p.delta / p.r0) * as3 / (rd * rd));

    for(unsigned int i=0;i<4;i++){
		for(unsigned int j=0;j<4;j++){
			p_m[0][i] += mat[i][j] * u[j];		
		}
	}

    //for polynome rho_i
    tmp0=rd*p.r0-1.;
    for(unsigned int nb=0;nb<4;nb++)
    {
      u[2]=exp(nbeta[nb]*tmp0);
      u[3]=(nbeta[nb]*p.r0)*u[2];
      for(unsigned int i=0;i<4;i++)
      {
      	for(unsigned int j=0;j<4;j++)
        {
      		p_m[nb+1][i] += mat[i][j] * u[j];		
      	}
      }
    }
    
    for(unsigned int i=0;i<5;i++)
		  for(unsigned int j=0;j<4;j++)
		   	p.PolynomeMeam[4*i+j]= p_m[i][j];
  }
  
  
  // NEW Functions implemented for the amr feature
  
  ONIKA_HOST_DEVICE_FUNC inline void screeningFunction(
      const double* __restrict__ drx_
    , const double* __restrict__ dry_
    , const double* __restrict__ drz_
    , double* S_
    , const size_t n) const
  {
    using onika::cuda::min;
    
    ONIKA_ASSUME_ALIGNED( drx_ );
    ONIKA_ASSUME_ALIGNED( dry_ );
    ONIKA_ASSUME_ALIGNED( drz_ );
    ONIKA_ASSUME_ALIGNED( S_ );
                   
    /* assumption : the number of neighbors are less than 10 000 */
    assert( n <= MeamParameters::MAX_PARTICLE_NEIGHBORS );

    for ( size_t j=0; j<n; j++)
    {
      //alignas(onika::memory::DEFAULT_ALIGNMENT) double CoeffC[ MeamParameters::MAX_PARTICLE_NEIGHBORS ];
      
      double C = 0.0;
      int UnderCmin = 0;

      // Load data vector ij
      const double drxj=drx_[j];  //rx(i,j)
      const double dryj=dry_[j];  //ry(i,j)
      const double drzj=drz_[j];  //rz(i,j)

      // could be declared locally, as const
      const double tmprij2    = drxj*drxj + dryj*dryj + drzj*drzj;    // compute rij  square  distance
      
      assert(tmprij2!=0.);
      
      const double invtmprij2 = 1./tmprij2;        
      //invtmprij  = -1.*std::sqrt(invtmprij2);
      
      /* if UnderCmin is equal to 1, Sij =0 */
      /* We check for each interaction(i,j) if an atom k is in the neighborhood of i and j can screen this iteraction */
      /* Estimated potential speedup: 1.710 with -mavx */
      #pragma omp simd reduction(max:UnderCmin) reduction(+:C)
      for ( size_t k=0; k<n; k++) 
      {
        double CoeffC_k = 0.0;
      
        /* Avoid */
        if(j==k) {CoeffC_k=0;}
        else{
        
          const double drxk = drx_[k];  //rx(i,k)
          const double dryk = dry_[k];  //ry(i,k)
          const double drzk = drz_[k];  //rz(i,k)
          
          const double tmprik2   = drxk*drxk + dryk*dryk + drzk*drzk;  // compute ik distance (squared)
          assert(tmprik2 != 0);  
          //invtmprik = 1./std::sqrt(tmprik2);

          /* compute kj distance by vector formula kj = ki + ij= ij - ki */
          const double tmprkj2 = (drxj-drxk)*(drxj-drxk) +
                                 (dryj-dryk)*(dryj-dryk) +
                                 (drzj-drzk)*(drzj-drzk);  
          
          /* In this case, atom k does not screen the interaction between atom i and atom j */
          if(tmprkj2>=p.Rcut*p.Rcut) {CoeffC_k=0.;}
          else {
          
            assert(tmprkj2 != 0);            
            
            /* Compute tmp2=(xik-xij)^2  */
            const double tmp5  = (tmprik2-tmprkj2)*invtmprij2;          // tmpd0=Xik-Xij    
            double tmp2   = tmp5*tmp5;             

            /* tmp0=2(Xik+Xij)-(xik-xij)^2-1 numerator */
            double tmp0   = 2.*(tmprik2+tmprkj2)*invtmprij2-tmp2-1.;  
            
            /*  tmp1=-(xik-xij)^2+1 denominator */
            double tmp1   = 1.-tmp2;                               
            
            const double tmp4   = tmp1;
            
            /* tmp2 = Cijk = ( 2(Xik+Xij)-(xik-xij)^2-1  )  / (  1-(xik-xij)^2  ) */
            tmp2   = tmp0/tmp1;                      
            
            /* We keep the minimum between Cikj and Cmax to avoid (C<Cmax) condition */
            /* Indeed, if Cikj >= Cmax -> Sikj = 1                                   */
            tmp0   = min<double>(tmp2,p.Cmax);                
            double tmp3   = (tmp0-p.Cmin);

            /* In this case, Cikj < Cmin -> Sijk = 0 -> Sij =0 */
            if(tmp3*tmp3<1e-15)
            { 
              UnderCmin = 1;
              CoeffC_k = 0.;
            }
            else {
              tmp3   = 1./tmp3;
              tmp1   = (p.Cmax-tmp0)*tmp3;  // (Cmax-C)/(C-Cmin)

              /* Exceptions */
              tmp1 = (tmp2> 1.e-5 && tmp4>=1.e-5? tmp1 : 0.);                     
              
              CoeffC_k=-tmp1*tmp1;
            }
          }
        }
        
        C += CoeffC_k;
      }
      
      if(UnderCmin==1) S_[j] = 0.;
      else S_[j] = exp(C);
#     if 0
      {
        /* Compute Sij = prod_k(Sikj) */
        C = 0.;
        /* prod exp(C1)*exp(C2) ... exp(Cn) = exp(Sum_k(Ck)) */
        /* estimated potential speedup: 3.880 with -mavx */
        #pragma omp simd reduction(+:C)
        for ( size_t k=0; k<n; k++) C+=CoeffC[k];
        
        S_[j] = exp(C);
      }
#     endif

    }
  }
  
  /// @brief Calculate term of derivative of screenterm
  /// @tparam T Type of the variables
  /// @param [in] dij Distances^2
  /// @param [in] dik Distances^2
  /// @param [in] djk Distances^2
  /// @param [out] Sx X component of the derivative
  /// @param [out] Sy Y component of the derivative
  /// @param [out] Sz Z component of the derivative
  ONIKA_HOST_DEVICE_FUNC inline void screeningDerived( double dij, double dik, double djk, double &Sx, double &Sy, double &Sz) const
  {
    using onika::cuda::min;

    double  invdij,invdik,invdkj,
      tmpd0, tmpd1, tmpd2, tmpd3,
       tmp0,  tmp1,  tmp2,  tmp3,
                     bool1,bool2;

    assert(dij !=0);
    assert(dik !=0);    
    assert(djk !=0);
    
    invdij  = 1./sqrt(dij);
    invdik  = 1./sqrt(dik);
    invdkj  = 1./sqrt(djk);
      
    tmpd0   = (dik-djk)*invdij*invdij;              // tmpd0=Xik-Xij    
    tmp2    = tmpd0*tmpd0;                          // tmp2=(xik-xij)^2    
      
    /* compute Cijk = ( 2(Xik+Xij)-(xik-xij)^2-1  )  / (  1-(xik-xij)^2  ) */  
    /* tmp 0 = numerator of Cikj   */
    /* tmp 1 = denominator of Cikj */
    tmp0    = 2.*(dik+djk)*invdij*invdij-tmp2-1.; 
    tmp1    = 1.-tmp2;                             
    
    assert(tmp1 !=0 );
    
    bool1   = sqrt(tmp1*tmp1);
    tmp3    = 1./tmp1;
    tmpd2   = tmp0*tmp3; /* tmpd2 = Cikj */                        
    bool2   = tmpd2;
    tmp3    = tmp3*tmp3;

    tmpd1   = ((2.-2.*tmpd0)*tmp1+2.*tmpd0*tmp0)*tmp3;         
    tmpd3   = ((2.+2.*tmpd0)*tmp1-2.*tmpd0*tmp0)*tmp3;

    tmp1    = -2.*(tmpd1*dik+tmpd3*djk)*invdij*invdij*invdij;
    tmp2    = 2.  *tmpd1*dik*invdik*invdij*invdij;
    tmp3    = 2.  *tmpd3*djk*invdkj*invdij*invdij;
        
    /* As in the function screeningFunction, we avoid the condition Cikj > Cmax -> Sijk = 1 */   
    tmpd0   = min(tmpd2,p.Cmax);                
    tmpd1   = (tmpd0-p.Cmin);                        //tmpd1<1e-5 déjà enlevé comme possibilité


    if(bool1<1.e-5||bool2<0.) { Sx=0; Sy=0; Sz=0; }
    else {
    
      tmpd1   = 1./tmpd1;                            //tmpd1=1/(C-Cmin)
      tmpd2   = (p.Cmax-p.Cmin)*tmpd1*tmpd1;             //tmpd2=(Cmax-Cmin)/(C-Cmin)^2 //ici j'ai inversé Cmin et Cmax
      tmpd3   = (p.Cmax-tmpd0)*tmpd1;                  //tmpd3=(Cmax-C)/(C-Cmin)
    
      tmpd1   = -2 * tmpd3 * tmpd2;                 // -2 cijk * Sijk * tmpd2 // on enlève Sijk car on le divise après
    
      Sx = tmp1 * tmpd1;
      Sy = tmp2 * tmpd1;
      Sz = tmp3 * tmpd1; 
    }
  }


  // PHI
  
	/// @brief C0+C1*x+C2*x^2+C3*x^3
	/// @param [in] x Variable
	/// @param [in] C0 Zero order coefficient
	/// @param [in] C1 First order coefficient
	/// @param [in] C2 Second order coefficient
	/// @param [in] C3 Third order coefficient
	/// @return Result
 ONIKA_HOST_DEVICE_FUNC static inline double UpTo3(const double x, const double C0, const double C1, const double C2, const double C3) {
		double x2 = x*x;
		return C0+C1*x+C2*x2+C3*x*x2;
		// <=> ( ( C3*x + C2 ) * x + C1 ) * x + C0
	}


	/// @brief Derivative of : C1*x+C2*x^2+C3*x^3
	/// @param [in] x Variable
	/// @param [in] C1 First order coefficient
	/// @param [in] C2 Second order coefficient
	/// @param [in] C3 Third order coefficient
	/// @return Result
	ONIKA_HOST_DEVICE_FUNC static inline double Derivative3(const double x,
			const double C1, const double C2, const double C3) {
		return C1+2.*C2*x+3.*C3*x*x;
		// <=> ( 3*C3*x + 2*C2 ) * x + C1
	}
		
  /// @brief Calculate term of the energy and force of pair interaction
  /// @tparam T Type of the variables
  /// @param [out] fx_ X components to the forces
  /// @param [out] fy_ Y components to the forces
  /// @param [out] fz_ Z components to the forces
  /// @param [out] fSx_ X contribution of screening term to the forces
  /// @param [out] fSy_ Y contribution of screening term to the forces
  /// @param [out] fSz_ Z contribution of screening term to the forces
  /// @param [out] en_ Potential energies
  /// @param [in] drx_ X components of the distances
  /// @param [in] dry_ Y components of the distances
  /// @param [in] drz_ Z components of the distances
  /// @param [in] S_ Screening function
  /// @param [in] n Number of elements
  ONIKA_HOST_DEVICE_FUNC inline void phi (
    double * __restrict__ fx_,
    double * __restrict__ fy_,
    double * __restrict__ fz_, 
    double * __restrict__ fSx_, 
    double * __restrict__ fSy_, 
    double * __restrict__ fSz_, 
    double * __restrict__ en_,
    const double * __restrict__ drx_, 
    const double * __restrict__ dry_, 
    const double * __restrict__ drz_, 
    const double * __restrict__ S_, 
    const size_t n) const
  {

    ONIKA_ASSUME_ALIGNED( fx_ );
    ONIKA_ASSUME_ALIGNED( fy_ );
    ONIKA_ASSUME_ALIGNED( fz_ );
    ONIKA_ASSUME_ALIGNED( fSx_ );
    ONIKA_ASSUME_ALIGNED( fSy_ );
    ONIKA_ASSUME_ALIGNED( fSz_ );
    ONIKA_ASSUME_ALIGNED( en_ );
    ONIKA_ASSUME_ALIGNED( drx_ );
    ONIKA_ASSUME_ALIGNED( dry_ );
    ONIKA_ASSUME_ALIGNED( drz_ );
    ONIKA_ASSUME_ALIGNED( S_ );

    double tmpr, 
          // tmp0, 
          tmp1,
          // tmp2,
          tmp3, tmp4, tmp5, tmp6,
          drx, dry, drz, 
          drxk, dryk, drzk,
          S, Sx, Sy, Sz, 
          // Sijk,
          F , dF, Emu, dEmu,
          rhoref, drhoref,
          // zero(0.),
          // v3(3.),
          tmprij2, tmprik2, tmprkj2,
          enj, enk,
          smoothing, nosmoothing;
          
    const double iZ(1./p.Z);                  
    const double iZ_neg(-iZ);

    /* estimated potential speedup: 2.960 with -mavx */
    #pragma omp simd
    for (size_t j=0; j<n; j++) 
    {    
      // Load data
      drx = drx_[j];  //rx(i,j)
      dry = dry_[j];  //ry(i,j)
      drz = drz_[j];  //rz(i,j)
      S   = S_  [j];  //Screening term

      /* compute r(i,j)² */
      tmpr = sqrt(drx*drx + dry*dry + drz*drz); 
      
      /* (r0*delta)/r */
      tmp5 = p.delta/(tmpr*p.r0);                       
      
      /* alpha(r*1./r0 - 1) */
      tmp1 = p.alpha*(tmpr*p.r0-1.);       
      
      /* ( alpha(r*1./r0 - 1) )² */
      tmp4 = tmp1*tmp1;        
                            
      tmp6 = 1./tmpr;                               
      
      /*  E0 * exp(-alpha(r*1./r0 - 1)) */
      tmp3 = p.E0*exp(-tmp1); 

      /* We use a smoothing polynome between rcut-rp and rcut. */
      smoothing   = UpTo3(tmpr,p.PolynomeMeam[0],p.PolynomeMeam[1],p.PolynomeMeam[2],p.PolynomeMeam[3]);
      
      /*-Emu = E0 * exp(-alpha(r*1./r0 - 1)) * (1+a+delta/a^3*rc/r) */
      nosmoothing = tmp3 * (1.+tmp1+(tmp1*tmp1*tmp1)*tmp5);

      Emu  =  tmpr>p.rd? (-1)*smoothing : nosmoothing; 

      smoothing   = Derivative3(tmpr,p.PolynomeMeam[1],p.PolynomeMeam[2],p.PolynomeMeam[3]);
      nosmoothing =  Emu*p.alpha*p.r0
                    -tmp3*( p.alpha*p.r0
                           +3.*p.delta*tmp4*p.alpha*tmp6
                           -tmp5*tmp4*tmp1*tmp6);

      dEmu =  tmpr>p.rd? smoothing : nosmoothing; 

      Rhoref(rhoref,drhoref,tmpr); 
      Fe(F,dF,rhoref);

      //Store energy
      en_[j]=iZ_neg*(Emu+F)*S;  // tmp4 = phi(r) 
      
      tmp5=2.*iZ*(dEmu+iZ_neg*drhoref*dF)*tmp6*S; 
    
      /* store dphi(rij)/dri^alpha * Sij */
      fx_[j]=drx*tmp5;      
      fy_[j]=dry*tmp5;
      fz_[j]=drz*tmp5;

    }

    /* estimated potential speedup: 5.140  with -mavx*/
    #pragma omp simd
    for(size_t j=0;j<n;j++)
    {
      fSx_[j]=0.;
      fSy_[j]=0.;
      fSz_[j]=0.;  
    }

    for (size_t j=0; j<n; j++)
    {
      //moche
/*
      alignas(onika::memory::DEFAULT_ALIGNMENT) double drxk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double dryk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double drzk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double en_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockfx[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockfy[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockfz[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockSx[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockSy[MeamParameters::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double stockSz[MeamParameters::MAX_PARTICLE_NEIGHBORS];
*/
      alignas(onika::memory::DEFAULT_ALIGNMENT) uint16_t inter_nbh[MeamParameters::MAX_PARTICLE_NEIGHBORS];

      drx=drx_[j];  //rx(i,j)
      dry=dry_[j];  //ry(i,j)
      drz=drz_[j];  //rz(i,j)

      /* Compute rij  square  distance */
      tmprij2=drx*drx+dry*dry+drz*drz;            
      enj=-2*en_[j];

      //si en_[k] est différent de 0, cela veut dire que les différents Sikj sont différents de 0
      //Donc on ne divise pas par 0. Sinon la dérivée est égale à 0. 
      int indx=0;

      /* We check for each interaction(i,j) if an atom k!=j is in the neighborhood of i and j */

      for (size_t k=0; k<n; k++)
      {
        if(j==k) continue;
        else if(  (  (drx_[j]-drx_[k])*(drx_[j]-drx_[k])+(dry_[j]-dry_[k])*(dry_[j]-dry_[k])+(drz_[j]-drz_[k])*(drz_[j]-drz_[k]))  <= p.Rcut*p.Rcut)
        {
//          drxk_bis[indx]=drx_[k];
//          dryk_bis[indx]=dry_[k];
//          drzk_bis[indx]=drz_[k];
//          en_bis[indx]=en_[k];
          inter_nbh[ indx++ ] = k;
        }
      }

      double stockfx=0.0, stockfy=0.0, stockfz=0.0, stockSx=0.0, stockSy=0.0, stockSz=0.0;

      #pragma omp simd reduction(+:stockfx,stockfy,stockfz,stockSx,stockSy,stockSz)
      for (int ki=0; ki<indx; ki++) 
//      for (size_t k=0; k<n; k++)
      {
        const int k = inter_nbh[ki];
        //if( j!=k && (  (drx_[j]-drx_[k])*(drx_[j]-drx_[k])+(dry_[j]-dry_[k])*(dry_[j]-dry_[k])+(drz_[j]-drz_[k])*(drz_[j]-drz_[k]))  <= p.Rcut*p.Rcut )
        //{
          drxk = drx_[k];	//rx(i,k)
          dryk = dry_[k];	//ry(i,k)
          drzk = drz_[k];	//rz(i,k)

          /* compute ik distance */
          tmprik2= drxk*drxk+dryk*dryk+drzk*drzk;          
          
          /* Compute kj distance by vector formula kj = ki + ij= ij - ki */
          tmprkj2= (drx-drxk)*(drx-drxk) + (dry-dryk)*(dry-dryk) + (drz-drzk)*(drz-drzk); 

          assert(tmprik2!=0);
          assert(tmprkj2!=0);
          
          /* get energy of the atom k */
          enk = en_[k]; 

          /* We compute first part of dSijk/drij et dSijk/dik */  
          screeningDerived(tmprij2,tmprik2,tmprkj2, Sx, Sy, Sz);      

          Sy/=sqrt(tmprik2);
          Sx/=sqrt(tmprij2);

          stockfx  +=  enj*(drxk*Sy  + drx*Sx);
          stockfy  +=  enj*(dryk*Sy  + dry*Sx);
          stockfz  +=  enj*(drzk*Sy  + drz*Sx);

          /* We compute second part part of dSijk/drij et dSijk/dik */  
          screeningDerived(tmprik2,tmprkj2,tmprij2, Sx, Sy, Sz); 

          Sy/=sqrt(tmprkj2);
          Sz/=sqrt(tmprij2);

          stockSx  += enk*((drx-drxk)*Sy + drx*Sz);
          stockSy  += enk*((dry-dryk)*Sy + dry*Sz);
          stockSz  += enk*((drz-drzk)*Sy + drz*Sz);
        //}
        
        /*else
        {
          stockfx[k] = 0.0;
          stockfy[k] = 0.0;
          stockfz[k] = 0.0;
          stockSx[k] = 0.0;
          stockSy[k] = 0.0;
          stockSz[k] = 0.0;
        }*/
      }

      /* Sum term of temporary arrays in derived of phi*S */
      //for(size_t i=0;i<n;i++)
      //{
        fx_[j]+=stockfx;
        fy_[j]+=stockfy;
        fz_[j]+=stockfz;

        fSx_[j]-=stockSx;
        fSy_[j]-=stockSy;
        fSz_[j]-=stockSz;
      //}

    }
  }

    // Rho
    
    /// @file
  /// @brief Vectorization of MEAM potential force computation : functions related to rho
  /// @brief Compute rho_i=exp(-beta_i(rij/r0-1))
  /// @tparam T Type of the variables
  /// @param [in] coeff coeff=-beta_i(rij/r0-1)
  /// @param [in] r Distance of rij to check if it's inferior to rd=rcut-rp
  /// @param [in] n Allow to determinate which rho_i it is (to use the good smoothing polynomial)
  /// @return rho_i=exp(-beta_i(rij/r0-1))
  ONIKA_HOST_DEVICE_FUNC inline double rhoi(const double  coeff ,const double  r, size_t n) const
  {
     assert(r>0.);
     
     /* We use a smoothing polynome between rcut-rp and rcut. */
     /* else rho_i = exp(-beta_i(rij/r0-1))*/
     double smoothing = UpTo3(r,p.PolynomeMeam[(n+1)*4],p.PolynomeMeam[(n+1)*4+1],p.PolynomeMeam[(n+1)*4+2],p.PolynomeMeam[(n+1)*4+3]);
     return r>p.rd ? smoothing : exp(coeff);
  }

  /// @brief Compute generic derivative of rho_i=exp(-beta_i(rij/r0-1))
  /// @param [in] coeff  coeff = rij/r0-1
  /// @param [in] nbeta  nbeta = -beta_n
  /// @param [in] r Distance of rij to check if it's inferior to rd=rcut-rp
  /// @param [in] n Allow to determinate which rho_i it is (to use the good smoothing polynomial)
  /// @return Derivative of rho_i=exp(-beta_i(rij/r0-1))
  ONIKA_HOST_DEVICE_FUNC inline double Drhoi(const double  coeff, const double  nbeta,const double  r, size_t n) const
  {
      assert(r>0);

      //coeff= Beta(n) (r/r0-1)
      double smoothing = Derivative3(r,p.PolynomeMeam[(n+1)*4+1], p.PolynomeMeam[(n+1)*4+2], p.PolynomeMeam[(n+1)*4+3]);
      return r>p.rd ? smoothing : nbeta*p.r0*exp(nbeta*coeff);
  }



  /// @brief Compute rho^0_i term which used to compute rho terms
  /// @tparam T Type of the variables
  /// @param [in] S_ Screening function
  /// @param [in] drx_ X components of the distances
  /// @param [in] dry_ Y components of the distances
  /// @param [in] drz_ Z components of the distances
  /// @param [in] n Number of elements
  /// @return rho^0_i
  ONIKA_HOST_DEVICE_FUNC double compute_rho_0(const double* S_, const double* drx_, const double* dry_, const double* drz_, const size_t n) const
  {
    double tmp0, tmp1; 

    double resultat=0.;

    /* estimated potential speedup: 3.690 -mavx without aligned*/
    /* estimated potential speedup: 4.230 -mavx with aligned*/
    #pragma omp simd reduction(+:resultat)
    for (size_t i=0; i<n; i++)
    {
          /* compute r(i,j) */
          tmp0 = sqrt(drx_[i]*drx_[i] + dry_[i]*dry_[i] + drz_[i]*drz_[i]);  
          
          /* compure r/r0-1 */
          tmp1 = tmp0*p.r0-1.;
          
          /* get sum of  exp(-beta_i(rij/r0-1)) * Sij */
          resultat += rhoi(-p.beta0*tmp1,tmp0,0)*S_[i];
    }

    return resultat;
  }

  /// @brief Compute gamma term
  /// @tparam T Type of the variables
  /// @param [in] rho0 rho0
  /// @param [in] rho1_2 rho1^2
  /// @param [in] rho2_2 rho2^2
  /// @param [in] rho3_2 rho3^2
  /// @return Gamma
  ONIKA_HOST_DEVICE_FUNC inline double Compute_Gamma(const double rho0, const double rho1_2, const double rho2_2, const double rho3_2) const
  {  
    assert(rho0!=0);
    /* compute (t1*rho 1 + t2*rho_2 + t3*rho_3 )/rho_0² */ 
    return (p.t1*rho1_2+p.t2*rho2_2+p.t3*rho3_2)/(rho0*rho0);  
  }

  /// @brief Calculate the density contribution of a neighbor atom
  /// @tparam T Type of the variables
  /// @param [out] out_ Density contribution
  /// @param [in] rho0 rho0
  /// @param [in] rho1_2 rho1^2
  /// @param [in] rho2_2 rho2^2
  /// @param [in] rho3_2 rho3^2
  ONIKA_HOST_DEVICE_FUNC inline void rho (double &out_, const double rho0, const double rho1_2, const double rho2_2, const double rho3_2) const
  {
    /* compute (t1*rho 1 + t2*rho_2 + t3*rho_3 )/rho_0² */
    double Gamma=Compute_Gamma(rho0,rho1_2,rho2_2,rho3_2);
    
    /*                    2*rho_0                        */
    /* ------------------------------------------------- */
    /* 1 + e^{-(t1*rho 1 + t2*rho_2 + t3*rho_3 )/rho_0²} */
    out_=(2.*rho0)/(1.+exp(-Gamma));        
  }

  /// @brief Calculate the density contribution of a neighbor atom reference, used in pair term
  /// @tparam T Type of the variables
  /// @param [out] out Density contribution
  /// @param [out] outf Derivative of density contribution
  /// @param [in] dr Distance
  ONIKA_HOST_DEVICE_FUNC inline void Rhoref (double &out,double &outf, double dr) const
  {
    /* (r*1./r0 - 1) */
    double r = dr*p.r0-1.;         
    
    /* rho_0 */     
    double temp0=rhoi(-p.beta0*r,dr,0);  
    
    /* rho_1 */
    double temp1=rhoi(-p.beta1*r,dr,1);  
    
    /* rho_2 */
    double temp2=rhoi(-p.beta2*r,dr,2);  
    
    /* rho_3 */
    double temp3=rhoi(-p.beta3*r,dr,3);  

    /* rhoref = sqrt( Sum_i=0,1,2,3(s_i*rho_i²)) */
    out = sqrt( p.s0*temp0*temp0
                    +p.s1*p.t1*temp1*temp1
                    +p.s2*p.t2*temp2*temp2
                    +p.s3*p.t3*temp3*temp3);    

    /* drho_0 */
    double tmp0=Drhoi(r,-p.beta0,dr,0);  
    
    /* drho_1 */
    double tmp1=Drhoi(r,-p.beta1,dr,1);  
    
    /* drho_2 */
    double tmp2=Drhoi(r,-p.beta2,dr,2);  
    
    /* drho_3 */
    double tmp3=Drhoi(r,-p.beta3,dr,3);  
     
    /* 1/rhoref*r0*(sum i (t^i . s^i . rho_i . -beta^i . rho_i )=sum( dRhoref_i . dRho_i ) */              
    outf= out>0. ? ( p.s0*temp0*tmp0
                    +p.s1*p.t1*tmp1*temp1
                    +p.s2*p.t2*tmp2*temp2
                    +p.s3*p.t3*tmp3*temp3) /out : 0.;
  }

  /// @brief F (vectorised version)
  /// @tparam T Type of the variables
  /// @param [out] out F
  /// @param [out] outf Derivate of F
  /// @param [in] rho Density contribution
  ONIKA_HOST_DEVICE_FUNC inline void Fe (double &out,double &outf, double rho) const
  {
    rho /= p.Z;
    
    /* log */
    assert(rho>=0);
    
    /* A*Ecoh*log(sqrt( Sum_i=0,1,2,3(s_i*rho_i²))/Z)*/
    double tmp0=p.A*p.Ecoh*log(rho);

    out  = rho>0.? rho*tmp0   : 0.;
    outf = rho>0.? p.A*p.Ecoh+tmp0 : 0.; 
  }


  // DRHO

  /// @brief Compute rho0, rho1, rho2 ,rho3 and their derived functions
  /// @tparam T Type of the variables
  /// @param [in] S_ Screening function
  /// @param [in] drx_ X components of the distances
  /// @param [in] dry_ Y components of the distances
  /// @param [in] drz_ Z components of the distances
  /// @param [out] rho0_x_ X components of derivate of rho0
  /// @param [out] rho0_y_ Y components of derivate of rho0
  /// @param [out] rho0_z_ Z components of derivate of rho0
  /// @param [out] rhod0 Reduction of rho0 derivatives on all the neighbors (vec3)
  /// @param [out] rho1 Reduction of rho1 on all the neighbors (vec3)
  /// @param [out] rho1_x_ X components of derivate of rho1
  /// @param [out] rho1_y_ Y components of derivate of rho1
  /// @param [out] rho1_z_ Z components of derivate of rho1
  /// @param [out] rhod1 Reduction of rho1 derivatives on all the neighbors (vec3)
  /// @param [out] rho2 Reduction of rho2 on all the neighbors (vec3)
  /// @param [out] rho2_x_ X components of derivate of rho2
  /// @param [out] rho2_y_ Y components of derivate of rho2
  /// @param [out] rho2_z_ Z components of derivate of rho2
  /// @param [out] rhod2 Reduction of rho2 derivatives on all the neighbors (vec3)
  /// @param [out] rho3 Reduction of rho3 on all the neighbors (vec3)
  /// @param [out] rho3_x_ X components of derivate of rho3
  /// @param [out] rho3_y_ Y components of derivate of rho3
  /// @param [out] rho3_z_ Z components of derivate of rho3
  /// @param [out] rhod3 Reduction of rho3 derivatives on all the neighbors (vec3)
  /// @param [in] n Number of elements
  ONIKA_HOST_DEVICE_FUNC inline void derivedRhoS(
    const double* __restrict__ S_,
    const double* __restrict__ drx_, 
    const double* __restrict__ dry_,
    const double* __restrict__ drz_,
    double* __restrict__ rho0_x_,
    double* __restrict__ rho0_y_,  
    double* __restrict__ rho0_z_, 
    Vec3d &rhod0,
		double &rho1, 
		double* __restrict__ rho1_x_, 
		double* __restrict__ rho1_y_,  
		double* __restrict__ rho1_z_, 
		Vec3d &rhod1,
		double &rho2, 
		double* __restrict__ rho2_x_,  
		double* __restrict__ rho2_y_,  
		double* __restrict__ rho2_z_, 
		Vec3d &rhod2,
		double &rho3, 
		double* __restrict__ rho3_x_,  
		double* __restrict__ rho3_y_,  
		double* __restrict__ rho3_z_, 
		Vec3d &rhod3,
		const size_t n) const
  {
    using onika::cuda::max;
  
    ONIKA_ASSUME_ALIGNED( S_ );
    ONIKA_ASSUME_ALIGNED( drx_ );
    ONIKA_ASSUME_ALIGNED( dry_ );
    ONIKA_ASSUME_ALIGNED( drz_ );
    ONIKA_ASSUME_ALIGNED( rho0_x_ );
    ONIKA_ASSUME_ALIGNED( rho0_y_ );
    ONIKA_ASSUME_ALIGNED( rho0_z_ );
    ONIKA_ASSUME_ALIGNED( rho1_x_ );
    ONIKA_ASSUME_ALIGNED( rho1_y_ );
    ONIKA_ASSUME_ALIGNED( rho1_z_ );    
    ONIKA_ASSUME_ALIGNED( rho2_x_ );
    ONIKA_ASSUME_ALIGNED( rho2_y_ );
    ONIKA_ASSUME_ALIGNED( rho2_z_ );
    ONIKA_ASSUME_ALIGNED( rho3_x_ );
    ONIKA_ASSUME_ALIGNED( rho3_y_ );
    ONIKA_ASSUME_ALIGNED( rho3_z_ );    

    double tmp0,tmp1,tmp2,tmp3,drx,dry,drz,S,Rho0ij;
          //, v0(0.), smoothing;
    //int j;

    double x,y,z,coeff,rho1I,rho2I,rho3I;

    /* to store rho1, rho2, rho3 terms */
    double  Ta[3],Tab[7],Tabc[13];                         

    for(int a=0;a<3;a++)   Ta[a]=0.;
    for(int a=0;a<7;a++)   Tab[a]=0.;
    for(int a=0;a<13;a++)  Tabc[a]=0.;

    double /*Vecrho0,*/ xVec,yVec,zVec,rI;

    /* estimated potential speedup: 3.080 -mavx */
    #pragma omp simd
    for (size_t j=0; j<n; j++)
    {   
      drx = drx_[j];  //rx(i,j)
      dry = dry_[j];  //ry(i,j)
      drz = drz_[j];  //rz(i,j)
      S   = S_[j];    //Sij

      /* Comute r(i,j) */
      tmp3 = sqrt(drx*drx + dry*dry + drz*drz);	
      
      assert(tmp3 !=0);
      
      /* rI = 1/r */
      rI = 1./tmp3; 
      
      /* r/r0-1 */                        
      tmp1 = tmp3*p.r0-1;		  

      /* r(i,j)_{x,y,z}/r */
      xVec=drx*rI;
      yVec=dry*rI;
      zVec=drz*rI;

      /* Sij*rho_1_ij */
      tmp2=S*rhoi(-p.beta1*tmp1,tmp3,1);  

      Ta[0]  += xVec*tmp2;
      Ta[1]  += yVec*tmp2;
      Ta[2]  += zVec*tmp2;

      /* Sij*rho_2_ij */
      tmp2=S*rhoi(-p.beta2*tmp1,tmp3,2);  
   
      Tab[0] += xVec*xVec*tmp2;
      Tab[1] += xVec*yVec*tmp2;
      Tab[2] += xVec*zVec*tmp2;
      Tab[3] += yVec*yVec*tmp2;
      Tab[4] += yVec*zVec*tmp2;
      Tab[5] += zVec*zVec*tmp2;

      /* correction term rho_2 */
      Tab[6] += tmp2;

      /* Sij*rho_3_ij */  
      tmp2=S*rhoi(-p.beta3*tmp1,tmp3,3); 

      Tabc[0]  += xVec*xVec*xVec*tmp2;
      Tabc[1]  += xVec*xVec*yVec*tmp2;
      Tabc[2]  += xVec*xVec*zVec*tmp2;
      Tabc[3]  += xVec*yVec*yVec*tmp2;
      Tabc[4]  += xVec*yVec*zVec*tmp2;
      Tabc[5]  += xVec*zVec*zVec*tmp2;
      Tabc[6]  += yVec*yVec*yVec*tmp2;
      Tabc[7]  += yVec*yVec*zVec*tmp2;
      Tabc[8]  += yVec*zVec*zVec*tmp2;
      Tabc[9]  += zVec*zVec*zVec*tmp2;

      /* correction term rho_3 */
      Tabc[10] += xVec*tmp2;
      Tabc[11] += yVec*tmp2;
      Tabc[12] += zVec*tmp2;

    }

    /* compute rho_1, rho_2, rho_3 thank to previous compute */
    rho1 = Ta[0]*Ta[0]+Ta[1]*Ta[1]+Ta[2]*Ta[2];

    rho2 = Tab[0]*Tab[0]+2.*Tab[1]*Tab[1] +2.*Tab[2]*Tab[2]   
          +Tab[3]*Tab[3]+2.*Tab[4]*Tab[4 ]+   Tab[5]*Tab[5]   
          -Tab[6]*Tab[6]/3.;

    rho3 = Tabc[0]*Tabc[0]+Tabc[6]*Tabc[6]+Tabc[9]*Tabc[9]
         +3.*(Tabc[1]*Tabc[1]+Tabc[2]*Tabc[2]+Tabc[3]*Tabc[3]
         +  Tabc[5]*Tabc[5]+Tabc[7]*Tabc[7]+Tabc[8]*Tabc[8])
         +6.* Tabc[4]*Tabc[4]-3.*(Tabc[10]*Tabc[10]+Tabc[11]*Tabc[11]+Tabc[12]*Tabc[12])/5.;

    /* Stamp condition because if rho=0, ln(rho) = bug */
    rho1=max(rho1,1e-30);
    rho2=max(rho2,1e-30);
    rho3=max(rho3,1e-30);

    rho1I=1./sqrt(rho1);
    rho2I=1./sqrt(rho2);
    rho3I=1./sqrt(rho3);

    /* We obtain _Tbeta,gamma=sum_j!=i ( Sij T_ij^(3betagammamu) ) with mu,gamma,beta=x,y,z */
    #pragma omp simd
    for(int a=0;a<3;a++)   Ta[a]*=rho1I;
    
    #pragma omp simd
    for(int a=0;a<7;a++)   Tab[a]*=rho2I;
    
    #pragma omp simd
    for(int a=0;a<13;a++)  Tabc[a]*=rho3I;

    /* Will be used during the sum of Sij T^(3betagammamu)ij */
    double CoeffXXX(-Tabc[0]),CoeffXXY(-3.*Tabc[1]),CoeffXXZ(-3.*Tabc[2]),CoeffXYY(-3.*Tabc[3])
           ,CoeffXYZ(-6.*Tabc[4]),CoeffXZZ(-3.*Tabc[5]),CoeffYYY(-Tabc[6]),CoeffYYZ(-3.*Tabc[7])
           ,CoeffYZZ(-3.*Tabc[8]),CoeffZZZ(-Tabc[9]),CoeffX1(Tabc[10]*0.6),CoeffY1(Tabc[11]*0.6)
           ,CoeffZ1(Tabc[12]*0.6);       

    /* sum Sij T^(2betagamma)ij */
    double CoeffXX(-Tab[0]),CoeffXY(-2.*Tab[1]),CoeffXZ(-2.*Tab[2]),CoeffYY(-Tab[3])
           ,CoeffYZ(-2.*Tab[4]),CoeffZZ(-Tab[5]),CoeffunTier(-1.*Tab[6]/3.);  

    /* sum Sij T^(1beta)ij */
    double CoeffX(-Ta[0]),CoeffY(-Ta[1]),CoeffZ(-Ta[2]);     

    double tmpx,tmpy,tmpz,Rij,Coeff;

    /* estimated potential speedup: 3.610 -mavx */
    #pragma omp simd
    for (size_t i=0; i<n; i++)
    {
      drx = drx_[i];	//rx(i,j)
      dry = dry_[i];	//ry(i,j)
      drz = drz_[i];	//rz(i,j)
      S   = S_[i];	 //Sij

      /* compute r(i,j) */
      Rij=sqrt(drx*drx+dry*dry+drz*drz);  
      
      /* r/r0 - 1 */
      Coeff=p.r0*Rij -1.;      // Coeff 
      
      /* 1/r */
      tmp0=1./Rij;                 

      /* Begin by first derivative of rho0 */
      /* drhoij/drij * 1/rij */
      Rho0ij = S*tmp0*Drhoi(Coeff,-p.beta0,Rij,0);    

      /* store -comp_rin[j] * drhoa0  * Screening term[i][n] == r^alpha/rij* drho0_ij*Sij */
      rho0_x_[i] = drx*Rho0ij; 
      rho0_y_[i] = dry*Rho0ij; 
      rho0_z_[i] = drz*Rho0ij;         

       /* tmp0 = 1/r */
      drx=tmp0*drx;   //r^x_ij
      dry=tmp0*dry;   //r^y_ij
      drz=tmp0*drz;   //r^z_ij
      
      /* derivative of rho1 */
      /* rho_1/r */
      tmp2=rhoi(-p.beta1*Coeff,Rij,1)*tmp0 ;       

      /* rho1ij/rij - d(rho)/drin -> -beta1/r0*exp(-beta1 (rij/r0-1) */
      tmp1=tmp2-Drhoi(Coeff,-p.beta1,Rij,1); 

      /* warning d(r^beta(in))/dr^alpha_i == (comp_rin[j]*comp_rin[k] - 1)* distance_inI; */
      /* if alpha=beta, else comp_rin[j]*comp_rin[k]* distance_inI */
      
      tmpx=CoeffX*S*(drx*drx*tmp1-tmp2);
      tmpy=CoeffX*S*dry*drx*tmp1;
      tmpz=CoeffX*S*drz*drx*tmp1;

      tmpx+=CoeffY*S*drx*dry*tmp1;
      tmpy+=CoeffY*S*(dry*dry*tmp1-tmp2);
      tmpz+=CoeffY*S*drz*dry*tmp1;

      tmpx+=CoeffZ*S*drx*drz*tmp1;
      tmpy+=CoeffZ*S*dry*drz*tmp1;
      tmpz+=CoeffZ*S*(drz*drz*tmp1-tmp2);

      //store 
      rho1_x_[i] = tmpx;
      rho1_y_[i] = tmpy;
      rho1_z_[i] = tmpz;        

      /* rho_2 derived */
      /* rho_2_ij/r */
      tmp2=rhoi(-p.beta2*Coeff,Rij,2)*tmp0;    
      
      /*    d(rho_2_ij)/drin */
      tmp3=Drhoi(Coeff,-p.beta2,Rij,2); 
      
      /* 2*rho_2_ij - d(rho_2_ij)/drin*/
      tmp1=2.*tmp2 - tmp3;                   

      tmp3*=CoeffunTier;

      /* d(r^beta(in))/dr^alpha_i == 1 if alpha=beta, else 0 */
      /* rho2.cpp + smoothing */
      tmpx=CoeffXX*(drx*drx*drx*tmp1-drx*tmp2-drx*tmp2);
      tmpy=CoeffXX*dry*drx*drx*tmp1;
      tmpz=CoeffXX*drz*drx*drx*tmp1;

      tmpx+=CoeffXY*(drx*drx*dry*tmp1-dry*tmp2);
      tmpy+=CoeffXY*(dry*drx*dry*tmp1-drx*tmp2);
      tmpz+=CoeffXY*drz*drx*dry*tmp1;

      tmpx+=CoeffXZ*(drx*drx*drz*tmp1-drz*tmp2);
      tmpy+=CoeffXZ*dry*drx*drz*tmp1;
      tmpz+=CoeffXZ*(drz*drx*drz*tmp1-drx*tmp2);

      tmpx+=CoeffYY*drx*dry*dry*tmp1;
      tmpy+=CoeffYY*(dry*dry*dry*tmp1-dry*tmp2-dry*tmp2);
      tmpz+=CoeffYY*drz*dry*dry*tmp1;

      tmpx+=CoeffYZ*drx*dry*drz*tmp1;
      tmpy+=CoeffYZ*(dry*dry*drz*tmp1-drz*tmp2);
      tmpz+=CoeffYZ*(drz*dry*drz*tmp1-dry*tmp2);

      tmpx+=CoeffZZ*drx*drz*drz*tmp1;
      tmpy+=CoeffZZ*dry*drz*drz*tmp1;
      tmpz+=CoeffZZ*(drz*drz*drz*tmp1-drz*tmp2-drz*tmp2);

      tmpx+=drx*tmp3;
      tmpy+=dry*tmp3;
      tmpz+=drz*tmp3;

      /* store part of derivative of rho2 */
      rho2_x_[i] = S*tmpx;
      rho2_y_[i] = S*tmpy;
      rho2_z_[i] = S*tmpz;


      /* finaly rho3 */
      tmp3=Drhoi(Coeff,-p.beta3,Rij,3);
      
      /* rho_3_ij/r */
      tmp2=rhoi(-p.beta3*Coeff,Rij,3)*tmp0; 
      
      /* d(rho3)/drin */
      tmp1=3.*tmp2-tmp3;                     

      tmp3=(tmp2-tmp3);

      /* code rho3.cpp */
      tmpx=CoeffXXX*(drx*drx*drx*drx*(tmp1)-drx*drx*tmp2-drx*drx*tmp2-drx*drx*tmp2);
      tmpy=CoeffXXX*(dry*drx*drx*drx*(tmp1));
      tmpz=CoeffXXX*(drz*drx*drx*drx*(tmp1));

      tmpx+=CoeffXXY*(drx*drx*drx*dry*(tmp1)-drx*dry*tmp2-drx*dry*tmp2);
      tmpy+=CoeffXXY*(dry*drx*drx*dry*(tmp1)-drx*drx*tmp2);
      tmpz+=CoeffXXY*(drz*drx*drx*dry*(tmp1));

      tmpx+=CoeffXXZ*(drx*drx*drx*drz*(tmp1)-drx*drz*tmp2-drx*drz*tmp2);
      tmpy+=CoeffXXZ*(dry*drx*drx*drz*(tmp1));
      tmpz+=CoeffXXZ*(drz*drx*drx*drz*(tmp1)-drx*drx*tmp2);

      tmpx+=CoeffXYY*(drx*drx*dry*dry*(tmp1)-dry*dry*tmp2);
      tmpy+=CoeffXYY*(dry*drx*dry*dry*(tmp1)-drx*dry*tmp2-drx*dry*tmp2);
      tmpz+=CoeffXYY*(drz*drx*dry*dry*(tmp1));

      tmpx+=CoeffXYZ*(drx*drx*dry*drz*(tmp1)-dry*drz*tmp2);
      tmpy+=CoeffXYZ*(dry*drx*dry*drz*(tmp1)-drx*drz*tmp2);
      tmpz+=CoeffXYZ*(drz*drx*dry*drz*(tmp1)-drx*dry*tmp2);

      tmpx+=CoeffXZZ*(drx*drx*drz*drz*(tmp1)-drz*drz*tmp2);
      tmpy+=CoeffXZZ*(dry*drx*drz*drz*(tmp1));
      tmpz+=CoeffXZZ*(drz*drx*drz*drz*(tmp1)-drx*drz*tmp2-drx*drz*tmp2);

      tmpx+=CoeffYYY*(drx*dry*dry*dry*(tmp1));
      tmpy+=CoeffYYY*(dry*dry*dry*dry*(tmp1)-dry*dry*tmp2-dry*dry*tmp2-dry*dry*tmp2);
      tmpz+=CoeffYYY*(drz*dry*dry*dry*(tmp1));

      tmpx+=CoeffYYZ*(drx*dry*dry*drz*(tmp1));
      tmpy+=CoeffYYZ*(dry*dry*dry*drz*(tmp1)-dry*drz*tmp2-dry*drz*tmp2);
      tmpz+=CoeffYYZ*(drz*dry*dry*drz*(tmp1)-dry*dry*tmp2);

      tmpx+=CoeffYZZ*(drx*dry*drz*drz*(tmp1));
      tmpy+=CoeffYZZ*(dry*dry*drz*drz*(tmp1)-drz*drz*tmp2);
      tmpz+=CoeffYZZ*(drz*dry*drz*drz*(tmp1)-dry*drz*tmp2-dry*drz*tmp2);

      tmpx+=CoeffZZZ*(drx*drz*drz*drz*(tmp1));
      tmpy+=CoeffZZZ*(dry*drz*drz*drz*(tmp1));
      tmpz+=CoeffZZZ*(drz*drz*drz*drz*(tmp1)-drz*drz*tmp2-drz*drz*tmp2-drz*drz*tmp2);

      //Correctif terms
      tmpx+=CoeffX1*(drx*drx*tmp3-tmp2);
      tmpy+=CoeffX1*(dry*drx*tmp3);
      tmpz+=CoeffX1*(drz*drx*tmp3);

      tmpx+=CoeffY1*(drx*dry*tmp3);
      tmpy+=CoeffY1*(dry*dry*tmp3-tmp2);
      tmpz+=CoeffY1*(drz*dry*tmp3);

      tmpx+=CoeffZ1*(drx*drz*tmp3);
      tmpy+=CoeffZ1*(dry*drz*tmp3);
      tmpz+=CoeffZ1*(drz*drz*tmp3-tmp2);

      //store 
      rho3_x_[i]=S*tmpx;
      rho3_y_[i]=S*tmpy;
      rho3_z_[i]=S*tmpz;
    }

    for (size_t j=0; j<n; j++)
    {
      rhod0.x+=rho0_x_[j];
      rhod0.y+=rho0_y_[j];
      rhod0.z+=rho0_z_[j];

      rhod1.x+=rho1_x_[j];
      rhod1.y+=rho1_y_[j];
      rhod1.z+=rho1_z_[j];

      rhod2.x+=rho2_x_[j];
      rhod2.y+=rho2_y_[j];
      rhod2.z+=rho2_z_[j];

      rhod3.x+=rho3_x_[j];
      rhod3.y+=rho3_y_[j];
      rhod3.z+=rho3_z_[j];
    }

    double drxk, dryk, drzk, xk, yk, zk, temp;
    double drxj, dryj, drzj, xj, yj, zj;
    double rij,rik,rkj,rij2,rik2,rkj2/*,rijI*/,rkiI, rkjI;
    double rho0ij,rho1ij,rho2ij,rho3ij;
    double rho0ik,rho1ik,rho2ik,rho3ik;
    double Sk,Sij, Sik, Skj,v3(3.),v6(6.),v3inv5(0.6); //-> dS x,y,z 
    double dSx,dSy,dSz;
    double rij_seq,rij2_seq,rijI_seq;;

    alignas(onika::memory::DEFAULT_ALIGNMENT) uint16_t indx_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
/*
    alignas(onika::memory::DEFAULT_ALIGNMENT) double drzk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    alignas(onika::memory::DEFAULT_ALIGNMENT) double drxk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    alignas(onika::memory::DEFAULT_ALIGNMENT) double dryk_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS]; 
    alignas(onika::memory::DEFAULT_ALIGNMENT) double S_bis[MeamParameters::MAX_PARTICLE_NEIGHBORS];
*/    
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_0_0; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_0_1; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_0_2; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_1_0; //[MeamParameters::MAX_PARTICLE_NEIGHBORS]; 
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_1_1; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_1_2; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_2_0; //[MeamParameters::MAX_PARTICLE_NEIGHBORS]; 
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_2_1; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_2_2; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_3_0; //[MeamParameters::MAX_PARTICLE_NEIGHBORS]; 
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_3_1; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    /*alignas(onika::memory::DEFAULT_ALIGNMENT)*/ double stock_rho_3_2; //[MeamParameters::MAX_PARTICLE_NEIGHBORS];
    
    double tmpr0x(0.), tmpr0y(0.), tmpr0z(0.), 
           tmpr1x(0.), tmpr1y(0.), tmpr1z(0.), 
           tmpr2x(0.), tmpr2y(0.), tmpr2z(0.), 
           tmpr3x(0.), tmpr3y(0.), tmpr3z(0.); 
           
    /* Compute derivative of screening function part */
    for (size_t j=0; j<n; j++)
    {    
      /* rij */
      rij2_seq= drx_[j]*drx_[j]
               +dry_[j]*dry_[j]
               +drz_[j]*drz_[j];   
      
      assert(rij2_seq!=0);
      
      rij_seq=sqrt(rij2_seq);
      coeff=rij_seq*p.r0-1;
      
      rijI_seq=1./rij_seq;

      rij=rij_seq;
      rij2=rij2_seq;
      //rijI=rijI_seq;

      /* r(i,j)_{x,y,z}/rij */
      x=drx_[j]*rijI_seq;
      y=dry_[j]*rijI_seq;
      z=drz_[j]*rijI_seq;

      rho0ij=rhoi(-p.beta0*coeff,rij,0);          //rho0_ij
      rho0ij=rho0ij*S_[j];

      temp=(Ta[0]*x+Ta[1]*y+Ta[2]*z);
      rho1ij=temp*rhoi(-p.beta1*coeff,rij,1)*S_[j]; //rho1_ij/rho1

      temp=Tab[0]*x*x+2.*Tab[1]*x*y+2.*Tab[2]*x*z
        +Tab[3]*y*y+2.*Tab[4]*y*z
        +Tab[5]*z*z-Tab[6]/3.; 

      rho2ij=temp*rhoi(-p.beta2*coeff,rij,2)*S_[j]; //rho2_ij/rho2

      temp=Tabc[0]*x*x*x+3.*Tabc[1]*x*x*y+3.*Tabc[2]*x*x*z
        +3.*Tabc[3]*x*y*y+6.*Tabc[4]*x*y*z+3.*Tabc[5]*x*z*z+Tabc[6]*y*y*y
        +3.*Tabc[7]*y*y*z+3.*Tabc[8]*y*z*z+Tabc[9]*z*z*z
        -0.6*(Tabc[10]*x+Tabc[11]*y+Tabc[12]*z);

      rho3ij=temp*rhoi(-p.beta3*coeff,rij,3)*S_[j];           //rho3_ij/rho3

      xj=x;
      yj=y;
      zj=z;

      drxj=drx_[j];
      dryj=dry_[j];
      drzj=drz_[j];

      size_t indx=0;


      double tmpr0jx(0.), tmpr0jy(0.), tmpr0jz(0.), 
             tmpr1jx(0.), tmpr1jy(0.), tmpr1jz(0.), 
             tmpr2jx(0.), tmpr2jy(0.), tmpr2jz(0.), 
             tmpr3jx(0.), tmpr3jy(0.), tmpr3jz(0.);        

      //sort of neighbors. 
      for (size_t k=j+1; k<n; k++)
      {
        if(  (  (drx_[j]-drx_[k])*(drx_[j]-drx_[k])+(dry_[j]-dry_[k])*(dry_[j]-dry_[k])+(drz_[j]-drz_[k])*(drz_[j]-drz_[k]))  <= p.Rcut*p.Rcut)
        {
          //drxk_bis[indx]=drx_[k];
          //dryk_bis[indx]=dry_[k];
          //drzk_bis[indx]=drz_[k];
          indx_bis[indx]=k;
          //S_bis[indx]= S_[k];
          indx++;
        }

      }

      /* estimated potential speedup: 3.280 -mavx */
      #pragma omp simd reduction(+:tmpr0jx, tmpr0jy, tmpr0jz, tmpr1jx, tmpr1jy, tmpr1jz, tmpr2jx, tmpr2jy, tmpr2jz, tmpr3jx, tmpr3jy, tmpr3jz,  \
           tmpr0x, tmpr0y, tmpr0z, tmpr1x, tmpr1y, tmpr1z, tmpr2x, tmpr2y, tmpr2z, tmpr3x, tmpr3y, tmpr3z)     
      for (size_t k=0; k<indx; k++)
      {   
        /* load data about atom k */
        drxk = drx_[indx_bis[k]]; //drxk_bis[k];	
        dryk = dry_[indx_bis[k]]; //dryk_bis[k];	
        drzk = drz_[indx_bis[k]]; //drzk_bis[k];
        Sk = S_[indx_bis[k]]; //S_bis[k];

        /* compute distances r(i,k),r(i,k)²,1./r(i,k) */
        rik2=drxk*drxk+dryk*dryk+drzk*drzk;
        
        assert(rik2 != 0);
        
        rik=sqrt(rik2);
        rkiI=1./rik;

        Coeff=rik*p.r0-1.;

        rkj2=(drxj-drxk)*(drxj-drxk)+(dryj-dryk)*(dryj-dryk)+(drzj-drzk)*(drzj-drzk); 
        
        assert(rkj2 != 0);
        
        rkj=sqrt(rkj2); 
        rkjI=1./rkj;

        xk=drxk*rkiI;
        yk=dryk*rkiI;
        zk=drzk*rkiI;

        /* Compute rhoik, k=0,1,2,3 */

        /* rho0_ik */
        rho0ik=rhoi(-p.beta0*Coeff,rik,0)*Sk;                            

        //rho1_ik
        rho1ik= (Ta[0]*xk+Ta[1]*yk+Ta[2]*zk) *rhoi(-p.beta1*Coeff,rik,1)*Sk; 
        
        //rho2_ik/rho2 
        rho2ik=Tab[0]*xk*xk+2.*Tab[1]*xk*yk+2.*Tab[2]*xk*zk
            +Tab[3]*yk*yk+2.*Tab[4]*yk*zk+Tab[5]*zk*zk-(1./3.)*Tab[6];      

        rho2ik=rho2ik*rhoi(-p.beta2*Coeff,rik,2)*Sk;

        /* rho3_ik/rho3 */ 
        rho3ik=Tabc[0]*xk*xk*xk+v3*Tabc[1]*xk*xk*yk
            +v3*Tabc[2]*xk*xk*zk+v3*Tabc[3]*xk*yk*yk
            +v6*Tabc[4]*xk*yk*zk+v3*Tabc[5]*xk*zk*zk
            +Tabc[6]*yk*yk*yk+v3*Tabc[7]*yk*yk*zk+v3*Tabc[8]*yk*zk*zk
            +Tabc[9]*zk*zk*zk-v3inv5*(Tabc[10]*xk+Tabc[11]*yk+Tabc[12]*zk);

        rho3ik=rho3ik*rhoi(-p.beta3*Coeff,rik,3)*Sk;

        /* this function use distance² */
        screeningDerived(rij2,rik2,rkj2, Sij, Sik, Skj);			

        dSx=( xj * Sij + xk  * Sik );
        dSy=( yj * Sij + yk  * Sik );
        dSz=( zj * Sij + zk  * Sik );

        tmpr0x -= rho0ij * dSx;
        tmpr0y -= rho0ij * dSy;
        tmpr0z -= rho0ij * dSz;

        tmpr1x -= rho1ij * dSx;
        tmpr1y -= rho1ij * dSy;
        tmpr1z -= rho1ij * dSz;

        tmpr2x -= rho2ij * dSx;
        tmpr2y -= rho2ij * dSy;
        tmpr2z -= rho2ij * dSz;
        
        tmpr3x -= rho3ij * dSx;
        tmpr3y -= rho3ij * dSy;
        tmpr3z -= rho3ij * dSz;

        Skj=Skj*rkjI;

        dSx= ( xj * Sij - (drxk-drxj)  * Skj );     
        dSy= ( yj * Sij - (dryk-dryj)  * Skj );
        dSz= ( zj * Sij - (drzk-drzj)  * Skj );

        tmpr0jx -= rho0ij * dSx; 
        tmpr0jy -= rho0ij * dSy;
        tmpr0jz -= rho0ij * dSz;

        tmpr1jx -= rho1ij * dSx; 
        tmpr1jy -= rho1ij * dSy;
        tmpr1jz -= rho1ij * dSz;

        tmpr2jx -= rho2ij * dSx; 
        tmpr2jy -= rho2ij * dSy;
        tmpr2jz -= rho2ij * dSz;

        tmpr3jx -= rho3ij * dSx; 
        tmpr3jy -= rho3ij * dSy;
        tmpr3jz -= rho3ij * dSz;

        dSx= ( xk * Sik + (drxk-drxj)  * Skj);
        dSy= ( yk * Sik + (dryk-dryj)  * Skj);
        dSz= ( zk * Sik + (drzk-drzj)  * Skj);

        stock_rho_0_0/*[k]*/=rho0ij * dSx;
        stock_rho_1_0/*[k]*/=rho1ij * dSx;
        stock_rho_2_0/*[k]*/=rho2ij * dSx;
        stock_rho_3_0/*[k]*/=rho3ij * dSx;
        
        stock_rho_0_1/*[k]*/=rho0ij * dSy;
        stock_rho_1_1/*[k]*/=rho1ij * dSy;
        stock_rho_2_1/*[k]*/=rho2ij * dSy;
        stock_rho_3_1/*[k]*/=rho3ij * dSy;
        
        stock_rho_0_2/*[k]*/=rho0ij * dSz;
        stock_rho_1_2/*[k]*/=rho1ij * dSz;
        stock_rho_2_2/*[k]*/=rho2ij * dSz;
        stock_rho_3_2/*[k]*/=rho3ij * dSz;
 
        /* this function use with distance² */
        screeningDerived(rik2,rij2,rkj2, Sik, Sij, Skj);			              
        Skj=Skj*rkjI;
          
        dSx= ( xj * Sij - (drxk-drxj)  * Skj);
        dSy= ( yj * Sij - (dryk-dryj)  * Skj);
        dSz= ( zj * Sij - (drzk-drzj)  * Skj);

        tmpr0jx -= rho0ik * dSx;
        tmpr0jy -= rho0ik * dSy;
        tmpr0jz -= rho0ik * dSz;

        tmpr1jx -= rho1ik * dSx;
        tmpr1jy -= rho1ik * dSy;
        tmpr1jz -= rho1ik * dSz;

        tmpr2jx -= rho2ik * dSx;
        tmpr2jy -= rho2ik * dSy;
        tmpr2jz -= rho2ik * dSz;

        tmpr3jx -= rho3ik * dSx;
        tmpr3jy -= rho3ik * dSy;
        tmpr3jz -= rho3ik * dSz;

        dSx=( xj * Sij + xk  * Sik );
        dSy=( yj * Sij + yk  * Sik );
        dSz=( zj * Sij + zk  * Sik );

        tmpr0x -= rho0ik * dSx;
        tmpr0y -= rho0ik * dSy;
        tmpr0z -= rho0ik * dSz;

        tmpr1x -= rho1ik * dSx;
        tmpr1y -= rho1ik * dSy;
        tmpr1z -= rho1ik * dSz;

        tmpr2x -= rho2ik * dSx;
        tmpr2y -= rho2ik * dSy;
        tmpr2z -= rho2ik * dSz;
       
        tmpr3x -= rho3ik * dSx;
        tmpr3y -= rho3ik * dSy;
        tmpr3z -= rho3ik * dSz;

        dSx= ( xk * Sik + (drxk-drxj)  * Skj);
        dSy= ( yk * Sik + (dryk-dryj)  * Skj);
        dSz= ( zk * Sik + (drzk-drzj)  * Skj);
        
        stock_rho_0_0/*[k]*/+=rho0ik * dSx;
        stock_rho_1_0/*[k]*/+=rho1ik * dSx;
        stock_rho_2_0/*[k]*/+=rho2ik * dSx;
        stock_rho_3_0/*[k]*/+=rho3ik * dSx;
        
        stock_rho_0_1/*[k]*/+=rho0ik * dSy;
        stock_rho_1_1/*[k]*/+=rho1ik * dSy;
        stock_rho_2_1/*[k]*/+=rho2ik * dSy;
        stock_rho_3_1/*[k]*/+=rho3ik * dSy;
        
        stock_rho_0_2/*[k]*/+=rho0ik * dSz;
        stock_rho_1_2/*[k]*/+=rho1ik * dSz;
        stock_rho_2_2/*[k]*/+=rho2ik * dSz;
        stock_rho_3_2/*[k]*/+=rho3ik * dSz;
/*
      }
      
      for (size_t k=0; k<indx; k++)
      {
*/
        const int indxk = indx_bis[k];
        rho0_x_[indxk]-=stock_rho_0_0; //[k];
        rho0_y_[indxk]-=stock_rho_0_1; //[k];
        rho0_z_[indxk]-=stock_rho_0_2; //[k];
        rho1_x_[indxk]-=stock_rho_1_0; //[k];
        rho1_y_[indxk]-=stock_rho_1_1; //[k];
        rho1_z_[indxk]-=stock_rho_1_2; //[k];
        rho2_x_[indxk]-=stock_rho_2_0; //[k];
        rho2_y_[indxk]-=stock_rho_2_1; //[k];
        rho2_z_[indxk]-=stock_rho_2_2; //[k];
        rho3_x_[indxk]-=stock_rho_3_0; //[k];
        rho3_y_[indxk]-=stock_rho_3_1; //[k];
        rho3_z_[indxk]-=stock_rho_3_2; //[k];
     }
           
     rho0_x_[j] += tmpr0jx; rho0_y_[j]+= tmpr0jy; rho0_z_[j] += tmpr0jz; 
     rho1_x_[j] += tmpr1jx; rho1_y_[j]+= tmpr1jy; rho1_z_[j] += tmpr1jz; 
     rho2_x_[j] += tmpr2jx; rho2_y_[j]+= tmpr2jy; rho2_z_[j] += tmpr2jz; 
     rho3_x_[j] += tmpr3jx; rho3_y_[j]+= tmpr3jy; rho3_z_[j] += tmpr3jz;       
    }
    
     rhod0.x += tmpr0x; rhod0.y+= tmpr0y; rhod0.z += tmpr0z; 
     rhod1.x += tmpr1x; rhod1.y+= tmpr1y; rhod1.z += tmpr1z; 
     rhod2.x += tmpr2x; rhod2.y+= tmpr2y; rhod2.z += tmpr2z; 
     rhod3.x += tmpr3x; rhod3.y+= tmpr3y; rhod3.z += tmpr3z; 
  }

  /// @brief Compute derived function of rho_i=0..3 and Gamma
  /// @tparam T Type of the variables
  /// @param [in] rho0 rh0 value
  /// @param [in] rho1 rho1 value
  /// @param [in] rho2 rho2 value
  /// @param [in] rho3 rho3 value
  /// @param [out] tmp0 rho0 derived function
  /// @param [out] tmp1 rho1 derived function
  /// @param [out] tmp2 rho2 derived function
  /// @param [out] tmp3 rho3 derived function
  /// @param [out] tmpG Gamma
  /// @param [out] tmpdG Derived Gamma
  ONIKA_HOST_DEVICE_FUNC void getCoeffDerivateRho(const double rho0,const double rho1,const double rho2,const double rho3, double& tmp0, double& tmp1, double& tmp2, double& tmp3, double& tmpG, double& tmpdG ) const
  {
    double Gamma,rho0i,r1,r2,r3,r12,r22,r32;
    
    Gamma=Compute_Gamma(rho0,rho1,rho2,rho3);

    tmp1  =   exp(-Gamma);
    tmp2  =   1.+tmp1;
    tmp0  =   1./tmp2;
    tmpG  =   2.*tmp0;                //drho0* Gamma *drho0
    tmpdG   =   tmpG*tmp1*tmp0;       //d (G(gamma))/ d gamma

    rho0i   =   1./rho0;
    r1    =   sqrt(rho1) * rho0i;
    r2    =   sqrt(rho2) * rho0i;
    r3    =   sqrt(rho3) * rho0i;

    r12   =   r1 * r1;
    r22   =   r2 * r2;	
    r32   =   r3 * r3;
  	
    if(rho0i==0.) ONIKA_CU_ABORT(); //exit(0);

    tmp0  =   (-2. * rho0i * (p.t1 * r12 + p.t2 * r22 + p.t3 * r32)  );// *  drho0;	
    tmp1  =   (2 * p.t1 * rho0i * r1);    //  *   drho1;	
    tmp2  =   (2 * p.t2 * rho0i * r2);    //  *   drho2;	
    tmp3  =   (2 * p.t3 * rho0i * r3);    //  *   drho3;	
  }


  /// @brief Acessor to the parameters
  ONIKA_HOST_DEVICE_FUNC const MeamParameters& getParameters() { return p; }

  // the only member
  MeamParameters p; ///< Parameters
};


