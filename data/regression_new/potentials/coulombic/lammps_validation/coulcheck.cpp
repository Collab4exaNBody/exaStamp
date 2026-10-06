// Standalone re-implementation of src/potential/coulombic kernels EXACTLY as coded,
// in exaStamp internal units (Å, Da, ps, e), compared with LAMMPS metal dumps.
// usage: coulcheck <wolf|dsf|ewald> <lammps_dump> [kmax_user] [g_ewald_user] [fix_kmax(0/1)] [use_codata_const(0/1)]
// e.g.  g++ -O2 coulcheck.cpp -o coulcheck && ./coulcheck ewald dump_ewald_fixed.10.txt 8 0.3 1 1
// (kernels as coded before the 2026-10 fixes : fix_kmax=0 / use_codata_const=0 reproduce the old Ewald bugs)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <complex>

// ---- unit system (onika constants) ----
static const double e_C   = 1.602176634e-19;
static const double amu   = 1.66053906892e-27;
static const double eps0  = 8.8541878188e-12;
static const double E_int = amu * 1e-20 / 1e-24;       // J per internal energy unit
static const double eV    = e_C / E_int;               // eV -> internal
static const double bar_per_eVA3 = 1.602176634e6;      // 1 eV/Å^3 = 1.602e6 bar
// epsilon0 in internal units: C^2 s^2 m^-3 kg^-1 -> e^2 ps^2 Å^-3 Da^-1
static const double eps0_int = eps0 * (1e-30 * amu) / (e_C * e_C * 1e-24);
// hex constants from coulombic/ewald.h
static const double hex_epsilonZero = 0x1.337f14f782bb5p-21;
static const double hex_fpe0        = 0x1.e303a50b0d1ap-18;

struct Atom { long id; int type; double x[3]; double f[3]; double pe; double st[6]; };

static std::vector<Atom> read_dump(const char* fn, double L[3])
{
  std::ifstream in(fn); std::string line; std::vector<Atom> a; long n=0;
  while( std::getline(in,line) ) {
    if( line.rfind("ITEM: NUMBER",0)==0 ) { std::getline(in,line); n=atol(line.c_str()); }
    else if( line.rfind("ITEM: BOX",0)==0 ) { for(int d=0;d<3;d++){ std::getline(in,line); double lo,hi; sscanf(line.c_str(),"%lf %lf",&lo,&hi); L[d]=hi-lo; } }
    else if( line.rfind("ITEM: ATOMS",0)==0 ) {
      a.resize(n);
      for(long i=0;i<n;i++){ std::getline(in,line); std::istringstream ss(line); Atom& t=a[i];
        ss>>t.id>>t.type>>t.x[0]>>t.x[1]>>t.x[2]>>t.f[0]>>t.f[1]>>t.f[2]>>t.pe; for(int k=0;k<6;k++) ss>>t.st[k]; }
      break;
    }
  }
  return a;
}

// ---- kernels, copied from coulombic/*.h (return e, de in internal units) ----
struct WolfP { double alpha, rc, qqrd2e=14.399645, e_shift, f_shift; };
static void wolf_compute_energy(const WolfP& p, double c, double r, double& e, double& de)
{
  const double prefactor = p.qqrd2e * c / r;
  const double r2 = r*r; const double alpha2 = p.alpha*p.alpha;
  const double erfcc = erfc(p.alpha*r); const double erfcd = exp(-alpha2*r2);
  const double v_sh = (erfcc - p.e_shift*r) * prefactor;
  e = v_sh * eV;
  const double dvdrr = (erfcc/r2 + 2.0*p.alpha/sqrt(M_PI)*erfcd/r) + p.f_shift;
  const double forcecoul = dvdrr*r2*prefactor;
  const double fpair = -forcecoul/r;
  de = fpair * eV;   // eV/ang -> internal (length unit is ang)
}
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429
static void dsf_compute_energy(const WolfP& p, double c, double r, double& e, double& de)
{
  double MY_PIS = sqrt(M_PI); double rsq=r*r; double cut_coul=p.rc; double cut_coulsq=cut_coul*cut_coul;
  double erfcc = erfc(p.alpha*cut_coul); double erfcd = exp(-p.alpha*p.alpha*cut_coul*cut_coul);
  double f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*p.alpha*erfcd/cut_coul);
  double e_shift = erfcc/cut_coul - f_shift*cut_coul;
  double prefactor = p.qqrd2e*c/r;
  erfcd = exp(-p.alpha*p.alpha*rsq);
  double t = 1.0/(1.0+EWALD_P*p.alpha*r);
  erfcc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*erfcd;
  double forcecoul = prefactor*(erfcc/r + 2.0*p.alpha/MY_PIS*erfcd + r*f_shift)*r;
  double fpair = -forcecoul/r;
  double ecoul = prefactor*(erfcc - r*e_shift - rsq*f_shift);
  e = ecoul*eV; de = fpair*eV;
}
struct EwP { double g_ewald, gm_sr, bt_sr; };
static void ewald_compute_energy(const EwP& p, double c, double r, double& e, double& de)
{
  const double cf = p.gm_sr*c; const double grij = p.g_ewald*r; const double expm2 = exp(-grij*grij);
  const double t = 1.0/(1.0+EWALD_P*grij);
  const double erfcv = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2;
  e = cf*erfcv/r; de = -(cf*p.bt_sr*expm2 + e)/r;
}

int main(int argc, char** argv)
{
  if(argc<3){ fprintf(stderr,"usage\n"); return 1; }
  std::string mode = argv[1];
  long kmax_user = argc>3 ? atol(argv[3]) : 0;
  double g_user = argc>4 ? atof(argv[4]) : 0.0;
  bool fix_kmax = argc>5 && atoi(argv[5]);
  bool codata = argc>6 && atoi(argv[6]);
  double L[3]; auto at = read_dump(argv[2], L);
  const long N = at.size(); const double V = L[0]*L[1]*L[2];
  std::vector<double> q(N); double qsq=0, qsum=0;
  for(long i=0;i<N;i++){ q[i] = at[i].type==1 ? 2.2208 : -1.1104; qsq+=q[i]*q[i]; qsum+=q[i]; }

  printf("units: eV->int=%.15g  eps0_int=%.15g hex_eps0=%.15g (rel %.2e)  4pi eps0_int=%.15g hex_fpe0=%.15g (rel %.2e)\n",
    eV, eps0_int, hex_epsilonZero, hex_epsilonZero/eps0_int-1, 4*M_PI*eps0_int, hex_fpe0, hex_fpe0/(4*M_PI*eps0_int)-1);
  printf("1/(4pi eps0) in eV.A/e^2 from hex = %.10f (LAMMPS qqr2e 14.399645)\n", 1.0/(4*M_PI*hex_epsilonZero)/eV);
  printf("N=%ld V=%.6f qsum=%.3e qsq=%.6f\n",N,V,qsum,qsq);

  const double rc = 10.0;
  WolfP wp{0.2, rc}; wp.e_shift = erfc(wp.alpha*wp.rc)/wp.rc;
  wp.f_shift = -(wp.e_shift + 2.0*wp.alpha/sqrt(M_PI)*exp(-wp.alpha*wp.alpha*wp.rc*wp.rc))/wp.rc;
  EwP ep{};

  // ---------- Ewald init exactly as coulombic/ewald.h ----------
  long kxmax=0,kymax=0,kzmax=0,kmax=kmax_user; double g_ewald=0;
  std::vector<double> G; // Gx Gy Gz Gc
  if(mode=="ewald"){
    double acc=1e-5;
    double g = acc*sqrt(N*rc*V)/(2.0*qsq);
    if(g>=1.0) g=(1.35-0.15*log(acc))/rc; else g=sqrt(-log(g))/rc;
    g_ewald = g_user>0 ? g_user : g;
    auto err=[&](long km,double len){ return 2.0*qsq*g_ewald/len*sqrt(1.0/(M_PI*km*N))*exp(-M_PI*M_PI*km*km/(g_ewald*g_ewald*len*len)); };
    if(kmax<=0){
      kxmax=1; while(err(kxmax,L[0])>acc) kxmax++;
      kymax=1; while(err(kymax,L[1])>acc) kymax++;
      kzmax=1; while(err(kzmax,L[2])>acc) kzmax++;
      kmax=std::max(kxmax,std::max(kymax,kzmax));
    } else if(fix_kmax) { kxmax=kymax=kzmax=kmax; } // else: BUG reproduced: kxmax..kzmax stay 0
    double ux=2*M_PI/L[0], uy=2*M_PI/L[1], uz=2*M_PI/L[2];
    double GnMax=std::max(ux*ux*kxmax*kxmax,std::max(uy*uy*kymax*kymax,uz*uz*kzmax*kzmax));
    const double fpe0 = codata ? 1.0/(14.399645*eV) : hex_fpe0; const double epsZ = fpe0/(4*M_PI);
    double bt = 2.*M_PI/fpe0/V; ep.g_ewald=g_ewald; ep.bt_sr=2.*g_ewald/sqrt(M_PI);
    double gm=1./(4.*g_ewald*g_ewald); ep.gm_sr=0.25/M_PI/epsZ;
    long nexact=0;
    for(long kx=-kxmax;kx<=kxmax;kx++) for(long ky=-kymax;ky<=kymax;ky++) for(long kz=-kzmax;kz<=kzmax;kz++){
      if(kx*kx+ky*ky+kz*kz==0) continue;
      double gx=kx*ux,gy=ky*uy,gz=kz*uz, gn=gx*gx+gy*gy+gz*gz;
      if(gn<=GnMax*1.00001) nexact++;
      if(gn<=(fix_kmax?GnMax*1.00001:GnMax)){ G.push_back(gx);G.push_back(gy);G.push_back(gz);G.push_back(bt*exp(-gm*gn)/gn); }
    }
    printf("ewald: g_ewald=%.8f kmax=%ld k(x,y,z)max=%ld,%ld,%ld nknz(full sphere)=%zu  [with LAMMPS 1.00001 margin: %ld]\n",
      g_ewald,kmax,kxmax,kymax,kzmax,G.size()/4,nexact);
  }

  // ---------- short range pairs (cell list, min image) ----------
  std::vector<double> F(3*N,0.0), E(N,0.0); double W[6]={0,0,0,0,0,0};
  int nc[3]; for(int d=0;d<3;d++) nc[d]=std::max(3,(int)floor(L[d]/rc));
  std::vector<std::vector<long>> cells(nc[0]*nc[1]*nc[2]);
  auto cid=[&](int a,int b,int c){ return ((a+nc[0])%nc[0]) + nc[0]*(((b+nc[1])%nc[1]) + nc[1]*((c+nc[2])%nc[2])); };
  for(long i=0;i<N;i++){ int c[3]; for(int d=0;d<3;d++){ double s=at[i].x[d]/L[d]; s-=floor(s); c[d]=std::min(nc[d]-1,(int)(s*nc[d])); } cells[cid(c[0],c[1],c[2])].push_back(i); }
  #pragma omp parallel for schedule(dynamic) reduction(+:W[:6])
  for(long i=0;i<N;i++){
    int c[3]; for(int d=0;d<3;d++){ double s=at[i].x[d]/L[d]; s-=floor(s); c[d]=std::min(nc[d]-1,(int)(s*nc[d])); }
    for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++)for(int cc=-1;cc<=1;cc++)
      for(long j: cells[cid(c[0]+a,c[1]+b,c[2]+cc)]){
        if(j==i) continue;
        double dr[3]; for(int d=0;d<3;d++){ dr[d]=at[j].x[d]-at[i].x[d]; dr[d]-=L[d]*round(dr[d]/L[d]); }
        double d2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]; if(d2>rc*rc) continue;
        double r=sqrt(d2), e=0, de=0, cq=q[i]*q[j];
        if(mode=="wolf") wolf_compute_energy(wp,cq,r,e,de);
        else if(mode=="dsf") dsf_compute_energy(wp,cq,r,e,de);
        else ewald_compute_energy(ep,cq,r,e,de);
        de/=r;
        for(int d=0;d<3;d++) F[3*i+d]+=de*dr[d];   // same as ForceOp: f += de*dr, dr = rb - ra
        E[i]+=0.5*e;
        // LAMMPS-convention virial: -0.5 * dr (x) f_on_i ... sum r_ij F_ij /2 per atom
        W[0]+=-0.5*dr[0]*de*dr[0]; W[1]+=-0.5*dr[1]*de*dr[1]; W[2]+=-0.5*dr[2]*de*dr[2];
        W[3]+=-0.5*dr[0]*de*dr[1]; W[4]+=-0.5*dr[0]*de*dr[2]; W[5]+=-0.5*dr[1]*de*dr[2];
      }
  }
  // self terms
  if(mode=="wolf"||mode=="dsf"){
    double es = (mode=="wolf") ? wp.e_shift : 0;
    if(mode=="dsf"){ double ec=erfc(wp.alpha*wp.rc), ed=exp(-wp.alpha*wp.alpha*wp.rc*wp.rc);
      double fs=-(ec/(wp.rc*wp.rc)+2.0/sqrt(M_PI)*wp.alpha*ed/wp.rc); es = ec/wp.rc - fs*wp.rc; }
    for(long i=0;i<N;i++) E[i] += -(es/2.0 + wp.alpha/sqrt(M_PI))*q[i]*q[i]*wp.qqrd2e*eV;
  }
  double Erecip=0; std::vector<double> Erecip_atom(N,0.0); double Wrec[6]={0,0,0,0,0,0};
  if(mode=="ewald"){
    const double fpe0 = codata ? 1.0/(14.399645*eV) : hex_fpe0;
    for(long i=0;i<N;i++) E[i] -= 1./fpe0*g_ewald/sqrt(M_PI)*q[i]*q[i];
    size_t nk=G.size()/4; std::vector<std::complex<double>> rho(nk,0.0);
    #pragma omp parallel
    { std::vector<std::complex<double>> loc(nk,0.0);
      #pragma omp for
      for(long i=0;i<N;i++) for(size_t k=0;k<nk;k++){ double ps=at[i].x[0]*G[4*k]+at[i].x[1]*G[4*k+1]+at[i].x[2]*G[4*k+2]; loc[k]+=q[i]*std::complex<double>(cos(ps),sin(ps)); }
      #pragma omp critical
      for(size_t k=0;k<nk;k++) rho[k]+=loc[k]; }
    #pragma omp parallel for
    for(long i=0;i<N;i++){ double qq=2*q[i], l[3]={0,0,0};
      for(size_t k=0;k<nk;k++){ double ps=at[i].x[0]*G[4*k]+at[i].x[1]*G[4*k+1]+at[i].x[2]*G[4*k+2];
        double al=qq*G[4*k+3]*(rho[k].real()*sin(ps)-rho[k].imag()*cos(ps)); l[0]+=al*G[4*k]; l[1]+=al*G[4*k+1]; l[2]+=al*G[4*k+2]; }
      for(int d=0;d<3;d++) F[3*i+d]+=l[d]; }
    for(size_t k=0;k<nk;k++) Erecip += G[4*k+3]*std::norm(rho[k]);
    // Phase-1 additions (validated here): per-atom recip energy and recip virial.
    // E_i = q_i * sum_k Gc Re( conj(S_k) e^{i k.r_i} )  (sums to Erecip)
    // W_ab = sum_k Gc |S_k|^2 ( d_ab - 2 (1 + G^2/(4g^2)) G_a G_b / G^2 )
    double Esum=0;
    #pragma omp parallel for reduction(+:Esum)
    for(long i=0;i<N;i++){ double ei=0;
      for(size_t k=0;k<nk;k++){ double ps=at[i].x[0]*G[4*k]+at[i].x[1]*G[4*k+1]+at[i].x[2]*G[4*k+2];
        ei += G[4*k+3]*(rho[k].real()*cos(ps)+rho[k].imag()*sin(ps)); }
      Erecip_atom[i]=q[i]*ei; Esum+=q[i]*ei; }
    for(size_t k=0;k<nk;k++){ double gx=G[4*k],gy=G[4*k+1],gz=G[4*k+2], g2=gx*gx+gy*gy+gz*gz;
      double u=G[4*k+3]*std::norm(rho[k]), f=2.0*(1.0+g2/(4.0*g_ewald*g_ewald))/g2;
      Wrec[0]+=u*(1-f*gx*gx); Wrec[1]+=u*(1-f*gy*gy); Wrec[2]+=u*(1-f*gz*gz);
      Wrec[3]+=u*(-f*gx*gy); Wrec[4]+=u*(-f*gx*gz); Wrec[5]+=u*(-f*gy*gz); }
    printf("per-atom recip sum = %.12e eV vs Erecip %.12e eV\n", Esum/eV, Erecip/eV);
  }

  // ---------- compare ----------
  double Etot=0; for(long i=0;i<N;i++) Etot+=E[i]; Etot+=Erecip;
  double Ltot=0, maxdf=0, maxf=0, maxde=0; long imax=0;
  for(long i=0;i<N;i++){ Ltot+=at[i].pe;
    for(int d=0;d<3;d++){ double df=fabs(F[3*i+d]/eV - at[i].f[d]); if(df>maxdf){maxdf=df; imax=i;} maxf=std::max(maxf,fabs(at[i].f[d])); }
    maxde=std::max(maxde,fabs(E[i]/eV-at[i].pe)); }
  printf("[%s] E_total exaStamp-formula = %.12e eV   LAMMPS sum pe/atom = %.12e eV   diff=%.3e  (recip part %.6e eV)\n",
    mode.c_str(), Etot/eV, Ltot, Etot/eV-Ltot, Erecip/eV);
  printf("[%s] max|dF| = %.3e eV/A (max|F| %.3e, atom id %ld)   max|dE_atom| = %.3e eV%s\n", mode.c_str(), maxdf, maxf, at[imax].id, maxde,
    mode=="ewald"?" (per-atom recip missing in exaStamp: expected large)":"");
  double P[6]; for(int k=0;k<6;k++) P[k] = (W[k]+Wrec[k])/eV/V*1.6021765e6;   // LAMMPS nktv2p (metal)
  if(mode=="ewald"){ double mde=0; for(long i=0;i<N;i++) mde=std::max(mde,fabs((E[i]+Erecip_atom[i])/eV-at[i].pe));
    printf("[ewald] WITH per-atom recip: max|dE_atom| = %.3e eV\n",mde); }
  printf("[%s] pressure tensor SR+recip (bar) xx yy zz xy xz yz = %.10e %.10e %.10e %.10e %.10e %.10e\n",mode.c_str(),P[0],P[1],P[2],P[3],P[4],P[5]);
  FILE* fo=fopen(("exastamp_formula_"+mode+(kmax_user>0?"_kmax"+std::to_string(kmax_user):"")+".txt").c_str(),"w");
  fprintf(fo,"# id fx fy fz (eV/A) pe_atom (eV, no recip)  | E_total %.15e eV  Erecip %.15e eV\n",Etot/eV,Erecip/eV);
  for(long i=0;i<N;i++) fprintf(fo,"%ld %.15e %.15e %.15e %.15e\n",at[i].id,F[3*i]/eV,F[3*i+1]/eV,F[3*i+2]/eV,E[i]/eV);
  fclose(fo);
  return 0;
}
