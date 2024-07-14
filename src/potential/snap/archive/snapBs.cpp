
#include "snapBs.h"
#include <exanb/core/physics_constants.h>

#include <iostream>
#include <fstream>
#include <iomanip> // setw.
#include <cmath> // sqrt.
#include <strings.h>
#include <algorithm> 
#include <set> 
#include <map>
#include <list>
#include <cassert>
// #define SNAP_VERBOSE_DBG 1

#include <onika/force_assert.h>
#include <onika/integral_constant.h>

#include "snap_idx_cache.h"
#include "snap_gsh_opt.h"

namespace SnapExt
{

std::map<int,snapBsIdxCache*> snapBsIdxCache::s_cache_map;
std::mutex snapBsIdxCache::s_mutex;

// =======================================================
// =======================================================
// =======================================================
template<class MapT>
static int snap_gen_gsh_recursion_depth(int blockStart, int i, MapT & gsh_depends)
{
  if( gsh_depends.find(i) == gsh_depends.end() ) return 0;
  if(i<blockStart) return 0;
  int a = gsh_depends[i].first_idx;
  int b = gsh_depends[i].second_idx;
  return std::max( snap_gen_gsh_recursion_depth(blockStart,a,gsh_depends) , snap_gen_gsh_recursion_depth(blockStart,b,gsh_depends) ) + 1;
}

/*
template<class MapT>
static void snap_gen_gsh_expr( int i, MapT & gsh_depends)
{
  if( i==0 ) std::cout<<"0";
  else if( i==1 ) std::cout<<"1";
  else
  {
    const char* ft[2] = { "Zm", "Zp" };
    const char* st[2] = { "Xp", "Xm" };  
    auto e = gsh_depends[i];
  }
}
*/

static void snap_gen_gsh_recursion_dependencies(
        const SnapExt::SnapBsIndexT* main_data
      , const double* cgprod_data
      , const SnapExt::SubIterT* subit_data
      , int n_bs
      , int m_two_jmax
      , int KernelBlockSize
      )
{
  static constexpr bool verbose = false;

  struct BSWorkItem { double cg_prod; int bs; int idx_a; int idx_b; };  
  std::map<int , std::vector<BSWorkItem> > progressive_bs_compute;
  int main_cursor = 0;
  int sub_cursor = 0;
  //int main_idx = main_idx_start; // bs_steps_data[bs_idx].start_main_idx;
  for(int bs_step=0;bs_step<n_bs;bs_step++) // constant number of steps
  {
    int nsubiters = main_data[main_cursor]; //.nsubiters;
    ++ main_cursor;
    const auto* lcgprod = cgprod_data + sub_cursor;
    const auto* lsubit = subit_data + sub_cursor;
    sub_cursor += nsubiters;
    ONIKA_FORCE_ASSERT( nsubiters > 0 );
    ONIKA_FORCE_ASSERT( lsubit[0].idx_a == 0 && lsubit[0].idx_b == bs_step );
    for(int si=0;si<nsubiters;si++)
    {
      auto subit = lsubit[si];
      int a=subit.idx_a , b=subit.idx_b;
      if( a > b ) std::swap(a,b);
      progressive_bs_compute[b].push_back( {lcgprod[si],bs_step,a,b} );
    }
  }
  
  if(verbose)
  {
    for(const auto& pbs:progressive_bs_compute)
    {
      std::cout<<"CMM step "<<pbs.first<<" ("<<pbs.second.size()<<") :";
      for(const auto& wi:pbs.second) std::cout<<" "<<wi.bs;
      std::cout<<"\n";      
    }
  }

  struct GshRecurs
  {
    int first_a;
    int first_b;
    int first_idx;
    int second_a;
    int second_b;
    int second_idx;
    bool m2eqmj;
  };

  std::map<int,int> last_used;
  std::map<int,int> first_used;
  std::map< int , GshRecurs > gsh_depends;
  
  auto add_recursion_step = [&last_used,&first_used](int r, int a, int b)
  {
    if(first_used.find(a)==first_used.end()) first_used[a]=r;
    if(first_used.find(b)==first_used.end()) first_used[b]=r;
    last_used[a] = r;
    last_used[b] = r;
  };

  for (int j = 1 /*0 OK: recursion initialized*/; j <= m_two_jmax; ++j) 
  {
    for (int m2 = -j; m2 <= j; m2 += 2) 
    {
      for (int m1 = -j; m1 <= j; m1 += 2) 
      {
        // Compute gsh for j, m1, m2.
        const int idx_j_m1_m2 = snap_gsh_idx(j,m1,m2);

        if (m2 == -j) { // Bartok PhD (B.8)
          const int idx_jm1_m1p1_m2p1 = snap_gsh_idx(j-1,m1+1,m2+1);
          const int idx_jm1_m1m1_m2p1 = snap_gsh_idx(j-1,m1-1,m2+1);
          //double first  = sqrt((double)(j-m1)/(double)(j-m2));
          //double second = sqrt((double)(j+m1)/(double)(j-m2));
          // z_plus , x_minus
          gsh_depends[idx_j_m1_m2] = { j-m1, j-m2, idx_jm1_m1p1_m2p1, j+m1, j-m2, idx_jm1_m1m1_m2p1, true };
          add_recursion_step(idx_j_m1_m2,idx_jm1_m1p1_m2p1,idx_jm1_m1m1_m2p1);
        }
        else { // Bartok PhD (B.8)
          const int idx_jm1_m1m1_m2m1 = snap_gsh_idx(j-1,m1-1,m2-1);
          const int idx_jm1_m1p1_m2m1 = snap_gsh_idx(j-1,m1+1,m2-1);
          //double first  = sqrt((double)(j+m1)/(double)(j+m2));
          //double second = sqrt((double)(j-m1)/(double)(j+m2));
          // z_minus, x_plus
          gsh_depends[idx_j_m1_m2] = { j+m1, j+m2, idx_jm1_m1m1_m2m1, j-m1, j+m2, idx_jm1_m1p1_m2m1, false };
          add_recursion_step(idx_j_m1_m2,idx_jm1_m1m1_m2m1,idx_jm1_m1p1_m2m1);
        }
      }
    }
  }

  if(verbose)
  {
    int i=3;
    std::cout<<"GSH start "<<i<<"\n";
    while( i <= n_bs )
    {
      int j=i;
      while( j<=n_bs && ( gsh_depends.find(j)==gsh_depends.end() || ( gsh_depends[j].first_idx<i && gsh_depends[j].second_idx<i ) ) ) j++;
      std::cout<<"GSH round ["<<i<<";"<<j<<"]\n";
      i = j;
    }
  }

/*
  if(verbose)
  {
    for(const auto& gd:gsh_depends)
    {
      std::cout<<"GSH #"<<gd.first<<" depth=" << snap_gen_gsh_recursion_depth(0,gd.first,gsh_depends) <<"\n";
    }
  }
*/

  if(verbose)
  {
    for(const auto& d:gsh_depends)
    {
      std::cout<<d.first<<" <- "<<d.second.first_idx<<" , "<<d.second.second_idx<<" : used in ["<<first_used[d.first]<<";"<<last_used[d.first]<<"] : range="<<last_used[d.first]-d.first+1<<"\n";
    }
  }
}

void snap_gen_gsh_sqrt_table(std::ofstream& fout, int m_two_jmax)
{
  std::map< std::pair<int,int> , double > coef;
  for (int j = 1 /*0 OK: recursion initialized*/; j <= m_two_jmax; ++j) 
  {
    for (int m2 = -j; m2 <= j; m2 += 2) 
    {
      for (int m1 = -j; m1 <= j; m1 += 2) 
      {
        if (m2 == -j) { // Bartok PhD (B.8)
          coef[ {j-m1,j-m2} ] = sqrt((double)(j-m1)/(double)(j-m2));
          coef[ {j+m1,j-m2} ] = sqrt((double)(j+m1)/(double)(j-m2));
        }
        else
        {
          coef[ {j+m1,j+m2} ] = sqrt((double)(j+m1)/(double)(j+m2));
          coef[ {j-m1,j+m2} ] = sqrt((double)(j-m1)/(double)(j+m2));
        }
      }
    }
  }

  int maxn=0, maxq=0;    
  for(const auto & p : coef)
  {
    assert( p.first.first%2 == 0 && p.first.second%2 == 0 );
    if(p.first.first>maxn) maxn=p.first.first;
    if(p.first.second>maxn) maxq=p.first.second;
  }
  
  fout << "static constexpr double gsh_sqrt_frac[" <<maxn/2+1<<"]["<<maxq/2+1<<"] = {\n";
  for(int n=0;n<=maxn/2;n++)
  {
    fout << (n==0? " {" : ",{");
    for(int q=0;q<=maxq/2;q++)
    {
      fout<< (q==0?' ':',') << std::hexfloat << coef[{n*2,q*2}] << std::defaultfloat;
    }
    fout << " }\n";      
  }
  fout << "};\n";
}


static inline void snap_gen_bulk_lbs_tables(
        const SnapExt::SnapBsIndexT* main_data
      , const double* cgprod_data
      , const SnapExt::SubIterT* subit_data
      , int n_bs
      , std::vector< SnapExt::BSFullBlockWorkItem >& bs_fblock
      )
{
  //struct BSMultiBlockWorkItem{ double cg_prod; uint16_t bs; uint16_t wo; uint16_t idx_a; uint16_t idx_b; };
  struct BSFullBlockWork{ SnapExt::BSFullBlockWorkItem work_item[SnapExt::CUDA_BLOCK_SIZE]; };
  struct BSMultiBlockWork{ int n_work_items; int n_write_passes; SnapExt::BSFullBlockWorkItem work_item[SnapExt::CUDA_BLOCK_SIZE]; };

  static constexpr int blockSize = SnapExt::CUDA_BLOCK_SIZE;
  static constexpr bool verbose = false;

  struct WorkItem { double cg_prod; int idx_a; int idx_b; };
  struct BSWorkBatch 
  {
    int bs;
    std::list< WorkItem > workitem;
  };
  
  std::vector< BSWorkBatch > bs_work_batches( n_bs );
  int total_subit = 0;
  int main_cursor = 0;
  int sub_cursor = 0;
  //int main_idx = main_idx_start; // bs_steps_data[bs_idx].start_main_idx;
  for(int bs_step=0;bs_step<n_bs;bs_step++) // constant number of steps
  {
    int nsubiters = main_data[main_cursor]; //.nsubiters;
    ++ main_cursor;
    const auto* __restrict__ lcgprod = cgprod_data + sub_cursor;
    const auto* __restrict__ lsubit = subit_data + sub_cursor;
    sub_cursor += nsubiters;
    bs_work_batches[bs_step].bs = bs_step;
    std::vector< WorkItem > tmp;
    for(int si=0;si<nsubiters;si++)
    {
      auto subit = lsubit[si];
      int a=subit.idx_a , b=subit.idx_b;
      if( a > b ) std::swap(a,b);
      tmp.push_back( { lcgprod[si], a , b } );
    }
    //std::vector< WorkItem > tmp ( bs_work_batches[bs_step].workitem.begin() , bs_work_batches[bs_step].workitem.end() );
    std::sort( tmp.begin() , tmp.end() , [](const WorkItem& a, const WorkItem& b)->bool{ return a.idx_a<b.idx_a || (a.idx_a==b.idx_a && a.idx_b<b.idx_b); } );
    bs_work_batches[bs_step].workitem.assign( tmp.begin() , tmp.end() );
    total_subit += nsubiters;
  }
  std::sort( bs_work_batches.begin() , bs_work_batches.end() , [](const BSWorkBatch& a, const BSWorkBatch& b)->bool{return a.workitem.size()>b.workitem.size();} );

  struct BSWorkItem { double cg_prod; int bs; int idx_a; int idx_b; };
  std::vector< std::vector<BSWorkItem> > bs_work_blocks;
  bool multiblock_strategy_switch = false;
  while( total_subit > 0 )
  {
    std::vector<BSWorkItem> work_block;
    if( multiblock_strategy_switch || bs_work_batches.size() < blockSize || bs_work_batches[blockSize-1].workitem.empty() )
    {
      std::sort( bs_work_batches.begin() , bs_work_batches.end() , [](const BSWorkBatch& a, const BSWorkBatch& b)->bool{return a.workitem.size()>b.workitem.size();} );
    }
    int i;
    for(i=0;i<blockSize && !bs_work_batches[i].workitem.empty() ; i++ )
    {
      auto wi = bs_work_batches[i].workitem.front();
      work_block.push_back( {wi.cg_prod, bs_work_batches[i].bs, wi.idx_a, wi.idx_b} );
      bs_work_batches[i].workitem.pop_front();
      -- total_subit;
    }
    std::sort( work_block.begin(), work_block.end(), [](const BSWorkItem& a,const BSWorkItem& b)->bool{return a.idx_a<b.idx_a || ( a.idx_a==b.idx_a && a.idx_b<b.idx_b); } );
//    std::cout<<"add work block size = "<<work_block.size()<<"\n";
    if( work_block.size() < blockSize ) multiblock_strategy_switch = true;
    bs_work_blocks.push_back( std::move(work_block) );
  }
  for(const auto& bswb:bs_work_batches) { ONIKA_FORCE_ASSERT( bswb.workitem.empty() ); }
  
  int tfb_count=0;
  int tmb_count=0;
  int mbsz=0;
  for(size_t i=0;i<bs_work_blocks.size();i++)
  {
    if( bs_work_blocks[i].size() == blockSize ) ++ tfb_count;
    else
    {
      mbsz+=bs_work_blocks[i].size();
      if( mbsz>=blockSize ) { ++tmb_count; mbsz-=blockSize; }
    }
  }
  if(mbsz>0) ++tmb_count;


  // === generate full blocks ===

  bs_fblock.clear(); bs_fblock.reserve(tfb_count*blockSize);
  //int fbcount = 0;
  int flushcount = 0;
  std::set<int> bs_prev_access;
  for(size_t i=0;i<bs_work_blocks.size();i++)
  {
    if( bs_work_blocks[i].size() == blockSize )
    {
      BSFullBlockWork fblock;
      int nprevcommon=0;
      std::set<int> nextaccess;
      for(int j=0;j<blockSize;j++)
      {
        int bs = bs_work_blocks[i][j].bs;
        fblock.work_item[j] = SnapExt::BSFullBlockWorkItem { bs_work_blocks[i][j].cg_prod , uint16_t(bs) , uint16_t(bs_work_blocks[i][j].idx_a) , uint16_t(bs_work_blocks[i][j].idx_b) , 0 };
        if( bs_prev_access.empty() || bs_prev_access.find(bs)!=bs_prev_access.end() ) ++ nprevcommon;
        nextaccess.insert(bs);
      }
      std::sort( fblock.work_item , fblock.work_item+blockSize , [](const SnapExt::BSFullBlockWorkItem& a, const SnapExt::BSFullBlockWorkItem& b)->bool{return a.bs<b.bs;} );
      if( nprevcommon < blockSize ) { ++flushcount; }
      bs_prev_access = std::move( nextaccess );
      // std::cout<<"FB "<<(fbcount++)<<" : ncommon="<<nprevcommon<<"\n";
      for(int k=0;k<blockSize;k++) bs_fblock.push_back( fblock.work_item[k] );
    }
  }

  if( verbose )
  {
    std::cout << "Block LBS tables for blockSize="<<blockSize<<" :\n";
    std::cout << bs_fblock.size()/blockSize <<" full blocks , "<<(flushcount+1)<<" flushes\n";
  }

  // === generate partial blocks ===
  if( verbose )  std::cout << tmb_count <<" multi-pass blocks\n";

  int cursor = bs_fblock.size();
  for(size_t i=0;i<bs_work_blocks.size();i++)
  {
    if( bs_work_blocks[i].size() < blockSize )
    {
      for(size_t j=0;j<bs_work_blocks[i].size();j++)
      {
        auto wi = bs_work_blocks[i][j];
        bs_fblock.push_back( SnapExt::BSFullBlockWorkItem{ wi.cg_prod , uint16_t(wi.bs) , uint16_t(wi.idx_a), uint16_t(wi.idx_b), 0 } );
      }
    }
  }
  std::sort( bs_fblock.begin()+cursor , bs_fblock.end() , [](const SnapExt::BSFullBlockWorkItem& a, const SnapExt::BSFullBlockWorkItem& b)->bool{return a.bs<b.bs;} );

  if( verbose )
  {
    for(size_t i=0;i<bs_fblock.size();i++)
    {
      if(i%blockSize==0) std::cout<<"\n";
      std::cout<<" "<<bs_fblock[i].bs;
    }
    std::cout<<"\n";
  }
  
}

using namespace std;


void snapBs::generate_constant_tables(std::ofstream& fout)
{
//  const auto* main_data = m_idx_cache->m_main.data();
//  const auto* cgprod_data = m_idx_cache->m_cgprod.data();
//  const auto* subit_data = m_idx_cache->m_subiter.data();
  const int main_size = m_idx_cache->m_main.size();
  const int cgprod_size = m_idx_cache->m_cgprod.size();
  const int subit_size = m_idx_cache->m_subiter.size();
  const int n_bs = m_idx_cache->m_n_bs;
  const unsigned int jmi = static_cast<unsigned int>(m_jmax);
//  const char jm[2] = { static_cast<char>('0'+jmi) , '\0' };

  fout<<"template<> struct SnapConstants<"<<jmi<<">\n{\n";
  fout<<"static constexpr int compact_gsh_size = "<<m_idx_cache->compact_size()<<";\n";
  fout<<"static constexpr int gsh_start = 1;\n"; // the one here is different from m_idx_cache->gsh_start(), because we use snap_gsh_opt function
  fout<<"static constexpr int m_jmax = "<<m_jmax<<";\n";
  fout<<"static constexpr int m_two_jmax = "<<m_two_jmax<<";\n";
  fout<<"static constexpr int gsh_data_size = "<<m_gsh.size()<<";\n";
  fout<<"static constexpr int n_bs = "<<n_bs<<";\n";
  
  fout<<"static constexpr size_t main_data_size = "<<main_size<<";\n";
  fout<<"static constexpr size_t cgprod_data_size = "<<cgprod_size<<";\n";
  fout<<"static constexpr size_t subit_data_size = "<<subit_size<<";\n";

  // cmm initialization with ones at some places only
  {
    const int compact_gsh_size = m_idx_cache->compact_size();
    const int gsh_start = m_idx_cache->gsh_start();
    std::cout<<"gsh_start="<<gsh_start<<", compact_gsh_size="<<compact_gsh_size<<std::endl;

    int n=0;
    for (int j = 0; j <= int(m_two_jmax); ++j) for (int m1 = -j; m1 <= j; m1 += 2) ++n;
    char sep=' ';
    fout<<"static constexpr size_t n_cmm_ones ="<<n<<";\nstatic constexpr ConstOneIndices<n_cmm_ones> cmm_init_ones = { {";
    for (int j = 0; j <= int(m_two_jmax); ++j)
    {
      for (int m1 = -j; m1 <= j; m1 += 2)
      {
        int k = m_idx_cache->compact_idx(j,m1,m1);
        assert( k >= 0 && k < compact_gsh_size );
        fout<<sep<<k; sep=',';
      }
    }
    fout<<" } };\n\n";
  }

  snap_gen_gsh_sqrt_table(fout, m_two_jmax);

  fout<<"};\n\n";
}


void snapBs::generate_lbs_compute_blocks(std::vector< SnapExt::BSFullBlockWorkItem >&  bs_fblock )
{
  snap_gen_bulk_lbs_tables( m_idx_cache->m_main.data(), m_idx_cache->m_cgprod.data(), m_idx_cache->m_subiter.data(), m_idx_cache->m_n_bs, bs_fblock );  
}

snapBs::snapBs(double const jmax, const snapCg &cg, double const *coefs, double factor)
  : m_nidx(n_idx_bs(floor(2*jmax)))
  , m_two_jmax(floor(2*jmax))
  , m_gsh_size(pow(floor(2*jmax)+1,3))
  , m_jmax(jmax)
  , m_factor(factor)
  , m_coefs(coefs)
  , m_cg(cg)
{
  m_idx_cache = snapBsIdxCache::getIdxCache(m_two_jmax,m_cg,coefs,m_gsh_size);  
  m_gsh.set_jmax(m_jmax);
  compute_bs0();

  snap_gen_gsh_recursion_dependencies(m_idx_cache->m_main.data(), m_idx_cache->m_cgprod.data(), m_idx_cache->m_subiter.data(), m_idx_cache->m_n_bs, m_two_jmax, 32 );
  
#if 0
  std::ofstream fout(std::string("snap_constants") + std::to_string(int(m_jmax))+ ".h");
  generate_constant_tables(fout);
#endif
}

//Compute the number of bispectrum components
int snapBs::n_idx_bs(int m_two_jmax)
{
  int n_idx=0;
  for (int j1=0; j1 <=m_two_jmax; j1+=1)
  {
    for (int j2=0; j2 <= j1; j2+=1)
    { //Thompson-Lammps
      for (int j=abs(j1-j2); j<=min(m_two_jmax,j1+j2); j+=2)
      {
        if ((j+j1+j2)%2==1)  continue;
        if (j<j1) continue;
        n_idx+=1;
      }
    } //Thompson-Lammps
  }
  return n_idx;
}

//Update the list of neighbours
int snapBs::set_neighbours(double const *rx, double const *ry, double const *rz, const double * rd, double rcut, size_t N_atom)
{
  m_N_atom = N_atom;
  m_rx = rx;
  m_ry = ry;
  m_rz = rz;
  m_radius = rd;
  
//  std::cout << "N_atom = "<<N_atom << std::endl;
  
  for (int i=0; i<m_N_atom; ++i)
  {
    if ( (m_radius[i]<1.e-12) || (m_radius[i]>rcut) )
    {
      std::cerr << "SnapV2: Illegal radius for neighbor "<<i<<" r="<<m_radius[i]<<", rcut="<<rcut<<", min_r="<<1.e-12<<std::endl;
      std::abort();
    }
#if !defined(NDEBUG) && defined(SNAP_VERBOSE_DBG)
    if (m_radius[i]<rcut) { cout << "neighbour: " << i << ", x,y,z,r : " << rx[i] << ", " << ry[i] << ", " << rz[i] << ", " << m_radius[i] << endl; }
#endif
  }
  return 0;
}

//Compute the C_{mm] coefficients (linear combination of GSH)
void snapBs::compute_cmm( double rcut, double rfac0, double rmin0 )
{

  //Initialization
  //const CMMData cmm_zero = { {0.,0.} , {0.,0.} , {0.,0.} , {0.,0.} };
  const int compact_gsh_size = m_idx_cache->compact_size();
  const int gsh_start = m_idx_cache->gsh_start();

  if( m_N_atom == 0 ) { return; }
  
  //m_cmm.assign( m_N_atom * compact_gsh_size , cmm_zero );
  m_cmm.resize( m_N_atom * compact_gsh_size );
  for(int k=0;k<compact_gsh_size;k++)
  {
    m_cmm[k].cmm = 0.;
  }
  for (int j = 0; j <= m_two_jmax; ++j)
  {
    for (int m1 = -j; m1 <= j; m1 += 2)
    {
      int k = m_idx_cache->compact_idx(j,m1,m1);
      assert( k >= 0 && k < compact_gsh_size );
      m_cmm[k].cmm = 1.;
    }
  }

  for (int i=0; i<m_N_atom; ++i)
  { //loop on neighbours
    double3d r_vec = { m_rx[i] , m_ry[i] , m_rz[i] };
    m_gsh.compute_gsh(r_vec,rcut,rfac0,rmin0); //Compute the gsh for atom i
    const auto* __restrict__ gsh = m_gsh.gsh_data();
    const auto* __restrict__ dgsh = m_gsh.dgsh_data();
    double fcut=0.5*(cos(M_PI*(m_radius[i]-rmin0)/(rcut-rmin0))+1.)*m_factor;
    double dfcut=-0.5*M_PI/(rcut-rmin0)*sin(M_PI*(m_radius[i]-rmin0)/(rcut-rmin0))*m_factor;
    for(int cmm_k=0;cmm_k<compact_gsh_size;cmm_k++)
    {
      int gsh_k = cmm_k + gsh_start;
      m_cmm[cmm_k].cmm += gsh[gsh_k] * fcut;
      complex3d dc = dfcut * r_vec * gsh[gsh_k] / m_radius[i] + fcut * dgsh[gsh_k];
      m_cmm[cmm_k+i*compact_gsh_size].dcmm_x = dc.x;
      m_cmm[cmm_k+i*compact_gsh_size].dcmm_y = dc.y;
      m_cmm[cmm_k+i*compact_gsh_size].dcmm_z = dc.z;
    }    
  }

  for (int i=1; i<m_N_atom; ++i)
  {
    for(int k=0;k<compact_gsh_size;k++)
    {
      m_cmm[k+i*compact_gsh_size].cmm = m_cmm[k].cmm;
    }
  }

}

void snapBs::compute_bs()
{
  m_energy = 0.;

  if( m_N_atom == 0 )
  {
    return;
  }

  const int compact_gsh_size = m_idx_cache->compact_size();
  const auto* __restrict__ main_data = m_idx_cache->m_main.data(); KMP_ASSUME_ALIGNED(main_data);
  const auto* __restrict__ cgprod_data = m_idx_cache->m_cgprod.data(); KMP_ASSUME_ALIGNED(cgprod_data);
  const auto* __restrict__ subit_data = m_idx_cache->m_subiter.data(); KMP_ASSUME_ALIGNED(subit_data);

  //initialize bs and dbs to zero
  if( m_N_atom > static_cast<int>(m_force.size()) )
  {
    m_force.resize(m_N_atom);
  }
  //int bs_didx = 0;
  const int n_bs = m_idx_cache->m_n_bs; // long term constant
  const int main_idx_start = m_idx_cache->m_main_idx_start;
  //std::cout << "main_idx_start = " << main_idx_start << std::endl;
    
  for (int i=0; i < m_N_atom ; i++) // N_atom can be bound to a maximum constant value
  {
    int main_cursor = 0;
    int sub_cursor = 0;
    double3d f_val = { 0., 0., 0. };

    int main_idx = main_idx_start; // bs_steps_data[bs_idx].start_main_idx;
    for(int bs_step=0;bs_step<n_bs;bs_step++) // constant number of steps
    {
      int nsubiters = main_data[main_cursor]; //.nsubiters;
      ++ main_cursor;

      const auto* __restrict__ lcgprod = cgprod_data + sub_cursor;
      const auto* __restrict__ lsubit = subit_data + sub_cursor;
      sub_cursor += nsubiters;

      Complexd lbs = { 0. , 0. };
      Complexd dlbs_x = {0.,0.};
      Complexd dlbs_y = {0.,0.};
      Complexd dlbs_z = {0.,0.};
      for(int si=0;si<nsubiters;si++)
      {
        const double cg_prod = lcgprod[si];
        auto subit = lsubit[si];

        // load cmm data        
        auto cmm_a = m_cmm[ i*compact_gsh_size + subit.idx_a ];
        auto cmm_b = m_cmm[ i*compact_gsh_size + subit.idx_b ];
        
        // compute
        lbs += cg_prod * cmm_a.cmm * cmm_b.cmm ;
        dlbs_x += cg_prod * ( cmm_a.dcmm_x * cmm_b.cmm + cmm_a.cmm * cmm_b.dcmm_x );
        dlbs_y += cg_prod * ( cmm_a.dcmm_y * cmm_b.cmm + cmm_a.cmm * cmm_b.dcmm_y );
        dlbs_z += cg_prod * ( cmm_a.dcmm_z * cmm_b.cmm + cmm_a.cmm * cmm_b.dcmm_z );
      }
      
      // load main cmm data
      auto main_cmm = m_cmm[ i*compact_gsh_size + main_idx + bs_step ];
      Complexd conj_cmm = conj( main_cmm.cmm );
      Complexd conj_dcmm_x = conj( main_cmm.dcmm_x );
      Complexd conj_dcmm_y = conj( main_cmm.dcmm_y );
      Complexd conj_dcmm_z = conj( main_cmm.dcmm_z );
      
      // compute
      Complexd c_dbs_x = conj_dcmm_x * lbs + conj_cmm * dlbs_x ;
      Complexd c_dbs_y = conj_dcmm_y * lbs + conj_cmm * dlbs_y ;
      Complexd c_dbs_z = conj_dcmm_z * lbs + conj_cmm * dlbs_z ;
      m_energy += (conj_cmm * lbs ).real();
      f_val.x -= c_dbs_x.real();
      f_val.y -= c_dbs_y.real();
      f_val.z -= c_dbs_z.real();
    }
      
    m_force[i] = f_val;
  }		//atoms

  m_energy *= m_factor / m_N_atom;
  m_energy += m_coefs[0];
}

void snapBs::compute_bs0()
{
  // static const double conv_energy_inv =  1e-4 * exanb::legacy_constant::elementaryCharge / exanb::legacy_constant::atomicMass;
  
  const unsigned int nbs_comp = n_idx_bs(m_two_jmax);
  double bzero[m_two_jmax+1]; 
  for(int j=0; j< (m_two_jmax+1); j++) bzero[j]=0.;

  // vector<double> bzero = {0.};
  // bzero.resize(m_two_jmax+1);
  double www = 1.;
  bool bnormflag = false;
  for(int j = 0; j <= m_two_jmax; j++)
    if (bnormflag)
      bzero[j] = www;
    else{
      bzero[j] = www*(j+1);  
    }

  int idxb[nbs_comp];
  for(unsigned int j=0; j<nbs_comp; j++) idxb[j]=0;
//  vector<int> idxb = {0};
//  idxb.resize(nbs_comp);
  
  int idxb_count = 0;
  for(int j1 = 0; j1 <= m_two_jmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= min(m_two_jmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count] = j;
          idxb_count++;
        }

  double e0 = 0.;  
  for (unsigned int jjb = 0; jjb < nbs_comp; jjb++) {
    const int j = idxb[jjb];
    e0 += m_coefs[jjb+1] * bzero[j];
  } 

  m_energy_zero = e0;
  //m_energy_zero = 11.557310000000001 * conv_energy_inv;    
}

}

