#include "snap_idx_cache.h"

#include <algorithm>
#include <iostream>

// =======================================================
// =============== Index cache ===========================
// =======================================================

namespace SnapExt
{


snapBsIdxCache::snapBsIdxCache( int twojmax, const snapCg &cg, const double* bs_coefs, int gsh_size)
{
  m_indices_two_jmax = twojmax;
  m_main.clear();
  m_cgprod.clear();
  m_subiter.clear();

  const auto* __restrict__ cg_tab = cg.cg_tab_data();
  
  std::vector<OrderedSubIter> tmp_subiter;
  std::vector<int> cmm_use_count( gsh_size , 0 );

//    std::set< std::pair<int,int> > subidx_pairs_set;

  int total_raw_subiters = 0;
  int total_nz_subiters = 0;
  int max_subiters = 0;
  int bs_idx = 0;
  std::vector<BSStepT> bs_steps;

  for (int j1=0; j1 <= m_indices_two_jmax; j1++)
  {
    // Thompson - Lammps
    for (int j2=0; j2 <= j1; j2++)
    { 
      for (int j=abs(j1-j2); j<=std::min(m_indices_two_jmax,j1+j2); j+=2)
      {
        if ((j+j1+j2)%2==1) continue;
        if (j<j1) continue; // Thompson

        // here, one new bispectrum starts
        //int bs_idx = m_bs.size();
        double bscoef = bs_coefs[bs_idx+1];

        for (int m1=-j; m1 <=j; m1+=2)
        {
          for (int m2=-j; m2 <=j; m2+=2)
          {
            //std::cout <<"j1="<<j1<<", j2="<<j2<<", j="<<j<<", m1="<<m1<<", m2="<<m2<<", idx_j_m1_m2="<<idx_j_m1_m2<<", didx_i_j_m1_m2="<<didx_i_j_m1_m2<<std::endl;
            
            SnapBsIndexT nsubiters = 0;
            tmp_subiter.clear();
            
            for (int m11=std::max(-j1,m1-j2); m11<=std::min(j1,m1+j2); m11+=2)
            {
              for (int m12=std::max(-j1,m2-j2); m12<=std::min(j1,m2+j2); m12+=2)
              {
                int m21 = m1-m11;
                int m22 = m2-m12;
                ++ nsubiters;

                SnapBsIndexT cg_idx_j_j1_j2_m1_m11_m21 = cg.index(j,j1,j2,m1,m11,m21);
                SnapBsIndexT cg_idx_j_j1_j2_m2_m12_m22 = cg.index(j,j1,j2,m2,m12,m22);
                SnapBsIndexT idx_j1_m11_m12 = snapGsh::idx(j1,m11,m12);
                SnapBsIndexT idx_j2_m21_m22 = snapGsh::idx(j2,m21,m22);

                ++ cmm_use_count[idx_j1_m11_m12];
                ++ cmm_use_count[idx_j2_m21_m22];
  
                if( idx_j1_m11_m12 > idx_j2_m21_m22 )
                {
                  std::swap( idx_j1_m11_m12, idx_j2_m21_m22 );
                }

                double cg_prod = cg_tab[cg_idx_j_j1_j2_m1_m11_m21] * cg_tab[cg_idx_j_j1_j2_m2_m12_m22] * bscoef;
                if( cg_prod != 0. )
                {
                  tmp_subiter.push_back( OrderedSubIter{cg_prod, idx_j1_m11_m12, idx_j2_m21_m22} );
                }
              }
            }
            total_raw_subiters += nsubiters;
            total_nz_subiters += tmp_subiter.size();


            BSStepT bs;
            bs.idx = snapGsh::idx(j,m1,m2);
            for(const auto& si:tmp_subiter)
            {
              if( si.cgprod!=0. )
              {
                bs.osi.push_back( si );
              }
            }
                        
//              std::cout << "BS : IDX="<<bs.idx<<" : NSI="<<bs.osi.size()<<std::endl;
            bs_steps.push_back( bs );              
          }   //m2
        }	//m1
        
        ++ bs_idx;
      }	//j
    }		//j2
  }		//j1  // Thompson

  if( bs_steps.empty() )
  {
    std::cerr<<"internal error: empty bs steps"<<std::endl;
    std::abort();
  }

  std::sort( bs_steps.begin() , bs_steps.end() );
  
  // merge bs steps from different bispectra
  {
    int last_main_idx = -1;
    int n_bs_steps = bs_steps.size();
    std::vector<BSStepT> tmp_bs_steps;
    for(int i=0;i<n_bs_steps;i++)
    {
      if( bs_steps[i].idx != last_main_idx )
      {
//          std::cout << "add BS : IDX "<<bs_steps[i].idx<<" : NSI="<<bs_steps[i].osi.size()<<std::endl;
        tmp_bs_steps.push_back( bs_steps[i] );
        last_main_idx = bs_steps[i].idx;
      }
      else
      {
//          std::cout << "merge BS : IDX "<<bs_steps[i].idx<<" : NSI="<<bs_steps[i].osi.size()<<std::endl;
        for(const auto& si:bs_steps[i].osi) { tmp_bs_steps.back().osi.push_back(si); }
      }
    }
    bs_steps = tmp_bs_steps;
  }

  // remove duplicates and accumulate corresponding coeficients
  for(auto& bss : bs_steps)
  {
    std::sort( bss.osi.begin() , bss.osi.end() );
    auto tmp_subiter = bss.osi;
    bss.osi.clear();
    OrderedSubIter accum { 0. , 0 , 0 };
    for(const auto& si:tmp_subiter)
    {
      if( si.idx_a==accum.idx_a && si.idx_b==accum.idx_b )
      {
        accum.cgprod += si.cgprod;
      }
      else
      {
        if( accum.cgprod!=0. )
        {
          bss.osi.push_back( accum );
        }
        accum = si;
      }
    }
    if( accum.cgprod!=0. )
    {
      bss.osi.push_back( accum );
    }
  }

  int last_main_idx = -1;
  if( ! bs_steps.empty() )
  {
    m_main_idx_start = bs_steps.front().idx;
    last_main_idx = m_main_idx_start - 1;
  }
  
  int bsscount=0;
  for( const auto& bss : bs_steps )
  {
    if( bss.idx != (last_main_idx+1) )
    {
      std::cerr << "Internal error: non monotonic main idx : "<<last_main_idx<<"+1 does not follow "<<bss.idx<<std::endl;
      std::abort();
    }
    last_main_idx = bss.idx;
  
    if( static_cast<int>(bss.osi.size()) > max_subiters ) { max_subiters = bss.osi.size(); }
    m_main.push_back( /*MainInfoT{*/bss.osi.size()/*,bss.idx}*/ );

    /*
    std::cout<<"BS : "<<bsscount<<" : MI="<< bss.idx<<" : N="<<bss.osi.size() ;//<<" :";
    int nsi = bss.osi.size();
    int last_main_idx_occurence = -1;
    for(int k=0;k<nsi;k++)
    {
      //std::cout <<" "<<bss.osi[k].cgprod;
      std::cout <<" ("<<bss.osi[k].idx_a<<","<<bss.osi[k].idx_b<<")";
    }
    std::cout << std::endl;
    */

    for(const auto& si : bss.osi)
    {
      m_cgprod.push_back( si.cgprod );
      m_subiter.push_back( SubIterT{ static_cast<SnapBsIndexT>(si.idx_a), static_cast<SnapBsIndexT>(si.idx_b) } );
      //subidx_pairs_set.insert( std::pair<int,int>(si.idx_a,si.idx_b) );
    }
    ++ bsscount;
  }
  
  m_n_bs = m_main.size();

  {
    int i=0;
    while(i<gsh_size && cmm_use_count[i]==0) ++i;
    m_gsh_start = i;
    m_gsh_end = i;
    for(;i<gsh_size;i++)
    {
      if( cmm_use_count[i]!=0 ) { m_gsh_end = i+1; }
    }
  }

/*
  for(int i=m_gsh_start;i<m_gsh_end;i++)
  {
    if( cmm_use_count[i]==0 ) { std::cout<<"cmm["<<i<<"] unused"<<std::endl; }
  }
  std::vector<int> subitcount( max_subiters , 0 );
  int total_si = 0;
  for(int nsi : m_main)
  {
    subitcount[ nsi ] ++;
    total_si += nsi;
  }
  for(int i=0;i<max_subiters;i++)
  {
    if( subitcount[i] != 0 ) { std::cout << "subit x"<<i<<" : "<<subitcount[i]<<" "<<i*subitcount[i]<<" "<< (i*subitcount[i])*100.0/total_si <<" %" <<std::endl; }
  }
*/

  // adjust indices according to m_gsh_start
  m_main_idx_start -= m_gsh_start;
  for(auto& si : m_subiter)
  {
    si.idx_a -= m_gsh_start;
    si.idx_b -= m_gsh_start;
  }

/*
  std::set<double> cgpset;
  for(double x : m_cgprod) { cgpset.insert(x); }
  std::cout<<"main="<<m_main.size()<<", subiter="<<m_subiter.size()<<", nz="<<total_nz_subiters<<", raw="<<total_raw_subiters
           <<", gshrange=["<<m_gsh_start<<','<<m_gsh_end<<"], ncgprod="<<cgpset.size() <<std::endl;
*/

  // CMM idx test this ensures that we can use a 0 to compact_size()-1 loop instead for the 3 nested loops
  size_t cmmidx = 0;
  for (int j = 0; j <= m_indices_two_jmax; ++j)
  {
    for (int m2 = -j; m2 <= j; m2 += 2)
    {
      //std::cout << "GSH diag "<<compact_idx(j,m2,m2) << std::endl;
      for (int m1 = -j; m1 <= j; m1 += 2)
      {
        int idx1 = snapGsh::idx(j,m1,m2) - m_gsh_start;
        int idx2 = compact_idx(j,m1,m2);
        if( idx1!=idx2 ) { std::abort(); }
        if( idx1 != ssize_t(cmmidx) ) { std::abort(); }
        ++ cmmidx;
      }
    }
  }
  // std::cout << "contiguous CMM indices in [0;"<<cmmidx<<"[ compact size = "<<compact_size()<<std::endl;
}

// Create a new Index cache or retreive an existing one
snapBsIdxCache* snapBsIdxCache::getIdxCache(int twojmax, const snapCg &cg, const double* bs_coefs, int gsh_size)
{
  const std::lock_guard<std::mutex> lock(s_mutex);  
  auto it = s_cache_map.find( twojmax );
  if( it == s_cache_map.end() )
  {
    snapBsIdxCache* ic = new snapBsIdxCache(twojmax,cg,bs_coefs,gsh_size);
    //std::cout << "gsh compact size = " << ic->compact_size() << " , gsh start = "<< ic->gsh_start()<<" , m_n_bs = "<<ic->m_n_bs<<std::endl;
    s_cache_map[twojmax] = ic;
    return ic;
  }
  else
  {
    return it->second;
  }
}


}

