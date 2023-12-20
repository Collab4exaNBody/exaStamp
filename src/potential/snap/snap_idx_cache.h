#pragma once

#include "snap_vector.h"
#include "snapCg.h"
#include "snapGsh.h"
#include "snap_ext.h"
#include <vector>
#include <map>
#include <map>
#include <mutex>
		
// =======================================================
// =============== Index cache ===========================
// =======================================================

namespace SnapExt
{

  struct snapBsIdxCache
  {
    using SnapBsIndexT = SnapExt::SnapBsIndexT;
    using SubIterT = SnapExt::SubIterT;

    struct OrderedSubIter
    {
      double cgprod;
      int idx_a;
      int idx_b;
      inline bool operator < (const OrderedSubIter& osi) const { return idx_a < osi.idx_a; }
    };

    struct BSStepT
    {
      int idx;
      std::vector<OrderedSubIter> osi;
      inline bool operator < (const BSStepT& bss) const { return idx < bss.idx; }
    };

    snapBsIdxCache( int twojmax, const snapCg &cg, const double* bs_coefs, int gsh_size);
    static snapBsIdxCache* getIdxCache(int twojmax, const snapCg &cg, const double* bs_coefs, int gsh_size);

    inline int gsh_start() const { return m_gsh_start; }
    inline int compact_size() const { return m_gsh_end - m_gsh_start; }
    inline int compact_idx(int j, int m1, int m2) const { return snapGsh::idx(j,m1,m2) - m_gsh_start; }

    SnapVectorT<SnapBsIndexT> m_main;
    SnapVectorT<double> m_cgprod;
    SnapVectorT<SubIterT> m_subiter;

    int m_main_idx_start = 0;
    int m_n_bs = 0;
    int m_indices_two_jmax = 0;
    int m_gsh_start = 0;
    int m_gsh_end = 0;

    static std::map<int,snapBsIdxCache*> s_cache_map;
    static std::mutex s_mutex;
  };

}

