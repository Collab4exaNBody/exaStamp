#include <exanb/core/log.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <memory>

#undef _CLASS_NAME
#undef _CLASS_BASE

#ifdef USTAMP_POTENTIAL_WITH_VIRIAL
#define _CLASS_NAME PairPotentialCompboundVirialOperator
#define _CLASS_BASE PairPotentialComputeVirialOperator
#define _VIRIAL_PARAM ,Mat3d& virial
#define _VIRIAL_ARG ,virial
#else
#define _CLASS_NAME PairPotentialCompboundOperator
#define _CLASS_BASE PairPotentialComputeOperator
#define _VIRIAL_PARAM
#define _VIRIAL_ARG
#endif


namespace exaStamp
{
  using namespace exanb;

  class _CLASS_NAME final : public _CLASS_BASE
  {
  public:

    inline _CLASS_NAME( const std::vector< std::shared_ptr<_CLASS_BASE> >& ops, const std::vector<double>& rcuts )
    {
      m_ops.clear();
      if( rcuts.size() != ops.size() )
      {
        lerr << "rcut list size inconsitent with operators list size" << std::endl;
        std::abort();
      }
      m_name.clear();
      size_t n = ops.size();
      for(size_t i=0;i<n;i++)
      {
        m_ops.push_back( { rcuts[i]*rcuts[i] , ops[i] } );
        if(i>0) m_name += "+";
        m_name += ops[i]->name();
      }
      std::sort( m_ops.begin(), m_ops.end() , std::greater<>() );
    }

    ~_CLASS_NAME() override final = default;

    inline void set_parameters(const PairPotentialParameters& params) override final
    {
      for(auto& op:m_ops) op.second->set_parameters( params );
    }

    inline void set_rcut( double rcut ) override final
    {
      m_rcut2 = rcut*rcut;
    }

    inline uint64_t signature() const override final
    {
      std::vector<uint64_t> hdata;
      for(const auto& op:m_ops)
      {
        hdata.push_back( op.second->signature() );
      }
      const char* s = reinterpret_cast<const char*>( hdata.data() );
      std::string hstr( s , hdata.size()*sizeof(uint64_t) );
      return std::hash<std::string>{}(hstr);
    }

    inline const std::string& name() const override final
    {
      return m_name;
    }

    template<bool UseWeights>
    inline void doit (ComputePairBuffer2<UseWeights,false>& tab,double& ep, double& ax, double& ay, double& az _VIRIAL_PARAM) const noexcept
    {
      if( m_ops.empty() ) return;
      assert( m_rcut2 >= m_ops.front().first );
      double last_rcut2 = m_rcut2;
      for(auto& op:m_ops)
      {
        if( op.first < last_rcut2 )
        {
          last_rcut2 = op.first;
          for(int i=0;i<tab.count;)
          {
            if( tab.d2[i] > last_rcut2 ) tab.copy( -- tab.count , i );
            else ++i;
          }
        }
        // std::cout << "rcut2="<<last_rcut2<<" count="<<tab.count<<std::endl<<std::flush;
        (*op.second) (tab,ep,ax,ay,az _VIRIAL_ARG );
      }
    }

    inline void operator() (ComputePairBuffer2<false,false>& tab,double& ep, double& ax, double& ay, double& az _VIRIAL_PARAM ) const noexcept override final
    {
      doit(tab,ep,ax,ay,az _VIRIAL_ARG);
    }
    inline void operator() (ComputePairBuffer2<true,false>& tab,double& ep, double& ax, double& ay, double& az _VIRIAL_PARAM ) const noexcept override final
    {
      doit(tab,ep,ax,ay,az _VIRIAL_ARG);
    }

  private:
    std::vector< std::pair< double , std::shared_ptr<_CLASS_BASE> > > m_ops;
    double m_rcut2 = 0.0;
    std::string m_name;
  };

}

#undef _CLASS_NAME
#undef _CLASS_BASE
#undef _VIRIAL_PARAM
#undef _VIRIAL_ARG

