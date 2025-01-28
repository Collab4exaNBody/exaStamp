#pragma once

#include <memory>
#include <cstdint>
#include <map>
#include <iostream>
#include <algorithm>
#include <onika/log.h>
#include <onika/cuda/cuda_math.h>
#include <onika/memory/allocator.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

namespace SnapInternal
{
  using bigint = long long;

  struct Size3 { size_t i=0,j=0,k=0; };

  class ArrayDescriptor
  {
  public:
    virtual ~ArrayDescriptor() {}
    virtual size_t element_bytes() const =0;
    virtual size_t block_size() const =0;
    virtual Size3 dim() const =0;
    virtual const char* name() const =0;

    inline size_t allocation_bytes() const
    {
      const Size3 size = dim();
      const size_t ebytes = element_bytes();
      const size_t elements = size.i * size.j * size.k;
      const size_t bs = block_size();
      const size_t data_bytes = ( elements * ebytes * bs + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( size.k * sizeof(void**) ) + ( size.k * size.j * sizeof(void*) );
      const size_t alloc_elements = ( data_bytes + pointer_bytes + ebytes-1 ) / ebytes;
      return alloc_elements * ebytes;
    }

    template<class StreamT>
    inline StreamT& print(StreamT& out) const
    {
      const Size3 size = dim();
      const size_t ebytes = element_bytes();
      const size_t elements = size.i * size.j * size.k;
      const size_t bs = block_size();
      const size_t data_bytes = ( elements * ebytes * bs + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( size.k * sizeof(void**) ) + ( size.k * size.j * sizeof(void*) );
      out << "alloc "<<name()<<" , bs="<<bs<<" , elements=";
      if(size.k>1) out<<size.k<<"x";
      if(size.j>1) out<<size.j<<"x";
      out << size.i;
      if( size.k>1 || size.j>1 ) out <<"="<<elements;
      out <<" , data_bytes="<<data_bytes<<" , pointer_bytes="<<pointer_bytes << " , total="<<allocation_bytes() <<std::endl;
      return out;
    }
  };

  template<class T>  
  class ArrayDescriptorImpl : public ArrayDescriptor
  {
    size_t m_idim = 0;
    size_t m_jdim = 0;
    size_t m_kdim = 0;
    size_t m_block_size = 0;
    //T * m_iptr = nullptr; // pointer to flat array valid for ptr[idx] with idx>=0 && idx<(dim[0]*dim[1]*dim[2])
    T * m_iptr = nullptr;
    T ** m_jptr = nullptr; // arrays of kdim*jdim pointers to i rows
    T *** m_kptr = nullptr; // array of kdim pointers to jdim pointers
    char m_name[32] = { '\0' };

  public:
    inline size_t element_bytes() const override final { return sizeof(T); }
    inline Size3 dim() const override final { return { m_idim , m_jdim , m_kdim }; }
    inline const char* name() const override final { return m_name; }
    inline size_t block_size() const override final { return m_block_size; }
       
    inline T* ptr(size_t i=0, size_t j=0, size_t k=0)
    {
      assert( i < m_idim && j < m_jdim && k < m_kdim );
      return m_iptr + ( ( ( k * m_jdim + j ) * m_idim + i ) * m_block_size );
    }
    inline T** jptr(size_t j=0, size_t k=0)
    {
      assert( j < m_jdim && k < m_kdim );
      return m_jptr + ( k * m_jdim + j );
    }
    inline T***  kptr(size_t k=0)
    {
      assert( k < m_kdim );
      return m_kptr + k;
    }

    inline ArrayDescriptorImpl( const char* s, size_t bs, size_t idim, size_t jdim=1, size_t kdim=1 )
    : m_idim(idim), m_jdim(jdim), m_kdim(kdim), m_block_size(bs)
    {
      std::strncpy(m_name,s,32); m_name[31]='\0';
      const size_t elements = m_idim * m_jdim * m_kdim;
      const size_t data_bytes = ( elements * sizeof(T) * m_block_size + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( m_kdim * sizeof(T**) ) + ( m_kdim * m_jdim * sizeof(T*) );
      const size_t alloc_elements = ( data_bytes + pointer_bytes + sizeof(T)-1 ) / sizeof(T);
      assert( alloc_elements*sizeof(T) == allocation_bytes() );
      
      m_iptr = onika::memory::CudaManagedAllocator<T>{}.allocate( alloc_elements );
      std::memset( m_iptr, 0, alloc_elements * sizeof(T) );
      
      m_jptr = reinterpret_cast<T**>( reinterpret_cast<uint8_t*>(m_iptr) + data_bytes );
      m_kptr = reinterpret_cast<T***>( m_jptr + ( m_kdim * m_jdim ) );

      for(size_t k=0;k<m_kdim;k++)
      for(size_t j=0;j<m_jdim;j++)
      {
        m_jptr[ ( k * m_jdim ) + j ] = ptr(0,j,k);
      }
      for(size_t k=0;k<m_kdim;k++)
      {
        m_kptr[k] = m_jptr + ( k * m_jdim );
      }
      
      for(size_t k=0;k<m_kdim;k++)
      for(size_t j=0;j<m_jdim;j++)
      for(size_t i=0;i<m_idim;i++)
      {
        assert( ptr(i,j,k) == & m_kptr[k][j][i* m_block_size] );
      }
    }
    
    inline ~ArrayDescriptorImpl()
    {
      const size_t elements = m_idim * m_jdim * m_kdim;
      const size_t data_bytes = ( elements * sizeof(T) * m_block_size + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( m_kdim * sizeof(T**) ) + ( m_kdim * m_jdim * sizeof(T*) );
      const size_t alloc_elements = ( data_bytes + pointer_bytes + sizeof(T)-1 ) / sizeof(T);
      assert( alloc_elements*sizeof(T) == allocation_bytes() );
      onika::memory::CudaManagedAllocator<T>{}.deallocate( m_iptr , alloc_elements );
      m_iptr = nullptr;
      m_jptr = nullptr;
      m_kptr = nullptr;
      m_idim = 0;
      m_jdim = 0;
      m_kdim = 0;
      m_block_size = 0;
    }
  };

  struct Memory
  {
    size_t m_block_size = 1;
    std::map< void* , std::shared_ptr<ArrayDescriptor> > m_allocs;

    template<class T>
    inline void destroy(T* p)
    {
      m_allocs.erase( (void*) p );
    }
    
    template<class T>
    inline void create(T * __restrict__ &p, size_t idim, const char* name )
    {
      m_allocs.erase( (void*) p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, m_block_size, idim );
      p = array->ptr();
      m_allocs[(void*)p] = array;
    }

    template<class T>
    inline void create(T *  __restrict__ *  __restrict__ &p, size_t jdim, size_t idim, const char* name )
    {
      m_allocs.erase( (void*) p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, m_block_size, idim, jdim );
      p = ( T* __restrict__ * ) array->jptr();
      m_allocs[(void*)p] = array;
    }

    template<class T>
    inline void create(T *  __restrict__ *  __restrict__ *  __restrict__ &p, size_t kdim, size_t jdim, size_t idim, const char* name )
    {
      m_allocs.erase( (void*) p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, m_block_size, idim, jdim, kdim );
      p = ( T* __restrict__ * __restrict__ * ) array->kptr();
      m_allocs[(void*)p] = array;
    }
    
    template<class StreamT>
    inline StreamT& print(StreamT& out) const
    {
      std::vector< std::shared_ptr<ArrayDescriptor> > desc;
      for(const auto& ap : m_allocs) { desc.push_back( ap.second ); }
      std::sort( desc.begin() , desc.end() , []( const std::shared_ptr<ArrayDescriptor>& a, const std::shared_ptr<ArrayDescriptor>& b ) -> bool { return a->allocation_bytes() > b->allocation_bytes(); } );
      out << "Allocation summary (block_size="<<m_block_size<<")"<<std::endl;
      size_t total = 0;
      for(const auto& array : desc)
      {
        total += array->allocation_bytes();
        array->print( out );
      }
      out << "Total="<<total<<std::endl;
      return out;
    }
  };


}

