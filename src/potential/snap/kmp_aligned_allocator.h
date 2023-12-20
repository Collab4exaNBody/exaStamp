#pragma once

#include <cstdlib>
#include <cstdint>

#ifdef ENABLE_INTEL_KMP_MALLOC
#include <omp.h>
#endif

// alocator for multi-thread friendly heap allocation
template <class T>
struct kmp_aligned_allocator
{
  static constexpr unsigned int Alignment = 64;

  typedef T value_type;
  
  kmp_aligned_allocator() noexcept {}
  template <class U> kmp_aligned_allocator (const kmp_aligned_allocator<U>&) noexcept {}
  
  inline T* allocate (std::size_t n)
  {
    size_t bytes = n*sizeof(T) + Alignment;
#   ifdef ENABLE_INTEL_KMP_MALLOC
    uint8_t* alloc_addr = (uint8_t*) kmp_malloc( bytes );
#   else
    uint8_t* alloc_addr = new uint8_t[ bytes ];
#   endif
    size_t addri = reinterpret_cast<size_t>( alloc_addr );
    size_t offset = Alignment - ( addri % Alignment );
    //std::cout<<"NEW: bytes="<<bytes<<", offset="<<offset<<" @"<< ((void*)alloc_addr) << std::endl;
    if( offset==0 || offset>255 )
    {
      std::abort();
    }
    *alloc_addr = static_cast<uint8_t>(offset);
    for(size_t i=1;i<offset;i++) { alloc_addr[i] = 0; }
    return reinterpret_cast<T*>( alloc_addr + offset );
  }
  
  inline void deallocate (T* p, std::size_t n)
  {
//    size_t bytes = n*sizeof(T) + Alignment;
    int offset = 1;
    uint8_t* alloc_addr = reinterpret_cast<uint8_t*>( p );
    while( *(alloc_addr-offset) != offset && offset<255) { ++offset; }
    if( offset > 255 )
    {
      std::abort();
    }
    //std::cout<<"FRE: bytes="<< bytes <<", offset="<<offset<<" @"<< ((void*)(alloc_addr-offset)) << std::endl;
#   ifdef ENABLE_INTEL_KMP_MALLOC
    kmp_free( alloc_addr - offset );
#   else
    delete [] ( alloc_addr - offset );
#   endif
  }
  
  template<typename U>
  inline bool operator != (const U&) const { return true; }

  inline bool operator != (const kmp_aligned_allocator&) const { return false; }
};

