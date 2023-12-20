#pragma once

#include "kmp_aligned_allocator.h"
#include <vector>

namespace SnapExt
{

template<typename T> using SnapVectorT = std::vector< T , kmp_aligned_allocator<T> >;

}

#define KMP_ASSUME_ALIGNED(x) x = ( decltype(x) __restrict__ ) __builtin_assume_aligned( x , ::kmp_aligned_allocator<int>::Alignment )

