/// @file
/// @brief Testing compilers and versions, set flags

#ifndef __COMPILE_HPP_INCLUDED
#define __COMPILE_HPP_INCLUDED


#if defined __INTEL_COMPILER


#define __INTEL_CXX_VERSION (__INTEL_COMPILER*100 + __INTEL_COMPILER_UPDATE)

#if __INTEL_CXX_VERSION < 140002
#error This version of the code requires icc 14.0.2 or higher
#endif

#define __USE_SVML


#elif defined __GNUG__


#define __GNU_CXX_VERSION (__GNUC__*10000 + __GNUC_MINOR__*100 + __GNUC_PATCHLEVEL__)

#if __GNU_CXX_VERSION < 40800
#error This version of the code requires g++ 4.8.0 or higher
#endif


#endif

#endif // __COMPILE_HPP_INCLUDED
