#pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

// main Onika Cuda abstraction layer :
// contains function attributes such as ONIKA_HOST_DEVICE 
// and other hybride functions working on both sides
#include <onika/cuda/cuda.h> 

// unified memory handling : defines a onika::CudaMMVector template,
// which is a std::vector with a unified memory allocator
#include <onika/memory/allocator.h>

namespace exaStamp
{
  using namespace exanb;

  // same as std::vector<double>, but allocated in unified memory
  using MyDoubleVector = onika::CudaMMVector<double>;

  // does x += y
  // function can be called from host or GPU code
  ONIKA_HOST_DEVICE_FUNC
  void array_add_value(double & x , double y)
  {
    x += y;
  }
  
  // GPU entry point function, a.k.a GPU kernel.
  // the only kind of GPU function callable from host
  ONIKA_DEVICE_KERNEL_FUNC
  void array_add_kernel(double * x_array , long n, double y)
  {
    long TotalThreads = ONIKA_CU_GRID_SIZE * ONIKA_CU_BLOCK_SIZE;
    long i = ONIKA_CU_BLOCK_IDX * ONIKA_CU_BLOCK_SIZE + ONIKA_CU_THREAD_IDX;
    // each thread in each block processes a different value
    // and progress advances by TotalThreads elements
    for( ; i < n ; i += TotalThreads )
    {
      array_add_value( x_array[i] , y );  // same function can be used on both sides
    }
  }

  // Our test component, simply add a contant value to all elements in my_array
  template< class GridT >
  class MyGpuTest : public OperatorNode
  {
    ADD_SLOT( double , value , INPUT , REQUIRED );
    ADD_SLOT( MyDoubleVector , my_array  , INPUT_OUTPUT);

  public:
    inline void execute () override final
    {
      my_array->assign(10000,1.0);

      long N = my_array->size();
      double * X_array = my_array->data();
      double Y = *value;

      auto exec_ctx = parallel_execution_context();
      bool gpu_present = exec_ctx != nullptr
                      && exec_ctx->has_gpu_context()
                      && exec_ctx->gpu_context()->has_devices();
      if( gpu_present )
      {
        lout << "GPU execution enabled" << std::endl;
        int NbBlocks = 128;   // number of thread blocks, choosen upon number of SMs available on current GPU
        int BlockSize = 32;   // number of threads per blocks, should be a multiple of 32
        ONIKA_CU_LAUNCH_KERNEL( NbBlocks, BlockSize, 0, exec_ctx->gpu_stream()   // kernel launch parameters
                              , array_add_kernel                                // kernel to launch
                              , X_array , N , Y ); // kernel function parameters
      }
      else
      {
        lout << "No GPU available, executing with OpenMP" << std::endl;
        #pragma omp parallel for
        for(long i=0;i<N;i++)
        {
          array_add_value( X_array[i] , Y ); // same function can be used on both sides
        }
      }
    }
    
  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("my_gpu_test" , make_grid_variant_operator< MyGpuTest > );
  }

}

