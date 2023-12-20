#pragma once


template<int N> struct SetArraysHelper;

template<int N>
struct SetArraysHelper
{

	template<typename... T>
	static inline void set_arrays( size_t i, const std::tuple<T...>& value, const std::tuple<T*...>& arrays )
	{
		std::get<N-1>( arrays ) [i] = std::get<N-1>( value );
		SetArraysHelper<N-1>::set_arrays( i , value , arrays );
	}
};

template<>
struct SetArraysHelper<0>
{
	template<typename... T>
	static inline void set_arrays( size_t i, const std::tuple<T...>& value, const std::tuple<T*...>& arrays ) {}
};

template<typename... T>
static inline void set_arrays( size_t i, const std::tuple<T...>& value, const std::tuple<T*...>& arrays )
{
	SetArraysHelper< sizeof...(T) >::set_arrays( i, value, arrays );
};

// 
template<class T> inline void swapData (size_t i, size_t j, T* arr)
{
  T tmp  = arr[j];
  arr[j] = arr[i];
  arr[i] = tmp;
};


//
template<class Elem, class... T>
inline void swapData(size_t i, size_t j, Elem* arr, T*... data)
{
  swapData(i,j, arr);
	swapData(i,j,data...);
};


/* 
* Quicksort algorithm modified for Introsort 
*/  

template<class SORTED, class... T> inline size_t partition(size_t lo, size_t hi, SORTED pivot, SORTED* toSort, T*... data)  
{  
	size_t i=lo;
	size_t j=hi;  

	while (1)  
	{  
		while (toSort[i] < pivot) 
			i++;

		j=j-1;  
		while (pivot < toSort[j]) 
			--j;  

		if(!(i < j))  
			return i;

		swapData(i,j,toSort,data...);  
		i++;  
	}  
}  
      
template<class SORTED> SORTED medianof3(size_t lo, size_t mid, size_t hi, SORTED* toSort)  
{  

	assert(lo>=0);
	assert(mid>=0);
	assert(hi>=0);

	if (toSort[mid] < toSort[lo])  
	{  
		if (toSort[hi]<toSort[mid])  
			return toSort[mid];  
		else  
		{  
			if (toSort[hi] < toSort[lo])  
				return toSort[hi];  
			else  
				return toSort[lo];  
		}  
	}  
	else  
	{  
		if (toSort[hi]< toSort[mid])  
		{  
			if (toSort[hi]<toSort[lo])  
				return toSort[lo];  
			else  
				return toSort[hi];  
		}  
		else  
			return toSort[mid];  
	}  
}


template<class SORTED, class... T> inline void downheap(size_t i, size_t n, size_t lo, SORTED* toSort, T*... data)  
{  

	assert(lo+i-1>=0);

	SORTED pivot(toSort[lo+i-1]);
	std::tuple<SORTED,T...> tmp = std::make_tuple(toSort[lo+i-1],data[lo+i-1]...);

	size_t child;

	while (i<=n*0.5)  
	{  
		child = 2*i;  
		if (child < n && (toSort[lo+child-1]< toSort[lo+child]) )  child++;  

		if (!(pivot < toSort[lo+child-1])) break;  
    
		swapData(lo+i-1,lo+child-1,toSort,data...);  
		i = child;  
	}  
      
	set_arrays(lo+i-1, tmp, std::make_tuple(toSort,data...));

} 


template<class SORTED, class... T>  inline void heapsort(size_t lo, size_t hi, SORTED* toSort, T*... data)  
{  
	assert(hi>lo);
	size_t n = hi-lo;  

	for (size_t i=n*0.5; i>=1; i=i-1)  
	{  
		downheap(i,n,lo,toSort,data...);  
	}  

	for (size_t i=n; i>1; i=i-1)  
	{  
	    	swapData(lo,lo+i-1,toSort,data...);  
		downheap(1,i-1,lo,toSort,data...);  
	}  
} 

template<class... T, class SORTED> inline void introsort_loop (size_t lo, size_t hi, size_t depth_limit, SORTED* toSort, T*... data)  
{  
	assert(hi>=lo);
	while (hi-lo > 16)  
	{  
		if (depth_limit == 0)  
		{  
			heapsort(lo,hi, toSort, data...);
			return;  
		}  
    
		depth_limit=depth_limit-1;  
		size_t p=partition(lo, hi, medianof3(lo, lo+((hi-lo)*0.5)+1, hi-1, toSort),  toSort, data...);  
		introsort_loop(p, hi, depth_limit, toSort, data...);  
		hi=p;  
	}  
}  






template<class SORTED, class... T> inline void insertSort(size_t n, SORTED* toSort, T*... data)  
{  
	for (size_t i=1 ; i < n ; ++i)  
	{  
		size_t j;
		SORTED pivot(toSort[i]);
		std::tuple<SORTED,T...> tmp = std::make_tuple(toSort[i],data[i]...);
    
		for (j = i; j > 0 && toSort[j-1] >pivot; --j)
			swapData(j,j-1, toSort, data...);

		set_arrays(j, tmp, std::make_tuple(toSort,data...));
	} 
}


 
  




template <class SORTED, class... T>
inline void introsort( const size_t n, SORTED* toSort, T*... data)  
{  
 	/* Check that the number of atoms is not inferior to 0 */
	assert(n>=0);

	/* More costly if atoms are already sorted or almost */
	/* This algo complexity is in O(n log(n)) */
	/* Intro sort is a mixed between heap, quick and insert sort */

	if(n>50) 
		introsort_loop(0, n, 2*floor(n), toSort, data...);
        
	/* Complexity O(n²) */
	insertSort(n, toSort , data...); 
}  

