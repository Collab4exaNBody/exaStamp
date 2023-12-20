#pragma once
 
void defineRecDirectionInOctree(std::vector < std::set < size_t> > & inout)
{
	// To find recursively the edge leaf cells in octrees
	// Indexing of cells in an octree
	// For example Xm corresponds to the underside according to x in an octree
	size_t Xm[4] = {0,1,2,3}; 
	size_t Ym[4] = {0,1,4,5};
	size_t Zm[4] = {0,2,4,6};
	size_t Xp[4];
	size_t Yp[4];
	size_t Zp[4];

	for(size_t i = 0; i < 4; ++i)
	{
		Xp[i] = Xm[i]+4;
		Yp[i] = Ym[i]+2;
		Zp[i] = Zm[i]+1;
	} 
  
	// rename latter
	inout.resize(27);

	auto linear = [] (int x, int y, int z) -> size_t { return x+1 + 3 * (y+1) + 9 * (z+1); };


	auto fill = [] (std::set<size_t>& inout, size_t p[4]) -> void { 
		for(size_t i = 0; i< 4; i++)
			inout.insert(p[i]);
	};

	for(int z=-1; z<=1; ++z)
		for(int y=-1; y<=1; ++y)
			for(int x=-1; x<=1; ++x)
			{
				std::set<size_t> & refSet = inout[linear(x,y,z)];

				if(x == -1) fill(refSet,Xm);
				if(x == 1)  fill(refSet,Xp);
				if(y == -1) fill(refSet,Ym);
				if(y == 1)  fill(refSet,Yp);
				if(z == -1) fill(refSet,Zm);
				if(z == 1)  fill(refSet,Zp);
			}	
}

template<size_t N> struct pushNOrderHelper;

template<> 
struct pushNOrderHelper<1>
{
	template<class U, class... T> 
	static inline void pushNOrder(double time, size_t numberOfElements, U* arr1, U* arr2,  T*... arrays)
	{
		pushNOrderHelper<1>::pushNOrder(time, numberOfElements, arr1, arr2);
		pushNOrderHelper<1>::pushNOrder(time, numberOfElements, arrays...);
	}

	template<class U> static inline void pushNOrder(double time, size_t numberOfElements, U* arr1, U* arr2)
	{
		simdAMR::kernels::push1stOrder(time, arr1, arr2, numberOfElements);
	}  
};

template<> 
struct pushNOrderHelper<2>
{
	template<class U, class... T> 
	static inline void pushNOrder(double time, size_t numberOfElements, U* arr1, U* arr2, U* arr3,  T*... arrays)
	{
		pushNOrderHelper<2>::pushNOrder(time, numberOfElements, arr1, arr2, arr3);
		pushNOrderHelper<2>::pushNOrder(time, numberOfElements, arrays...);
	}
  
	template<class U> 
	static inline void pushNOrder(double time, size_t numberOfElements, U* arr1, U* arr2, U* arr3)
	{
		simdAMR::kernels::push2ndOrder(time, arr1, arr2, arr3, numberOfElements);
	}  
};

template<size_t N>
struct pushNOrderHelperAssert
{
	template<class... T> static inline void pushNOrderAssert(double time, size_t numberOfElements, T*... arrays)
	{
		assert(N>0);  // N=0 impossible
		assert(N<=2); // not implemented
		assert((sizeof...(T) % (N+1)) == 0); // logique ^^
		assert(numberOfElements < 1000000); // no-sens
		assert(time>0.0); // non-sens
		pushNOrderHelper<N>::pushNOrder(time, numberOfElements, arrays...);
	}
};




