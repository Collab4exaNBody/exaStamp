#pragma once

#include<functional>

class coloringOctree {

	private :

	/* uint -> index of the octree */
	/* size_t -> corresponding wave number */
	std::vector<std::pair< uint,size_t>> m_info;
	std::vector<std::pair<size_t,std::vector<size_t> >> m_dep;
	vec3<size_t> m_wave[8];
	vec3<int> m_unlock[8];

	public :

	coloringOctree() : m_info(), m_dep()
	{
		m_wave[0] = vec3<size_t> (0,0,0);
		m_wave[1] = vec3<size_t> (1,0,0);
		m_wave[2] = vec3<size_t> (1,1,0);
		m_wave[3] = vec3<size_t> (0,1,0);
		m_wave[4] = vec3<size_t> (0,1,1);
		m_wave[5] = vec3<size_t> (0,0,1);
		m_wave[6] = vec3<size_t> (1,0,1);
		m_wave[7] = vec3<size_t> (1,1,1);

		m_unlock[0] = vec3<std::size_t> (0,0,0);
		m_unlock[1] = vec3<std::size_t> (1,0,0);
		m_unlock[2] = vec3<std::size_t> (0,1,0);
		m_unlock[3] = vec3<std::size_t> (1,0,0);
		m_unlock[4] = vec3<std::size_t> (0,0,1);
		m_unlock[5] = vec3<std::size_t> (0,1,0);
		m_unlock[6] = vec3<std::size_t> (1,0,0);
		m_unlock[7] = vec3<std::size_t> (0,1,0);
	}

	size_t size();
	size_t getWaveNumber		(size_t idx);
	uint getIndex			(size_t idx);
	size_t numberOfDependencies	(int i); 
	size_t valueOfDependency	(int i ); 
	size_t*arrayOfDependencies	(int i); 

	std::vector<std::pair<uint,size_t>>::iterator begin();
	std::vector<std::pair<uint,size_t>>::iterator end();

	std::pair<uint,size_t>& operator[](size_t i);
	template<class O, class P> void resort(TraversalManager::CellTraversal T, TraversalManager& tManager, O* octree, P &info,  vec3<int> numberOfCellsPerDim, bool hybrid);
	template<class O, class P> void resortAllTraversal(TraversalManager& tManager, O* octree, P &info,  vec3<int> numberOfCellsPerDim, bool hybrid);
 
};

/* accesor */
inline std::vector<std::pair<uint, size_t>>::iterator coloringOctree::begin()
{
	return m_info.begin();
}

inline std::vector<std::pair<uint, size_t>>::iterator coloringOctree::end()
{
	return m_info.end();
}

inline std::pair< uint,size_t>& coloringOctree::operator[](size_t i)
{
	assert(i<m_info.size());
	return m_info[i];
}

inline size_t coloringOctree::getWaveNumber(size_t idx)
{
	assert( m_info[idx].second >= 0  && m_info[idx].second < 8) ;
	return m_info[idx].second;
}

inline uint coloringOctree::getIndex(size_t idx)
{
	assert(m_info[idx].first>=0);
	return m_info[idx].first;
}

inline size_t coloringOctree::numberOfDependencies(int i) 
{
	return m_dep[i].second.size();
}

inline size_t coloringOctree::valueOfDependency(int i ) 
{
	return m_dep[i].first;
}  

inline size_t * coloringOctree::arrayOfDependencies(int i)
{
	return m_dep[i].second.data();
}  

/* Return the number of element */
inline size_t coloringOctree::size()
{
	return m_info.size();
}


class isEqualTo /* used for std::find_if */
{
	public :
	uint m_tmp;
	isEqualTo(int i) : m_tmp(uint(i)) {}
	bool operator () (std::pair<uint, size_t> const& p)
	{
		return (p.first == m_tmp);
	}
};


#include <parallel/algorithm>

template<class OCTREE, class Info> 
inline void coloringOctree::resort(TraversalManager::CellTraversal T, TraversalManager& tManager, OCTREE *octree, Info &info, vec3<int> numberOfCellsPerDim, bool hybrid)
{

	auto& cells = tManager.getTraversal(T);
	m_info.clear(); /* linear in size */
	m_info.reserve(cells.size());
	assert(m_info.empty());

	for(size_t i = 0; i < cells.size(); i++) 
	{
		/* skip empty tasks */
		if(octree[cells[i]].size == 0) continue ;

		vec3<size_t> paternPosition = auxMod(octree[cells[i]].position,2);

		for(size_t nWave = 0; nWave<8 ; nWave++)
			if(paternPosition == m_wave[nWave])
				m_info.push_back(std::make_pair(cells[i], nWave));
	}

	__gnu_parallel::stable_sort(m_info.begin(), m_info.end(), 
		[] (std::pair<uint,size_t> elem_a, std::pair<uint,size_t> elem_b) -> bool {return elem_a.second< elem_b.second; }
	);

	assert(m_info.size() <= cells.size() );

	m_dep.resize(m_info.size());

	auto linear = [&] (vec3<int> o)->size_t { 
		auto tmp_o = o - info.getCellInf();
		return tmp_o.x + (tmp_o.y) * (numberOfCellsPerDim.x+2) + (tmp_o.z) * (numberOfCellsPerDim.x+2) *(numberOfCellsPerDim.y+2) ; 
	};

	/* Declare here for recursive call */
	std::function<void(vec3<int> ,int ,vec3<int>, int it_octree)> recSearch;

	recSearch = [this, &info, &linear, &recSearch] ( vec3<int> posNeighOctree, int wave, vec3<int> posOctree, int it_octree)->void
	{
		if(wave<0) return; /* no following task */

		/* Check distance : an octree is in the neighborhood if position.{x,y,z}+-1 */
		vec3<bool> diff = auxAbs(posNeighOctree-posOctree) > 1;
		if( diff.x || diff.y || diff.z) return;

		uint idx_cell = info.indexForAMR(posNeighOctree);  

		auto it = std::find_if(m_info.begin(), m_info.end(), isEqualTo(idx_cell)); /* linear complexity (distance (begin(),end())) */

		/* The distance check implies that idx_cell is a index of real or ghost octree */  
		/* Reason : all neighbor octrees of a real/inside/edge are a real/inside/edge/ghost octree*/
		/* Note that idx_cell = -1 if idx_cell correspond to a ghost octree */
		if( it != m_info.end()) /* skip empty octree and octree which are not included in the travesal */ 
		{
			this->m_dep[it_octree].second.push_back(linear(posNeighOctree));	
		}
		else
		{
			recSearch(posNeighOctree+m_unlock[wave], wave-1, posOctree, it_octree); 
			recSearch(posNeighOctree-m_unlock[wave], wave-1, posOctree, it_octree); 
		}
	};


	#pragma omp parallel for schedule(runtime)
	for(int i = 0; i<m_info.size() ; i++)
	{
		this->m_dep[i].second.clear();

		vec3<int> tmp_posOctree = octree[getIndex(i)].position ; /* Local position */
		int nWave2 = getWaveNumber(i); 
		this->m_dep[i].first = linear(tmp_posOctree); /* unlocked Task */

		nWave2--;

		recSearch(tmp_posOctree+m_unlock[getWaveNumber(i)], nWave2, tmp_posOctree, i); 
		recSearch(tmp_posOctree-m_unlock[getWaveNumber(i)], nWave2, tmp_posOctree, i); 
	}
}


template<class OCTREE, class Info> /*for MEAM*/
inline void coloringOctree::resortAllTraversal( TraversalManager& tManager, OCTREE *octree, Info &info, vec3<int> numberOfCellsPerDim, bool hybrid)
{
	auto& cells = tManager.getTraversal(TraversalManager::ALL);

	m_info.clear(); /* linear in size */
	assert(m_info.empty());

	for(size_t i = 0; i < cells.size(); i++) 
	{
		/* skip empty tasks */
		if(octree[cells[i]].size == 0) continue ;

		vec3<size_t> paternPosition = auxMod(octree[cells[i]].position,2);

		for(size_t nWave = 0; nWave<8 ; nWave++) // check wave
			if(paternPosition == m_wave[nWave])
				m_info.push_back(std::make_pair(cells[i], nWave));
	}

	__gnu_parallel::stable_sort(m_info.begin(), m_info.end(), 
		[] (std::pair<uint,size_t> elem_a, std::pair<uint,size_t> elem_b) -> bool {return elem_a.second< elem_b.second; }
	);

	assert(m_info.size() <= cells.size() );

	m_dep.resize(m_info.size());

	auto linear = [&] (vec3<int> o)->size_t { 
		auto tmp_o = o - info.getCellInf();
		return tmp_o.x + (tmp_o.y) * (numberOfCellsPerDim.x+2) + (tmp_o.z) * (numberOfCellsPerDim.x+2) *(numberOfCellsPerDim.y+2) ; 
	};

	/* Declare here for recursive call */
	std::function<void(vec3<int> ,int ,vec3<int>, int it_octree)> recSearch;

	recSearch = [this, &info, &linear, &recSearch] ( vec3<int> posNeighOctree, int wave, vec3<int> posOctree, int it_octree)->void
	{
		if(wave<0) return; /* no following task */

		/* Check distance : an octree is in the neighborhood if position.{x,y,z}+-1 */
		vec3<bool> diff = auxAbs(posNeighOctree-posOctree) > 1;
		if( diff.x || diff.y || diff.z) return;

		uint idx_cell = info.indexForAMRWithGhost(posNeighOctree);  

		auto it = std::find_if(m_info.begin(), m_info.end(), isEqualTo(idx_cell)); /* linear complexity (distance (begin(),end())) */

		/* The distance check implies that idx_cell is a index of real or ghost octree */  
		/* Reason : all neighbor octrees of a real/inside/edge are a real/inside/edge/ghost octree*/
		/* Note that idx_cell = -1 if idx_cell correspond to a ghost octree */
		if( it != m_info.end()) /* skip empty octree and octree which are not included in the travesal */ 
		{
			this->m_dep[it_octree].second.push_back(linear(posNeighOctree));	
		}
		else
		{
			recSearch(posNeighOctree+m_unlock[wave], wave-1, posOctree, it_octree); 
			recSearch(posNeighOctree-m_unlock[wave], wave-1, posOctree, it_octree); 
		}
	};


	#pragma omp parallel for schedule(runtime)
	for(int i = 0; i<m_info.size() ; i++)
	{
		this->m_dep[i].second.clear();

		vec3<int> tmp_posOctree = octree[getIndex(i)].position; /* Local position */
		int nWave2 = getWaveNumber(i); 
		this->m_dep[i].first = linear(tmp_posOctree); /* unlocked Task */

		nWave2--;

		recSearch(tmp_posOctree+m_unlock[getWaveNumber(i)], nWave2, tmp_posOctree, i); 
		recSearch(tmp_posOctree-m_unlock[getWaveNumber(i)], nWave2, tmp_posOctree, i); 
	}

}

