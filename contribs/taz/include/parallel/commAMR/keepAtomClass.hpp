#pragma once
#include <parallel/commAMR/commUtils.hpp>

class keepAtom{

  public :

  size_t m_nbElem;
  size_t m_sizeOfData;
  std::vector<infoCellGhost> m_info;
  std::vector<char> m_data;
  
  public :
  
  keepAtom() {}
  ~keepAtom() {}

  inline std::vector<infoCellGhost>& getInfo()
  {
    return m_info;
  }
  
  /// @brief Define the size of memory per particle to fit buffer size
  /// @tparam T variadic template of arrays that will be copied into the buffer
  /// @param [in] data use to get information on data that will be copied
  template<class... T> inline void defineSizeOfData(T*... data)
  {
    m_sizeOfData = sizeOfData(data...);

    assert(m_sizeOfData > 0); // be > 0

    if(m_sizeOfData == 0) // only defined during the init
    {
      std::cerr << " The description of the particles for the cell update was not performed correctly " << std::endl;
      abort();
    }
  }
  
  /// @brief  Clear members
  inline void clear()
  {
    m_info.clear();
    m_data.clear();
    m_nbElem = 0;
  }
    
  inline void computeShift()
  {
    for(int i = 1; i < m_info.size(); i++)
    {  
      m_info[i].shift = m_info[i-1].nbElem + m_info[i-1].shift;
    }
    
    if(m_info.size() > 0)
      m_nbElem =  m_info[m_info.size()-1].nbElem + m_info[m_info.size()-1].shift;
    else
      m_nbElem = 0;
  }
  
  inline void allocateBuffer()
  {
    size_t totalSizeOfData = m_nbElem * m_sizeOfData;
    m_data.resize(totalSizeOfData);
  }
  
  template<class ... T>
  inline void fillBuffer( size_t shift, size_t nbElem, T... data)
  {
    assert(shift >= 0);
    assert(nbElem != 0);
    assert(shift + nbElem <= m_nbElem );
    packCopy(m_data.data(), shift, 0, nbElem, m_nbElem, data...);
  }
  
    template<class ... T>
  inline void decodeBuffer( size_t shift, size_t nbElem, T... data)
  {
    assert(shift >= 0);
    assert(nbElem != 0);
    assert(shift + nbElem <= m_nbElem);
    unpackCopy(m_data.data(), shift, nbElem, m_nbElem, data...);
  }

};
