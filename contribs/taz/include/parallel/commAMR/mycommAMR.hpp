#pragma once

#include <parallel/commAMR/commUtils.hpp>


/// @brief Intermediate class to fill in the inCellGhost class from the AMR grid
class infoOctreeGhost{

  public :
  int nbElem;         ///< number of element per cell
  int cellIndex;      ///< conserve the octree index for pack
  vec3<int> destCell; ///< cartesian position on the grid
  int shift;          ///< displacement on the octree storage 

  /// @brief constructor
  infoOctreeGhost() {}

  /// @brief destructor
  ~infoOctreeGhost() {}
};


/// @brief Class used to describe the memory sections to be copied into the buffers when transferring phantom particles
class infoCellGhost {

  public :
  int nbElem;         ///< number of element per cell
  vec3<int> destCell; ///< cartesian position on the grid
  int shift;          ///< displacement on the buffer storage

  /// @brief constructor
  infoCellGhost() {}

  /// @brief destructor
  ~infoCellGhost() {}

  // Use to define the mpi data structure
  // Should not be modified

  /// @brief get the number of type
  int numberOfType () { 
    return 1; // only integer datatype
  } 

  /// @brief get the mpi data type
  MPI_Datatype dataType() {return MPI_INT;}

  /// @brief get the number of element for the mpi data type 
  int numberOfElem(){ 
    return 5; // size of the infCellGhost structure 
  } 
};


typedef std::map<size_t,std::vector<infoCellGhost>> mapInfo;
typedef std::map<size_t,std::vector<infoOctreeGhost>> mapInfoOctree;
typedef std::vector<std::pair<size_t , size_t>> vectorOfPair;
typedef std::vector<vectorOfPair> vector2DofPair;
typedef std::vector<infoCellGhost> recvInfo;


/// @brief Function allowing to switch from octrees to classical cells while keeping the information in the copy vector
/// @param [out] mapInfoCellGhost structure to use to pack the information
/// @param [in]  mapInfoOctreeGhost contains information of memory section to be copied (leafcells)
/// @param [out] copy vector of vector of pair that contain the first and last index of infoOctreeGhost corresponding of infoCellGhost (used for pack)
/// @param [in]  domain contains the numbers of the neighboring domains
void octreeToCell(mapInfo& mapInfoCellGhost,  const mapInfoOctree& mapInfoOctreeGhost, vector2DofPair& copy, const std::set<size_t>& domain)
{
  copy.clear();
  mapInfoCellGhost.clear();
  copy.resize(mapInfoOctreeGhost.size());
  size_t tmp_it_data [domain.size()];
  size_t tmp_it_size = 0;

  for(size_t node:domain)
  {
    auto sendIt = mapInfoOctreeGhost.find(node);
    if(sendIt != mapInfoOctreeGhost.end())
    {
       
      if( sendIt->second.size() > 0) {
        tmp_it_data[tmp_it_size] = node;
        tmp_it_size++;
      }
    }    
  }

  spin_mutex mutex;

  // parallelisation on mpi message 
  parallel_region(0, tmp_it_size , [&](const size_t begin, const size_t end) 
  { 

  for(size_t it = begin ; it < end; it++)
  {
    // temporary structs to fill vectors
    infoCellGhost tmp_sendInfo;
    std::pair<int,int> tmp_copy; 
    std::vector<infoCellGhost> tmp_sendInfo_vec;
    std::vector<std::pair<int,int>> tmp_copy_vec;
    size_t node = tmp_it_data[it];
    auto sendIt = mapInfoOctreeGhost.find(node);
    const infoOctreeGhost * ptr_m_send = sendIt->second.data();

    assert(sendIt->second.size() > 0);

    // fill first value;
    tmp_sendInfo. destCell = ptr_m_send[0]. destCell;
    tmp_sendInfo. nbElem   = ptr_m_send[0]. nbElem;
    tmp_sendInfo. shift    = 0;

    tmp_copy.first=0;
    size_t size = sendIt->second.size();

    for(size_t i = 1; i < size; i++)
    {
      if(tmp_sendInfo.destCell == ptr_m_send[i].destCell)
      {
        tmp_sendInfo.nbElem += ptr_m_send[i].nbElem;
      }
      else
      {
        tmp_copy.second=i;
        tmp_copy_vec.push_back(tmp_copy);
        tmp_copy.first=i;

        tmp_sendInfo_vec.push_back(tmp_sendInfo);
        tmp_sendInfo. shift   += tmp_sendInfo.  nbElem;
        tmp_sendInfo. nbElem   = ptr_m_send[i]. nbElem;
        tmp_sendInfo. destCell = ptr_m_send[i]. destCell;
      }
    }
    tmp_copy.second=size;

    // fill last value;
    tmp_copy_vec.push_back(tmp_copy);

    tmp_sendInfo_vec.push_back(tmp_sendInfo);

    copy[it].assign(tmp_copy_vec.begin(),tmp_copy_vec.end());

    mutex.lock();
      mapInfoCellGhost[node].assign(tmp_sendInfo_vec.begin(),tmp_sendInfo_vec.end());
    mutex.unlock();
  }

  });
}


/// @brief Class used to manage ghost information (pack/unpack data and mpi communications)
class exchangeGhost{

  public:

  MPI_Comm comm;                              ///< Communicator for all nodes
  int numberOfNodes;                          ///< Number of nodes on the communicator
  int rank;                                   ///< Rank of the current node
  mapInfo m_sendInfo;                         ///< information used to pack data into the buffer
  recvInfo m_recvInfo;                        ///< information used to unpack data of the buffer 
  std::vector<size_t> sendCommCount;          ///< used to count the number of infoCellGhost to send
  std::vector<size_t> recvCommCount;          ///< used to count the number of infoCellGhost to recv

  //découpage des messages.
  size_t m_BytesPerMsg;
  std::vector<size_t> m_nbOfMsgPerMsgSend;    
  std::vector<size_t> m_nbOfMsgPerMsgRecv;
  
  // Buffer information
  size_t nbMsgSend;                           ///< number of mpi message sent 
  std::vector<size_t> m_nbAtomSend;           ///< number of atom sent per message
  std::vector<size_t> m_destNode;             ///< domain number 
  std::vector<std::vector<char>> m_buffSend;  ///< send buffer
  std::vector<size_t> m_shift;                ///< displacement on the send buffer toward the domain number where the buffer is sent
  size_t nbMsgRecv;                           ///< number of mpi message received
  size_t m_nbAtomRecvTotal;                   ///< total number of atom received per message
  std::vector<size_t>  m_nbAtomRecv;          ///< number of atom received per message
  std::vector<size_t>  m_srcNode;             ///< domain number 
  std::vector<char> m_buffRecv;               ///< receive buffer
  std::vector<size_t>  m_shiftRecvAtom;       ///< displacement on the receive buffer toward the domain number where the buffer is received
  size_t m_sizeData;                          ///< size of data transfert in Byte

  // MPI
  MPI_Datatype MPI_structInfoCellGhost;       ///< MPI data structure of infoCellGhost = 5 int
  std::vector<MPI_Request> reqSend;
  std::vector<MPI_Request> reqRecv;
  std::vector<MPI_Status> staSend;
  std::vector<MPI_Status> staRecv;

	/// @brief constructor need a MPI_Comm
	exchangeGhost(MPI_Comm Comm) :  comm(Comm), m_sendInfo(), m_recvInfo(), nbMsgSend(0), m_nbAtomSend(0), m_destNode(), m_buffSend(), m_shift(), nbMsgRecv(0), m_nbAtomRecvTotal(0), m_nbAtomRecv(), m_srcNode(), m_buffRecv(), m_shiftRecvAtom(), m_sizeData(0), reqSend(), reqRecv(), staSend(), staRecv()
	{
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &numberOfNodes);
		sendCommCount = std::vector<size_t>(numberOfNodes,0); 
		recvCommCount = std::vector<size_t>(numberOfNodes,0); 

	//unlimited sized
	m_BytesPerMsg = (size_t(1)<<31);

		//assume that infoCellGhost is only composed of integers
		infoCellGhost getInfomation;

		assert(getInfomation.numberOfType()==1); 

		int numberOfElem = getInfomation.numberOfElem();//assumption getInfomation.numberOfType()==1
		MPI_Datatype dataType = getInfomation.dataType(); //assumption
		MPI_Aint displ = 0;

		// Declare the MPI data structure used to send/recv m_sendInfo/m_recvInfo -> adjustComm
		MPI_Type_create_struct( getInfomation.numberOfType(), &numberOfElem, &displ, &dataType, &MPI_structInfoCellGhost);
		MPI_Type_commit(&MPI_structInfoCellGhost);
	}

  /// @brief destructor
  ~exchangeGhost() {}


  uint64_t getNumberOfParticlesSend(){

	uint64_t tmp=0;
	for(size_t numNode = 0 ; numNode < nbMsgSend ; numNode++)
		tmp+=m_nbAtomSend[numNode];
	return tmp;
  }

  uint64_t getNumberOfParticlesRecv(){
	return uint64_t(m_nbAtomRecvTotal);
  }

  /// @brief return information on cell that are copied
  mapInfo& getMapData() { 
    return m_sendInfo;
  }

  /// @brief returns the number of mpi messages receive
  inline size_t getNbMsgRecv() {
    return nbMsgRecv;
  } 

  /// @brief returns the number of mpi messages send
  inline size_t getNbMsgSend() {
    return nbMsgSend;
  } 

  /// @brief returns the vector of inCellGhost 
  inline recvInfo & getRecvInfo() {
    return m_recvInfo;
  }

  /// @brief get the first element on the send buffer for the message number msgId
  /// @param [in] msgId message number
  inline size_t getBeginIndexBuffer(size_t msgId) {
    return m_shift[msgId];
  }

  /// @brief get the last element on the send buffer for the message number msgId
  /// @param [in] msgId message number
  inline size_t getEndIndexBuffer(size_t msgId) {
    return m_shift[msgId] + recvCommCount[m_srcNode[msgId]];
  }

  /// @brief get the domain number of the message number msgId
  /// @param [in] msgId message number
  inline size_t getDestNode(size_t msgId) { 
    return m_destNode[msgId];
  }

  /// @brief Define the size of memory per particle to fit buffer size
  /// @tparam T variadic template of arrays that will be copied into the buffer
  /// @param [in] data use to get information on data that will be copied
  template<class... T> inline void defineSizeOfData(T*... data)
  {
    m_sizeData = sizeOfData(data...);

    assert(m_sizeData > 0); // be > 0

    if(m_sizeData == 0) // only defined during the init
    {
      std::cerr << " The description of the particles for the cell update was not performed correctly " << std::endl;
      abort();
    }
  }

  /// @brief defineSizeOfDataPerMsg has to be used after defineSizeOfData
  /// @param [in] numberOfElemPerMsg define the maximum number of elem send by mpi message
  inline void defineSizeOfDataPerMsg(size_t numberOfElemPerMsg)
  {
    assert(numberOfElemPerMsg > 1);
    m_BytesPerMsg = numberOfElemPerMsg * m_sizeData;
  }

  /// @brief Clear members
  inline void clear()
  {
    m_sendInfo. clear();
    m_recvInfo. clear();
    m_srcNode.  clear();
    m_destNode. clear();
    m_nbOfMsgPerMsgRecv.clear();
    m_nbOfMsgPerMsgSend.clear();

		sendCommCount.assign(numberOfNodes,0); 
		recvCommCount.assign(numberOfNodes,0); 

    nbMsgSend = 0;
    nbMsgRecv = 0;
  }

  /// @brief In this function we propagate the information of the number of infoCellGhost to send and to whom we send it. 
  /// In addition, the send buffers are allocated.
  /// @param [in] domain contains the numbers of the neighboring domains
  inline void defineCommunication(const std::set<size_t>& domain)
  {
    // get the number of elements to send toward the node
    for(auto numNode: domain)
    {
      mapInfo::iterator it = m_sendInfo.find(numNode);

      if(it != m_sendInfo.end())
      {
        sendCommCount[numNode] = it->second.size();
        if(sendCommCount[numNode] > 0)
        {
          nbMsgSend += 1;
          m_destNode. push_back(it->first);
        }
      }
      else
        sendCommCount[numNode] = 0;
    }

    MPI_Status  sta[domain.size()*2];
    MPI_Request req[domain.size()*2];

    size_t i=0;
    
    const int tagDefine=667;

    for(size_t numNode: domain)
    {
      MPI_Isend(
        &sendCommCount[numNode], 
        1, 
        MPI_INT, 
        numNode, // send node
        tagDefine, // tag
        comm, 
        &req[2*i]
      );  

      MPI_Irecv(
        &recvCommCount[numNode], 
        1,
        MPI_INT, 
        numNode,  // src node
        tagDefine, // tag
        comm, 
        &req[2*i+1]
      );

      i++;
    }

    allocationBufferSend();

    MPI_Waitall(domain.size()*2, req, sta);

    size_t countNbOfElement = 0;

    // If their are not any elemen -> no communication
    for(size_t numNode: domain)
    {
      countNbOfElement += recvCommCount[numNode];   
      if(recvCommCount[numNode] > 0) 
      {
        nbMsgRecv+=1;
        m_srcNode.push_back(numNode);
      }
    }

    m_shift.resize(nbMsgRecv);
    m_recvInfo.resize(countNbOfElement);

    // Compute the cumsum of the number of elements (infoCellsGhost) per node (for send)
    if(nbMsgRecv > 0) 
    {
      m_shift[0] = 0;
      for(size_t i=1; i<nbMsgRecv; i++)
        m_shift[i] = m_shift[i-1] + recvCommCount[m_srcNode[i-1]];

      assert(countNbOfElement == m_shift[nbMsgRecv-1]+recvCommCount[m_srcNode[nbMsgRecv-1]]);
    }
  }

  /// @brief Sending information on memory sections to copy/extract from buffers and allocations of reception buffers
  inline void adjustCommunication()
  {
    MPI_Status  sta[nbMsgRecv+nbMsgSend];
    MPI_Request req[nbMsgRecv+nbMsgSend];
    
    const int tagAdjust=666;

    for(size_t numNode=0; numNode<nbMsgRecv; numNode++)
    {
        MPI_Irecv(
          &m_recvInfo[m_shift[numNode]], // data 
          recvCommCount[m_srcNode[numNode]], // nb elem
          MPI_structInfoCellGhost, 
          m_srcNode[numNode], // type,src node
          tagAdjust, // tag
          comm, 
          &req[numNode]
        );
    }
 
    for(size_t numNode=0; numNode<nbMsgSend; numNode++)
    {
      mapInfo::iterator it = m_sendInfo.find(m_destNode[numNode]);
 
      MPI_Isend(  
        it->second.data(), 
        it->second.size(),
        MPI_structInfoCellGhost, 
        m_destNode[numNode], // send node
        tagAdjust, // tag
        comm, 
        &req[nbMsgRecv+numNode]
      );
    }

    MPI_Waitall(nbMsgRecv+nbMsgSend, req, sta);
    
    allocationBufferRecv();
    
    computeSplitMsgIndex();    
  } 

  template<class T> //--> T = AnyGridInfo ou RectilinearGridInfo
  inline bool assertRecvInfo(T& info)
  {
	  for(int i=0; i<m_recvInfo.size(); i++)
	  {
		  if(info.index (m_recvInfo[i].destCell) < 0 ) 
		  {	
			  std::cout << " bug : l'octree " << info.index (m_recvInfo[i].destCell) << " dont la position est "<< m_recvInfo[i].destCell << " n'appartient pas a la grille " <<std::endl; 
			  return false;
		  }
	  }
	  return true;
	}
	
  inline bool assertSendInfo()
  {
    for(size_t numNode=0; numNode<nbMsgSend; numNode++)
    {
      mapInfo::iterator it = m_sendInfo.find(m_destNode[numNode]);
      
      std::vector<infoCellGhost>& vecSend = it->second;
      if(vecSend.size() > 0)
      {
        if( (vecSend[vecSend.size()-1] . shift + vecSend[vecSend.size()-1] . nbElem) !=  m_nbAtomSend[numNode] )
        {
          std::cout << " shift = " << vecSend[vecSend.size()-1] . shift 
                    << " nbElem = " << vecSend[vecSend.size()-1] . nbElem 
                    << " != m_nbAtomSend["<<numNode<<"] " << m_nbAtomSend[numNode] << std::endl;
          return false;
        }
      }  
    }
    
  	return true;
  };

  /// @brief Allocation of the send buffers
  inline void allocationBufferSend()
  {
    m_buffSend.resize(nbMsgSend);
    m_nbAtomSend.resize(nbMsgSend);

    #pragma omp parallel for
    for(size_t numNode = 0 ; numNode < nbMsgSend ; numNode++)
    {
      auto it = m_sendInfo.find(m_destNode[numNode]); 
      m_nbAtomSend[numNode] = 0;

      for(int i=0; i<it->second.size() ; i++) 
        m_nbAtomSend[numNode] += it->second[i].nbElem;

      assert( m_nbAtomSend[numNode] * m_sizeData < size_t(-1) );
      m_buffSend[numNode].resize(m_nbAtomSend[numNode] * m_sizeData);
    }
  }

  
  /// @brief Allocation of the recieving buffer
  inline void allocationBufferRecv()
  {
    m_shiftRecvAtom.resize(nbMsgRecv); 
    m_nbAtomRecv.resize(nbMsgRecv);

    int nbAtomRecvTotal  = 0;
 
    #pragma omp simd reduction(+:nbAtomRecvTotal)
    for( int it = 0 ; it < m_recvInfo.size() ; it++)
      nbAtomRecvTotal += m_recvInfo[it]. nbElem;     

    m_nbAtomRecvTotal = nbAtomRecvTotal;

    #pragma omp parallel for
    for(size_t i = 0 ; i < nbMsgRecv ; i++)
    {
      m_nbAtomRecv[i]=0;
      for(int j = m_shift[i] ; j  < m_shift[i] + recvCommCount[m_srcNode[i]] ; j++)
        m_nbAtomRecv[i] += m_recvInfo[j]. nbElem;    
    }

     
    // Compute the cumsum of the number of elements (infoCellsGhost) per node (for recv)
    if(nbMsgRecv > 0)
    {
      m_shiftRecvAtom[0] = 0;
      for(int i=1; i<nbMsgRecv ; i++)
         m_shiftRecvAtom[i] =  m_shiftRecvAtom[i-1] + m_nbAtomRecv[i-1];
    }

    assert(m_nbAtomRecvTotal * m_sizeData < size_t(-1));

    m_buffRecv.resize(m_nbAtomRecvTotal * m_sizeData); 
  }
  
  inline void computeSplitMsgIndex()
  {
    //begin by send
    m_nbOfMsgPerMsgSend.resize(nbMsgSend);
    
    for(size_t numNode=0; numNode < nbMsgSend; numNode++)
    {
      //get information
      size_t sizeBuffer = m_buffSend[numNode].size(); // array of char: size -> nb of Bytes
      
      const size_t sizeMsgMax = m_BytesPerMsg;
      if(sizeBuffer % sizeMsgMax == 0) m_nbOfMsgPerMsgSend[numNode] = sizeBuffer/sizeMsgMax;
      else m_nbOfMsgPerMsgSend[numNode] = sizeBuffer/sizeMsgMax + 1;
    }
    
    //end by recv
    m_nbOfMsgPerMsgRecv.resize(nbMsgRecv);
    
    for(size_t numNode=0; numNode < nbMsgRecv; numNode++)
    {
      //get information
      size_t sizeBuffer = m_nbAtomRecv[numNode] * m_sizeData;
      const size_t sizeMsgMax = m_BytesPerMsg;
      
      if(sizeBuffer % sizeMsgMax == 0) m_nbOfMsgPerMsgRecv[numNode] = sizeBuffer/sizeMsgMax;
      else m_nbOfMsgPerMsgRecv[numNode] = sizeBuffer/sizeMsgMax + 1;
    }
    
    size_t nbSend = 0;
    size_t nbRecv = 0;
    
    for(size_t i = 0 ; i < nbMsgSend ; i++)
      nbSend += m_nbOfMsgPerMsgSend[i];

    for(size_t i = 0 ; i < nbMsgRecv ; i++)
      nbRecv += m_nbOfMsgPerMsgRecv[i];
      
    reqRecv.resize(nbRecv);
    staRecv.resize(nbRecv);  
    
    reqSend.resize(nbSend);
    staSend.resize(nbSend);  
  } 


  /// @brief copy arrays data from index begin to begin+nbElem in the buffer msgId 
  /// @tparam T variadic template of arrays  that will be copied into the buffer
  /// @param [in] msgId MPI message number 
  /// @param [in] displ displacement on the buffer
  /// @param [in] begin index of the first element to be copied
  /// @param [in] nbElem number of element to be copied
  /// @param [in] data sources to be copied in the buffer
  template <class... T> inline void pack(const size_t msgId, const size_t displ, const size_t begin, const size_t end, T*... data)
  {
    assert(msgId < m_nbAtomSend.size());
    assert(msgId < nbMsgSend);
    assert(displ < m_nbAtomSend[msgId]); 
    assert(begin <= end);
    packCopy(m_buffSend[msgId].data(), displ, begin, end, m_nbAtomSend[msgId], data...); 
  }



  /// @brief copy information in the buffer corresponding to the msgId from shift to shift+nbElem modulo type in the data arrays
  /// @tparam T variadic template of arrays  that will be copied into the buffer
  /// @param [in] msgId MPI message number 
  /// @param [in] shift displacement on cell/octree storage (=0 for cell)
  /// @param [in] nbElem number of element to be copied
  /// @param [in] data sources to be copied in the buffer
  template <class... T> inline void unpack(const size_t msgId, const size_t shift, const size_t nbElem, T*... data)
  {
      const int shiftMsgRecv = m_shiftRecvAtom[msgId] *m_sizeData;
      int nbElemPerMsg = m_nbAtomRecv[msgId];
      const char* buffRecv = &m_buffRecv[shiftMsgRecv];

      unpackCopy(buffRecv, shift // leaf cell 
                           , nbElem //elem per octree
                           , nbElemPerMsg
                           , data...);

  }

  /// @brief send buffers
  inline void send()
  {
    size_t nN = numberOfNodes;
    size_t inc = 0;

    for(size_t i = 0; i<nbMsgSend ; i++)
    {
    
      size_t sizeOfThisMsg = m_nbAtomSend[i] * m_sizeData; 
      size_t shift = 0;  
      
      for(size_t j = 0; j < m_nbOfMsgPerMsgSend[i] ; j++)
      {
  
        size_t sizeDataSend = std::min<size_t> ( sizeOfThisMsg-shift, m_BytesPerMsg);
        
        assert(shift < sizeOfThisMsg || sizeDataSend == 0); 
        
        if(sizeDataSend == 0) continue;
        
        char* tmpBuffSend = m_buffSend[i].data() + shift; 
        
        MPI_Isend(     
          tmpBuffSend, 
          sizeDataSend, 
          MPI_CHAR,           
          m_destNode[i],
          j, 
          comm,      
          &reqSend[inc]
        );
        
        inc++;
        shift += sizeDataSend;
      }
    } 
    assert(inc == reqSend.size());
  }

  /// @brief receive buffers
  inline void recv()
  {
    size_t nN = numberOfNodes;
    size_t inc = 0;

    for(size_t i = 0 ; i<nbMsgRecv ; i++)
    {
      size_t sizeOfThisMsg = m_nbAtomRecv[i] * m_sizeData; 
      size_t shift = 0;
      
      for(size_t j = 0 ; j < m_nbOfMsgPerMsgRecv[i] ; j++)
      {
        size_t sizeDataRecv = std::min<size_t> ( sizeOfThisMsg-shift, m_BytesPerMsg);
        
        assert(shift < sizeOfThisMsg || sizeDataRecv == 0);
        
        if(sizeDataRecv == 0) continue;
        
        char* tmpBuffRecv = &m_buffRecv[m_shiftRecvAtom[i] * m_sizeData + shift];

        MPI_Irecv(tmpBuffRecv, 
          sizeDataRecv,
          MPI_CHAR, 
          m_srcNode[i],
          j,
          comm, 
          &reqRecv[inc]
        );
        
        inc++;
        shift += sizeDataRecv;
      }
    } 
    
    assert(inc == reqRecv.size());
  }
  
  /// @brief count the number of element received
  inline size_t numberOfElementsRecv()
  {
    size_t compteur=0;
    for(size_t i = 0 ; i<nbMsgRecv ; i++)
      compteur += m_nbAtomRecv[i];
    
    return compteur;
  } 

  /// @brief wait for all ghost communications
  inline void WaitAll()
  {
    MPI_Waitall(reqSend.size(), reqSend.data(), staSend.data());
    MPI_Waitall(reqRecv.size(), reqRecv.data(), staRecv.data());
  }
};

