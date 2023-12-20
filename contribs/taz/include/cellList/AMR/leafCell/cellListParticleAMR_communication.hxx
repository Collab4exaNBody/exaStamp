/// @brief Copy the embedding terms to another cell of the grid
/// @param [in] to Recipient cell
template <class CList> inline void leafCell::embCopyAMR(CList& to) const {

  assert(shift + size <= granny->size);
  granny->m_eamStorage.embCopy_additional_storage(to.granny->m_eamStorage, shift, size);
}

inline void leafCell::m_eamStorageReset()
{
  granny->m_eamStorage.reset();
}


/// @brief Copy embedding terms of the particles into an array to send to other domains
/// @param [in] destDomain Recipient domain
/// @param [in] destCell Recipient cell
/// @param [in,out] leaving Array to send
inline void leafCell::embCopyAMR(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, ExchangeEAM> >& leaving) const {

  // Set the recipient domain in the object to send
  auto tmp = std::make_tuple(destDomain, ExchangeEAM());

  // Set the recipient cell in the object to send
  std::get<1>(tmp).cell = destCell;

  auto& id_ = std::get<1>(tmp).id;
  auto& emb = std::get<1>(tmp).emb;

  // Loop on the particles of the cell
  for (uint i=0; i<size; ++i) {
  
    // Set the index and embedding term of the particle in the object to send
    id_ = getId(i);
    emb = granny->m_eamStorage.emb(i);

    // Add the object to the array
    leaving.push_back(tmp);

  }

}




