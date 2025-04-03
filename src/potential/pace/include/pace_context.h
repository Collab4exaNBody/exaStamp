#pragma once

namespace exaStamp
{

  using namespace exanb;

  struct PaceThreadContext
  {
    ACEImpl * aceimpl = nullptr;
  };

  struct PaceContext
  {
    //    LAMMPS_NS::LAMMPS * ptr = nullptr;
    std::vector<PaceThreadContext> m_thread_ctx;
    PaceConfig m_config;
  };

}
