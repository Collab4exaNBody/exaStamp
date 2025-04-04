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
    std::vector<PaceThreadContext> m_thread_ctx;
    ACEImpl* aceimpl = nullptr;
  };

}
