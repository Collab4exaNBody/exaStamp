#pragma once

#include <exaStamp/potential/snap/snap_config.h>
#include <string>

namespace SnapExt
{
  void snap_read_lammps(const std::string& paramFileName, const std::string& coefFileName, SnapConfig& config, bool conv_units = true);
}

