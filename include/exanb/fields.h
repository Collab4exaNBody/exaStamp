#pragma once

#include <exaStamp/atom_bond_connectivity.h>
#include <exaStamp/sliplink_field.h>
#include <exanb/core/basic_types_def.h>
#include <exanb/core/quaternion.h>
#include <exanb/core/declare_field.h>
#include <cstdint>

// exaStamp fields are defined in namespace xstamp::field
// for rx, use xstamp::field::rx as a field descriptor
XSTAMP_DECLARE_FIELD(uint64_t        ,id                ,"particle id");
XSTAMP_DECLARE_FIELD(uint8_t         ,type              ,"particle type");
XSTAMP_DECLARE_FIELD(uint64_t        ,id_type           ,"particle encoded id and type" );
XSTAMP_DECLARE_FIELD(double          ,ep                ,"particle potential energy");
XSTAMP_DECLARE_FIELD(double          ,rx                ,"particle position X");
XSTAMP_DECLARE_FIELD(double          ,ry                ,"particle position Y");
XSTAMP_DECLARE_FIELD(double          ,rz                ,"particle position Z");
XSTAMP_DECLARE_FIELD(double          ,rxf               ,"filtered particle position X");
XSTAMP_DECLARE_FIELD(double          ,ryf               ,"filtered particle position Y");
XSTAMP_DECLARE_FIELD(double          ,rzf               ,"filtered particle position Z");
XSTAMP_DECLARE_FIELD(double          ,vx                ,"particle velocity X");
XSTAMP_DECLARE_FIELD(double          ,vy                ,"particle velocity Y");
XSTAMP_DECLARE_FIELD(double          ,vz                ,"particle velocity Z");
XSTAMP_DECLARE_FIELD(double          ,ax                ,"particle acceleration X");
XSTAMP_DECLARE_FIELD(double          ,ay                ,"particle acceleration Y");
XSTAMP_DECLARE_FIELD(double          ,az                ,"particle acceleration Z");
XSTAMP_DECLARE_FIELD(double          ,charge            ,"particle electric charge");
XSTAMP_DECLARE_FIELD(::exanb::Mat3d ,virial            ,"particle pressure virial tensor");

// fields for molecule extension
XSTAMP_DECLARE_FIELD(uint64_t                       ,idmol  ,"molecule Id");
XSTAMP_DECLARE_FIELD(::exaStamp::AtomBondConnectivity ,cmol   ,"molecule bonds");

// specific for sliplinks extension
XSTAMP_DECLARE_FIELD(::exaStamp::SlipLinkField        ,sl     ,"sliplink info");

// for rigid bodies (rigid molecules)
XSTAMP_DECLARE_FIELD(::exanb::Quaternion      ,orient   ,"rigid molecule orientation");
XSTAMP_DECLARE_FIELD(::exanb::Vec3d           ,angmom   ,"rigid molecule angular momentum");
XSTAMP_DECLARE_FIELD(::exanb::Vec3d           ,couple   ,"rigid molecule couple");

// aliases
XSTAMP_DECLARE_ALIAS( fx, ax )
XSTAMP_DECLARE_ALIAS( fy, ay )
XSTAMP_DECLARE_ALIAS( fz, az )

namespace exanb
{
  static inline constexpr bool HAS_POSITION_BACKUP_FIELDS = true;
  static constexpr auto PositionBackupFieldX = ::exanb::field::rxf;
  static constexpr auto PositionBackupFieldY = ::exanb::field::ryf;
  static constexpr auto PositionBackupFieldZ = ::exanb::field::rzf;
}

// should disapear ASAP, old depenncy to exaStamp in grid.h => particle_type method should not exist
#define EXANB_PARTICLE_TYPE_FIELD ::exanb::field::type

