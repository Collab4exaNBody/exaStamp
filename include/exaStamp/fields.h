/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

//#include <exanb/core/grid_fields.h>
#ifndef XNB_HAS_GRID_FIELDS_DEFINTIONS
#error Cannot be included outside of exanb/core/grid_fields.h
#endif

#include <onika/math/basic_types_def.h>
#include <onika/math/quaternion.h>
#include <exaStamp/atom_bond_connectivity.h>
#include <cstdint>


// exaStamp fields are defined in namespace exanb::field
XNB_DECLARE_FIELD(double          ,ep                ,"particle potential energy");
XNB_DECLARE_FIELD(double          ,rxf               ,"filtered particle position X");
XNB_DECLARE_FIELD(double          ,ryf               ,"filtered particle position Y");
XNB_DECLARE_FIELD(double          ,rzf               ,"filtered particle position Z");
XNB_DECLARE_FIELD(double          ,charge            ,"particle electric charge");
XNB_DECLARE_FIELD(::exanb::Mat3d  ,virial            ,"particle pressure virial tensor");

// fields for molecule extension
XNB_DECLARE_FIELD(uint64_t                         ,idmol  ,"molecule Id");
XNB_DECLARE_FIELD(::exaStamp::AtomBondConnectivity ,cmol   ,"molecule bonds");

// for rigid bodies (rigid molecules)
XNB_DECLARE_FIELD(::exanb::Quaternion      ,orient   ,"rigid molecule orientation");
XNB_DECLARE_FIELD(::exanb::Vec3d           ,angmom   ,"rigid molecule angular momentum");
XNB_DECLARE_FIELD(::exanb::Vec3d           ,couple   ,"rigid molecule couple");

// common flat arrays
XNB_DECLARE_FIELD(double          ,flat_rx           ,"coord X (flat)");
XNB_DECLARE_FIELD(double          ,flat_ry           ,"coord Y (flat)");
XNB_DECLARE_FIELD(double          ,flat_rz           ,"coord Z (flat)");
XNB_DECLARE_FIELD(double          ,flat_fx           ,"force X (flat)");
XNB_DECLARE_FIELD(double          ,flat_fy           ,"force Y (flat)");
XNB_DECLARE_FIELD(double          ,flat_fz           ,"force Z (flat)");
XNB_DECLARE_FIELD(double          ,flat_ep           ,"force X (flat)");
XNB_DECLARE_FIELD(uint8_t         ,flat_type         ,"atom type (flat)");
XNB_DECLARE_FIELD(::exanb::Mat3d  ,flat_virial       ,"virial tensor (flat)");


namespace exaStamp
{
  using namespace ::exanb;

  using MultiMatFieldSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_id,field::_type>;

  using MoleculeFieldSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_charge,field::_virial,field::_id,field::_idmol,field::_cmol,field::_type>;
  using RigidMoleculeFieldSet = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_id,field::_orient, field::_angmom, field::_couple,field::_type >;

  using FullFieldMechSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_rxf,field::_ryf,field::_rzf,field::_virial,field::_id,field::_type>;
  using MultimatMechFieldSet = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_virial,field::_id,field::_type>;
  
  static inline constexpr exanb::FieldSets< MultiMatFieldSet
# ifdef EXASTAMP_ENABLE_MOLECULE
    , MoleculeFieldSet, RigidMoleculeFieldSet
# endif
# ifdef EXASTAMP_ENABLE_MECHANICAL
    , FullFieldMechSet, MultimatMechFieldSet
# endif
    > available_field_sets_v = {};
}

#define HAS_POSITION_BACKUP_FIELDS true
#define PositionBackupFieldX ::exanb::field::rxf
#define PositionBackupFieldY ::exanb::field::ryf
#define PositionBackupFieldZ ::exanb::field::rzf

#define XNB_AVAILABLE_FIELD_SETS ::exaStamp::available_field_sets_v

