#pragma once

#include <exanb/fields.h>
#include <exanb/core/field_set_proto.h>

namespace exanb
{
  // rx, ry and rz are added implicitly

  // different assemblies of fields that will be used
  using MinimalFieldSet      = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz>;
  using MultiMatFieldSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_id,field::_type>;

  using MoleculeFieldSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_charge,field::_virial,field::_id,field::_idmol,field::_cmol,field::_type>;

  using FullFieldMechSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_rxf,field::_ryf,field::_rzf,field::_virial,field::_id,field::_type>;

  using MultimatMechFieldSet = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_virial,field::_id,field::_type>;
  
  // specific field set for grids holding chains' anchors in sliplink model
  using SlipLinkFieldSet     = FieldSet<field::_id, field::_sl>;

  // specific field set for rigid molecules
  using RigidMoleculeFieldSet = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_id,field::_orient, field::_angmom, field::_couple,field::_type >;

  // using UniversalFieldSet     = FieldSet<field::_ep,field::_fx,field::_fy,field::_fz,field::_vx,field::_vy,field::_vz,field::_id,field::_type,field::_idmol,field::_charge,field::_virial,field::_rxf,field::_ryf,field::_rzf>;

  // the standard set of FieldSet
  using StandardFieldSets = FieldSets< MinimalFieldSet, MultiMatFieldSet, MoleculeFieldSet, FullFieldMechSet, MultimatMechFieldSet, RigidMoleculeFieldSet >;
}

#include <exanb/core/field_set_utils.h>

