/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

/*
 * _exastamp_data — registers exaStamp-specific slot extractors into pyonika's
 * slot_as_array registry (imported by pyexastamp/__init__.py):
 *
 *   ThermodynamicState     →  dict (internal units + values as printed in the log)
 *   ParticleSpecies        →  list of dict
 *   MoleculeSpeciesVector  →  dict { "molecules": [...], "bridge_molecules": [...] }
 *
 * Every returned object is a copy and stays valid after pyonika.end().
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <onika/scg/operator_slot.h>
#include <onika/python/slot_extractor_api.h>
#include <onika/physics/constants.h>
#include <exaStamp/thermo_state/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>

#include <functional>
#include <typeindex>

namespace py = pybind11;
using OSB = onika::scg::OperatorSlotBase;

using RegisterFn = onika::python::RegisterSlotExtractorFn;

// versioned registration function published by pyonika (ImportError if incompatible)
static RegisterFn get_register_fn()
{
  return onika::python::slot_extractor_registry();
}

template<class T> static const T* slot_value(OSB& slot)
{
  auto* typed = static_cast<onika::scg::OperatorSlot<T>*>(&slot);
  if (!typed->has_value()) return nullptr;
  return &(**typed);
}

static py::array_t<double> vec3(const onika::math::Vec3d& v)
{
  py::array_t<double> a(3);
  auto r = a.mutable_unchecked<1>();
  r(0) = v.x; r(1) = v.y; r(2) = v.z;
  return a;
}

static py::array_t<double> mat3(const onika::math::Mat3d& m)
{
  py::array_t<double> a({3, 3});
  auto r = a.mutable_unchecked<2>();
  r(0,0) = m.m11; r(0,1) = m.m12; r(0,2) = m.m13;
  r(1,0) = m.m21; r(1,1) = m.m22; r(1,2) = m.m23;
  r(2,0) = m.m31; r(2,1) = m.m32; r(2,2) = m.m33;
  return a;
}

// ---------------------------------------------------------------------------
// ThermodynamicState. "internal" holds raw internal-unit totals; the top level
// holds the values print_thermodynamic_state logs (same conversion factors):
// per-atom energies in eV, temperature in K, pressure and stress in Pa.
// ---------------------------------------------------------------------------
static py::object thermodynamic_state_to_python(const exaStamp::ThermodynamicState& s)
{
  using namespace onika::physics;
  const double conv_temperature = 1.e4 * atomicMass / boltzmann;        // internal units to Kelvin
  const double conv_energy      = 1.e4 * atomicMass / elementaryCharge;  // internal units to eV
  const double conv_pressure    = 1.e4 * atomicMass * 1e30;              // internal units to Pascal
  const double conv_density     = atomicMass * 1e3 * 1e24;               // internal units to g/cm^3

  const double n = static_cast<double>(s.particle_count());
  const double per_atom = n > 0 ? conv_energy / n : 0.;
  const double dof = s.temperature_dof_count();

  py::dict internal;
  internal["kinetic_energy"]   = vec3(s.kinetic_energy());
  internal["potential_energy"] = s.potential_energy();
  internal["internal_energy"]  = s.internal_energy();
  internal["chemical_energy"]  = s.chemical_energy();
  internal["total_energy"]     = s.total_energy();
  internal["temperature"]      = vec3(s.temperature());
  internal["virial"]           = mat3(s.virial());
  internal["ke_tensor"]        = mat3(s.ke_tensor());
  internal["kinetic_momentum"] = vec3(s.kinetic_momentum());
  internal["mass"]             = s.mass();
  internal["volume"]           = s.volume();

  py::dict d;
  d["particle_count"]      = static_cast<unsigned long long>(s.particle_count());
  d["total_energy"]        = s.total_energy() * per_atom;
  d["kinetic_energy"]      = s.kinetic_energy_scal() * per_atom;
  d["potential_energy"]    = s.potential_energy() * per_atom;
  d["temperature"]         = s.temperature_scal() / dof * conv_temperature;
  d["temperature_xyz"]     = vec3((conv_temperature / dof) * s.temperature());
  d["pressure"]            = s.pressure_scal() * conv_pressure;
  d["stress"]              = mat3(conv_pressure * s.full_stress_tensor());
  d["volume"]              = s.volume();
  d["density"]             = s.volume() > 0 ? s.mass() / s.volume() * conv_density : 0.;
  d["units"]               = py::dict(py::arg("total_energy") = "eV/atom", py::arg("kinetic_energy") = "eV/atom",
                                      py::arg("potential_energy") = "eV/atom", py::arg("temperature") = "K",
                                      py::arg("pressure") = "Pa", py::arg("stress") = "Pa",
                                      py::arg("volume") = "ang^3", py::arg("density") = "g/cm^3");
  d["internal"]            = internal;
  return d;
}

// ---------------------------------------------------------------------------
// ParticleSpecies (masses in Dalton, charges in elementary charges)
// ---------------------------------------------------------------------------
static py::object species_to_python(const exaStamp::ParticleSpecies& species)
{
  py::list l;
  for (size_t t = 0; t < species.size(); ++t)
  {
    const auto& sp = species[t];
    py::dict d;
    d["type"]          = t;
    d["name"]          = sp.name();
    d["mass"]          = sp.m_mass;
    d["charge"]        = sp.m_charge;
    d["z"]             = sp.m_z;
    d["molecule_name"] = sp.molecule_name();
    d["ff_name"]       = sp.FFname();
    if (sp.m_rigid_atom_count > 0)
    {
      py::list atoms;
      for (unsigned int i = 0; i < sp.m_rigid_atom_count; ++i)
      {
        atoms.append(py::dict(py::arg("name") = sp.rigid_atom_name(i),
                              py::arg("type") = sp.m_rigid_atoms[i].m_atom_type,
                              py::arg("position") = vec3(sp.m_rigid_atoms[i].m_pos)));
      }
      d["rigid_atoms"] = atoms;
      d["inertia"] = vec3(sp.m_minert);
    }
    l.append(d);
  }
  return l;
}

// ---------------------------------------------------------------------------
// MoleculeSpeciesVector: connectivity as small integer arrays
// ---------------------------------------------------------------------------
template<class T, size_t N, size_t M>
static py::array_t<int> index_table(const T (&tab)[N][M], unsigned int count)
{
  py::array_t<int> a({static_cast<py::ssize_t>(count), static_cast<py::ssize_t>(M)});
  auto r = a.mutable_unchecked<2>();
  for (unsigned int i = 0; i < count; ++i)
    for (size_t j = 0; j < M; ++j) r(i, j) = static_cast<int>(tab[i][j]);
  return a;
}

static py::list molecules_to_python(const onika::memory::CudaMMVector<exaStamp::MoleculeSpecies>& mols)
{
  py::list l;
  for (const auto& m : mols)
  {
    py::dict d;
    d["name"] = m.name();
    py::array_t<int> types(m.m_nb_atoms);
    for (unsigned int i = 0; i < m.m_nb_atoms; ++i) types.mutable_at(i) = m.m_atom_type[i];
    d["atom_types"] = types;
    d["bonds"]      = index_table(m.m_bonds, m.m_nb_bonds);
    d["angles"]     = index_table(m.m_bends, m.m_nb_bends);
    d["torsions"]   = index_table(m.m_torsions, m.m_nb_torsions);
    d["impropers"]  = index_table(m.m_impropers, m.m_nb_impropers);
    l.append(d);
  }
  return l;
}

template<class T, class Fn> static void register_type(RegisterFn reg, Fn to_python)
{
  reg(std::type_index(typeid(T)), [to_python](OSB& slot) -> py::object {
    const T* v = slot_value<T>(slot);
    return v ? to_python(*v) : py::none();
  });
}

PYBIND11_MODULE(_exastamp_data, m)
{
  m.doc() = "exaStamp-specific slot extractors for pyonika's slot_as_array registry";

  RegisterFn reg = get_register_fn();
  register_type<exaStamp::ThermodynamicState>(reg, thermodynamic_state_to_python);
  register_type<exaStamp::ParticleSpecies>(reg, species_to_python);
  register_type<exaStamp::MoleculeSpeciesVector>(reg, [](const exaStamp::MoleculeSpeciesVector& v) -> py::object {
    py::dict d;
    d["molecules"] = molecules_to_python(v.m_molecules);
    d["bridge_molecules"] = molecules_to_python(v.m_bridge_molecules);
    return d;
  });
}
