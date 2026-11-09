"""
Structured description of a simulation, turned into an input file by
``render.render_input``.

The agent fills these fields instead of writing YAML: everything that must be
consistent (domain arithmetic, lattice types, grid flavor, block order) is then
derived by the renderer. Quantities are strings with units, e.g. "3.3 ang",
"300 K", "1.0e-3 ps", "0.0104 eV" (valid units: see onika_mcp.units).
"""

from typing import Annotated, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field

__all__ = ["SimulationSpec", "Species", "LatticeSystem", "XYZFileSystem", "EAMAlloyPotential",
           "LJPotential", "Ensemble", "Outputs"]

Quantity = Annotated[str, Field(description="number followed by a unit, e.g. '3.3 ang'")]


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Species(_Model):
    name: str = Field(description="species name; use the chemical symbol (e.g. 'Ta'), it must match the element names of an EAM file")
    mass: Optional[Quantity] = Field(None, description="e.g. '180.95 Da'; default: from the periodic table when name is a chemical symbol")
    z: Optional[int] = Field(None, description="atomic number; default: from the periodic table")
    charge: Quantity = "0 e-"


class LatticeSystem(_Model):
    kind: Literal["lattice"] = "lattice"
    structure: Literal["SC", "BCC", "FCC", "HCP", "2BCT", "c-DIA", "h-DIA", "graphite"]
    lattice_parameter: Union[Quantity, list[Quantity]] = Field(
        description="unit cell edge: one value for cubic cells, or [a, b, c]")
    repeats: list[int] = Field(min_length=3, max_length=3, description="number of unit cells along x, y, z")
    types: Optional[list[str]] = Field(
        None, description="species of each atom of the unit cell (SC 1, BCC 2, FCC/HCP/2BCT 4, diamond/graphite 8); default: the first species everywhere")
    periodic: list[bool] = Field([True, True, True], min_length=3, max_length=3)


class XYZFileSystem(_Model):
    kind: Literal["xyz_file"] = "xyz_file"
    filename: str = Field(description="extended XYZ file (Lattice= and Properties= on line 2), relative to the work directory")
    periodic: list[bool] = Field([True, True, True], min_length=3, max_length=3)


class EAMAlloyPotential(_Model):
    kind: Literal["eam_alloy"] = "eam_alloy"
    file: str = Field(description="LAMMPS eam/alloy (setfl) file name: use list_potential_files, never guess")
    rcut: Optional[Quantity] = Field(None, description="default: the cutoff stored in the file")


class LJPotential(_Model):
    kind: Literal["lj"] = "lj"
    epsilon: Quantity = Field(description="well depth, e.g. '0.0104 eV'")
    sigma: Quantity = Field(description="e.g. '3.4 ang'")
    rcut: Quantity = Field(description="cutoff, typically 2.5 sigma, e.g. '8.5 ang'")


class Ensemble(_Model):
    kind: Literal["nve", "langevin", "berendsen"] = "nve"
    temperature: Optional[Quantity] = Field(None, description="target temperature (langevin, berendsen)")
    temperature_start: Optional[Quantity] = Field(None, description="linear ramp start, instead of temperature")
    temperature_stop: Optional[Quantity] = Field(None, description="linear ramp end, instead of temperature")
    gamma: Quantity = Field("0.1 ps^-1", description="langevin friction")
    tau: Quantity = Field("0.1 ps", description="berendsen relaxation time")


class Outputs(_Model):
    thermo_every: int = Field(10, description="print and log the thermodynamic state every N steps")
    snapshot: Literal["none", "xyz", "paraview", "lammps_data"] = "none"
    snapshot_every: int = 100
    restart_every: Optional[int] = Field(None, description="write a restart (atoms dump) every N steps")


class SimulationSpec(_Model):
    title: Optional[str] = Field(None, description="one-line description, written as the first comment")
    species: list[Species] = Field(min_length=1)
    system: Annotated[Union[LatticeSystem, XYZFileSystem], Field(discriminator="kind")]
    potential: Annotated[Union[EAMAlloyPotential, LJPotential], Field(discriminator="kind")]
    ensemble: Ensemble = Field(default_factory=Ensemble)
    initial_temperature: Optional[Quantity] = Field(
        None, description="Maxwell-Boltzmann velocities at this temperature (lattice systems); default: ensemble temperature or start temperature, none for NVE")
    timestep: Quantity = "1.0e-3 ps"
    steps: int = Field(1000, description="number of time steps")
    outputs: Outputs = Field(default_factory=Outputs)
    rcut_inc: Quantity = Field("1.0 ang", description="neighbor list skin")
    grid_flavor: Optional[Literal["grid_flavor_full", "grid_flavor_multimat", "grid_flavor_minimal"]] = Field(
        None, description="default: chosen from the species and the potential")
    extra: Optional[str] = Field(
        None, description="raw YAML appended to the input, for features the description does not cover; checked like the rest")
