#!/usr/bin/python3
"""
Writes test fields for the IGAR potential (igar_read_field / igar_force_interp) as legacy VTK STRUCTURED_POINTS files.

With no argument, writes igar_test.vtk, the field read by LJ_igar.msp and CH_rebo_read_xyz.msp: a 128^3 sinusoid,
18 periods along z, values in [0, 1], stored as an 8-bit BINARY image (levels 0-255). Read it with
igar_read_field { storage: uint8, value_scale: 0.00392156862745098 } (1/255).

Usage:
    generate_vtk_grid.py [--shape sinusoid|spheres] [-n N] [--output FILE] [--ascii]

--shape spheres writes 4 random periodic spheres (seed 165): voxels inside a sphere at level 255, the others at 0.
"""

import argparse
import numpy as np

def write_ovito_grid(filename, data_array, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering='F'):
    """
    Writes a 3D NumPy array to a LAMMPS-style grid dump (ITEM: GRID ..., one value per line), readable by OVITO.
    origin and spacing are unused: the box bounds are 0..N along each axis.

    Parameters:
        filename (str): The output file name.
        data_array (np.ndarray): A 3D NumPy array of shape (Nx, Ny, Nz).
        origin (tuple): The origin of the grid (default is (0.0, 0.0, 0.0)).
        spacing (tuple): The spacing between grid points (default is (1.0, 1.0, 1.0)).
    """
    if len(data_array.shape) != 3:
        raise ValueError("Input array must be 3D.")

    Nx, Ny, Nz = data_array.shape

    ofile = open(filename,'w')
    ofile.write('ITEM: TIMESTEP\n')
    ofile.write('0\n')
    ofile.write('ITEM: BOX BOUNDS pp pp pp\n')
    ofile.write('0 %d\n' %(Nx))
    ofile.write('0 %d\n' %(Ny))
    ofile.write('0 %d\n' %(Nz))
    ofile.write('ITEM: DIMENSION\n')
    ofile.write('3\n')
    ofile.write('ITEM: GRID SIZE nx ny nz\n')
    ofile.write('%d %d %d\n' %(Nx,Ny,Nz))
    ofile.write('ITEM: GRID CELLS MASK\n')
    data_array = data_array.ravel(order=ordering)
    for data in data_array:
        ofile.write('%5.8f\n' %(data))
    ofile.close()

def write_vtk_structured_points(filename, data_array, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering='F', binary=True):
    """
    Writes a 3D NumPy array to a legacy VTK structured points file.

    Parameters:
        filename (str): The output VTK file name (e.g., 'points.vtk').
        data_array (np.ndarray): A 3D NumPy array of shape (Nx, Ny, Nz). Its dtype is kept in the file:
            uint8 -> unsigned_char (1 byte per value, read by igar_read_field with storage: uint8), otherwise double.
        origin (tuple): The origin of the grid (default is (0.0, 0.0, 0.0)).
        spacing (tuple): The spacing between grid points (default is (1.0, 1.0, 1.0)).
        ordering (str): 'F' (x fastest, default) or 'C' (z fastest).
        binary (bool): BINARY file (default, much smaller and faster to read) or ASCII.
    """
    if len(data_array.shape) != 3:
        raise ValueError("Input array must be 3D.")

    Nx, Ny, Nz = data_array.shape
    is_uint8 = data_array.dtype == np.uint8
    values = data_array.ravel(order=ordering)

    # header written by hand: the VTK writer turns uint8 point data into COLOR_SCALARS (values / 255 in ASCII files)
    with open(filename, 'wb') as f:
        f.write(("# vtk DataFile Version 3.0\nvtk output\n%s\nDATASET STRUCTURED_POINTS\n"
                 "DIMENSIONS %d %d %d\nSPACING %g %g %g\nORIGIN %g %g %g\nPOINT_DATA %d\n"
                 "SCALARS ScalarData %s\nLOOKUP_TABLE default\n"
                 % ("BINARY" if binary else "ASCII", Nx, Ny, Nz, *spacing, *origin, values.size,
                    "unsigned_char" if is_uint8 else "double")).encode())
        if binary:
            # legacy VTK binary data is big-endian
            f.write(values.astype(">u1" if is_uint8 else ">f8").tobytes())
            f.write(b"\n")
        else:
            np.savetxt(f, values, fmt="%d" if is_uint8 else "%.17g")


def populate_3d_array_with_spheres(Nx, Ny, Nz, num_spheres, mean_radius, radius_variance, seed=None):
    """
    Create a 3D numpy array with dimensions Nx, Ny, Nz, populate it with 0s, and
    fill it with 1s in N spheres of random radii following a Gaussian distribution,
    ensuring periodicity.

    Parameters:
        Nx, Ny, Nz (int): Dimensions of the 3D array.
        num_spheres (int): Number of spheres to populate.
        mean_radius (float): Mean radius of the spheres.
        radius_variance (float): Variance of the sphere radii.
        seed (int, optional): Seed for random number generator.

    Returns:
        numpy.ndarray: The populated 3D array.
    """
    if seed is not None:
        np.random.seed(seed)

    # Initialize the array with zeros
    array = np.zeros((Nx, Ny, Nz), dtype=int)

    # Generate random sphere properties
    for _ in range(num_spheres):
        
        # Generate center in reduced coordinates [0, 1)
        cx_red = np.random.rand()
        cy_red = np.random.rand()
        cz_red = np.random.rand()
        
        # Scale to actual grid coordinates
        cx = int(np.floor(cx_red * Nx))
        cy = int(np.floor(cy_red * Ny))
        cz = int(np.floor(cz_red * Nz))
        
        # Random radius from Gaussian distribution
        radius = max(1, np.random.normal(mean_radius*min(Nx,Ny,Nz), np.sqrt(radius_variance*Nx)))

        # Fill the sphere into the array, using the periodic distance to its center
        dx = np.abs(np.arange(Nx) - cx); dx = np.minimum(dx, Nx - dx)
        dy = np.abs(np.arange(Ny) - cy); dy = np.minimum(dy, Ny - dy)
        dz = np.abs(np.arange(Nz) - cz); dz = np.minimum(dz, Nz - dz)
        distance = np.sqrt(dx[:, None, None]**2 + dy[None, :, None]**2 + dz[None, None, :]**2)
        array[distance <= radius] = 1

    return array


def to_uint8(data_array, val_min=0.0, val_max=1.0):
    """
    Maps [val_min, val_max] onto the 256 levels of an 8-bit image: level = rint(255 * (value - val_min) / (val_max - val_min)).
    With igar_read_field, use value_scale = (val_max - val_min) / 255 to get back values in [0, val_max - val_min].
    """
    levels = np.rint(255.0 * (data_array - val_min) / (val_max - val_min))
    if levels.min() < 0 or levels.max() > 255:
        raise ValueError("values outside [val_min, val_max]")
    return levels.astype(np.uint8)


def populate_3d_array_with_sinusoid(Nx, Ny, Nz, val_min, val_max, n_periods, direction='x'):
    """
    Create a 3D numpy array with a periodic sinusoidal pattern along one direction.

    The value at each grid point is:
        val_min + (val_max - val_min) * 0.5 * (1 + sin(2*pi*n_periods*i/N))
    where i is the index along `direction` and N is the grid size in that direction.

    Parameters:
        Nx, Ny, Nz (int): Dimensions of the 3D array.
        val_min (float): Minimum value of the sinusoid.
        val_max (float): Maximum value of the sinusoid.
        n_periods (float): Number of full periods across the domain in `direction`.
        direction (str or int): Direction of variation: 'x'/0, 'y'/1, or 'z'/2.

    Returns:
        numpy.ndarray: 3D array of shape (Nx, Ny, Nz) with sinusoidal values.
    """
    dir_map = {'x': 0, 'y': 1, 'z': 2}
    if isinstance(direction, str):
        direction = dir_map[direction.lower()]

    N = (Nx, Ny, Nz)[direction]
    coords = np.arange(N)
    sinusoid = val_min + (val_max - val_min) * 0.5 * (1.0 + np.sin(2.0 * np.pi * n_periods * coords / N))

    array = np.empty((Nx, Ny, Nz))
    if direction == 0:
        array[:, :, :] = sinusoid[:, np.newaxis, np.newaxis]
    elif direction == 1:
        array[:, :, :] = sinusoid[np.newaxis, :, np.newaxis]
    else:
        array[:, :, :] = sinusoid[np.newaxis, np.newaxis, :]

    return array


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shape", default="sinusoid", choices=["sinusoid", "spheres"])
    parser.add_argument("-n", type=int, default=128, help="voxels per axis (default 128)")
    parser.add_argument("--output", default="igar_test.vtk")
    parser.add_argument("--ascii", action="store_true", help="ASCII file instead of BINARY")
    args = parser.parse_args()

    Nx, Ny, Nz = args.n, args.n, args.n
    if args.shape == "sinusoid":
        data = populate_3d_array_with_sinusoid(Nx, Ny, Nz, val_min=0., val_max=1.0, n_periods=18, direction='z')
    else:
        data = populate_3d_array_with_spheres(Nx, Ny, Nz, num_spheres=4, mean_radius=0.25, radius_variance=0.1, seed=165)
    write_vtk_structured_points(args.output, to_uint8(data, 0.0, 1.0), origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0),
                                ordering='F', binary=not args.ascii)
    print("%s: %s, %d^3, 8-bit %s" % (args.output, args.shape, args.n, "ASCII" if args.ascii else "BINARY"))
