#!/usr/bin/python3

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def write_ovito_grid(filename, data_array, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering='F'):
    """
    Writes a 3D NumPy array to a VTK structured points file.

    Parameters:
        filename (str): The output VTK file name (e.g., 'points.vtk').
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

def write_vtk_structured_points(filename, data_array, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering='F'):
    """
    Writes a 3D NumPy array to a VTK structured points file.

    Parameters:
        filename (str): The output VTK file name (e.g., 'points.vtk').
        data_array (np.ndarray): A 3D NumPy array of shape (Nx, Ny, Nz).
        origin (tuple): The origin of the grid (default is (0.0, 0.0, 0.0)).
        spacing (tuple): The spacing between grid points (default is (1.0, 1.0, 1.0)).
    """
    if len(data_array.shape) != 3:
        raise ValueError("Input array must be 3D.")

    Nx, Ny, Nz = data_array.shape

    # Create the structured points
    structured_points = vtk.vtkImageData()
    structured_points.SetDimensions(Nx, Ny, Nz)
    structured_points.SetOrigin(origin)
    structured_points.SetSpacing(spacing)

    # Convert the NumPy array to a VTK array
    vtk_data_array = numpy_to_vtk(data_array.ravel(order=ordering), deep=True)
    vtk_data_array.SetName("ScalarData")

    # Add the data to the structured points
    structured_points.GetPointData().SetScalars(vtk_data_array)

    # Write the structured points to a file
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(filename)
    writer.SetInputData(structured_points)
    writer.Write()

import numpy as np

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
        
        print(cx,cy,cz)
        
        # Random radius from Gaussian distribution
        radius = max(1, np.random.normal(mean_radius*min(Nx,Ny,Nz), np.sqrt(radius_variance*Nx)))

        # Fill the sphere into the array
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    # Compute periodic distance to the sphere's center
                    dx = min(abs(x - cx), Nx - abs(x - cx))
                    dy = min(abs(y - cy), Ny - abs(y - cy))
                    dz = min(abs(z - cz), Nz - abs(z - cz))
                    distance = np.sqrt(dx**2 + dy**2 + dz**2)

                    if distance <= radius:
                        array[x, y, z] = 1

    return array


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


# Example usage
if __name__ == "__main__":
    Nx, Ny, Nz = 128,128,128
    ordering = 'F'

    # Sinusoidal grid: 3 periods along x, values in [0.5, 2.0]
    data_sin = populate_3d_array_with_sinusoid(Nx, Ny, Nz, val_min=0., val_max=1.0, n_periods=18, direction='z')
    write_vtk_structured_points(
#        "sinusoid_%dx%dx%d_%sorder.vtk" % (Nx, Ny, Nz, ordering),
        "igar_test.vtk",
        data_sin, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering=ordering)

    # # Sphere grid (original)
    # num_spheres = 4
    # mean_radius = 0.25
    # radius_variance = 0.1
    # data_sph = populate_3d_array_with_spheres(Nx, Ny, Nz, num_spheres, mean_radius, radius_variance, seed=165)
    # write_vtk_structured_points(
    #     "points_%dx%dx%d_%sorder.vtk" % (Nx, Ny, Nz, ordering),
    #     data_sph, origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), ordering=ordering)
