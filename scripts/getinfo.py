#!/usr/bin/env python3

# Simone Conti, 2019

import sys
import math
import random
import numpy as np
import mdtraj
from timeit import default_timer as timer

random.seed(1989)

# Print time info or not
_TIMING = False

# Van der Waals radii per element
_ATOMIC_RADII = {'H'   : 0.120, 'He'  : 0.140, 'Li'  : 0.076, 'Be' : 0.059,
                 'B'   : 0.192, 'C'   : 0.170, 'N'   : 0.155, 'O'  : 0.152,
                 'F'   : 0.147, 'Ne'  : 0.154, 'Na'  : 0.102, 'Mg' : 0.086,
                 'Al'  : 0.184, 'Si'  : 0.210, 'P'   : 0.180, 'S'  : 0.180,
                 'Cl'  : 0.181, 'Ar'  : 0.188, 'K'   : 0.138, 'Ca' : 0.114,
                 'Sc'  : 0.211, 'Ti'  : 0.200, 'V'   : 0.200, 'Cr' : 0.200,
                 'Mn'  : 0.200, 'Fe'  : 0.200, 'Co'  : 0.200, 'Ni' : 0.163,
                 'Cu'  : 0.140, 'Zn'  : 0.139, 'Ga'  : 0.187, 'Ge' : 0.211,
                 'As'  : 0.185, 'Se'  : 0.190, 'Br'  : 0.185, 'Kr' : 0.202,
                 'Rb'  : 0.303, 'Sr'  : 0.249, 'Y'   : 0.200, 'Zr' : 0.200,
                 'Nb'  : 0.200, 'Mo'  : 0.200, 'Tc'  : 0.200, 'Ru' : 0.200,
                 'Rh'  : 0.200, 'Pd'  : 0.163, 'Ag'  : 0.172, 'Cd' : 0.158,
                 'In'  : 0.193, 'Sn'  : 0.217, 'Sb'  : 0.206, 'Te' : 0.206,
                 'I'   : 0.198, 'Xe'  : 0.216, 'Cs'  : 0.167, 'Ba' : 0.149,
                 'La'  : 0.200, 'Ce'  : 0.200, 'Pr'  : 0.200, 'Nd' : 0.200,
                 'Pm'  : 0.200, 'Sm'  : 0.200, 'Eu'  : 0.200, 'Gd' : 0.200,
                 'Tb'  : 0.200, 'Dy'  : 0.200, 'Ho'  : 0.200, 'Er' : 0.200,
                 'Tm'  : 0.200, 'Yb'  : 0.200, 'Lu'  : 0.200, 'Hf' : 0.200,
                 'Ta'  : 0.200, 'W'   : 0.200, 'Re'  : 0.200, 'Os' : 0.200,
                 'Ir'  : 0.200, 'Pt'  : 0.175, 'Au'  : 0.166, 'Hg' : 0.155,
                 'Tl'  : 0.196, 'Pb'  : 0.202, 'Bi'  : 0.207, 'Po' : 0.197,
                 'At'  : 0.202, 'Rn'  : 0.220, 'Fr'  : 0.348, 'Ra' : 0.283,
                 'Ac'  : 0.200, 'Th'  : 0.200, 'Pa'  : 0.200, 'U'  : 0.186,
                 'Np'  : 0.200, 'Pu'  : 0.200, 'Am'  : 0.200, 'Cm' : 0.200,
                 'Bk'  : 0.200, 'Cf'  : 0.200, 'Es'  : 0.200, 'Fm' : 0.200,
                 'Md'  : 0.200, 'No'  : 0.200, 'Lr'  : 0.200, 'Rf' : 0.200,
                 'Db'  : 0.200, 'Sg'  : 0.200, 'Bh'  : 0.200, 'Hs' : 0.200,
                 'Mt'  : 0.200, 'Ds'  : 0.200, 'Rg'  : 0.200, 'Cn' : 0.200,
                 'Uut' : 0.200, 'Fl'  : 0.200, 'Uup' : 0.200, 'Lv' : 0.200,
                 'Uus' : 0.200, 'Uuo' : 0.200}

def rotation_matrix(axis, theta):
    """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta degrees.
    """
    theta *= np.pi/180.0
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def rotate_coor(coor, axis, theta):
    """
        Rotate the coordinates by theta degrees around axis.
        Return the new coordinates.
    """
    return np.matmul(coor, rotation_matrix(axis, theta)) 

def bounding_box(traj, frame=0, probe_radius=0.0):
    """
        Compute the surface area and the volume of the cuboid that completely
        enclose the molecule (considering van der Waals radii).
    """
    coor = traj.xyz[frame]
    atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in traj.topology.atoms]
    radii = np.array(atom_radii, np.float32) + probe_radius
    xmin = np.amin([x-r for x,r in zip(coor[:,0], radii)])
    xmax = np.amax([x+r for x,r in zip(coor[:,0], radii)])
    ymin = np.amin([y-r for y,r in zip(coor[:,1], radii)])
    ymax = np.amax([y+r for y,r in zip(coor[:,1], radii)])
    zmin = np.amin([z-r for z,r in zip(coor[:,2], radii)])
    zmax = np.amax([z+r for z,r in zip(coor[:,2], radii)])
    lx = xmax-xmin
    ly = ymax-ymin
    lz = zmax-zmin
    area = 2*(lx*ly + ly*lz + lz*lx)
    volume = lx*ly*lz
    return area, volume

def minimum_bounding_box(traj, frame=0, npoints=36):
    """
        Rotate the molecule in 3D on a grid of npoints (per axis) and find the 
        orientation that minimize the volume of the cuboid containing the 
        molecule (considering van der Waals radii).
        Return the oriented molecule, the area and the volume.
    """
    time_start = timer()
    drot = 360.0/npoints
    min_area = 1E200
    min_volume = 1E200
    xrot = yrot = zrot = 0.0
    for x in range(npoints):
        traj.xyz[frame] = rotate_coor(traj.xyz[frame], [1,0,0], drot)
        for y in range(npoints):
            traj.xyz[frame] = rotate_coor(traj.xyz[frame], [0,1,0], drot)
            for z in range(npoints):
                traj.xyz[frame] = rotate_coor(traj.xyz[frame], [0,0,1], drot)
                area, volume = bounding_box(traj)
                if volume<min_volume:
                    min_area = area
                    min_volume = volume
                    xrot, yrot, zrot = (x+1)*drot,(y+1)*drot,(z+1)*drot
                #print(x*drot, y*drot, z*drot, area*100.0, min_area*100.0)
                #traj.save('rotate_%g_%g_%g.pdb' % (x*drot,y*drot,z*drot))
    traj.xyz[frame] = rotate_coor(traj.xyz[frame], [1,0,0], xrot)
    traj.xyz[frame] = rotate_coor(traj.xyz[frame], [0,1,0], yrot)
    traj.xyz[frame] = rotate_coor(traj.xyz[frame], [0,0,1], zrot)
    area, volume = bounding_box(traj)
    time_end = timer()
    if _TIMING: print('Time in minimum_bounding_box:', time_end-time_start);
    return traj, area, volume


def volume_monte_carlo(traj, frame=0, probe_radius=0.14, ppa=50, nrepeat=10):
    """
        Compute the van der Waals volume using hit and miss Monte Carlo method.
        First, align the molecule, so that its bounding box has the minimum 
        volume. Then "shoot" volume_bbox*ppa points inside the bounding box
        and count how many of them are inside an atom. Get the volume as 
        ratio between hit and total points. To increase the accuracy, increase
        the number of points per angstrom cube (ppa).
    """
    time_start = timer()
    traj, surface_bbox, volume_bbox = minimum_bounding_box(traj, frame)
    coor = traj.xyz[frame]
    npoints = int(volume_bbox*ppa*1000.0)
    atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in traj.topology.atoms]
    radii = np.array(atom_radii, np.float32) + probe_radius
    xmin = np.amin([x-r for x,r in zip(coor[:,0], radii)])
    xmax = np.amax([x+r for x,r in zip(coor[:,0], radii)])
    ymin = np.amin([y-r for y,r in zip(coor[:,1], radii)])
    ymax = np.amax([y+r for y,r in zip(coor[:,1], radii)])
    zmin = np.amin([z-r for z,r in zip(coor[:,2], radii)])
    zmax = np.amax([z+r for z,r in zip(coor[:,2], radii)])
    volumes = list()
    for i in range(nrepeat):
        hit = 0
        for i in range(npoints):
            p = np.array([random.uniform(xmin,xmax), random.uniform(ymin, ymax), random.uniform(zmin, zmax)])
            for atom, radius in zip(coor, radii):
                diff = p-atom
                dist = math.sqrt(np.dot(diff, diff))
                if dist<radius:
                    hit += 1
                    break
        volume_vdw = (hit/npoints) * volume_bbox
        volumes.append(volume_vdw)
    volume_vdw = np.average(volumes)
    volume_vdw_error = np.std(volumes)/math.sqrt(nrepeat)
    time_end = timer()
    if _TIMING: print('Time in volume_monte_carlo:', time_end-time_start);
    return traj, surface_bbox, volume_bbox, volume_vdw, volume_vdw_error


def gyration_radius(traj, frame=0):
    """
        Compute the radius of gyration.
    """
    coor = traj.xyz[frame]
    avg = np.average(coor, axis=0)
    rgyr = 0.0
    for atom in coor:
        diff = atom-avg
        rgyr += np.dot(diff, diff)
    rgyr /= len(coor)
    rgyr = math.sqrt(rgyr)
    return rgyr

def get_all_info(molfile):
    frame = 0
    traj = mdtraj.load(molfile)
    traj.center_coordinates(mass_weighted=False)
    mass = 0.0
    for atom in traj.topology.atoms:
        mass += atom.element.mass
    traj, surface_bbox, volume_bbox, volume_vdw, volume_vdw_error = volume_monte_carlo(traj, frame=frame, probe_radius=0.0, ppa=1000, nrepeat=10)
    surface_bbox *= 100.0
    volume_vdw *= 1000.0
    rgyr = 10.0*gyration_radius(traj, frame=frame)
    inertia = mdtraj.compute_inertia_tensor(traj)[frame]
    eigvals, eigvect = np.linalg.eig(inertia)
    eigvals *= 100.0
    return mass, surface_bbox, volume_vdw, rgyr, eigvals

def print_info(molfile):
    mass, surface_bbox, volume_vdw, rgyr, inertia = get_all_info(molfile)
    print('mass = %.2f' % (mass))
    print('bbox = %.2f' % (surface_bbox))
    print('vvdw = %.2f' % (volume_vdw))
    print('rgyr = %.2f' % (rgyr))
    print('rotations = 3')
    print('%.3f' % (inertia[0]))
    print('%.3f' % (inertia[1]))
    print('%.3f' % (inertia[2]))


# Solvent properties
_SOLVENTS = [
    # Name           Density Acentric Permittivity Expansion
    ['ethylene',        0.57,  0.089,  1.0,  2.40],
    ['iodine',          3.96,  0.229,  4.0,  0.00],
    ['cyclohexane',     0.778, 0.212,  2.0,  1.21],
    ['benzene',         0.88,  0.212,  2.3,  1.25],
    ['toluene',         0.867, 0.263,  2.4,  1.08],
    ['m-xylene',        0.86,  0.325,  2.4,  0.99],
    ['o-xylene',        0.88,  0.31,   2.6,  0.00],
    ['pentane',         0.626, 0.251,  1.4,  1.58],
    ['isopentane',      0.616, 0.227,  1.8,  0.00],
    ['hexane',          0.66,  0.299,  1.9,  1.41],
    ['octane',          0.703, 0.398,  2.0,  1.14],
    ['chloroform',      1.49,  0.218,  4.8,  1.27],
    ['dioxane',         0.796, 0.307,  2.3,  1.12],
    ['acetaldehyde',    0.788, 0.303, 21.1,  1.69],
    ['acetone',         0.784, 0.304, 20.7,  1.43],
    ['ethyl_acetate',   0.81,  0.329,  6.0,  1.38],
    ['acetic_acid',     1.05,  0.447,  6.2,  1.10],
    ['acetonitrile',    0.786, 0.278, 37.5,  1.36],
    ['dimethyl_ether',  0.74,  0.2,    5.3,  0.00],
    ['diethyl_ether',   0.713, 0.281,  4.3,  1.60],
    ['helium',          0.13, -0.365,  1.1, -1.49],
    ['neon',            1.21, -0.029,  1.5, 15.40],
    ['argon',           1.40,  0.001,  1.5,  4.80],
    ['krypton',         2.41,  0.005,  1.7,  0.00],
    ['xenon',           2.94,  0.008,  1.9,  0.00],
    ['water',           1.00,  0.344, 78.5,  0.21],
    ['methanol',        0.796, 0.556, 32.6,  1.09],
    ['ethanol',         0.796, 0.644, 24.6,  1.09],
    ['propanol',        0.803, 0.623, 20.1,  0.79],
    ['isopropanol',     0.786, 0.665, 17.9,  0.00],
    ['butanol',         0.81,  0.593, 17.8,  0.75],
    ['isobutanol',      0.802, 0.592, 17.3,  0.94]
]

def make_solvent_file():
    """
        Calculate molecular properties for all solvents and print therm 
        formatted for inclusion in C file.
    """
    for name, density, acentricity, permittivity, expansion in _SOLVENTS:
        mass, bbox, vvdw, rgyr, inertia = get_all_info('../molfiles/%s.pdb' % (name))
        print('{ .name=%-15s, .density=%-8.3f, .acentricity=%-8.3f, .permittivity=%-8.3f, .expansion=%-8.3f, .mass=%-8.3f, .bbox=%-8.3f, .vvdw=%-8.3f, .rgyr=%-8.3f },' % 
            (name, density, acentricity, permittivity, expansion, mass, bbox, vvdw, rgyr))


if __name__=='__main__':

    if len(sys.argv)!=2:
        print('''
This script compute molecular properties of small molecules, in particular:
mass, van der Waals volume, minimum bounding box area, radius of gyration, and
inertia moments. This script is most useful for getting the necessary
properties for solvation entropies as computed by A.J. Garza.

Usage: %s molfile.mol2

molfile.mol2 is an optimized 3D model of the molecule of interest. You can
create it with Avogadro or similar software. Other file formats are supported,
like PDB.
''' % (sys.argv[0]))
        quit()
    print_info(sys.argv[1])

