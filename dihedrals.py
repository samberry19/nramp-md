'''
dihedrals.py
Sam Berry
27 May 2021

Calculate the dihedrals angles of a residue across an MD simulation.
'''

import numpy as np
import pandas as pd
import mdtraj as md
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sim", help='Path to simulation file (.dcd)')
parser.add_argument("top", help='Path to topology file (.pdb or .psf)')
parser.add_argument("res", help="Residue name w/ 3-letter code followed by number, e.g. 'ASP56'")
parser.add_argument("-u", "--u", dest="unit", default="deg", help="unit to report angles in. Options are 'deg' or 'rad' (default 'deg')", type=str)
parser.add_argument("-o", "--out", dest="outpath", default="Auto", help="custom output file path", type=str)
parser.add_argument("-m", "--min_jumps", dest="min_jumps", default=True, help='Whether to minimize jumps when reporting angles for nicer plotting', type=bool)
parser.add_argument("-d", "--dihtable", dest="dihtable", default="dihtable.csv", help='Custom dihedral angle table path', type=str)

args = parser.parse_args()

'''
Let's define some functions first!
'''

def change_boundary_condition(angles, shift, period=360):

    '''This is a helper function to go with minimize_jumps().
    It takes in some angles and simply adds or subtracts one period
    from some of them.'''

    shifted_angles = []
    for angle in angles:
        if angle + shift > period:
            shifted_angles.append(angle - period)
        else:
            shifted_angles.append(angle)

    shifted_angles  = np.array(shifted_angles)

    if np.all(shifted_angles < 0):
        shifted_angles = shifted_angles + period

    return shifted_angles


def minimize_jumps(angles, period, step=1):

    '''
    A little function that picks a quadrant to display a graph of angles in based on which one minimizes the summed difference
    between angles over time points.

    Basically, this is just my attempt to deal with the fact that we need to plot dihedrals continuously, but of course angles
    cycle back around after 360 degrees, and I don't want the graph to look like it has a massive jump in angle just because it went from 350 to 0
    (an actual difference of only 10 degrees.)
    '''

    a = []

    # Calculate the sum of the 'jumps' between points across the trajectory
    for i in np.arange(0, period, step):
        a.append(np.sum(np.abs(change_boundary_condition(angles,i)[1:] - change_boundary_condition(angles,i)[:-1])))

    opt_shift = np.argmin(a)

    return change_boundary_condition(angles, opt_shift)

def AtomsToIndices(traj):
    topology=traj.topology
    atom_to_index = {}
    for n,atom in enumerate(list(topology.atoms)):
        atom_to_index[str(atom)] = n
    return atom_to_index


def CalcDihedrals(traj, resi, dihtable, min_jumps=True, unit='deg'):

    '''
    Calculates the dihedrals for a given residue in a given trajectory.

    Inputs:
        traj: an mdtraj trajectory
        resi:
    '''

    atom_to_index = AtomsToIndices(traj)
    indices = []

    for atoms in dihtable.loc[resi[:3]].dropna():
        atom_list = atoms.split('-')
        indices.append([atom_to_index[resi+'-'+atom] for atom in atom_list])

    dihedrals = md.compute_dihedrals(traj, np.array(indices))

    if unit in ('deg', 'degree', 'degrees', 'd'):
        dihedrals = dihedrals*180/np.pi
        period = 360
        step = 1
    else:
        period = np.pi * 2
        step = np.pi/12

    if min_jumps:

        dihnew = np.zeros(dihedrals.shape)

        for k in range(dihedrals.shape[1]):
            dihnew[:,k] = minimize_jumps(dihedrals[:,k], period, step=step)

        dihedrals = dihnew

    return dihedrals


traj = md.load(args.sim, top = args.top)

dihtable=pd.read_csv(args.dihtable, index_col=0)

dihs = CalcDihedrals(traj, args.res, dihtable, unit=args.unit, min_jumps=args.min_jumps)

if args.outpath == 'Auto':
    outpath = args.sim.split('.')[0] + '_'+args.res+'_dihedrals.txt'
else:
    outpath = args.outpath

np.savetxt(outpath, dihs)
