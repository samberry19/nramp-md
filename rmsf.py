'''
rmsf.py
Sam Berry
27 May 2021

Calculate the root mean square fluctuation (RMSF) of each residue in a molecular dynamics
simulation against the first frame (or another reference) using MDtraj.

'''

import numpy as np
import mdtraj as md
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sim", help='Path to simulation file (.dcd)')
parser.add_argument("top", help='Path to topology file (.pdb or .psf)')
parser.add_argument("-s", "--sel", dest="sel", default="protein", help="selection string to pass mdtraj for RMSD calculation")
parser.add_argument("-o", "--out", dest="outpath", default="Auto", help="custom output file path")
parser.add_argument("-f", "--frame", dest="frame", default=0, help='Frame to index to (default 0)')

args = parser.parse_args()

traj = md.load(args.sim, top = args.top)

rmsd = md.rmsf(traj, traj, frame=args.frame, atom_indices = traj.topology.select(args.sel))

if args.outpath == 'Auto':
    outpath = args.sim.split('.')[0] + '_rsmf.txt'
else:
    outpath = args.outpath

np.savetxt(outpath, rmsd)
