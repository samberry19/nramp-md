import numpy as np
import mdtraj as md
import argparse
from mdtools import *

parser = argparse.ArgumentParser()
parser.add_argument("traj", help='Raw .dcd trajectory file')
parser.add_argument("top", help="Topology file")
parser.add_argument("--cfg", default='water_residues_cfg.csv')
parser.add_argument("--thresh", default=3.2, type=float)
parser.add_argument("--prefix", help="prefix for output files", default="./")
parser.add_argument("--batch", help="size of frame batches", type=int, default=100)
args = parser.parse_args()

traj = md.load(args.traj, top=args.top)
print('Loaded simulation')

resi_dict = read_water_residue_cfg(args.cfg)
print('Calculating occupancy at', len(resi_dict), 'sites with threshold of', str(args.thresh), 'Å')

if args.prefix[-1] != '/':
    args.prefix = args.prefix + '_'

for ns,site_residues in resi_dict.items():
    print('Site', ns, site_residues)
    wocc,nwats = water_occupancy(traj, site_residues, dist_thresh=args.thresh, frame_batch_size=int(args.batch))
    np.savetxt(args.prefix + '_'.join([str(res) for res in site_residues]) + '_closest_waters.txt', nwats, fmt='%s')
    np.savetxt(args.prefix + '_'.join([str(res) for res in site_residues]) + '_water_occupancy_t'+str(args.thresh)+'.txt', wocc)
