import numpy as np
import pandas as pd
import mdtraj as md

from sklearn.cluster import SpectralClustering

def selection_distance_matrix(traj, a, b):

    t1 = traj.xyz[:,a,None,:]
    t2 = traj.xyz[:,None,b,:]

    return np.sqrt(np.sum((t1 - t2)**2, axis=3))*10

def water_occupancy(traj, sites, dist_thresh=2.5, identify_waters=True, frame_batch_size=100):
    '''
    For each frame in a trajectory, returns whether there is a water within a distance threshold (default 3 Å) of
    a set of residues, passed as a list or tuple as the sites argument. Optionally, will also determine the identity of
    the closest water to the site and return it in a separate list. (This will slow it down)

    Arguments:
        traj:  an mdtraj trajectory (obviously must include protein and waters)
        sites: protein residues around which to calculate water occupancy. the water must be within DIST_THRESH Å of all of these residues.

        optional:
            dist_thresh: the distance threshold (defaults to 2.5 Å)
            identify_waters: return the identity of each nearest water (defaults to TRUE, turn to FALSE to speed up)
            frame_batch_size:
    '''

    if traj.n_frames > frame_batch_size * 1.5:

        occupancy = []; water_id_list = []

        for i in range(int(np.ceil(traj.n_frames / frame_batch_size))):
            output = water_occupancy(traj[frame_batch_size*i:frame_batch_size*(i+1)], sites, dist_thresh, identify_waters, frame_batch_size)
            if identify_waters:
                occupancy = occupancy + list(output[0])
                water_id_list = water_id_list + list(output[1])

            else:
                occupancy = occupancy + list(output)

        if identify_waters:
            return occupancy, water_id_list

        else:
            return occupancy

    else:

        X = []
        for n in range(len(sites)):
            X.append(selection_distance_matrix(traj,
                                               traj.topology.select("water"),
                                               traj.topology.select("protein and residue "+str(sites[n]))))

        pass_thresh = [np.min(X[i],axis=2)<dist_thresh for i in range(len(sites))]
        occupancy = np.any(np.all(np.array(pass_thresh),axis=0),axis=1)

        if identify_waters:
            residue_arr = np.array([atom.residue for atom in traj.topology.atoms])[traj.topology.select('water')]

            Y = np.array([np.min(resi_water_dist_mat,axis=2) for resi_water_dist_mat in X])
            water_id_list = residue_arr[np.argmin(np.sum(Y**2, axis=0),axis=1)]

            return occupancy, water_id_list

        else:
            return occupancy

def read_water_residue_cfg(filename):

    sets = {}

    with open(filename) as infile:
        for line in infile:
            L = line.strip().split(',')
            name = L[0]
            sets[name] = [int(i) for i in L[1:]]

    return sets

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
        dihedrals = dihedrals * 180 / np.pi
        period = 360
        shift = 0
    else:
        period = np.pi*2
        step = np.pi/12

    if min_jumps:

        dihnew = np.zeros(dihedrals.shape)

        for k in range(dihedrals.shape[1]):
            dihnew[:,k] = minimize_jumps(dihedrals[:,k], period)

        dihedrals = dihnew

    return dihedrals

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

def NSphereTransform(dihs):

    '''
    Converts a set of n angles to Cartesian points on an n-dimensional hypersphere for clustering.
    If two angles are input, you can plot the resulting sphere using spherical_plot().
    If you put in more than two angles, it will create a valid hypersphere, but you cannot
    visualize it.

    Works for up to 6 dimensions, to be useable with any amino acid sidechain. Would need to be extended
    for other applications.

    Takes a numpy array of dihedral angles
    '''

    # Dimensionality of the resulting hypersphere
    dim = dihs.shape[1] + 1

    # First coordinate will always be the same!
    s1 = np.cos(dihs[:,0])

    if dim==2:   # Makes a circle
        s2 = np.sin(dihs[:,0])
        return np.array([s1,s2])

    elif dim==3:  # Converts to spherical coordinates
        s2 = np.sin(dihs[:,0])*np.cos(dihs[:,1])
        s3 = np.sin(dihs[:,0])*np.sin(dihs[:,1])
        return np.array([s1,s2,s3])

    elif dim==4:  # Converts to 4-spherical coordinates
        s2 = np.sin(dihs[:,0])*np.cos(dihs[:,1])
        s3 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.cos(dihs[:,2])
        s4 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])
        return np.array([s1,s2,s3,s4])

    elif dim==5:  # Converts to 5-spherical coordinates
        s2 = np.sin(dihs[:,0])*np.cos(dihs[:,1])
        s3 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.cos(dihs[:,2])
        s4 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])*np.cos(dihs[:,3])
        s5 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])*np.sin(dihs[:,3])
        return np.array([s1,s2,s3,s4,s5])

    elif dim==6:  # Converts to 6-spherical coordinates
        s2 = np.sin(dihs[:,0])*np.cos(dihs[:,1])
        s3 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.cos(dihs[:,2])
        s4 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])*np.cos(dihs[:,3])
        s5 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])*np.sin(dihs[:,3])*np.cos(dihs[:,4])
        s6 = np.sin(dihs[:,0])*np.sin(dihs[:,1])*np.sin(dihs[:,2])*np.sin(dihs[:,3])*np.sin(dihs[:,4])
        return np.array([s1,s2,s3,s4,s5,s6])

    else:
        print('Please enter a set of at least one and no more than 5 angles')
        return None


def AcceleratedSpectralClustering(data, n_clusters, nsamples, gamma=0.01, weighting_constant=1):

    '''
    This function is an accelerated version of Spectral Clustering that is capable of performing similarly
    using much less memory. It is especially designed for datasets where clusterings are of very unequal
    sizes, such that some clusters have enormous numbers of data points while others do not; this kind
    of data structure can be difficult to get reliable clustering on.

    The major problem with Spectral Clustering is that it is approximately O(n^2) in both time and memory
    usage, which makes it intractable for large datasets. This function gets around this problem by subsampling
    the alignment, running spectral clustering on this subsampling, and then assigning the rest of the points to be
    the same as their nearest neighbor. Doing this naïvely would run the risk of picking too few points from clusters
    with relatively few points; the function corrects for this by first weighting each point based on how many
    points are in its local neighborhood. The extent of this weighting can be controlled with the weighting_constant
    parameter (>1 will lead to more neighborhood reweighting and <1 will lead to less). It will therefore sample
    close to every point in very low-density regions and a very small fraction of points in high-density regions.

    Inputs
    ------------------
    data: your dataset as a numpy array of shape N_points x N_features
    nclusters: how many clusters to select
    nsamples: how many samples to draw from the dataset for fitting the embedding model.
    gamma: gamma parameter for spectral clustering, defaults to 0.01.
    weighting constant: scaling factor for setting sequence weights
    '''

    print('Calculating weights per point...')
    x = np.array([1/np.sum(np.exp(-cdist([d], data)**2/2/weighting_constant**2)) for d in data])
    xnorm = x/np.sum(x)

    print('Performing spectral clustering...')
    sc = SpectralClustering(n_clusters=n_clusters, gamma=gamma)
    data_subset = data[np.random.choice(np.arange(len(data)), size=nsamples, replace=False, p=xnorm)]

    sc.fit(data_subset)

    final_labels = []

    print('Assigning nearest neighbors...')
    for d in data:
        final_labels.append(sc.labels_[np.argmin(cdist([d], data_subset))])

    return np.array(final_labels)
