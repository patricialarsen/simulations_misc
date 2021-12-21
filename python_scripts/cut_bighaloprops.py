
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ranks = comm.Get_size()


import numpy as np
import pygio as gio

# read in ids and steps from the lightcone for the most massive halos
if rank==0:
    ids_hm = np.loadtxt('/lus/eagle/projects/LastJourney/prlarsen/cluster_lensing/cosmo-cutout/halo_lists/ids_highmass.txt').astype(int)
    steps_hm = np.loadtxt('/lus/eagle/projects/LastJourney/prlarsen/cluster_lensing/cosmo-cutout/halo_lists/steps_highmass.txt').astype(int)
    len_hm = len(ids_hm)
else:
    len_hm = None
len_hm = comm.bcast(len_hm,root=0)
if (rank!=0):
    ids_hm = np.empty(len_hm,dtype=int)
    steps_hm = np.empty(len_hm,dtype=int)

comm.Bcast(ids_hm,root=0)
comm.Bcast(steps_hm,root=0)

stepsu = np.unique(steps_hm)
# cut out step 421 because this is for some reason non existent
for step in stepsu:
    if (rank==0):
        print(step)
    if (step!=421):
        ids_step = ids_hm[steps_hm==step]
        print('computed ids_step',rank)
        data = gio.read_genericio('/eagle/LastJourney/heitmann/OuterRim/M000/L4225/HACC000/analysis/Halos/b0168/Bighalop/STEP'+str(step)+'/m000.'+str(step)+'.bighaloparticles')
        num_elems = len(next(iter(data.values())))
        print('num_elems ',num_elems,rank)
        data_id = data['fof_halo_tag']
        print(len(data_id),rank)
        mask_data = np.isin(data_id,ids_step)
        num_masked = np.sum(mask_data)
        num_masked_total = comm.allreduce(num_masked)
        if num_masked_total>0:
            print('writing')
            # note when writing that the rotations mean axis flipping is required
            gio.write_genericio("bhp_cut_"+str(step)+".gio",
            variables = {"x": data['x'][mask_data], "y": data['y'][ mask_data], "z": data['z'][mask_data], "id": data['id'][mask_data], "tag": data['fof_halo_tag'][mask_data]},
            phys_scale = [1, 1, 1],
            phys_origin = [0, 0, 0],
            method = gio.pygio.PyGenericIO.FileIO.FileIOPOSIX
            )

