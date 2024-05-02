"""This function prepare tags and their positions
to add text (tags) on a ImageView and browse the stack.
"""


import numpy as np
from skimage.morphology import label
from skimage.measure import regionprops


class TagsAndPositions:
    """Only class, does all the job"""
    def __init__(self, mtx_tot_ok):

        steps     =  mtx_tot_ok.shape[0]                                                    # number of frames
        idxs_tot  =  mtx_tot_ok[mtx_tot_ok != 0]                                            # list of all the tags
        idxs_tot  =  np.unique(idxs_tot)
        poss      =  np.zeros((idxs_tot.size, 2, steps))                                    # initialize positions
        tags      =  np.zeros((idxs_tot.size))                                              # intialieze tags

        for k in range(idxs_tot.size):
            tags[k]  =  idxs_tot[k]                                                         # add tag
            bffr     =  (mtx_tot_ok == idxs_tot[k])                                         # isolate single tag object
            for t in range(steps):                                                          # in each time step
                bffr_t         =  label(bffr[t, :, :], connectivity=1)                      # label
                rgp_t          =  regionprops(bffr_t.astype(int))                           # measure regioinproperties
                poss[k, :, t]  =  rgp_t[0]['centroid'][0], rgp_t[0]['centroid'][1]          # add centroids

        self.tags  =  tags
        self.poss  =  poss

