"""This function fits the activation pattern in y coordinate with erf.

This can be considered as a DV coordinate system, if of course data shows a DV
activation pattern.
"""

import numpy as np
from scipy.signal import medfilt
from skimage.measure import regionprops_table
from sklearn.cluster import KMeans


class BorderFinder:
    """Find the activation border based on raw data"""
    def __init__(self, raw_spts):

        prof     =  np.sum(raw_spts, axis=(0, 1))           # sum over x axis to find the activation pattern
        prof_f   =  medfilt(prof, 151)                      # median filter to smooth the curve
        prof_f  -=  prof_f.min()                            # normalize and rescale
        prof_f  /=  prof_f.max()
        y_pos    =  np.argmin(np.abs(prof_f - 0.25))        # define the border of the activation pattern as the point in which the curve reaches the 25%

        self.y_pos  =  y_pos


class CellsStripesHor:
    """Organize nuclei in strips"""
    def __init__(self, pseudo_cells):

        rgp_psclls  =  regionprops_table(pseudo_cells, properties=["label", "centroid"])      # regionprops of the pseudo cells
        num_rows    =  int(np.round(np.sqrt(len(rgp_psclls["label"]))))                       # define the number of row of the organized cells as the square root of the tot number of cells
        lbls_ctrs   =  np.zeros((2, len(rgp_psclls["label"])))                                # 3D coordinates of the centroids
        for cntr, k in enumerate(rgp_psclls["label"]):
            lbls_ctrs[:, cntr]  =  k, rgp_psclls["centroid-2"][cntr]                          # matrix with labels and centroid y-coordinates

        kmeans         =  KMeans(n_clusters=num_rows, random_state=0).fit(lbls_ctrs[1, :].reshape((-1, 1)))         # y-coordinate clusters
        new_cells      =  np.zeros_like(pseudo_cells)                                                               # matrix with the new tag
        new_cells_bff  =  np.zeros_like(pseudo_cells)                                                               # matrix with the new tag
        old_tags       =  lbls_ctrs[0, :].reshape((-1, 1))
        for cc, ll in enumerate(kmeans.labels_):
            new_cells_bff  +=  (ll + 1) * (pseudo_cells == old_tags[cc])                           # pseudo cell matrix with the tags given by the stripe

        rgp_new_cells           =  regionprops_table(new_cells_bff, properties=["label", "centroid"])                   # regionprops: label and centroid to give organized tags to the stripes
        new_cells_lb_ctr        =  np.zeros((2, len(rgp_new_cells["label"])))                                           # store label and y centroid coordinate (stripes are horizontal)
        new_cells_lb_ctr[0, :]  =  rgp_new_cells["label"]                                                               # add labels
        new_cells_lb_ctr[1, :]  =  rgp_new_cells["centroid-2"]                                                          # add centroid y coordinates in corrispondent positions
        new_cells_lb_ctr        =  new_cells_lb_ctr[:, new_cells_lb_ctr.argsort()[1][::-1]]                             # order the matrix with respect to the y coordinate (now we have a way to map labels)

        for cnt, lb in enumerate(new_cells_lb_ctr[0, :]):                                                               # thanks to 'new_cells_lb_ctr' we can give label in geometrical order -increasing in y- to final strips
            new_cells  +=  (new_cells_bff == lb) * (cnt + 1)

        self.new_cells  =  new_cells


class CellsStripesVer:
    """Organize nuclei in strips"""
    def __init__(self, pseudo_cells):

        rgp_psclls  =  regionprops_table(pseudo_cells, properties=["label", "centroid"])      # regionprops of the pseudo cells
        num_clmns   =  int(np.round(np.sqrt(len(rgp_psclls["label"]))))                       # define the number of row of the organized cells as the square root of the tot number of cells
        lbls_ctrs   =  np.zeros((2, len(rgp_psclls["label"])))                                # 3D coordinates of the centroids
        for cntr, k in enumerate(rgp_psclls["label"]):
            lbls_ctrs[:, cntr]  =  k, rgp_psclls["centroid-1"][cntr]                          # matrix with labels and centroid x-coordinates

        kmeans         =  KMeans(n_clusters=num_clmns, random_state=0).fit(lbls_ctrs[1, :].reshape((-1, 1)))        # x-coordinate clusters
        new_cells      =  np.zeros_like(pseudo_cells)                                                               # matrix with the new tag
        new_cells_bff  =  np.zeros_like(pseudo_cells)                                                               # matrix with the new tag
        old_tags       =  lbls_ctrs[0, :].reshape((-1, 1))
        for cc, ll in enumerate(kmeans.labels_):
            new_cells_bff  +=  (ll + 1) * (pseudo_cells == old_tags[cc])                           # pseudo cell matrix with the tags given by the stripe

        rgp_new_cells           =  regionprops_table(new_cells_bff, properties=["label", "centroid"])                   # regionprops: label and centroid to give organized tags to the stripes
        new_cells_lb_ctr        =  np.zeros((2, len(rgp_new_cells["label"])))                                           # store label and y centroid coordinate (stripes are horizontal)
        new_cells_lb_ctr[0, :]  =  rgp_new_cells["label"]                                                               # add labels
        new_cells_lb_ctr[1, :]  =  rgp_new_cells["centroid-1"]                                                          # add centroid y coordinates in corrispondent positions
        new_cells_lb_ctr        =  new_cells_lb_ctr[:, new_cells_lb_ctr.argsort()[1]]                                   # order the matrix with respect to the y coordinate (now we have a way to map labels)

        for cnt, lb in enumerate(new_cells_lb_ctr[0, :]):                                                               # thanks to 'new_cells_lb_ctr' we can give label in geometrical order -increasing in y- to final strips
            new_cells  +=  (new_cells_bff == lb) * (cnt + 1)

        self.new_cells  =  new_cells


# class CellsStripes:
#     """Organize nuclei in strips"""
#     def __init__(self, pseudo_cells):
#
#         rgp_psclls  =  regionprops_table(pseudo_cells, properties=["label", "centroid"])      # regionprops of the pseudo cells
#         num_rows    =  int(np.round(np.sqrt(len(rgp_psclls["label"]))))                       # define the number of row of the organized cells as the square root of the tot number of cells
#         lbls_ctrs   =  np.zeros((2, len(rgp_psclls["label"])))                                # 3D coordinates of the centroids
#         for cntr, k in enumerate(rgp_psclls["label"]):
#             lbls_ctrs[:, cntr]  =  k, rgp_psclls["centroid-2"][cntr]                          # matrix with labels and centroid y-coordinates
#
#         kmeans     =  KMeans(n_clusters=num_rows, random_state=0).fit(lbls_ctrs[1, :].reshape((-1, 1)))         # y-coordinate clusters
#         new_cells  =  np.zeros_like(pseudo_cells)                                                               # matrix with the new tag
#         old_tags   =  lbls_ctrs[0, :].reshape((-1, 1))
#         for cc, ll in enumerate(kmeans.labels_):
#             new_cells  +=  (ll + 1) * (pseudo_cells == old_tags[cc])                           # pseudo cell matrix with the tags given by the stripe
#
#         self.new_cells  =  new_cells

