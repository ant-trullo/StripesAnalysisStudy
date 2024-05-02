"""This function loads raw data for smiFish analysis.

Raw data are .czi, .tif, .lsm, .lif. The input is the file name and the channel numbers,
the output are the matrices for the multi z nuclei and the multi z
for spots channels, plus the mip of the nuclei channel. If the nuclei channel is
not asked, the nucs_mip is going to be a zeros array.
chs_spts_nucs is a list that has 3 channel numbers for spots_a, spots_b and nuclei
respectively. If one of this number is -1, the relative channel is absent.
"""


import numpy as np
import czifile
import tifffile
from tifffile import TiffFile
import read_lif
from PyQt5 import QtWidgets  # , Qt
from PyQt5.QtCore import Qt


class RawDataLoader:
    """This function loads image files in several formats and organize it in a np array."""
    def __init__(self, raw_data_fname, chs_spts_nucs):

        if raw_data_fname[-4:] == ".czi" or raw_data_fname[-4:] == ".CZI":
            a      =  czifile.CziFile(raw_data_fname)                                           # read info about pixel size
            b      =  a.metadata()
            start  =  b.find("ScalingZ")
            end    =  b[start + 9:].find("ScalingZ")
            if start != -1:
                self.pix_sizeZ  =  float(b[start + 9:start + 7 + end]) * 1000000
                start           =  b.find("ScalingX")
                end             =  b[start + 9:].find("ScalingX")
                self.pix_sizeX  =  float(b[start + 9:start + 7 + end]) * 1000000
            elif start == -1:                                                                   # in case metadata file are not easy to read
                start2  =  b.find('<Distance Id="X"')
                if start2 != -1:
                    kk              =  b[start2:].find('</Value>')
                    self.pix_sizeX  =  float(b[start2 + 35:start2 + kk]) * 1000000
                    cc              =  b.find('<Distance Id="Z">')
                    cc1             =  b[cc:].find('<Value>')
                    cc2             =  b[cc:].find('</Value>')
                    self.pix_sizeZ  =  float(b[cc + cc1 + 7:cc + cc2]) * 1000000
                elif start2 == -1:                                                               # in case metadata file are not easy to read
                    [self.pix_sizeX, self.pix_sizeZ]  =  SetPixelSize().getPixelsValues()

            filedata  =  np.squeeze(czifile.imread(raw_data_fname))

            if len(filedata.shape) > 3:                                                                                 # filedata can be even single channel
                if chs_spts_nucs[0] != -1:
                    spts_a_bff   =  filedata[chs_spts_nucs[0]]
                    self.spts_a  =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_a.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.spts_a[z]  =  np.rot90(spts_a_bff[z, :, ::-1])

                if chs_spts_nucs[1] != -1:
                    spts_b_bff   =  filedata[chs_spts_nucs[1]]
                    self.spts_b  =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_b.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.spts_b[z]  =  np.rot90(spts_b_bff[z, :, ::-1])

                if chs_spts_nucs[2] != -1:
                    nucs_bff   =  filedata[chs_spts_nucs[2]]
                    self.nucs  =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.nucs.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.nucs[z]  =  np.rot90(nucs_bff[z, :, ::-1])

                if chs_spts_nucs[3] != -1:
                    memb_bff   =  filedata[chs_spts_nucs[3]]
                    self.memb  =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.memb.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.memb[z]  =  np.rot90(memb_bff[z, :, ::-1])

            else:
                if chs_spts_nucs[0] == 0:
                    spts_a_bff   =  filedata[chs_spts_nucs[0]]
                    self.spts_a  =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)      # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_a.shape[0]):                                                                      # rotated and mirrored to match ImageJ standards
                        self.spts_a[z] = np.rot90(spts_a_bff[z, :, ::-1])

                if chs_spts_nucs[1] == 0:
                    spts_b_bff   =  filedata[chs_spts_nucs[1]]
                    self.spts_b  =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)      # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_b.shape[0]):                                                                      # rotated and mirrored to match ImageJ standards
                        self.spts_b[z]  =  np.rot90(spts_b_bff[z, :, ::-1])

                if chs_spts_nucs[2] == 0:
                    nucs_bff   =  filedata[chs_spts_nucs[2]]
                    self.nucs  =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)             # if images are not square, initialization of the final matrix is needed
                    for z in range(self.nucs.shape[0]):                                                                         # rotated and mirrored to match ImageJ standards
                        self.nucs[z]  =  np.rot90(nucs_bff[z, :, ::-1])

                if chs_spts_nucs[3] == 0:
                    memb_bff   =  filedata[chs_spts_nucs[3]]
                    self.memb  =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.memb.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.memb[z]  =  np.rot90(memb_bff[z, :, ::-1])

        if raw_data_fname[-4:] == ".tif" or raw_data_fname[-4:] == ".TIF" or raw_data_fname[-4:] == ".TIFF" or raw_data_fname[-4:] == ".lsm":
            aa      =  tifffile.TiffFile(raw_data_fname)                                                                                   # read info about pixel size

            if raw_data_fname[-4:] == ".tif" or raw_data_fname[-4:] == ".TIF":
                try:
                    start           =  aa.imagej_metadata['Info'].find('ScalingZ')
                    if start != -1:
                        self.pix_sizeZ  =  np.float32(aa.imagej_metadata['Info'][start:start + aa.imagej_metadata['Info'][start:].find('\n')][14:]) * 1000000             # result in µm
                        start           =  aa.imagej_metadata['Info'].find('ScalingX')
                        self.pix_sizeX  =  np.float32(aa.imagej_metadata['Info'][start:start + aa.imagej_metadata['Info'][start:].find('\n')][14:]) * 1000000             # result in µm
                    elif start == -1:                                                                   # in case metadata file are not easy to read
                        [self.pix_sizeX, self.pix_sizeZ]  =  SetPixelSize().getPixelsValues()
                except KeyError:
                    self.pix_sizeX, self.pix_sizeZ  =  SetPixelSize().getPixelsValues()


            if raw_data_fname[-4:] == ".lsm":
                self.pix_sizeZ  =  TiffFile(str(raw_data_fname)).lsm_metadata["VoxelSizeZ"] * 1000000             # result in µm
                self.pix_sizeX  =  TiffFile(str(raw_data_fname)).lsm_metadata["VoxelSizeX"] * 1000000             # result in µm

            filedata  =  np.squeeze(tifffile.imread(raw_data_fname))

            if len(filedata.shape) > 3:                                                                                 # filedata can be even single channel
                if chs_spts_nucs[0] != -1:
                    spts_a_bff   =  filedata[:, chs_spts_nucs[0]]
                    self.spts_a  =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)     # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_a.shape[0]):                                                                    # rotated and mirrored to match ImageJ standards
                        self.spts_a[z, :, :]  =  np.rot90(spts_a_bff[z, :, ::-1])

                if chs_spts_nucs[1] != -1:
                    spts_b_bff   =  filedata[:, chs_spts_nucs[1]]
                    self.spts_b  =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)     # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_b.shape[0]):                                                                    # rotated and mirrored to match ImageJ standards
                        self.spts_b[z]  =  np.rot90(spts_b_bff[z, :, ::-1])

                if chs_spts_nucs[2] != -1:
                    nucs_bff   =  filedata[:, chs_spts_nucs[2]]
                    self.nucs  =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)           # if images are not square, initialization of the final matrix is needed
                    for z in range(self.nucs.shape[0]):                                                                      # rotated and mirrored to match ImageJ standards
                        self.nucs[z]  =  np.rot90(nucs_bff[z, :, ::-1])

                if chs_spts_nucs[3] != -1:
                    memb_bff   =  filedata[:, chs_spts_nucs[3]]
                    self.memb  =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.memb.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.memb[z]  =  np.rot90(memb_bff[z, :, ::-1])

            else:                                                                                                           # filedata has only one channel, just check its name
                if chs_spts_nucs[0] == 0:
                    spts_a_bff  =  filedata
                    self.spts_a  =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)   # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_a.shape[0]):                                                                   # rotated and mirrored to match ImageJ standards
                        self.spts_a[z]  =  np.rot90(spts_a_bff[z, :, ::-1])

                if chs_spts_nucs[1] == 0:
                    spts_b_bff   =  filedata[:, chs_spts_nucs[1]]
                    self.spts_b  =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)   # if images are not square, initialization of the final matrix is needed
                    for z in range(self.spts_b.shape[0]):                                                                   # rotated and mirrored to match ImageJ standards
                        self.spts_b[z] = np.rot90(spts_b_bff[z, :, ::-1])

                if chs_spts_nucs[2] == 0:
                    nucs_bff   =  filedata[:, chs_spts_nucs[2]]
                    self.nucs  =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)         # if images are not square, initialization of the final matrix is needed
                    for z in range(self.nucs.shape[0]):                                                                     # rotated and mirrored to match ImageJ standards
                        self.nucs[z]  =  np.rot90(nucs_bff[z, :, ::-1])

                if chs_spts_nucs[3] == 0:
                    memb_bff   =  filedata[chs_spts_nucs[3]]
                    self.memb  =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                    for z in range(self.memb.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                        self.memb[z]  =  np.rot90(memb_bff[z, :, ::-1])

        if raw_data_fname[-4:] == ".lif":
            aa              =  read_lif.Reader(raw_data_fname)                                                                             # read info about pixel size
            series          =  aa.getSeries()[0]
            self.pix_sizeZ  =  series.getMetadata()['voxel_size_z']
            self.pix_sizeX  =  series.getMetadata()['voxel_size_x']

            if chs_spts_nucs[0] != -1:
                spts_a_bff   =  series.getFrame(T=0, channel=chs_spts_nucs[0], dtype=np.int16)
                self.spts_a  =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(self.spts_a.shape[0]):                                                                                       # rotated and mirrored to match ImageJ standards
                    self.spts_a[z, :, :]  =  np.rot90(spts_a_bff[z, :, ::-1])

            if chs_spts_nucs[1] != -1:
                spts_b_bff   =  series.getFrame(T=0, channel=chs_spts_nucs[1], dtype=np.int16)
                self.spts_b  =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(self.spts_b.shape[0]):                                                                                        # rotated and mirrored to match ImageJ standards
                    self.spts_b[z, :, :]  =  np.rot90(spts_b_bff[z, :, ::-1])

            if chs_spts_nucs[2] != -1:
                nucs_bff   =  series.getFrame(T=0, channel=chs_spts_nucs[2], dtype=np.int16)
                self.nucs  =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(self.nucs.shape[0]):                                                                                           # rotated and mirrored to match ImageJ standards
                    self.nucs[z, :, :]  =  np.rot90(nucs_bff[z, :, ::-1])

            if chs_spts_nucs[3] != -1:
                memb_bff   =  series.getFrame(T=0, channel=chs_spts_nucs[3], dtype=np.int16)
                self.memb  =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(self.spts_b.shape[0]):                                                                                        # rotated and mirrored to match ImageJ standards
                    self.memb[z, :, :]  =  np.rot90(memb_bff[z, :, ::-1])


class MosaicRawDataLoader:
    """Load and combine 3 different images as a mosaic."""
    def __init__(self, raw_data_fname, chs_spts_nucs):
        """This function loads only .czi file organized as 3 mosaic images.
        This class mounts them removing the black separation. It works with only
        this setting."""

        a      =  czifile.CziFile(raw_data_fname)                                           # read info about pixel size
        b      =  a.metadata()
        start  =  b.find("ScalingZ")
        end    =  b[start + 9:].find("ScalingZ")
        if start != -1:
            self.pix_sizeZ  =  float(b[start + 9:start + 7 + end]) * 1000000
            start           =  b.find("ScalingX")
            end             =  b[start + 9:].find("ScalingX")
            self.pix_sizeX  =  float(b[start + 9:start + 7 + end]) * 1000000
        elif start == -1:                                                                   # in case metadata file are not easy to read
            start2  =  b.find('<Distance Id="X"')
            if start2 != -1:
                kk              =  b[start2:].find('</Value>')
                self.pix_sizeX  =  float(b[start2 + 35:start2 + kk]) * 1000000
                cc              =  b.find('<Distance Id="Z">')
                cc1             =  b[cc:].find('<Value>')
                cc2             =  b[cc:].find('</Value>')
                self.pix_sizeZ  =  float(b[cc + cc1 + 7:cc + cc2]) * 1000000
            elif start2 == -1:                                                               # in case metadata file are not easy to read
                [self.pix_sizeX, self.pix_sizeZ]  =  SetPixelSize().getPixelsValues()

        filedata  =  np.squeeze(czifile.imread(raw_data_fname))

        if len(filedata.shape) > 3:                                                                                 # filedata can be even single channel
            if chs_spts_nucs[0] != -1:
                spts_a_bff  =  filedata[chs_spts_nucs[0]]
                spts_a2     =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(spts_a2.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                    spts_a2[z]  =  np.rot90(spts_a_bff[z, :, ::-1])

                numb_zeros                                       =  np.where(spts_a2[spts_a2.shape[0] // 2, :, spts_a2.shape[2] // 2] == 0)[0]
                shift_numb                                       =  np.argmax(np.diff(numb_zeros)) + 1
                self.spts_a                                      =  np.zeros((spts_a2.shape[0], spts_a2.shape[1] - numb_zeros.size, spts_a2.shape[2]), dtype=spts_a2.dtype)
                self.spts_a[:, :numb_zeros[0]]                   =  spts_a2[:, :numb_zeros[0]]
                self.spts_a[:, numb_zeros[0]:2 * numb_zeros[0]]  =  spts_a2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.spts_a[:, 2 * numb_zeros[0]:]               =  spts_a2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[1] != -1:
                spts_b_bff  =  filedata[chs_spts_nucs[1]]
                spts_b2     =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(spts_b2.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                    spts_b2[z]  =  np.rot90(spts_b_bff[z, :, ::-1])

                numb_zeros                                       =  np.where(spts_b2[spts_b2.shape[0] // 2, :, spts_b2.shape[2] // 2] == 0)[0]
                shift_numb                                       =  np.argmax(np.diff(numb_zeros)) + 1
                self.spts_b                                      =  np.zeros((spts_b2.shape[0], spts_b2.shape[1] - numb_zeros.size, spts_b2.shape[2]), dtype=spts_b2.dtype)
                self.spts_b[:, :numb_zeros[0]]                   =  spts_b2[:, :numb_zeros[0]]
                self.spts_b[:, numb_zeros[0]:2 * numb_zeros[0]]  =  spts_b2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.spts_b[:, 2 * numb_zeros[0]:]               =  spts_b2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[2] != -1:
                nucs_bff  =  filedata[chs_spts_nucs[2]]
                nucs2     =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(nucs2.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                    nucs2[z]  =  np.rot90(nucs_bff[z, :, ::-1])

                numb_zeros                                     =  np.where(nucs2[nucs2.shape[0] // 2, :, nucs2.shape[2] // 2] == 0)[0]
                shift_numb                                     =  np.argmax(np.diff(numb_zeros)) + 1
                self.nucs                                      =  np.zeros((nucs2.shape[0], nucs2.shape[1] - numb_zeros.size, nucs2.shape[2]), dtype=nucs2.dtype)
                self.nucs[:, :numb_zeros[0]]                   =  nucs2[:, :numb_zeros[0]]
                self.nucs[:, numb_zeros[0]:2 * numb_zeros[0]]  =  nucs2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.nucs[:, 2 * numb_zeros[0]:]               =  nucs2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[3] != -1:
                memb_bff  =  filedata[chs_spts_nucs[3]]
                memb2     =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)  # if images are not square, initialization of the final matrix is needed
                for z in range(memb2.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                    memb2[z]  =  np.rot90(memb_bff[z, :, ::-1])
                numb_zeros                                     =  np.where(memb2[memb2.shape[0] // 2, :, memb2.shape[2] // 2] == 0)[0]
                shift_numb                                     =  np.argmax(np.diff(numb_zeros)) + 1
                self.memb                                      =  np.zeros((memb2.shape[0], memb2.shape[1] - numb_zeros.size, memb2.shape[2]), dtype=memb2.dtype)
                self.memb[:, :numb_zeros[0]]                   =  memb2[:, :numb_zeros[0]]
                self.memb[:, numb_zeros[0]:2 * numb_zeros[0]]  =  memb2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.memb[:, 2 * numb_zeros[0]:]               =  memb2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

        else:
            if chs_spts_nucs[0] == 0:
                spts_a_bff  =  filedata[chs_spts_nucs[0]]
                spts_a2     =  np.zeros(np.rot90(spts_a_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_a_bff.dtype)      # if images are not square, initialization of the final matrix is needed
                for z in range(spts_a2.shape[0]):                                                                      # rotated and mirrored to match ImageJ standards
                    spts_a2[z] = np.rot90(spts_a_bff[z, :, ::-1])
                numb_zeros                                       =  np.where(spts_a2[spts_a2.shape[0] // 2, :, spts_a2.shape[2] // 2] == 0)[0]
                shift_numb                                       =  np.argmax(np.diff(numb_zeros))
                self.spts_a                                      =  np.zeros((spts_a2.shape[0], spts_a2.shape[1] - numb_zeros.size, spts_a2.shape[2]), dtype=spts_a2.dtype)
                self.spts_a[:, :numb_zeros[0]]                   =  spts_a2[:, :numb_zeros[0]]
                self.spts_a[:, numb_zeros[0]:2 * numb_zeros[0]]  =  spts_a2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.spts_a[:, 2 * numb_zeros[0]:]               =  spts_a2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[1] == 0:
                spts_b_bff  =  filedata[chs_spts_nucs[1]]
                spts_b2     =  np.zeros(np.rot90(spts_b_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=spts_b_bff.dtype)      # if images are not square, initialization of the final matrix is needed
                for z in range(spts_b2.shape[0]):                                                                      # rotated and mirrored to match ImageJ standards
                    spts_b2[z]  =  np.rot90(spts_b_bff[z, :, ::-1])
                numb_zeros                                       =  np.where(spts_b2[spts_b2.shape[0] // 2, :, spts_b2.shape[2] // 2] == 0)[0]
                shift_numb                                       =  np.argmax(np.diff(numb_zeros))
                self.spts_b                                      =  np.zeros((spts_b2.shape[0], spts_b2.shape[1] - numb_zeros.size, spts_b2.shape[2]), dtype=spts_b2.dtype)
                self.spts_b[:, :numb_zeros[0]]                   =  spts_b2[:, :numb_zeros[0]]
                self.spts_b[:, numb_zeros[0]:2 * numb_zeros[0]]  =  spts_b2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.spts_b[:, 2 * numb_zeros[0]:]               =  spts_b2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[2] == 0:
                nucs_bff  =  filedata[chs_spts_nucs[2]]
                nucs2     =  np.zeros(np.rot90(nucs_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=nucs_bff.dtype)             # if images are not square, initialization of the final matrix is needed
                for z in range(nucs2.shape[0]):                                                                         # rotated and mirrored to match ImageJ standards
                    nucs2[z]  =  np.rot90(nucs_bff[z, :, ::-1])
                numb_zeros                                     =  np.where(nucs2[nucs2.shape[0] // 2, :, nucs2.shape[2] // 2] == 0)[0]
                shift_numb                                     =  np.argmax(np.diff(numb_zeros)) + 1
                self.nucs                                      =  np.zeros((nucs2.shape[0], nucs2.shape[1] - numb_zeros.size, nucs2.shape[2]), dtype=nucs2.dtype)
                self.nucs[:, :numb_zeros[0]]                   =  nucs2[:, :numb_zeros[0]]
                self.nucs[:, numb_zeros[0]:2 * numb_zeros[0]]  =  nucs2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.nucs[:, 2 * numb_zeros[0]:]               =  nucs2[:, 2 * numb_zeros[0] + 2 * shift_numb:]

            if chs_spts_nucs[3] == 0:
                memb_bff  =  filedata[chs_spts_nucs[3]]
                memb2     =  np.zeros(np.rot90(memb_bff[:, :, ::-1], axes=(1, 2)).shape, dtype=memb_bff.dtype)             # if images are not square, initialization of the final matrix is needed
                for z in range(memb2.shape[0]):                                                                         # both are rotated and mirrored to match ImageJ standards
                    memb2[z]  =  np.rot90(memb_bff[z, :, ::-1])
                numb_zeros                                     =  np.where(memb2[memb2.shape[0] // 2, :, memb2.shape[2] // 2] == 0)[0]
                shift_numb                                     =  np.argmax(np.diff(numb_zeros))
                self.memb                                      =  np.zeros((memb2.shape[0], memb2.shape[1] - numb_zeros.size, memb2.shape[2]), dtype=memb2.dtype)
                self.memb[:, :numb_zeros[0]]                   =  memb2[:, :numb_zeros[0]]
                self.memb[:, numb_zeros[0]:2 * numb_zeros[0]]  =  memb2[:, numb_zeros[0] + shift_numb:2 * numb_zeros[0] + shift_numb]
                self.memb[:, 2 * numb_zeros[0]:]               =  memb2[:, 2 * numb_zeros[0] + 2 * shift_numb:]


class SetPixelSize(QtWidgets.QDialog):
    """Choose if cluster nuclei horizontally or vertically"""
    def __init__(self, parent=None):
        super().__init__(parent)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        size_xy_lbl  =  QtWidgets.QLabel("X-Y size (µ)", self)
        size_xy_lbl.setFixedSize(int(ksf_h * 100), int(ksf_w * 22))

        size_xy_edt  =  QtWidgets.QLineEdit()
        size_xy_edt.setFixedSize(int(ksf_h * 50), int(ksf_w * 22))
        size_xy_edt.textChanged[str].connect(self.size_xy_var)

        size_xy_box  =  QtWidgets.QHBoxLayout()
        size_xy_box.addWidget(size_xy_lbl)
        size_xy_box.addWidget(size_xy_edt)

        size_z_lbl  =  QtWidgets.QLabel("Z size (µ)", self)
        size_z_lbl.setFixedSize(int(ksf_h * 100), int(ksf_w * 22))

        size_z_edt  =  QtWidgets.QLineEdit()
        size_z_edt.setFixedSize(int(ksf_h * 50), int(ksf_w * 22))
        size_z_edt.textChanged[str].connect(self.size_z_var)

        size_z_box  =  QtWidgets.QHBoxLayout()
        size_z_box.addWidget(size_z_lbl)
        size_z_box.addWidget(size_z_edt)

        send_btn  =  QtWidgets.QPushButton("Ok")
        send_btn.clicked.connect(self.send)
        send_btn.setToolTip("Insert values")
        send_btn.setFixedSize(int(ksf_h * 60), int(ksf_w * 25))

        send_box  =  QtWidgets.QHBoxLayout()
        send_box.addStretch()
        send_box.addWidget(send_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(size_xy_box)
        layout.addLayout(size_z_box)
        layout.addLayout(send_box)

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, int(ksf_h * 20), int(ksf_w * 25))
        self.setWindowTitle("Set Pixel Size")

    def size_xy_var(self, text):
        """Set X-Y pixel size."""
        self.size_xy_value  =  float(text)

    def size_z_var(self, text):
        """Set Z pixel size."""
        self.size_z_value  =  float(text)

    def params(self):
        """Function to send choice."""
        return [self.size_xy_value, self.size_z_value]

    def send(self):
        """Input values."""
        self.close()

    @staticmethod
    def getPixelsValues(parent=None):
        """Send choice."""
        dialog  =  SetPixelSize(parent)
        result  =  dialog.exec_()
        sizess  =  dialog.params()
        return sizess






# from czifile import CziFile
#
# with CziFile(raw_data_fname) as czi:
#     for attachment in czi.attachments():
#         print(attachment)
#
#         if attachment.attachment_entry.name == 'Scaling':
#             timestamps = attachment.data()
#             break
#     else:
#         raise ValueError('TimeStamps not found')
#
# print(timestamps)
