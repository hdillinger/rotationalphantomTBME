# /home/dillin01/rotationalphantom/.venv/bin/python /home/dillin01/rotationalphantom/recon_all_new.py

import os
import string
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import mapvbvd

from scipy import ndimage
from tqdm import tqdm

from pint import Quantity

from natsort import natsorted
from pathlib import Path
import re

# from explore import Explore

def recon_dat(data_file_path=None):
    # Load and prepare raw data
    print(f'\nLoading `{data_file_path}` ...')
    twix_obj = mapvbvd.mapVBVD(data_file_path)

    if isinstance(twix_obj, list):
        if 'cmr' in data_file_path.absolute().as_posix().lower():
            twix_obj[-1].image.flagRemoveOS = False
            # with warnings.catch_warnings(record=True) as w:
            rawdata = np.array(twix_obj[-1].image.unsorted(), dtype=np.complex64)
                # if w:
                    # return None
            
        if 'BEAT' in data_file_path.absolute().as_posix() or 'cor' in data_file_path.absolute().as_posix():
            twix_obj[-1].image.flagRemoveOS = True
            # load only first heart phase due to memory issues
            # rawdata = np.array(twix_obj[-1].image[:, :, :, :, :, :, :, :, :, :, 0, :, :, :, :, :], dtype=np.complex64)
            
            rawdata = np.array(twix_obj[-1].image.unsorted(), dtype=np.complex64)
            
            print(twix_obj[-1]['hdr']['Meas']['aflMaxwellCoefficients'])
        
        else:
            twix_obj[-1].image.flagRemoveOS = False
            # rawdata = np.array(twix_obj[-1].image[''], dtype=np.complex64)
            
            # load only first heart phase due to memory issues
            rawdata = np.array(twix_obj[-1].image[:, :, :, :, :, :, :, :, :, :, 0, :, :, :, :, :], dtype=np.complex64)
        
        seq_hash_twix = twix_obj[-1].hdr['Dicom']['tSequenceVariant']
    else:
        # rawdata = np.array(twix_obj.image.unsorted(), dtype=np.complex64)
        rawdata = np.array(twix_obj[-1].image[''], dtype=np.complex64)
        seq_hash_twix = twix_obj.hdr['Dicom']['tSequenceVariant']

    if len(seq_hash_twix) == 32:
        print(f'Raw data contains Pulseq-file signature {seq_hash_twix}')

    print(rawdata.shape)
    
    # if 'avgs' in data_file_path.absolute().as_posix():
    #     try:
    #         rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], 16, rawdata.shape[0], 4))
            
    #         rawdata = np.mean(rawdata, axis=2)
    #     # incorrectly named avg scan?
    #     except: 
    #         rawdata = np.reshape(rawdata[:, :, :4*rawdata.shape[0]], (rawdata.shape[0], rawdata.shape[1], rawdata.shape[0], 4))
    # else:
        # rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], rawdata.shape[0], 4))
    
    if 'cmr' in data_file_path.absolute().as_posix().lower():
        if 'dyn' in data_file_path.absolute().as_posix() or 'avgs' in data_file_path.absolute().as_posix():
            try:      
                rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], -1, rawdata.shape[0], 4))
                rawdata = np.mean(rawdata, axis=2, keepdims=True)
            except:
                return None
        else:
            rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], rawdata.shape[0], 4))
    
    if 'BEAT' in data_file_path.absolute().as_posix() or 'FLASH' in data_file_path.absolute().as_posix() or 'cor' in data_file_path.absolute().as_posix():            
        if 'dyn' in data_file_path.absolute().as_posix():
            try:      
                rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], -1, rawdata.shape[0]+1, 4))
                
                rawdata = np.mean(rawdata, axis=2, keepdims=True)
            except:
                return None
        elif 'avgs' in data_file_path.absolute().as_posix():
            try:  
                print('hello')
                rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], -1, 10, 4))
                rawdata = np.mean(rawdata, axis=3, keepdims=True)
            except:
                return None
        else:
            rawdata = np.reshape(rawdata, (rawdata.shape[0], rawdata.shape[1], rawdata.shape[0], 4))
        
        
        #     # if rawdata.shape[10] > 1:
        #         # rawdata = rawdata[:, :, :, :, :, :, :, :, :, :, 0, ...]
            
        #     rawdata = np.mean(rawdata, axis=5, keepdims=True)
        # else:
        #     rawdata = np.mean(rawdata, axis=10, keepdims=True)
        
    
    rawdata = rawdata.squeeze()
    print(rawdata.shape)
    # Explore(rawdata, loop=[1, 3], lim=(0, np.max(np.abs(rawdata[:]))))
    
    # coil compression with svd
    shape_without_coils = ((rawdata.shape[0],) + rawdata.shape[2:])
    
    non_singleton_dims = [i for i, size in enumerate(rawdata.shape) if size > 1]
    print(non_singleton_dims)
    
    shape_list = (non_singleton_dims[0],) + (*non_singleton_dims[2:],) + (1,)
    print(shape_list)
    
    data_array_flattened = np.reshape(np.transpose(rawdata, shape_list), (-1, rawdata.shape[1]))

    u, s, vh = np.linalg.svd(data_array_flattened, full_matrices=False)
    data_array_single_coil = np.reshape(u[:, 0] * s[0], shape_without_coils)

    rawdata = data_array_single_coil
    
    data = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(rawdata, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))

    
    return data



if __name__ == '__main__':
    files_list = []

    # path = Path('./data/mr')
    path = Path(r'/echo/_allgemein/projects/Hannes/rotphantom_PUB_v2/raw')
    # path = Path(r'/echo/_allgemein/projects/Hannes/rotphantom_PUB_v2/raw/old_data')
    path_list = list(path.glob('**/*.dat'))
    
    # path_list = [p for p in path_list if "BEAT" not in p.name]
    # path_list = [p for p in path_list if "cmr" not in p.name]
    
    # path_list = [p for p in path_list if not "avgs" in p.name]
    
    path_list = [p for p in path_list if not "meas_MID00052_FID15911_BEAT_FQ_avgs_noRFspoil.dat" in p.name]
    
    
    # print(path_list)
    
    # path_list = [p for p in path_list if "venc50_" in p.name]
    
    path_list = natsorted(path_list, key=str)
    
    fig0, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 4)) # figsize: pixel size for latex
    
    
    random_matrix = np.random.random_sample((128, 128))
    venc = 100
    
    im0 = axes[0].imshow(random_matrix, cmap='gray', origin='lower')
    im1 = axes[1].imshow(random_matrix, cmap='twilight', origin='lower', vmin=-venc, vmax=venc)
    im2 = axes[2].imshow(random_matrix, cmap='twilight', origin='lower', vmin=-venc, vmax=venc)
    im3 = axes[3].imshow(random_matrix, cmap='twilight', origin='lower', vmin=-venc, vmax=venc)
    

    for file in tqdm(path_list, desc=f'Listing dat files in directory {path}'):
        
        file_posix = file.absolute().as_posix()
        # if 'cor' in file_posix:
        # if 'BEAT' in file_posix:
        if 'cmr' in file_posix.lower():
        
        # if 'localizer' not in file_posix and 'AdjCoilSens' not in file_posix and 'AdjFre' not in file_posix and 'AdjFieldMap' not in file_posix and 'meas_MID00165_FID16024_GOOD4_CMRBernstein_300rpm_venc200_BW500_16avgs_DIR2' not in file_posix:
            
            venc = re.search(r'venc(\d+)', file_posix)
            if venc is None:
                venc = 100
            else:
                venc = float(venc.group(1))
                
            rpm = re.search(r'rpm(\d+)', file_posix)
            if rpm is None:
                rpm = 0
            else:
                rpm = float(rpm.group(1))
                
            dir = re.search(r'dir(\d+)', file_posix)
            if dir is None:
                dir = 0
            else:
                dir = float(dir.group(1))
                
            # dyn = re.search(r'dyn(\d+)', file_posix)
            # if dyn is None:
            #     dyn = 0
            # else:
            #     dyn = float(dyn.group(1))
            
            file_name_png = Path(file.stem + '.png')
            print(file_name_png)
            
            if file_name_png.exists():
                continue

            data = recon_dat(file)
            
            
            # try:
            #     data = recon_dat(file)
            # except ValueError:
            #     data = None
            
            # Explore(data)
            
            if data is None:
                continue
            
            centre = np.int64(np.round(ndimage.center_of_mass(np.abs(np.mean(data, axis=2)))))

            # homogeneous coordinate matrix, no rotation, no shift
            rot = np.eye(4, 4)
            # shift the center by integer offset
            offset = centre - tuple(ti//2 for ti in data.shape[0:2])
            print(f'Shifting image by {offset[0]} and {offset[1]} pixels.')
            # adjust matrix according to offset
            rot[0:2, -1] = offset

            # for plotting hlines and vlines later
            center_new = tuple(ti//2 for ti in data.shape[0:2]) + offset

            # apply translation to image
            # data = ndimage.affine_transform(input=data, matrix=rot, output=np.complex64, mode='constant', cval=np.nan)
            
            # mag = np.abs(np.mean(data, axis=2))
            mag = np.mean(np.abs(data), axis=2)
            
            # mag = np.abs(data[:, :, 0])
            
            # q1, q3 = np.percentile(mag, [25, 75])  # Compute Q1 and Q3
            q1, q3 = np.percentile(mag, [5, 95])  # Compute Q1 and Q3
            im0.set_clim(vmin=q1, vmax=q3)
            im0.set_data(mag)
            
            im1.set_data(np.angle(data[:, :, 1] / data[:, :, 0]) / np.pi * venc)
            im2.set_data(np.angle(data[:, :, 2] / data[:, :, 0]) / np.pi * venc)
            im3.set_data(np.angle(data[:, :, 3] / data[:, :, 0]) / np.pi * venc)
            
            
            print(f'Saving file as {file_name_png}')
            plt.savefig(file_name_png)
            
            np.save(file.stem + '.npy', data)
            