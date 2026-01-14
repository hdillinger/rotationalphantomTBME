import os
import argparse
import sys
import re

from natsort import natsorted
from pathlib import Path

import tensorflow as tf

# 3rd Party dependencies 
from pint import Quantity
import pyvista
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from copy import deepcopy
from tqdm.notebook import tqdm
from time import perf_counter
from IPython.display import clear_output

sys.path.insert(0, "./cmrsim")
sys.path.insert(0, "./cmrseq")
sys.path.insert(0, "./cmrsim/notebooks")

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ['CUDA_VISIBLE_DEVICES'] = "1"

gpu = tf.config.get_visible_devices("GPU")
print(gpu)
if gpu:
    tf.config.experimental.set_memory_growth(gpu[0], True)
    
# Project library cmrsim
import local_functions
import cmrsim
import cmrseq
import cmrsim.utils.particle_properties as part_factory

from cmrseq.io import sequence_from_json
from cmrseq.io import PulseSeqFile


from scipy.spatial.transform import Rotation


# taken from https://gist.github.com/jakevdp/91077b0cae40f8f8244a
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    # return base.from_list(cmap_name, color_list, N)
    return plt.cm.colors.ListedColormap(color_list, color_list, N)


class ConstantRotationTrajectory(cmrsim.trajectory._base.BaseTrajectoryModule):
    """Implements an analytic trajectory of particles rotating about an arbitrary
    rotational axis in 3D with constant angular velocity.
    
    On instantiation, a 4x4 matrix :math:`T` is computed from the specified position and
    direction of the rotational axis, to transform the problem to a rotation
    about the z-axis in a new coordinate system.

    .. math::

        \\vec{r}\\prime = T^{-1}R_zT\\vec{r}

    
    :param angular_velocity:
    :param center_position: (3, ) 3D vector specifying the positional offset of 
            to rotational axis, in units of meter
    :param rotation_axis: (3, ) 3D vector specifying the direction of the rotational
            axis (this is normalized internally).
    """
    #: angular velocity in rotations per second (Hz)
    angular_velocity: float
    #: 3D position of the bottom end of the rotational axis
    center_position: tf.Tensor
    #: Rotational axis orientation
    rotation_axis: tf.Tensor
    #: (4, 4) Rotation matrix to transform the coordinate system in which 
    #: the rotation is performed about the z-axis
    _rot_mat: np.ndarray
    #: Inverse of the coordinate transform _rot_mat
    _rot_mat_inv: np.ndarray
    #: (4, 4) Translation in homogenous coordinates to account for positional
    #: offset of the rotation center
    _translate: np.ndarray
    #: inverse of _translate
    _translate_inv: np.ndarray
    
    def __init__(self, angular_velocity: Quantity,
                 center_position: Quantity,
                 rotation_axis: np.ndarray):
        # Save arguments:
        self.angular_velocity = angular_velocity.m_as("Hz")
        self.center_position = center_position.m_as("m")
        self.rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

        # As rotation about an arbitrary axes is treated as rotation about z-axis
        # where the coordinates are transformed beforehand, the corresponding 
        # transformation matrices need to be computed:
        # Rotation
        tmp, _ = Rotation.align_vectors(self.rotation_axis, np.array([0., 0., 1.]))
        self._rotmat = np.eye(4, 4)
        self._rotmat[:3, :3] = tmp.as_matrix()
        tmp_inv, _ = Rotation.align_vectors(np.array([0., 0., 1.]), self.rotation_axis)
        self._rotmat_inv = np.eye(4, 4)
        self._rotmat_inv[:3, :3] = tmp_inv.as_matrix()
        
        # Translation
        self._translate = np.eye(4, 4)
        self._translate[:3, 3] = -self.center_position
        self._translate_inv = np.eye(4, 4)
        self._translate_inv[:3, 3] = self.center_position
        
    def __call__(self, initial_positions: tf.Tensor, timing: tf.Tensor,
                 **kwargs) -> (tf.Tensor, dict):
        """
        :param initial_positions: (N, 3) tensor containing the initial
                    position of N particles
        :param timing: (t, ) 1D-tensor containing all time-points at which the
                    trajectory will be evaluated
        :return: (new positions (N, t, 3), additional fields)
        """
        delta_phi = 2 * np.pi * self.angular_velocity / 1000 * timing
        rotation_matrices = tf.stack([tf.stack([tf.math.cos(delta_phi), tf.math.sin(delta_phi), tf.zeros_like(delta_phi), tf.zeros_like(delta_phi)], axis=-1),
                                      tf.stack([-tf.math.sin(delta_phi), tf.math.cos(delta_phi), tf.zeros_like(delta_phi), tf.zeros_like(delta_phi)], axis=-1),
                                      tf.stack([tf.zeros_like(delta_phi), tf.zeros_like(delta_phi), tf.ones_like(delta_phi), tf.zeros_like(delta_phi)], axis=-1),
                                      tf.stack([tf.zeros_like(delta_phi), tf.zeros_like(delta_phi), tf.zeros_like(delta_phi), tf.ones_like(delta_phi)], axis=-1)],
                                     axis=-2)
        tmp = tf.concat([initial_positions, tf.ones_like(initial_positions[:, :1])], 1)
        tmp = tf.einsum("ij, nj -> ni", self._translate, tmp)
        tmp = tf.einsum("ij, tj -> ti", self._rotmat, tmp)
        tmp = tf.einsum("tij, nj -> nti", rotation_matrices, tmp)
        tmp = tf.einsum("ij, ntj -> nti", self._rotmat_inv, tmp)
        tmp = tf.einsum("ij, ntj -> nti", self._translate_inv, tmp)
        return tmp[..., :3], {}

    def increment_particles(self, particle_positions: tf.Tensor, dt: tf.Tensor, **kwargs):
        """
        :param particle_positions: (N, 3) tensor containing the initial
                    position of N particles
        :param dt: 0D-tensor containing the time-step for which the trajectory
                    will be evaluated
        :return: (new positions (N, 3), additional fields)
        """
        dt_expand = tf.reshape(dt, [1,])
        new_positions, _ = self(particle_positions, dt_expand)
        return new_positions[:, 0], {}


def main():
    parser = argparse.ArgumentParser(description="Python script for cmrsim to batch simulate seq files with rotational phantom trajectories.")
    
    parser.add_argument("--seq", type=str, help="Input sequence file", required=True)
    parser.add_argument("--dir", type=str, help="Directory of sequence file if different from seq file name", required=True)
    parser.add_argument("--rpm", type=int, help="Rotational frequency", required=False)

    args = parser.parse_args()

    # Output based on arguments
    if args.seq:
        seq_name = args.seq
        print(f"Sequence file {seq_name} will be used.")
    if args.dir:
        folder_name = args.dir
        print(f"Input directory set to {folder_name}.")
    if args.rpm:
        rpm = args.rpm
        print(f"Rotational frequency set to {rpm}.")
        rps = rpm / 60
    else:
        match = re.search(r"rpm(\d+)", args.seq)
        if match:
            rpm = int(match.group(1))
            print(f"Extracted RPM: {rpm}")
            rps = rpm / 60
        else:
            print("RPM value not found.")
        
    sequence_list = []
    print("Exists:", Path(folder_name, seq_name).exists())
    
    pathlist = list(Path(folder_name, seq_name).rglob('*.json'))
    pathlist = natsorted(pathlist, key=str)
        

    for path in tqdm(pathlist, desc='Loading sequence from json'):
        sequence_list.append(sequence_from_json(str(path)))

    f, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [5, 2]})
    cmrseq.plotting.plot_sequence(sequence_list[0], axes=a1, adc_yoffset=15.7)
    
    for seq in tqdm(sequence_list[1:], desc='Plotting sequence'):
        cmrseq.plotting.plot_sequence(seq, axes=[f.axes[-1], a1, a1, a1], format_axes=False, adc_yoffset=15.7)
    for s in tqdm(sequence_list[1:], desc='Plotting k-space'):
        cmrseq.plotting.plot_kspace_2d(s, ax=a2, k_axes=[0, 1])
        
    system_specs = sequence_list[0]._system_specs

    plt.savefig(os.path.join(folder_name, seq_name, f'{seq_name}_{rpm:d}Hz_sequenceJSON.png'))
    
    spacing = (0.0005, 0.0005, 0.0005)
    # spacing = (0.01, 0.01, 0.01)

    dimensions = (300, 300, 5)

    origin = tuple((-np.array(spacing) * np.array(dimensions) / 2).tolist())
    print(spacing, dimensions, origin)

    uniform_mesh = pyvista.ImageData(spacing=spacing, dimensions=dimensions, origin=origin)
    inner_radius = 0.0
    outer_radius = 0.0647

    # initial_positions = np.array(uniform_mesh.points, dtype=np.float32)
    initial_positions = np.array(uniform_mesh.points + np.random.standard_normal(uniform_mesh.points.shape) * np.max(spacing) / 2, dtype=np.float32)
    # n_seeds = int(10)
    # r = dimensions[0] * np.sqrt(np.random.uniform(0, 1, n_seeds))
    # theta = np.random.uniform(0, 1, n_seeds)
    # initial_positions = np.array([r*np.cos(theta), r*np.sin(theta)])

    n_mesh = 150
    extent = 0.1
    x_mesh, y_mesh, z_mesh = np.meshgrid(np.linspace(-extent, extent, n_mesh), np.linspace(-extent, extent, n_mesh), np.linspace(-extent, extent, n_mesh))

    # initial_positions = np.float32(np.column_stack([x_mesh.ravel(), y_mesh.ravel(), z_mesh.ravel()]))

    # print(np.max(initial_positions[:]))

    initial_positions = np.float32(initial_positions + 1e-6*np.random.standard_normal(initial_positions.shape))
    # initial_positions = np.float32(initial_positions)

    print(initial_positions.shape)

    initial_positions = initial_positions[np.where(np.logical_and(np.linalg.norm(initial_positions[:, :2], axis=1) > inner_radius,
                                                                np.linalg.norm(initial_positions[:, :2], axis=1) < outer_radius))]
    n_parts = initial_positions.shape[0]
    print(n_parts)
    properties = {"M0":  np.ones(n_parts, np.float32),
                "T1":  np.ones(n_parts, np.float32) * 400,
                "T2":  np.ones(n_parts, np.float32) * 300,
                "magnetization": part_factory.norm_magnetization()(n_parts),
                "initial_position": initial_positions,
                #   "off_res": np.ones(n_parts, np.float32).reshape(-1, 1) * initial_positions[:, 1].reshape(-1, 1)*0.1
                }
    
    trajectory_module = ConstantRotationTrajectory(angular_velocity=Quantity(rps, "Hz"),
                                               center_position=Quantity([0., 0., 0.], "m"),
                                               rotation_axis=np.array([0., 0., 1.]))
    
    for idx_dummy, seq in enumerate(sequence_list):
        if len(seq.partial_sequence(partial_string_match=('adc_0'), copy_blocks=False).blocks) > 0:
            break

    print(f'Found {idx_dummy} dummy shots in sequence list.')


    # Define off-resonsance submodule
    # offres = cmrsim.bloch.submodules.OffResonance(gamma=system_specs.gamma.m_as("MHz/T"), device="GPU:0")


    time_rf_dummy, rf_grid_dummy, grad_grid_dummy, adc_on_grid_dummy = [np.stack(v) for v in cmrseq.utils.grid_sequence_list(sequence_list[:idx_dummy], force_uniform_grid=False)]
    print([v.shape for v in [time_rf_dummy, rf_grid_dummy, grad_grid_dummy, adc_on_grid_dummy]])

    print("Input shapes for dummy Bloch Module:\n  ", "\n   ".join([f"{time_rf_dummy.shape=}", f"{rf_grid_dummy.shape=}",
                                                            f"{grad_grid_dummy.shape=}", f"{adc_on_grid_dummy.shape=}"]))

    time_rf, rf_grid, grad_grid, adc_on_grid = [np.stack(v) for v in cmrseq.utils.grid_sequence_list(sequence_list[idx_dummy:], force_uniform_grid=False)]
    print([v.shape for v in [time_rf, rf_grid, grad_grid, adc_on_grid]])


    print("Input shapes for Bloch Module:\n  ", "\n   ".join([f"{time_rf.shape=}", f"{rf_grid.shape=}",
                                                            f"{grad_grid.shape=}", f"{adc_on_grid.shape=}"]))

    module = cmrsim.bloch.GeneralBlochOperator(name="rf", gamma=system_specs.gamma_rad.m_as("rad/mT/ms"),
                                            time_grid=time_rf[0], gradient_waveforms=grad_grid,
                                            rf_waveforms=rf_grid, adc_events=adc_on_grid, device="GPU:0",
                                            submodules=[])

    spoiling_mod = cmrsim.bloch.PerfectSpoiling("spoiling_op")

    n_samples = tf.shape(module.time_signal_acc[0])[0]
    n_repetitions = tf.shape(module.gradient_waveforms)[0]
    
    r = properties.pop("initial_position")
    m = properties.pop("magnetization")

    for rbatch, mbatch, dict_batch in tqdm(tf.data.Dataset.from_tensor_slices((r, m, properties)).batch(int(1e5)), desc="Phantom Batches", leave=False):
        _ = module(trajectory_module=trajectory_module, initial_position=rbatch, magnetization=mbatch, run_parallel=True, **dict_batch)

    print(tf.stack(module.time_signal_acc, axis=0).numpy().shape)

    # 4d flow
    time_signal = tf.stack(module.time_signal_acc, axis=0).numpy().reshape((224, 4, 224))

    time_signal = tf.transpose(time_signal, perm=[1, 0, 2])

    np.save(os.path.join(folder_name, seq_name, f'{seq_name}_{rpm:d}Hz_time_signal.npy'), time_signal.numpy())

    image = tf.signal.fftshift(tf.signal.ifft2d(tf.signal.ifftshift(time_signal, axes=(1, 2))), axes=(1, 2)).numpy()

    image_mag = np.abs(image[0, :, :])

    image = image[1, :, :] / image[0, :, :]
    print(image.shape)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)

    abs_plot = axes[0].imshow((np.abs(np.squeeze(image_mag))), cmap="gray")
    fig.colorbar(abs_plot, ax=axes[0], fraction=0.05, pad=0.005,    orientation="horizontal")

    phase_plot = axes[1].imshow(np.angle(np.squeeze(image)), cmap=discrete_cmap(10, "twilight"), vmin=-np.pi, vmax=np.pi)

    # masked
    # phase_plot = axes[1].imshow(np.angle(np.squeeze(image))* (image_mag / np.max(image_mag[:])), cmap="twilight", vmin=-np.pi/2, vmax=np.pi/2)

    # phase_plot = axes[1].imshow(np.angle(np.squeeze(image))* (image_mag / np.max(image_mag[:])), cmap="twilight", vmin=-np.pi, vmax=np.pi)
    fig.colorbar(phase_plot, ax=axes[1], fraction=0.05, pad=0.005, orientation="horizontal")
    # [(_.set_yticks([]), _.set_xticks([])) for _ in axes]
    [_.set_title(t, fontsize=14) for _, t in zip(axes, ["Magnitude", "Phase"])]

    plt.savefig(os.path.join(folder_name, seq_name, f'{seq_name}_{rpm:d}Hz_result.png'))
    # os.path.join(folder_name, f'sequence{idx:04d}')

if __name__ == "__main__":
    main()
