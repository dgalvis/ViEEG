# -*- coding: utf-8 - *-
# Volume based ViEEG signal reconstruction of ictal epochs using beamforming algorithsms
# tags: seizure, ictal, MEG, volume-based, source reconstruction, beamformers
# Author: Miao Cao
# Date: 28 Nov 2018
# Email: miao.cao@unimelb.edu.au


# import libraries
import os
import sys
import hdf5storage
import pprint
import json
import datetime
import platform
import socket

# libraries essential to source signal reconstruction
import numpy as np
import mne
from mne import io

print(__doc__)


def beamforming_time_series(epoch, fwd, beamforming_method, raw_cov=None):
    """
    The function to reconstruct source time-series using beamformer algorithms.
    Two algorithsms are implemented: LCMV and DICS

    Parameters
    ----------
        epoch: object
            mne epoch object

        fwd: object
            mne forward model object

        beamforming_method: string
            beamformer algorithsm to use, values 'lcmv' or 'dics'

        raw_cov: object
            mne covariance matrix object

    Returns
    -------
        list_src_ts: list of objects
            A list of reconstructed source-space activity object (mne)

        beamforming_parameter_dict: dictionary
            a dictionary that contains parameter information of source reconstruction methods
    """
    list_src_ts = None
    filters = None
    beamforming_parameter_dict = None
    # make filters according to different beamforming methods
    if beamforming_method is "lcmv":
        # 5.2.a.1 define parameters
        beamforming_parameter_dict = {
            "beamforming_name": beamforming_method,
            "info": epoch.info,
            "forward": fwd,
            "data_cov": None,
            "reg": 0.01,
            "noise_cov": None,
            "label": None,
            "pick_ori": "max-power",
            "rank": None,
            "weight_norm": "unit-noise-gain",
            # 'weight_norm': 'nai',
            "reduce_rank": False,
            "depth": None,  # not use depth correction for now
            "verbose": True,
        }

        # 5.2.a.2 generate data covariance matrix from epoch
        data_cov_method = "auto"
        data_tmin, data_tmax = None, None
        beamforming_parameter_dict["data_cov"] = mne.compute_covariance(
            epochs=epoch,
            keep_sample_mean=True,
            tmin=data_tmin,
            tmax=data_tmax,
            projs=None,
            method=data_cov_method,
            method_params=None,
            cv=3,
            scalings=None,
            n_jobs=12,
            return_estimators=False,
            on_mismatch="raise",
            rank=None,
            verbose=True,
        )
        # 5.2.a.3 generate noise covariance matrix from either pre-event data or raw_cov (pre-computed)
        noise_cov_method = "auto"
        noise_tmin, noise_tmax = None, None
        if raw_cov is None:
            noise_tmin, noise_tmax = None, 0
            beamforming_parameter_dict["noise_cov"] = mne.compute_covariance(
                epochs=epoch,
                keep_sample_mean=True,
                tmin=noise_tmin,
                tmax=noise_tmax,
                projs=None,
                method=noise_cov_method,
                method_params=None,
                cv=3,
                scalings=None,
                n_jobs=12,
                return_estimators=False,
                on_mismatch="raise",
                rank=None,
                verbose=True,
            )
        else:
            noise_tmin, noise_tmax = "start raw data", "end raw data"
            beamforming_parameter_dict["noise_cov"] = raw_cov

        # 5.2.a.4 make inverse operate using beamforming
        filters = mne.beamformer.make_lcmv(
            info=beamforming_parameter_dict["info"],
            forward=beamforming_parameter_dict["forward"],
            data_cov=beamforming_parameter_dict["data_cov"],
            reg=beamforming_parameter_dict["reg"],
            noise_cov=beamforming_parameter_dict["noise_cov"],
            label=beamforming_parameter_dict["label"],
            pick_ori=beamforming_parameter_dict["pick_ori"],
            rank=beamforming_parameter_dict["rank"],
            weight_norm=beamforming_parameter_dict["weight_norm"],
            reduce_rank=beamforming_parameter_dict["reduce_rank"],
            depth=beamforming_parameter_dict["depth"],
            verbose=beamforming_parameter_dict["verbose"],
        )

        # 5.3.1 apply inverse operator to a seizure epoch
        # 5.3.1.a use lcmv
        evoked = epoch.average()
        list_src_ts = mne.beamformer.apply_lcmv(
            evoked=evoked,
            filters=filters,
            # max_ori_out='signed',
            verbose=True,
        )

        # list_src_ts = mne.beamformer.apply_lcmv_epochs(epochs=epoch,
        #                                                filters=filters,
        #                                                max_ori_out='signed',
        #                                                return_generator=False,
        #                                                verbose=True)
        beamforming_parameter_dict["info"] = "raw_meg_data_info is not writable here."
        beamforming_parameter_dict["forward"] = pprint.pformat(
            list(fwd.values()))
        beamforming_parameter_dict["data_cov"] = (
            "covariance matrix estimated from seizure epoch from "
            + str(data_tmin)
            + "s to "
            + str(data_tmax)
            + "s"
            + " using method "
            + data_cov_method
        )
        beamforming_parameter_dict["noise_cov"] = (
            "covariance matrix estimated from seizure epoch from "
            + str(noise_tmin)
            + "s to "
            + str(noise_tmax)
            + "s"
            + " using method "
            + noise_cov_method
        )
    elif beamforming_method is "dics":
        # 5.2.b.1 define parameters
        data_tmin, data_tmax = None, None
        noise_tmin, noise_tmax = None, None
        beamforming_parameter_dict = {
            "beamforming_name": beamforming_method,
            "info": epoch.info,
            "forward": fwd,
            "csd": None,
            "reg": 0.01,
            "label": None,
            "pick_ori": "max-power",
            "rank": None,
            "inversion": "single",
            "weight_norm": "unit-noise-gain",
            "normalize_fwd": True,
            "real_filter": True,
            "reduce_rank": True,
            "verbose": True,
        }

        # 5.2.b.2 compute cross spectrum density matrix from data

        # 5.2.b.3 make inverse operator using beamforming
        filters = mne.beamformer.make_dics(
            info=beamforming_parameter_dict["info"],
            forward=beamforming_parameter_dict["forward"],
            csd=beamforming_parameter_dict["csd"],
            reg=beamforming_parameter_dict["reg"],
            label=beamforming_parameter_dict["label"],
            pick_ori=beamforming_parameter_dict["pick_ori"],
            rank=beamforming_parameter_dict["rank"],
            inversion=beamforming_parameter_dict["inversion"],
            weight_norm=beamforming_parameter_dict["weight_norm"],
            normalize_fwd=beamforming_parameter_dict["normalize_fwd"],
            real_filter=beamforming_parameter_dict["real_filter"],
            reduce_rank=beamforming_parameter_dict["reduce_rank"],
            verbose=beamforming_parameter_dict["verbose"],
        )

        # 5.3.1 apply inverse operator to a seizure epoch
        # 5.3.1.b use dics method
        list_src_ts = mne.beamformer.apply_dics_epochs(
            epochs=epoch, filters=filters, return_generator=False, verbose=True
        )
        beamforming_parameter_dict["info"] = "raw_meg_data_info is not writable here."
        beamforming_parameter_dict["forward"] = pprint.pformat(
            list(fwd.values()))
        beamforming_parameter_dict["data_cov"] = (
            "covariance matrix estimated from seizure epoch from "
            + str(data_tmin)
            + "s to "
            + str(data_tmax)
            + "s"
        )
        beamforming_parameter_dict["noise_cov"] = (
            "covariance matrix estimated from seizure epoch from "
            + str(noise_tmin)
            + "s to "
            + str(noise_tmax)
            + "s"
        )
    else:
        print("Error: it needs to use a valid beamforming method!")

    return list_src_ts, beamforming_parameter_dict


# Raw data path. This path is pointing to EEG-MEG data repository. It's
# normally a different location from analysis scripts.
# By default, the sample rate is 1000 Hz. Sample rates can be 1000 Hz or 5000 Hz across this dataset
MEG_sample_rate = 1000
# Data paths
raw_data_path = "path to raw MEG data folder"
preprocessed_data_path = "path to preprocessed data folder (save)"
post_analysis_data_path = "create a folder and path to post analysis data folder (save)"
freesurfer_subjects_dir = "path to Freesurfer folder"
hostname = socket.gethostname()
IPAddr = socket.gethostbyname(hostname)
# define switches for modules
# switch for resample block
switch_downsample_epoch_signals = False # do not downsample signals
epoch_resample_freq = 500
# switch for data crop in time after source recon
switch_crop_time_recon_src_data = True
recon_src_data_crop_time_tmin = -30
recon_src_data_crop_time_tmax = 30
# switch for visualising virtual grid coregistration
switch_check_virtual_gride_coordinate_locations = True
if IPAddr == "172.25.134.219": # if it is run at the computing server (a static ip address)
    freesurfer_subjects_dir = "/data/mcao/projects/freesurfer"
    switch_check_virtual_gride_coordinate_locations = False

# Apply band-pass filter to MEG sensor signals
# use FIR or IIR filter
use_FIR_filter = False
use_IIR_filter = True
l_freq = 0.1
h_freq = 200

# Patient list, a list of patient IDs.
patient_list = [
    "P012_EM",
    "P016_SD",
    "P035_IA",
    "P037_AB",
    "P040_EW",
    "P042_KD",
    "P046_KX",
    "P049_JC2",
    "P061_AE",
    "P070_PM",
    "P072_TB",
    "P074_AF",
]

patient_vg_config_analysis = {
    "P042_KD": {
        "vg_ver": "0.32",
        "vg_config": "_1x8x8_grid_2x1x6_depth_FS-RAS_ViEEG_v",
    },
    "P016_SD": {"vg_ver": "0.24", "vg_config": "_1x9x9_1x8x8_grid_FS-RAS_ViEEG_v"},
    "P061_AE": {"vg_ver": "0.22", "vg_config": "_2x10x5_grid_FS-RAS_ViEEG_v"},
    "P046_KX": {"vg_ver": "0.12", "vg_config": "_1x8x8_grid_FS-RAS_ViEEG_v"},
    "P040_EW": {
        "vg_ver": "0.24",
        "vg_config": "_1x8x8_1x8x4_1x8x8_projected_grid_FS-RAS_ViEEG_v",
    },
    "P012_EM": {"vg_ver": "0.14", "vg_config": "_1x9x9_grid_FS-RAS_ViEEG_v"},
    "P049_JC2": {"vg_ver": "0.22", "vg_config": "_1x10x10_1x4x4_grid_FS-RAS_ViEEG_v"},
    "P037_AB": {"vg_ver": "0.23", "vg_config": "_1x10x10_1x8x4_grid_FS-RAS_ViEEG_v"},
    "P072_TB": {"vg_ver": "0.21", "vg_config": "_1x10x10_1x8x8_grid_FS-RAS_ViEEG_v"},
    "P035_IA": {
        "vg_ver": "0.32",
        "vg_config": "_1x4x10_1x6x6_1x10x10_grid_FS-RAS_ViEEG_v",
    },
    "P074_AF": {"vg_ver": "0.15", "vg_config": "_1x8x12_grid_FS-RAS_ViEEG_v"},
    "P070_PM": {
        "vg_ver": "0.51",
        "vg_config": "_3x8x8_grid_2x1x8_depth_FS-RAS_ViEEG_v",
    },
}

# patient seizure profile. This includes seizure data file and seizure ids.
patient_seizure_profile = {
    "P012_EM": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2, 3, 4, 5, 6, 7],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P016_SD": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2, 3, 4, 5, 6],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P035_IA": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["STI002"],
        "trigger_thres": 4.8,
    },
    "P037_AB": {
        "seizure_data_file": "_eyesclosed_2_tsss.fif",
        "seizure_id": [1],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P040_EW": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2, 3, 4, 5, 6, 7, 8],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P042_KD": {
        "seizure_data_file": "_eyesclosed_MEGonly_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P046_KX": {
        "seizure_data_file": "_eyesclosed_MEGandEEG_FIXED_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P049_JC2": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P061_AE": {
        "seizure_data_file": "_MEGonly_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["NO_CHAN"],
        "trigger_thres": 1.8,
    },
    "P070_PM": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1, 2],
        "trigger_ch": ["NO_CHAN"],
        "trigger_thres": 1.8,
    },
    "P072_TB": {
        "seizure_data_file": "_MEGonly_tsss.fif",
        "seizure_id": [1],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
    "P074_AF": {
        "seizure_data_file": "_eyesclosed_tsss.fif",
        "seizure_id": [1],
        "trigger_ch": ["STI101"],
        "trigger_thres": 1.8,
    },
}

# Select a patient and get its patient ID
patient_ID = patient_list[5]

# Get patient number (patient number start with P and followed by a
# 3-digit number)
patient_number = patient_ID[0:4]
# Get patient initials
patient_initials = patient_ID[5:6]

# MEG data file name
if patient_ID in patient_seizure_profile:
    MEG_data_file_suffix_name = patient_seizure_profile[patient_ID]["seizure_data_file"]
    timing_trigger_ch_name = patient_seizure_profile[patient_ID]["trigger_ch"]
    timing_trigger_thres = patient_seizure_profile[patient_ID]["trigger_thres"]
else:
    sys.exit(patient_ID + " does not have any seizures in MEG recording.")

MEG_data_file = patient_ID + MEG_data_file_suffix_name

raw_MEG_data_path = raw_data_path + patient_ID + "/MEG/"

# --------------------------------------
# Important parameters for beamforming
# --------------------------------------
pre_seizure_time_length_sec = 600
post_seizure_time_length_sec = 600
seizure_epoch_baseline = (None, -10)
# ViEEG configuration identifier (internal)
vg_config = patient_vg_config_analysis[patient_ID]["vg_config"]
# ViEEG configuration version identifier (internal)
vg_config_ver = patient_vg_config_analysis[patient_ID]["vg_ver"]
# ViEEG electrode coodinates
fname_vg_coords = (
    raw_data_path
    + patient_ID
    + "/Virtual_Grid/"
    + patient_ID
    + vg_config
    + vg_config_ver
    + ".pom"
)
gen_time = datetime.datetime.now().strftime("%H-%M-%S_%d-%m-%Y")
# keep a record of system and toolbox information for each data file generated
str_system_info = (
    "mne-python version "
    + mne.__version__
    + ". Python version "
    + platform.python_version()
    + ". OS "
    + platform.platform()
    + ". Machine "
    + platform.machine()
    + "."
)
# keep a record of parameters used in the pipeline for each data file generated
result_parameter_list = {
    "matching_result_file_name": None,
    "patient_ID": patient_ID,
    "patient_seizure_profile": patient_seizure_profile[patient_ID],
    "MEG_raw_data_file": raw_MEG_data_path + MEG_data_file,
    "result_generation_time": gen_time,
    "BEM_parameter_list": None,
    "source_recon_family": "beamforming",
    "source_recon_parameter_list": None,
    "virtual_grid_configuration": vg_config,
    "virtual_grid_configuration_version": vg_config_ver,
    "virtual_grid_configuration_file": fname_vg_coords,
    "generating_script": __file__,
    "system_info": str_system_info,
}

# BEM model parameter list
# ico4 6.2mm, ico5 4.9mm 4098 sources per hemisphere
source_spacing = "N/A virtual grid defined"
bem_spacing = None
BEM_method = "curry"
src_surface_name = "N/A virtual grid defined"
# path to meFLASH images. Only used when meFLASH sequence is used.
meFLASH_image_path = ""

BEM_parameter_list = {
    "source_spacing": source_spacing,
    "bem_spacing": bem_spacing,
    "BEM_method": BEM_method,
    "src_surface_name": src_surface_name,
    "meFLASH_image_path": meFLASH_image_path,
}

result_parameter_list["BEM_parameter_list"] = BEM_parameter_list

# beamforming method
source_recon_family = "beamforming"
beamforming_method = "lcmv"
# beamforming_method = 'dics'


# --------------------------------------
# Virtual grid configuration
# --------------------------------------
# virtual grid version
vg_config = result_parameter_list["virtual_grid_configuration"]
vg_config_ver = result_parameter_list["virtual_grid_configuration_version"]
# --------------------------------------
# Read virtual grid coordinates from pom file generated from Curry
# --------------------------------------
fname_trans = (
    freesurfer_subjects_dir + "/" + patient_ID + "/bem/" + patient_ID + "-trans.fif"
)
print("Reading virtual grid configuration from " + fname_vg_coords)

# Read virtual grid coordinates from pom file generated from Curry
vg_file = open(fname_vg_coords, "r", encoding="utf-8-sig")
vg_file_lines = vg_file.readlines()
vg_coordinates = []
vg_elect_labels = []
num_vg_elects = int(vg_file_lines[0])
for read_cursor in np.arange(1, num_vg_elects + 1):
    line = vg_file_lines[read_cursor]
    line_list = line.strip("\n").split("\t")
    vg_coordinates.append([float(num) for num in line_list])

for read_cursor in np.arange(num_vg_elects + 1, num_vg_elects * 2 + 1):
    line = vg_file_lines[read_cursor]
    line.strip("\n").split("\t")
    vg_elect_labels.append(line)
vg_elect_labels = np.array(vg_elect_labels)

# Check the virtual electrode locations.
# Coregistered figure will be plotted when true
if switch_check_virtual_gride_coordinate_locations is True:
    raw_MEG_data = mne.io.read_raw_fif(
        raw_MEG_data_path + MEG_data_file, allow_maxshield=True, preload=False
    )
    # get measurement information from raw MEG data
    raw_MEG_data_info = raw_MEG_data.info
    np_vg_coordinates = np.array(vg_coordinates) * 0.001  # convert to meters
    fig1 = mne.viz.plot_alignment(
        info=raw_MEG_data_info,
        trans=fname_trans,
        subject=patient_ID,
        subjects_dir=freesurfer_subjects_dir,
        meg=True,
        dig=True,
        mri_fiducials=True,
        surfaces="pial",
        coord_frame="mri",
    )
    from mayavi import mlab

    mlab.points3d(
        np_vg_coordinates[:, 0],
        np_vg_coordinates[:, 1],
        np_vg_coordinates[:, 2],
        mode="2darrow",
        line_width=0.5,
        figure=fig1,
    )
# ---------------------------------------

# --------------------------------------
# read raw data
# Read and preprocess raw MEG data
# --------------------------------------

# --------------------------------------
# 2.1 Read bad channel names from an external file, badchannels
# This file is located in RawData, patient folder, Empty_Room_Recording
# Bad channels are labelled by the radiographer in Swinburne MEG centre.
bad_ch_names = []
bad_channel_file_path = raw_data_path + \
    patient_ID + "/Empty_Room_Recording/badchannels"
bad_ch_file = open(bad_channel_file_path, "r")
bad_ch_names = bad_ch_file.read().rstrip().split(",")
if bad_ch_names[0] == "":
    bad_ch_names = []
# --------------------------------------

# --------------------------------------
# 2.2 Read raw MEG data
# read MEG recording which has at lease one seizure/ictal events.
# --------------------------------------
raw_MEG_data = mne.io.read_raw_fif(
    raw_MEG_data_path + MEG_data_file, allow_maxshield=True, preload=True
)

# add bad channels to data.info structure
if len(bad_ch_names) != 0:
    raw_MEG_data.info["bads"] += bad_ch_names
# get MEG sample rate
MEG_sample_rate = int(raw_MEG_data.info["sfreq"])
# get time of end of recording
time_end_recording = int(raw_MEG_data.times[-1])
# recording acquisition time offset
start_sample_recording = raw_MEG_data.first_samp
# get measurement information from raw MEG data
raw_MEG_data_info = raw_MEG_data.info
# ---------------------------------------

# ---------------------------------------
# 2.3 Preprocess sensor signals
# ---------------------------------------
# only pick MEG channels (Magnetometres and gradiometres)
picks = mne.pick_types(raw_MEG_data.info, meg=True, eeg=False, eog=False)
# Apply notch filter to remove line noise - 50Hz and its harmonics
raw_MEG_data = raw_MEG_data.notch_filter(
    np.arange(50, 300, 50),
    picks=picks,
    filter_length="auto",
    n_jobs=12,
    phase="zero",
    verbose=True,
)

# Apply band-pass filter, using FIR or IIR filter

if use_FIR_filter is True:
    print("Band-pass filter raw MEG signals using a FIR filter")
    raw_MEG_data = raw_MEG_data.filter(
        l_freq=l_freq,
        h_freq=h_freq,
        picks=picks,
        l_trans_bandwidth="auto",
        h_trans_bandwidth="auto",
        n_jobs=12,
        method="fir",
        fir_window="hann",
        filter_length="auto",
        phase="zero",
        verbose=True,
    )
elif use_IIR_filter is True:
    print("Band-pass filter raw MEG signals using a IIR filter")
    iir_params = dict(order=4, ftype="butter", output="sos")
    iir_params = mne.filter.construct_iir_filter(
        iir_params=iir_params,
        f_pass=[l_freq, h_freq],
        f_stop=None,
        sfreq=MEG_sample_rate,
        btype="bandpass",
        return_copy=True,
        verbose=True,
    )
    raw_MEG_data = raw_MEG_data.filter(
        l_freq=l_freq,
        h_freq=h_freq,
        picks=picks,
        l_trans_bandwidth="auto",
        h_trans_bandwidth="auto",
        n_jobs=12,
        method="iir",
        iir_params=iir_params,
        fir_window="hann",
        fir_design="firwin",
        verbose=True,
    )
# ---------------------------------------

# ---------------------------------------
# 3. Compute noise covariance matrix
# from empty room recording
# ---------------------------------------
# Comments: empty room recording has to be
# IAS on and Maxwell-filtered.
raw_empty_room_fname = (
    raw_data_path
    + patient_ID
    + "/Empty_Room_Recording/"
    + patient_ID
    + "_ERM_IASon_tsss.fif"
)

# ---------------------------------------
# If empty room recording is not Maxfiltered, use MNE-Python built-in
# Maxwell_filter to filter signals.
switch_use_Empty_Room_Recording = False
if switch_use_Empty_Room_Recording:
    if os.path.isfile(raw_empty_room_fname):
        # load empty room recording
        raw_empty_room = mne.io.read_raw_fif(
            raw_empty_room_fname, allow_maxshield=True, preload=True
        )

        # add bad channels to empty room recording
        raw_empty_room.info["bads"] = [
            bb for bb in raw_MEG_data.info["bads"] if "EEG" not in bb
        ]

        # add projections to empty room recording
        # raw_empty_room.add_proj(
        #     [pp.copy() for pp in raw_MEG_data.info['projs'] if 'EEG' not in pp['desc']])

        # Compute noise covariance matrix
        ERM_noise_cov = mne.compute_raw_covariance(
            raw_empty_room, tmin=0, tmax=None)
    else:  # empty room recording is not maxwell-filtered, use MNE-Python built-in function
        raw_empty_room_fname = (
            raw_data_path
            + patient_ID
            + "/Empty_Room_Recording/"
            + patient_ID
            + "_ERM_IASon.fif"
        )
        # load empty room recording
        raw_empty_room = mne.io.read_raw_fif(
            raw_empty_room_fname, allow_maxshield=True, preload=True
        )

        # add bad channels to empty room recording
        raw_empty_room.info["bads"] = [
            bb for bb in raw_MEG_data.info["bads"] if "EEG" not in bb
        ]

        raw_empty_room_sss = mne.preprocessing.maxwell_filter(
            raw_empty_room, coord_frame="meg"
        )

        # add projections to empty room recording
        # raw_empty_room.add_proj(
        #     [pp.copy() for pp in raw_MEG_data.info['projs'] if 'EEG' not in pp['desc']])

        # Compute noise covariance matrix
        ERM_noise_cov = mne.compute_raw_covariance(
            raw_empty_room_sss, tmin=0, tmax=None
        )

        del raw_empty_room_sss, raw_empty_room
# ---------------------------------------


# ---------------------------------------
# 4. Epoch seizure events
# ---------------------------------------
epochs = []
seizure_range_list = []

# ---------------------------------------
# 4.1 Read event information
# ---------------------------------------
# read seizure/ictal events file
seizure_events_fname = (
    patient_ID + "_seizure" + MEG_data_file_suffix_name[0:-4] + ".eve"
)

events = mne.read_events(
    filename=preprocessed_data_path + patient_ID + "/MEG/" + seizure_events_fname
)

# calibrate timing to the first sample of acquisition
events[:, 0] = events[:, 0] + start_sample_recording
# number of events
num_events = int(events.shape[0] / 2)
# event ID
event_id = {"ictal event": 1}
# ---------------------------------------
# 4.2 Epoch the sensor signals
# ---------------------------------------
for i_event in np.arange(0, num_events):
    # epoch the signals according to seizure/ictal events
    seizure_length = (
        events[2 * i_event + 1][0] - events[2 * i_event][0]
    ) / MEG_sample_rate
    # tmin, tmax = -1 * seizure_length, seizure_length * 2
    tmin, tmax = (
        -1 * pre_seizure_time_length_sec,
        seizure_length + post_seizure_time_length_sec,
    )
    event_id = {"ictal event": 1}
    current_event = events[2 * i_event: 2 * i_event + 2, :]

    # check if this exceeds beginning or end of the recording.
    event_start_time = current_event[0, 0] / MEG_sample_rate

    # event start time relative to start of the recording
    event_recording_start_time = event_start_time - (
        start_sample_recording / MEG_sample_rate
    )
    # calibrate start of epoch and end of epoch if they exceed start or end of recording
    if tmin + event_recording_start_time <= 0:
        tmin = -1 * (event_recording_start_time - 0)
    if event_recording_start_time + tmax >= time_end_recording:
        tmax = time_end_recording - event_recording_start_time
    seizure_range_list.append([0, seizure_length])

    if seizure_epoch_baseline[1] + event_recording_start_time <= 0:
        # tim_seizure_epoch_baseline = -1 * \
        #     (event_recording_start_time - 0)
        tim_seizure_epoch_baseline = None
        tmax_seizure_epoch_baseline = -0.5 * (event_recording_start_time - 0)
        seizure_epoch_baseline = (
            tim_seizure_epoch_baseline,
            tmax_seizure_epoch_baseline,
        )

    epoch = mne.Epochs(
        raw_MEG_data,
        events=current_event,
        event_id=event_id,
        tmin=tmin,
        tmax=tmax,
        baseline=seizure_epoch_baseline,
        picks=picks,
        preload=True,
        detrend=1,
        verbose=True,
    )
    # if we downsample epoch signals
    if switch_downsample_epoch_signals is True:
        epoch = epoch.resample(
            sfreq=epoch_resample_freq,
            npad="auto",
            window="boxcar",
            n_jobs=12,
            pad="edge",
            verbose=True,
        )
    epochs.append(epoch)

# save some memory to avoid memory issues when processing using a workstaion (instead of a high performance computing server)
del raw_MEG_data

# ---------------------------------------
# 4.3 Save epoch sensor signals to a fif
# ---------------------------------------
# switch of this block
switch_save_epo_fif = False
if switch_save_epo_fif:
    for sz_n, epoch in enumerate(epochs):
        # crop_time_min = -1 * seizure_range_list[sz_n][1]
        # crop_time_max = 2 * seizure_range_list[sz_n][1]
        crop_time_min, crop_time_max = -15, 25
        crop_epoch = epoch.copy()
        # if this epoch is used to source recon we should not crop it
        crop_epoch = crop_epoch.crop(tmin=crop_time_min, tmax=crop_time_max)
        croped_seizure_epoch_fname = (
            preprocessed_data_path
            + patient_ID
            + "/MEG/MEG_maggrad_seizure_epoch_sz"
            + str(sz_n + 1)
            + "-epo.fif"
        )
        crop_epoch.save(croped_seizure_epoch_fname, verbose=True)

# ---------------------------------------
# 5. MEG source reconstruction
# using beamforming with
# method and parameters
# specified above
# ---------------------------------------

# ---------------------------------------
# 5.0 prepare variables and fif files
# ---------------------------------------
# create source estimate folder path
src_estimate_file_path = (
    post_analysis_data_path
    + patient_ID
    + "/MEG/"
    + "source_space/seizure/"
    + MEG_data_file_suffix_name[1:-4]
    + "/source_estimate/"
)
if not os.path.exists(src_estimate_file_path):
    os.makedirs(src_estimate_file_path)

# ---------------------------------------

# ---------------------------------------
# 5.1 make forward solution
# ---------------------------------------
# -bem-sol.fif BEM solution

# if forward model has been made by this or other pipelines, switch this off
# True when need to make a forward model
# read in virtual grid configurations
fwd = None
fname_bem = (
    freesurfer_subjects_dir
    + "/"
    + patient_ID
    + "/bem/"
    + patient_ID
    + "_"
    + BEM_method
    + "_ico"
    + str(bem_spacing)
    + "_bem.fif"
)
fname_bem_sol = (
    freesurfer_subjects_dir
    + "/"
    + patient_ID
    + "/bem/"
    + patient_ID
    + "_"
    + BEM_method
    + "_ico"
    + str(bem_spacing)
    + "_bem-sol.fif"
)
pos_rr = np.array([[1, 0, 0]])
pos_nn = np.array([[1, 0, 0]])
vg_coordinates = np.array(vg_coordinates) * 0.001  # convert to meters
for vg_coord in vg_coordinates:
    pos_rr = np.append(pos_rr, [vg_coord], axis=0)
    # pos_rr = np.append(pos_rr, [vg_coord], axis=0)
    # pos_rr = np.append(pos_rr, [vg_coord], axis=0)
    # pos_nn = np.append(pos_nn, [[1, 0, 0]], axis=0)
    # pos_nn = np.append(pos_nn, [[0, 1, 0]], axis=0)
    pos_nn = np.append(pos_nn, [[0, 0, 1]], axis=0)
pos_rr = np.delete(pos_rr, [0], axis=0)
pos_nn = np.delete(pos_nn, [0], axis=0)
pos = {"rr": pos_rr, "nn": pos_nn}
vol_src = mne.setup_volume_source_space(
    subject=patient_ID,
    pos=pos,
    mri=None,
    bem=fname_bem,
    surface=None,
    mindist=0.0,
    exclude=0.0,
    subjects_dir=freesurfer_subjects_dir,
    volume_label=None,
    add_interpolator=False,
    verbose=True,
)
# set up the forward model
fwd = mne.make_forward_solution(
    info=raw_MEG_data_info,
    trans=fname_trans,
    src=vol_src,
    bem=fname_bem_sol,
    meg=True,
    eeg=False,
    mindist=0.0,
    ignore_ref=False,
    n_jobs=12,
    verbose=True,
)
# Here we need to check if number of sources in volume source space
# match number of points pre_defined by virtual grid
# IMPORTANT NOTE: Forward model does not change the order of
# virtual grid points in volume source space.
switch_exit_when_num_vol_src_not_match = True
if fwd.get("src")[0]["nuse"] != num_vg_elects:
    inuse_vg_elects_indice = fwd.get("src")[0]["inuse"]
    inuse_vg_elects_indice_bool = np.array(
        [bool(x) for x in inuse_vg_elects_indice])
    vg_coordinates = vg_coordinates[inuse_vg_elects_indice_bool]
    vg_elect_labels = vg_elect_labels[inuse_vg_elects_indice_bool]
    if switch_exit_when_num_vol_src_not_match == True:
        sys.exit(
            "Number of volume sources does not match with number of electrodes in virtual grid that has been pre-defined."
            + " # of volume source in fwd is "
            + str(fwd.get("src")[0]["nuse"])
            + ", while # of virtual grid points is "
            + str(num_vg_elects)
            + ". Exit with this error!"
        )
    else:
        num_vg_elects = fwd.get("src")[0]["nuse"]
# ---------------------------------------


# ---------------------------------------
# 5.2 make inverse operator using beamforming
# and apply inverse operators to seizure epochs
# ---------------------------------------

# cycle through seizure/ictal events
for i_event in np.arange(0, num_events):

    # get this epoch
    epoch = epochs[i_event]
    # Whether apply a band-pass filter to this epoch/sensor signals
    # before source reconstruction.
    apply_bpf_to_epoch = False
    if apply_bpf_to_epoch == True:
        epoch = epoch.filter(
            l_freq=7,
            h_freq=25,
            picks="all",
            l_trans_bandwidth="auto",
            h_trans_bandwidth="auto",
            n_jobs=12,
            filter_length="auto",
            phase="zero",
            verbose=True,
        )

    # get the epoch time range in second
    seizure_range = seizure_range_list[i_event]
    # pre-define source estimate
    vol_src_est = None

    # ---------------------------------------
    # 5.2 make inverse operators using beamformers
    # ---------------------------------------

    # pre-define the filters and beamforming parameter list
    beamforming_parameter_dict = None
    vol_src_est, beamforming_parameter_dict = beamforming_time_series(
        epoch, fwd, beamforming_method, raw_cov=None
    )

    # IMPORTANT: If # of reconstructed signals does not match # of VG channels
    # then exit with errors.
    # if need to crop data in time, as parameters set
    if switch_crop_time_recon_src_data is True:
        recon_src_data_crop_time_tmax = recon_src_data_crop_time_tmax + \
            seizure_range[1]
        vol_src_est.crop(
            tmin=recon_src_data_crop_time_tmin,
            tmax=recon_src_data_crop_time_tmax,
            include_tmax=True,
        )
    src_sz_ts = vol_src_est.data
    if num_vg_elects != src_sz_ts.shape[0]:
        sys.exit(
            "Number of reconstructed signals does not match with number of electrodes in virtual grid that has been pre-defined."
            + " # of reconstructed signals is "
            + str(src_sz_ts.shape[0])
            + ", while # of virtual grid points is "
            + str(num_vg_elects)
            + ". Exit with this error!"
        )
    # ---------------------------------------
    # 5.3.2 save reconstructed source space to .stc files
    # ---------------------------------------
    # save reconstructed source signals to .stc file
    src_file_name = (
        src_estimate_file_path
        + patient_ID
        + "_seizure_sz"
        + str(i_event + 1)
        + "_VG"
        + vg_config_ver
        + "_vol-src-recon_"
        + source_recon_family
        + "_"
        + beamforming_method
        + "_time_"
        + gen_time
    )

    # save the parameters to a parameter file here
    # Here need to be done.
    # ---------------------------------------
    result_parameter_list["source_recon_parameter_list"] = beamforming_parameter_dict
    result_parameter_list["seizure_number"] = int(i_event + 1)
    result_parameter_list["seizure_range_second"] = seizure_range
    result_parameter_list["MEG_src_recon_sfreq_HZ"] = int(vol_src_est.sfreq)
    result_parameter_list["matching_result_file_name"] = src_file_name
    result_parameter_list_fname = src_file_name + "_parameter_list.json"
    with open(result_parameter_list_fname, "w") as file_pointer:
        json.dump(result_parameter_list, file_pointer)

    # ---------------------------------------
    # 5.3.3 save reconstructed source signals to a .mat file
    # ---------------------------------------
    src_sz_times = vol_src_est.times
    src_sz_ts = vol_src_est.data
    src_coords = fwd.get("source_rr")  # still in head space
    mat_file_name = src_file_name + ".mat"
    hdf5storage.savemat(
        mat_file_name,
        mdict={
            "patient_id": patient_ID,
            "results_generate_time": gen_time,
            "seizure_range": seizure_range,
            "seizure_number": i_event + 1,
            "sfreq": int(vol_src_est.sfreq),
            "vg_electrode_coordinates": src_coords,
            "vg_electrode_label": vg_elect_labels,
            "times": src_sz_times,
            "virtual_grid_seizure_time_series": src_sz_ts,
        },
    )
