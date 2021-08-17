"""
Simulation for forward and inverse modelling of
sythetic signals using beamforming.
Prepared for ViEEG manuscript (2021) submitted to
Nature Communications for revision.

This script uses the following Python packages:
Numpy (v1.20.2), scipy (v.1.6.2), mne (v0.22.1),
hdf5storage (v0.1.17), matplotlib.pylab (v1.20.2)

# tags: beamformers, forward models, inverse models
# Author: Miao Cao
# Date: 2 Apr. 2021
# Email: caomiao89@gmail.com
"""

import numpy as np
import scipy.signal
from scipy.io import loadmat, savemat
# import hdf5storage

import mne

import matplotlib.pyplot as plt
import sys

def read_ViEEG_coordinates(fname_vg_coords):
    '''
    This function reads the text file of ViEEG coordinates

    Parameters
    ----------
        fname_vg_coords: String
            A string, specifying path and file name of the text file containing coordinate
            and labels of ViEEG electrodes
    Returns
    -------
        no return
    '''
    print('Reading virtual grid configuration from ' + fname_vg_coords)

    vg_file = open(fname_vg_coords, 'r')
    vg_file_lines = vg_file.readlines()
    vg_coordinates = []
    vg_elect_labels = []
    num_vg_elects = int(vg_file_lines[0])
    for read_cursor in np.arange(1, num_vg_elects + 1):
        line = vg_file_lines[read_cursor]
        line_list = line.strip('\n').split('\t')
        # Read virtual grid coordinates from pom file generated from Curry
        vg_coordinates.append([float(num) for num in line_list])

    for read_cursor in np.arange(num_vg_elects + 1, num_vg_elects * 2 + 1):
        line = vg_file_lines[read_cursor]
        line.strip('\n').split('\t')
        vg_elect_labels.append(line)
    return np.array(vg_coordinates), np.array(vg_elect_labels)


def beamforming_lcmv_time_series(epoch, fwd):
    '''
    This function applies inverse modelling to
    sensor signals using beamformers implemented 
    by MNE-Python (Version 0.22.1)
    https://mne.tools/stable/index.html

    Parameters
    ----------
        epoch: mne.Epochs object
            Epochs (mne.Epochs object), containing segmented signals
        fwd: mne.Forward object
            Forward (mne.Forward object) model, containing all information
            including head model geometries to construct a forward solution

    Returns
    -------
        list_src_ts: list
            List of reconstructed source spaces
    '''
    beamforming_method = 'lcmv'
    filters = None
    beamforming_parameter_dict = None
    # make filters according to different beamforming methods
    if beamforming_method == 'lcmv':
        # 5.2.a.1 define parameters
        beamforming_parameter_dict = {'beamforming_name': beamforming_method,
                                      'info': epoch.info,
                                      'forward': fwd,
                                      'data_cov': None,
                                      'reg': 0.01,
                                      'noise_cov': None,
                                      'label': None,
                                      'pick_ori': 'max-power',
                                      'rank': 'info',
                                      'weight_norm': None,
                                      'reduce_rank': True,
                                      'depth': None,
                                      'inversion': 'single',
                                      'verbose': True}

        # 5.2.a.2 generate data covariance matrix from data
        data_tmin, data_tmax = None, None
        data_cov = mne.compute_covariance(epochs=epoch,
                                          keep_sample_mean=True,
                                          tmin=data_tmin,
                                          tmax=data_tmax,
                                          projs=None,
                                          method='shrunk',
                                          method_params=None,
                                          cv=3,
                                          scalings=None,
                                          n_jobs=1,
                                          return_estimators=False,
                                          on_mismatch='raise',
                                          rank=None,
                                          verbose=True)
        beamforming_parameter_dict['data_cov'] = data_cov
        # 5.2.a.3 generate noise covariance matrix from data
        noise_tmin, noise_tmax = None, 0
        noise_cov = mne.compute_covariance(epochs=epoch,
                                           keep_sample_mean=True,
                                           tmin=noise_tmin,
                                           tmax=noise_tmax,
                                           projs=None,
                                           method='shrunk',
                                           method_params=None,
                                           cv=3,
                                           scalings=None,
                                           n_jobs=1,
                                           return_estimators=False,
                                           on_mismatch='raise',
                                           rank=None,
                                           verbose=True)
        beamforming_parameter_dict['noise_cov'] = noise_cov

        # 5.2.a.4 make inverse operate using beamforming
        filters = mne.beamformer.make_lcmv(info=beamforming_parameter_dict['info'],
                                           forward=beamforming_parameter_dict['forward'],
                                           data_cov=beamforming_parameter_dict['data_cov'],
                                           reg=beamforming_parameter_dict['reg'],
                                           noise_cov=beamforming_parameter_dict['noise_cov'],
                                           label=beamforming_parameter_dict['label'],
                                           pick_ori=beamforming_parameter_dict['pick_ori'],
                                           rank=beamforming_parameter_dict['rank'],
                                           weight_norm=beamforming_parameter_dict['weight_norm'],
                                           reduce_rank=beamforming_parameter_dict['reduce_rank'],
                                           verbose=beamforming_parameter_dict['verbose'])

        # 5.3.1 apply inverse operator to a seizure epoch
        # 5.3.1.a use lcmv\
        evoked = epoch.average()
        list_src_ts = mne.beamformer.apply_lcmv(evoked=evoked,
                                                filters=filters,
                                                max_ori_out='signed',
                                                verbose=True)
        # list_src_ts = mne.beamformer.apply_lcmv_epochs(epochs=epoch,
        #                                                filters=filters,
        #                                                max_ori_out='signed',
        #                                                return_generator=False,
        #                                                verbose=True)
        beamforming_parameter_dict['info'] = 'raw_meg_data_info is not writable here.'
        beamforming_parameter_dict['forward'] = fwd.values()
        beamforming_parameter_dict['data_cov'] = 'covariance matrix estimated from seizure epoch from ' + str(
            data_tmin) + 's to ' + str(data_tmax) + 's'
        beamforming_parameter_dict['noise_cov'] = 'covariance matrix estimated from seizure epoch from ' + str(
            noise_tmin) + 's to ' + str(noise_tmax) + 's'

        return list_src_ts


def forward_inverse_modelling(time_series, vg_coordinates, base_folder, SNR=10, tmin_epoch=-5, tmax_epoch=35):
    '''
    This function implements forward and inverse modelling 
    using MNE-Python (Version 0.22.1) and beamformers
    https://mne.tools/stable/index.html

    Parameters
    ----------
        time_series: numpy.array
            Matrix with shape (n_channels, n_times) specifying values of n_channels signals over time
        vg_coordinates: numpy.array
            Matrix with shape (n_channels, x, y, z) specifying the coordinates (x, y, z) of n_channels
        base_folder: string
            String, specifying the folder path to where this script and data are placed
        SNR: float
            Signal to noise ratio (SNR), used to add the level of noise to forward modelled sensor signals

    Returns
    -------
        d_z_n: numpy.array
            Vector with shape (N,), the differential at time t
    '''
    # Select a patient and get its patient ID
    FreeSurfer_subjects_dir = base_folder
    patient_ID = 'Subject1'

    # Get patient number (patient number start with P and followed by a
    # 3-digit number)
    patient_number = patient_ID[0:4]
    # Get patient initials
    patient_initials = patient_ID[5:6]

    MEG_data_file = patient_ID + '_eyesclosed_MEGonly_tsss.fif'

    raw_MEG_data_path = base_folder

    src_estimate_file_path = base_folder

    RS_noise_cov_fname = base_folder + \
        'RS_Noise-cov.fif'

    # random seed to generate Gaussian noise added to sensor level
    # seed = np.random.randint(1, 10000)
    # scaler of Gaussian noise assigned to original source space
    scalor_Gaussian_noise_orig_src = np.sqrt(1 / SNR)

    # --------------------------------------
    # Important parameters
    # --------------------------------------
    # Source space parameters for MNE (not beamformer family)
    source_spacing = 'ico5'  # ico4 6.2mm, ico5 4.9mm 4098 sources per hemisphere
    bem_spacing = None
    BEM_method = 'Curry'
    src_surface_name = 'pial'
    # path to meFLASH images. Only used when meFLASH sequence is used.
    meFLASH_image_path = ''
    source_recon_family = 'beamforming'
    beamforming_method = 'lcmv'
    # --------------------------------------

    # --------------------------------------
    # 1. Read patient's raw MEG signals
    # and read in empty room recording
    # and compute noise covariance matrix
    raw_MEG_data = mne.io.read_raw_fif(
        raw_MEG_data_path + MEG_data_file, allow_maxshield=True, preload=False, on_split_missing='warn')

    MEG_sample_rate = int(raw_MEG_data.info['sfreq'])
    # Read noise covariance matrix from -cov.fif file
    RS_Noise_cov = mne.read_cov(RS_noise_cov_fname, verbose=True)
    # --------------------------------------

    # --------------------------------------
    # 2. Read forward model (pre-defined)
    # fwd_model_fname = src_estimate_file_path + patient_ID + \
    #     '_source_seizure_' + source_spacing + '_BEM-' + BEM_method + '_fwd.fif'
    # fwd = mne.read_forward_solution(
    #     fname=fwd_model_fname, verbose=True)
    # fwd = mne.convert_forward_solution(fwd=fwd,
    #                                    surf_ori=True,
    #                                    force_fixed=True,
    #                                    use_cps=True,
    #                                    verbose=True)

    # I think here I have to re-define volume-based source space and
    # generate another forward model???
    fname_bem = FreeSurfer_subjects_dir + patient_ID + '_' + BEM_method + \
        '_ico' + str(bem_spacing) + '_bem.fif'
    fname_trans = FreeSurfer_subjects_dir + patient_ID + '-trans.fif'
    fname_bem_sol = FreeSurfer_subjects_dir + patient_ID + '_' + BEM_method + \
        '_ico' + str(bem_spacing) + '_bem-sol.fif'
    num_ictal_vg_channels = time_series.shape[0]
    vg_coordinates = vg_coordinates[0: num_ictal_vg_channels, :]
    pos_rr = np.array([[1, 0, 0]])
    pos_nn = np.array([[1, 0, 0]])
    vg_coordinates = vg_coordinates * 0.001
    for vg_coord in vg_coordinates:
        pos_rr = np.append(pos_rr, [vg_coord], axis=0)
        # pos_rr = np.append(pos_rr, [vg_coord], axis=0)
        # pos_rr = np.append(pos_rr, [vg_coord], axis=0)
        pos_nn = np.append(pos_nn, [[1, 0, 0]], axis=0)
        # pos_nn = np.append(pos_nn, [[0, 1, 0]], axis=0)
        # pos_nn = np.append(pos_nn, [[0, 0, 1]], axis=0)
    pos_rr = np.delete(pos_rr, [0], axis=0)
    pos_nn = np.delete(pos_nn, [0], axis=0)
    pos = {'rr': pos_rr,
           'nn': pos_nn}

    # Check the virtual electrode locations.
    # Coregistered figure will be plotted when true
    switch_check_virtual_gride_coordinate_locations = False
    if switch_check_virtual_gride_coordinate_locations == True:
        np_vg_coordinates = np.array(
            vg_coordinates) * 0.001  # convert to meters
        fig1 = mne.viz.plot_alignment(
            info=raw_MEG_data.info,
            trans=fname_trans,
            subject=patient_ID,
            subjects_dir=FreeSurfer_subjects_dir,
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

    vol_src = mne.setup_volume_source_space(subject=patient_ID,
                                            pos=pos,
                                            mri=None,
                                            bem=fname_bem,
                                            surface=None,
                                            mindist=0.0,
                                            exclude=0.0,
                                            subjects_dir=FreeSurfer_subjects_dir,
                                            volume_label=None,
                                            add_interpolator=True,
                                            verbose=True)

    fwd = mne.make_forward_solution(info=raw_MEG_data.info,
                                    trans=fname_trans,
                                    src=vol_src,
                                    bem=fname_bem_sol,
                                    meg=True,
                                    eeg=False,
                                    mindist=0.0,
                                    ignore_ref=False,
                                    n_jobs=4,
                                    verbose=True)
    inuse_vertno = fwd['src'][0]['vertno']
    # --------------------------------------

    # --------------------------------------

    # --------------------------------------
    # 3. Set up simulated source space
    # a). read source space which has been set up.
    # b) set up simulated source space.
    # c) read simulated data from a Matlab file.
    # create volume based source estimate
    print('Now set up volume-based source estimate...')
    # downsample ictal spike time-series to 1000Hz
    # rs_time_series = mne.filter.resample(time_series, down=5, npad='auto')
    # channel_selection = np.arange(0, 32) # Only grid channels without removing bad channels
    # channel_selection = np.append(np.arange(0, 14), np.arange(15, 32)) # Only grid channels with bad channels removed.
    # grid+depth channels with bad channels removed
    # channel_selection = np.append(np.arange(0, 14), np.arange(15, 40))
    # channel_selection = np.arange(0, 40) # grid+depth channels without removing bad channels

    # ictal_spikes = rs_time_series[:, :]
    # Always normalise ictal spikes to between 0 nAm and 100 nAm. This is important for simulation!!!
    time_series = 1e-7 * \
        (time_series / np.max(np.abs(time_series)))
    time_series = time_series[-1 * len(inuse_vertno):, :]
    sim_stc_est = mne.VolSourceEstimate(data=time_series,
                                        vertices=[inuse_vertno],
                                        tmin=0.0,
                                        tstep=1 / MEG_sample_rate,
                                        subject=patient_ID,
                                        verbose=True)
    # sim_stc_est_data = np.zeros(sim_stc_est.data.shape)
    sim_stc_est_data = sim_stc_est.data
    sim_stc_est_vertices = sim_stc_est.vertices

    # # generate uncorrelated Gaussian random noise for original source space
    # num_ori_vertice = sim_stc_est_data.shape[0]
    # num_ori_time_samples = sim_stc_est_data.shape[1]
    # # to save time here. all zeros
    # # ori_src_Gaussian_bk_act = np.zeros(sim_stc_est_data.shape)
    # ori_src_Gaussian_bk_act = np.random.multivariate_normal(
    #     np.zeros(num_ori_vertice), np.identity(num_ori_vertice), num_ori_time_samples).T
    # # normalise the noise signals to the range of 0nAm ~ 100nAm or 1e-9Am ~ 100*1e-9Am
    # ori_src_Gaussian_bk_act = scalor_Gaussian_noise_orig_src * 1e-7 * \
    #     (ori_src_Gaussian_bk_act / np.max(np.abs(ori_src_Gaussian_bk_act)))
    #
    # # assign Gaussian background noise to original source space
    # sim_stc_est_data = ori_src_Gaussian_bk_act

    # src4vg_indices_inuse = src4vg_indices[0][0:ictal_spikes.shape[0]]
    # vg_lh_fwd_src_vertno = fwd['src'][0]['vertno'][src4vg_indices_inuse]
    #
    # index_sim_stc_est_data = np.array([], dtype=int)
    # for vert_vg in vg_lh_fwd_src_vertno:
    #     index_sim_stc_est_data = np.append(
    #         index_sim_stc_est_data, np.where(sim_stc_est_vertices[0] == vert_vg)[0])
    #
    # index_sim_stc_est_data = index_sim_stc_est_data
    # sim_stc_est_data[index_sim_stc_est_data,
    #                  :] = ictal_spikes
    # sim_stc_est.data = sim_stc_est_data

    # kwargs = dict(subjects_dir=FreeSurfer_subjects_dir, hemi='both', smoothing_steps=4,
    #               time_unit='s', initial_time=0.05, size=1200,
    #               views='lat', surface='inflated')
    # clim = dict(kind='value', pos_lims=[1e-9, 1e-8, 1e-7])
    # fig = mlab.figure(1)
    # brain_sim = sim_stc_est.plot(clim=clim, figure=fig, **kwargs)

    # --------------------------------------

    # --------------------------------------
    # 4. Apply forward modelling

    # sim_sensor_signals = mne.apply_forward(fwd=fwd,
    #                                        stc=sim_stc_est,
    #                                        info=raw_MEG_data.info,
    #                                        use_cps=True,
    #                                        verbose=True)
    # get data from sim_sensor_signals
    # transform into raw
    # sim_raw_MEG = mne.io.RawArray(sim_sensor_signals.data, sim_sensor_signals.info)
    # segment into epoch
    # epoch = mne.Epochs(sim_raw_MEG, events, event_id=event_id, tmin=-5, tmax=35, preload=True)

    # sim_sensor_signals = mne.apply_forward_raw(fwd=fwd,
    #                                            stc=sim_stc_est,
    #                                            info=raw_MEG_data.info,
    #                                            verbose=True)

    # here, add Gaussian noise to sensor signals
    sim_sensor_signals_evoked = mne.apply_forward(fwd=fwd,
                                                  stc=sim_stc_est,
                                                  info=raw_MEG_data.info,
                                                  use_cps=True,
                                                  verbose=True)
    # generate Gaussian noise based on covariance matrix
    Gaussian_noise_added_sensors = np.random.multivariate_normal(np.zeros(sim_sensor_signals_evoked.data.shape[0]),
                                                                 RS_Noise_cov.data,
                                                                 sim_sensor_signals_evoked.data.shape[1]).T
    Gaussian_noise_added_sensors = scalor_Gaussian_noise_orig_src * \
        np.abs(sim_sensor_signals_evoked.data).max(
        ) * Gaussian_noise_added_sensors / (np.abs(Gaussian_noise_added_sensors).max())
    sim_sensor_signals_evoked.data += Gaussian_noise_added_sensors
    sim_sensor_signals = mne.io.RawArray(
        sim_sensor_signals_evoked.data, sim_sensor_signals_evoked.info)
    events = np.array([[int(sim_sensor_signals.get_data().shape[1] / 3), 0, 1],
                       [int(sim_sensor_signals.get_data().shape[1] / 3) * 2, 1, 0]])
    event_id = {'ts': 1}
    epoch = mne.Epochs(sim_sensor_signals, events,
                       event_id=event_id, tmin=tmin_epoch, tmax=tmax_epoch, preload=True)
    # pass it to forward and inverse modelling
    # --------------------------------------

    # --------------------------------------
    # 5. Apply inverse operator
    recon_sim_src_est = beamforming_lcmv_time_series(epoch, fwd)
    # kwargs = dict(subjects_dir=FreeSurfer_subjects_dir, hemi='both', smoothing_steps=4,
    #               time_unit='s', initial_time=0.05, size=1200,
    #               views='lat')
    # clim = dict(kind='value', pos_lims=[1e-9, 1e-8, 1e-7])
    # figs = [mlab.figure(2)]
    # brain_sim_recon = recon_sim_src_est.plot(clim='auto', figure=figs, **kwargs)

    # --------------------------------------

    # --------------------------------------
    # 6. pull reconstructed signals
    recon_time_series = recon_sim_src_est.data.T
    time_series = np.transpose(time_series[:, int(sim_sensor_signals.get_data().shape[1] / 3) + tmin_epoch * MEG_sample_rate: int(
        sim_sensor_signals.get_data().shape[1] / 3) + tmax_epoch * MEG_sample_rate])
    # --------------------------------------

    return recon_time_series, time_series


if __name__ == '__main__':
    print(__doc__)
    # the folder where this script and data are placed
    base_folder_path = "./"
    # read/generate time-series

    # DANNY COMMENTS OUT
    #n_channels = 64  # no more than 64 channels
    #MEG_sfreq = 1000  # this is fixed in this simulation
    #time_sec_sim = 90  # simulate at least 60 seconds
    # here, you can replace your simulated time-series
    # just generated some Gaussian noise as an example
    #sim_time_series = np.random.randn(n_channels, MEG_sfreq*time_sec_sim)
    # here you can use the real iEEG seizure signals
    # ictal_data_dict = loadmat(
    #     base_folder_path + 'Subject1_iEEG_Raw_SZ-Epo_preprocessed_sz1.mat')
    # sim_time_series = mne.filter.resample(
    #     ictal_data_dict['sz_ts'][:40, :], down=5, npad='auto')

    #DANNY ADDITIONS
    ictal_data_dict = loadmat(base_folder_path + 'model_output' + sys.argv[1] + '.mat')
    n_channels = ictal_data_dict['sz_ts'].shape[0]
    MEG_sfreq = ictal_data_dict['sfreq'][-1, -1]
    time_sec_sim = ictal_data_dict['sz_times'][-1, -1]
    sim_time_series = ictal_data_dict['sz_ts'] # (I want to use the whole time series, with same sampling)
    # generate baseline Gaussian noise
    sigma, mu =0.05, 0
    baseline_sec = int(np.ceil(time_sec_sim / 2))
    baseline_time_series = sigma * np.random.randn(n_channels, MEG_sfreq*baseline_sec) + mu
    sim_time_series = np.concatenate((baseline_time_series, sim_time_series), axis=1)
    tmin_epoch, tmax_epoch = -10, 80 #(sz_times goes from 0 to 89.999 - 90000 points)
    # If possible can the line above be general to give 2/3 of the points back? Always with my times starting at t = 0 and ending at t = end
    # also can you include something that saves a new file with the output time series please? (as .mat file)



    # read virtual grid vg_coordinates here
    fname_vg_coords = base_folder_path + \
        'Subject1_ViEEG_8x8_grid_left_temporal_FS-RAS_Curry_COORDS_ONLY_v0.01.pom'
    vg_coordinates, vg_labels = read_ViEEG_coordinates(
        fname_vg_coords=fname_vg_coords)
    # Signal-to-Noise-Ratio (SNR) between simulated signals and noise signals
    SNR = 10
    # crop the data as outputs of forward and inverse modelling
    # tmin_epoch, tmax_epoch = -10, 30
    # apply forward and inverse modelling
    [recon_time_series, original_time_series] = forward_inverse_modelling(
        sim_time_series, vg_coordinates, base_folder_path, SNR=SNR, tmin_epoch=tmin_epoch, tmax_epoch=tmax_epoch)

    # # visualise original and reconstructed time-series
    fil_recon_time_series = np.transpose(mne.filter.filter_data(np.transpose(
        recon_time_series), sfreq=MEG_sfreq, l_freq=0.1, h_freq=70, n_jobs=5))
    
    savemat(base_folder_path + 'model_beamform' + sys.argv[1] + '.mat', {'data_beam': fil_recon_time_series})
    plt.figure()
    plt.subplot(3, 1, 1)
    offset_1 = np.max(np.abs(original_time_series))
    for m_row in np.arange(0, original_time_series.shape[1]):
       plt.plot(original_time_series[:, m_row] + offset_1 * m_row)
    plt.title('Ictal spikes assigned to original source space')

    plt.subplot(3, 1, 2)
    offset_2 = np.max(np.abs(recon_time_series))
    for m_row in np.arange(0, recon_time_series.shape[1]):
       plt.plot(recon_time_series[:, m_row] + offset_2 * m_row)
    plt.title('Source reconstructed (beamformers LCMV) ictal spikes from the same locations before applying any filtering')

    plt.subplot(3, 1, 3)
    offset_3 = np.max(np.abs(fil_recon_time_series))
    for m_row in np.arange(0, fil_recon_time_series.shape[1]):
       plt.plot(fil_recon_time_series[:, m_row] + offset_3 * m_row)
    plt.title('Source reconstructed (beamformers LCMV) ictal spikes from the same locations after applying low pass filtering (70Hz)')

    plt.show(block=False)
