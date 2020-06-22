#
# -*- coding: utf-8 -*-
#
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

from read_data import ReadClass

rc('text', usetex=True)
rc('font', family='serif')

Fr_str = '04'
# setting_name = "setting_Fr" + Fr_str
setting_name = "setting_restart_Fr" + Fr_str
# setting_name = "setting_restart"
figure_dir = './figure_Fr' + Fr_str + '/'
# setting_name = "setting"

# setting_name = "setting_Dossmann_Pr1_Re40000_restart"
# figure_dir = './figure_Dossmann_Pr1_Re40000/'
# setting_name = "setting_Dossmann_Pr1_Re40000"
# figure_dir = './figure_Dossmann_Pr1_Re40000/'

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)

K_axis = np.linspace(0.0, reader.N_truncate, reader.N_truncate+1)
M_axis = np.linspace(0.0, reader.N_truncate, reader.N_truncate+1)
K, M = np.meshgrid(K_axis, M_axis)

# line of sigma = omega / 2
M_dispersion = K_axis[:]  \
    * np.sqrt(reader.buoyancy_frequency**2 - reader.wave_frequency**2)  \
    / reader.wave_frequency
M_dispersion_2 = K_axis[:]  \
    * np.sqrt(4 * reader.buoyancy_frequency**2 - reader.wave_frequency**2)  \
    / reader.wave_frequency

c_level = 64
var_dict = {'E': 'Energy spectrum',
            'KE': 'Kinetic energy spectrum',
            'PE': 'Available potential energy spectrum'}
var = 'E'

for time in range(16, 26):
    it = time * 16
    str_it = '{0:04d}'.format(it)
    # record = int(time * reader.omg / (2 * np.pi) * 1.0001)
    # str_record = str(record)
    # str_record_pad = '{0:03d}'.format(record)
    spec_XY_Z = reader.read_spec_XY_Z(it, var)

    energy_all = np.sum(spec_XY_Z)
    spec_XY_Z[spec_XY_Z <= 0.0] = np.nan
    spec_XY_Z_log10 = np.log10(spec_XY_Z / energy_all)
    max_spec = np.nanmax(spec_XY_Z_log10)
    min_spec = np.nanmin(spec_XY_Z_log10)
    #
    print(it, max_spec, min_spec, energy_all)
    # print(spec_XY_Z_log10)
    # spec_range = np.linspace(min_spec, max_spec, c_level)
    spec_range = np.linspace(-9.0, 0.0, c_level)

    fig = plt.figure(figsize=[4, 3])
    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                        top=0.9, wspace=0.15, hspace=0.15)
    ax = fig.add_subplot(1, 1, 1)
    spec_color = ax.contourf(
        K, M, spec_XY_Z_log10, spec_range,
        cmap=cm.nipy_spectral, extend="both")
    spec_color_bar = plt.colorbar(spec_color)
    ax.plot(K_axis, M_dispersion, c='k')
    ax.plot(K_axis, M_dispersion_2, c='k')
    ax.patch.set_facecolor('black')
    ax.set_title(var_dict[var] + ' at t=' + str_it + 'T')
    #
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1, reader.N_truncate/2])
    ax.set_ylim([1, reader.N_truncate/2])
    #
    ax.set_xlabel("Horizontal wavenumber")
    ax.set_ylabel("Vertical wavenumber")
    # fileeps = figure_dir + 'spec_XY_Z_' + var + str_it + '.eps'
    fileeps = figure_dir + 'spec_XY_Z_log_' + var + str_it + '.eps'
    plt.savefig(fileeps)
    plt.close()
