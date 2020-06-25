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

setting_name = "setting"
figure_dir = './figure/'

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)

K_axis, L_axis = reader.set_axis()

c_level = 64

spec0 = reader.read_Q2_ave(0)

for time in range(1, 11):
    it = time * reader.interval_write_variable
    str_it = '{0:04d}'.format(it)
    spec = reader.read_Q2_ave(it) - spec0

    sum_spec = np.sum(spec)
    #
    print(it, sum_spec / np.sum(spec0))
    # spec_range = np.linspace(np.nanmin(spec), np.nanmax(spec), c_level)
    spec_range = np.linspace(0, np.nanmax(spec), c_level)

    fig = plt.figure(figsize=[4, 3])
    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                        top=0.9, wspace=0.15, hspace=0.15)
    ax = fig.add_subplot(1, 1, 1)
    # spec_color = ax.contourf(K_axis, L_axis, spec_log10, spec_range_log,
    #                          cmap=cm.nipy_spectral, extend="both")
    spec_color = ax.contourf(K_axis, L_axis, spec, spec_range,
                             cmap=cm.nipy_spectral, extend="both")
    spec_color_bar = plt.colorbar(spec_color)
    ax.patch.set_facecolor('black')
    ax.set_title('Averaged enstrophy' + ' at t=' + str_it + 'T')
    #
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_xlim([1, reader.NK_truncate])
     #ax.set_ylim([1, reader.NL_truncate])
    #
    ax.set_xlabel("Zonal wavenumber")
    ax.set_ylabel("Meridional wavenumber")
    fileeps = figure_dir + 'spec_Q2_diff_' + str_it + '.eps'
    plt.savefig(fileeps)
    plt.close()
