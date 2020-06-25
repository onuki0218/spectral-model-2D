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

for time in range(0, 11):
    it = time * reader.interval_write_variable
    str_it = '{0:04d}'.format(it)
    spec = reader.read_E_ave(it)

    energy_all = np.sum(spec)
    spec[spec <= 0.0] = np.nan
    spec_log10 = np.log10(spec / energy_all)
    max_spec = np.nanmax(spec_log10)
    min_spec = np.nanmin(spec_log10)
    #
    print(it, max_spec, min_spec, energy_all)
    # print(spec_XY_Z_log10)
    # spec_range = np.linspace(min_spec, max_spec, c_level)
    spec_range = np.linspace(-5.0, 0.0, c_level)

    fig = plt.figure(figsize=[4, 3])
    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                        top=0.9, wspace=0.15, hspace=0.15)
    ax = fig.add_subplot(1, 1, 1)
    spec_color = ax.contourf(K_axis, L_axis, spec_log10, spec_range,
                             cmap=cm.nipy_spectral, extend="both")
    spec_color_bar = plt.colorbar(spec_color)
    ax.patch.set_facecolor('black')
    ax.set_title('Averaged energy' + ' at t=' + str_it + 'T')
    #
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1, reader.NK_truncate])
    ax.set_ylim([1, reader.NL_truncate])
    #
    ax.set_xlabel("Zonal wavenumber")
    ax.set_ylabel("Meridional wavenumber")
    fileeps = figure_dir + 'spec_E' + str_it + '.eps'
    plt.savefig(fileeps)
    plt.close()
