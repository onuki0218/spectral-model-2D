#
# -*- coding: utf-8 -*-
#
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import fftpack

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

X, Y = reader.set_XY_axis()

c_level = 64
var_dict = {'E': 'Energy spectrum',
            'P': 'Enstrophy spectrum',
            'Q': 'Vorticity'}
var = 'Q'

ip = 1
str_ip = '{0:04d}'.format(ip)
Q0 = reader.read_real(0, ip, 'Q')
Q = np.zeros((reader.NK, reader.NL))
Q[:reader.NK_truncate-1, :reader.NK_truncate-1] = Q0[:, :]
Q_real = fftpack.dstn(Q, type=1, norm='ortho')
max_Q = Q_real.max()
min_Q = Q_real.min()

for time in range(0, 11):
    it = time * reader.interval_write_variable
    str_it = '{0:04d}'.format(it)
    Q0 = reader.read_real(it, 0, var)
    Q = np.zeros((reader.NK, reader.NL))
    Q[:reader.NK_truncate-1, :reader.NK_truncate-1] = Q0[:, :]

    enstrophy_K = np.sum(Q**2) / 2
    Q_real = fftpack.dstn(Q, type=1, norm='ortho')
    enstrophy_real = np.sum(Q_real**2) / 2
    #
    print(it, enstrophy_K, enstrophy_real)
    c_range = np.linspace(min_Q, max_Q, c_level)

    fig = plt.figure(figsize=[4, 3])
    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                        top=0.9, wspace=0.15, hspace=0.15)
    ax = fig.add_subplot(1, 1, 1)
    color = ax.contourf(X, Y, Q_real, c_range, extend="both")
    color_bar = plt.colorbar(color)
    # ax.patch.set_facecolor('black')
    ax.set_title(var_dict[var] + ' at t=' + str_it + 'T')
    #
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_xlim([1, reader.NK])
    # ax.set_ylim([1, reader.NL])
    #
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    fileeps = figure_dir + var + str_ip + '_' + str_it + '.eps'
    plt.savefig(fileeps)
    # plt.show()
    plt.close()
