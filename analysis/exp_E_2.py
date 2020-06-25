#
# -*- coding: utf-8 -*-
#
from matplotlib import rc
import matplotlib.pyplot as plt
# import pandas as pd
import numpy as np

from read_data import ReadClass

plt.style.use('ggplot')
rc('text', usetex=True)
rc('font', family='serif')
figure_dir = './figure/'

setting_name = "setting"

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)
aspect = reader.read_aspect_ratio(ssh_flag=True)
beta = reader.beta

free_energy = np.zeros(aspect.size)
for i in range(aspect.size):
    free_energy[i] = reader.set_free_energy(aspect=aspect[i], beta=beta)
enstrophy_constant = reader.set_enstrophy_constant(aspect=aspect[0], beta=beta)

exp_E = reader.read_exp_E(beta=beta)
exp_F = np.exp(- beta * (free_energy[:] - free_energy[0]) / enstrophy_constant)
e_array = exp_E - exp_F

fig = plt.figure(figsize=[9, 3])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)
ax.plot(exp_E)
ax.plot(exp_F)
ax.set_ylim(0, 2)
plt.savefig(figure_dir + 'exp_Energy' + setting_name + '.eps')

plt.show()
plt.close()
