#
# -*- coding: utf-8 -*-
#
from matplotlib import rc
import matplotlib.pyplot as plt
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
N = (reader.NK_truncate-1) * (reader.NL_truncate-1)

free_energy = np.zeros(aspect.size)
for i in range(aspect.size):
    free_energy[i] = reader.set_free_energy(aspect=aspect[i], beta=beta)
# enstrophy_constant = reader.set_enstrophy_constant(aspect=aspect[0],
#                                                    beta=beta)

exp_E = reader.read_exp_E(beta=beta, number_time=750, ssh_flag=False)
# exp_F = np.exp(- beta * (free_energy[:] - free_energy[0]) )
# e_array = exp_E - exp_F

delta_F = (free_energy[:] - free_energy[0]) / N
delta_F_experiment = - np.log(exp_E) / (beta * N)

fig = plt.figure(figsize=[9, 3])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)
ax.plot(delta_F, c='blue')
ax.plot(delta_F_experiment, 'rx')
# ax.set_ylim(0, 2)
plt.savefig(figure_dir + 'free_energy_' + setting_name + '.eps')

plt.show()
plt.close()
