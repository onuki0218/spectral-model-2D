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
beta = -19.56

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)
aspect = reader.read_aspect_ratio(ssh_flag=True)

N = (reader.NL_truncate-1) * (reader.NK_truncate-1)

energy = np.zeros(aspect.size)
Q2 = np.zeros([reader.NL_truncate, reader.NK_truncate, reader.N_process])
for j in range(reader.N_process):
    Q2[:, :, j] = reader.read_real(it=0, ip=j, var_name="Q")**2
    print('j =', j)

# t = 0
mu = reader.set_mu(aspect=aspect[0])
E0 = np.zeros([reader.NL_truncate-1, reader.NK_truncate-1, reader.N_process])
for j in range(reader.N_process):
    E0[:, :, j] = Q2[:-1, :-1, j] / mu[:, :]

E = np.zeros([reader.NL_truncate-1, reader.NK_truncate-1])
work = np.zeros(reader.N_process)
free_energy = np.zeros(aspect.size)
exp_E = np.zeros(aspect.size)
for i in range(aspect.size):
    free_energy[i] = reader.set_free_energy(aspect=aspect[i], beta=beta)
    mu = reader.set_mu(aspect=aspect[i])
    for j in range(reader.N_process):
        E[:, :] = Q2[:-1, :-1, j] / mu[:, :]
        work[j] = np.sum(E[:, :] - E0[:, :, j])
    exp_E[i] = np.mean(np.exp(-N * beta * work[:]))
    print('i =', i)

enstrophy_constant = reader.set_enstrophy_constant(aspect=aspect[0], beta=beta)

# exp_E = reader.read_exp_E(beta=beta)
exp_F = np.exp(- beta * (free_energy[:] - free_energy[0]) / enstrophy_constant)
e_array = exp_E - exp_F

fig = plt.figure(figsize=[9, 3])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)
ax.plot(exp_E)
ax.plot(exp_F)
# ax.set_ylim(0, 2)
plt.savefig(figure_dir + 'exp_Energy_RDT' + setting_name + '.eps')

plt.show()
plt.close()
