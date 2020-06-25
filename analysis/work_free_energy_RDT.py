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

setting_name = "setting_3"
# xlim = [0, 200]
# ylim = [0.0028, 0.003]

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)
aspect = reader.read_aspect_ratio(ssh_flag=True)
beta = reader.beta

work = reader.read_work(beta=beta, number_time=125, ssh_flag=False)
work_mean = work.mean(axis=1)

free_energy = np.zeros(aspect.size)
RDT_work = np.zeros(aspect.size)
for i in range(aspect.size):
    free_energy[i] = reader.set_free_energy(aspect=aspect[i], beta=beta)
    RDT_work[i] = reader.set_RDT_work(aspect_0=aspect[0], aspect=aspect[i],
                                      beta=beta)
delta_F = (free_energy[:] - free_energy[0]) / reader.N

fig = plt.figure(figsize=[9, 3])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)
ax.plot(work_mean, 'ro')
ax.plot(delta_F, c='blue')
ax.plot(RDT_work, c='black')
plt.savefig(figure_dir + 'work_free_energy_' + setting_name + '.eps')

plt.show()
plt.close()
