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

work = reader.read_work(beta=beta, number_time=201, ssh_flag=False)
work_mean = work.mean(axis=1)

time = np.zeros(work.shape)
process = np.arange(work.shape[1])
for i in range(time.shape[0]):
    time[i, :] = i * reader.time_write

it = 40
fig = plt.figure(figsize=[15, 6])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)
# ax.scatter(process, work[it, :], s=3)
# ax.set_ylim(-0.0002, 0.0002)
# plt.savefig(figure_dir + 'work_scatter_'
#             + str(it) + '_' + setting_name + '.eps')
ax.scatter(time, work, s=3)
ax.plot(time[:, 0], work_mean, 'k-')
ax.set_ylim(np.min(work) * 1.2, np.max(work) * 1.2)
plt.savefig(figure_dir + 'work_scatter_' + setting_name + '.eps')

plt.show()
plt.close()
