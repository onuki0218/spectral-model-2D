#
# -*- coding: utf-8 -*-
#
from matplotlib import rc
import matplotlib.pyplot as plt
import pandas as pd
# import numpy as np

from read_data import ReadClass

plt.style.use('ggplot')
rc('text', usetex=True)
rc('font', family='serif')
figure_dir = './figure/'

setting_name = "setting_3"
xlim = [0, 200]
ylim = [0.0028, 0.003]

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)

budget = reader.read_budget(number_time=118, ssh_flag=True)
columns = ['T', 'E', 'E-diag', 'E-lower', 'E-upper',
           'Q2', 'Q2-diag', 'Q2-lower', 'Q2-upper', 'PE']
df = pd.DataFrame(data=budget, columns=columns)
df = df.set_index('T')

for col in df.columns.tolist():
    df = df.rename(columns={col: reader.long_name(col)})

df.plot(y=[reader.long_name('E')],
           # reader.long_name('E-diag'),
           # reader.long_name('E-lower'),
           # reader.long_name('E-upper')],
        # logy=True,
        # xlim=xlim,
        # ylim=ylim,
        figsize=(9, 3)
)
plt.savefig(figure_dir + 'Energy' + setting_name + '.eps')

plt.show()
plt.close()
