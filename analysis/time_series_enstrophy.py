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

setting_name = "setting"
xlim = [0, 200]
ylim = [0.44, 0.48]

reader = ReadClass()
remote_data = reader.data_PID(setting_name)
reader.set_ssh(remote_data=remote_data)
reader.read_setting(setting_name)
reader.read_file_number(ssh_flag=True)

budget = reader.read_budget(number_time=420, ssh_flag=False)
columns = ['T', 'E', 'E-first', 'E-diag', 'E-lower', 'E-upper',
           'Q2', 'Q2-first', 'Q2-diag', 'Q2-lower', 'Q2-upper', 'PE']
df = pd.DataFrame(data=budget, columns=columns)
df = df.set_index('T')

for col in df.columns.tolist():
    df = df.rename(columns={col: reader.long_name(col)})
df['Q2-others'] = df[reader.long_name('Q2')] - df[reader.long_name('Q2-first')]
df['Q2-others-mean'] = df['Q2-others'] / (reader.N - 1)

df.plot(y=[reader.long_name('Q2-first'),
           reader.long_name('E-first'),
           'Q2-others'],
           # reader.long_name('Q2-lower'),
           # reader.long_name('Q2-upper')],
        logy=True,
        # xlim=xlim,
        # ylim=ylim
        figsize=(9, 3)
)
plt.savefig(figure_dir + 'Enstrophy_total' + setting_name + '.eps')

plt.show()
plt.close()
