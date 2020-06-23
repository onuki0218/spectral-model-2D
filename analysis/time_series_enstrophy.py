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

budget = reader.read_budget()
columns = ['T', 'E', 'E-diag', 'E-lower', 'E-upper',
           'Q2', 'Q2-diag', 'Q2-lower', 'Q2-upper', 'PE']
df = pd.DataFrame(data=budget, columns=columns)
df = df.set_index('T')

for col in df.columns.tolist():
    df = df.rename(columns={col: reader.long_name(col)})

df.plot(figsize=(9, 3),
        y=[reader.long_name('Q2')]
           # reader.long_name('Q2-diag'),
           # reader.long_name('Q2-lower'),
           # reader.long_name('Q2-upper')],
        # logy=True,
        # xlim=xlim,
        # ylim=ylim
)
plt.savefig(figure_dir + 'Enstrophy_total' + setting_name + '.eps')

plt.show()
plt.close()
