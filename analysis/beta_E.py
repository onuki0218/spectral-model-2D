#
# -*- coding: utf-8 -*-
#
# import sys
from matplotlib import rc
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import numpy as np

rc('text', usetex=True)
rc('font', family='serif')

figure_dir = './figure'

m_max = 128  # Number of zonal modes
n_max = 128  # Number of meridional modes
a = 1.0  # Aspect ratio
r_re = 0.0  # Inverse of Rossby radius
number_beta = 200

N = m_max * n_max

m_array = np.linspace(1, m_max, m_max)
n_array = np.linspace(1, n_max, n_max)
m_mesh, n_mesh = np.meshgrid(m_array, n_array)
lambda_array = (np.pi**2 * (m_mesh**2 / a**2 + n_mesh**2 * a**2) + r_re**2)
# print(lambda_array)

epsilon_array = np.zeros(number_beta)
beta_array = np.zeros(number_beta)
energy_array = np.zeros(number_beta)
enstrophy_array = np.zeros(number_beta)
energy_normal_array = np.zeros(number_beta)

# sys.exit()


def energy_epsilon(epsilon_N):
    epsilon = epsilon_N / N  # Energy of the first mode
    beta = - lambda_array[0, 0] + 1.0 / (2 * epsilon * N)
    energy_sum = 1.0 / (2 * N) * np.sum(1 / (beta + lambda_array[:, :]))
    enstrophy_sum = 1.0 / (2 * N)  \
        * np.sum(lambda_array[:, :] / (beta + lambda_array[:, :]))
    return epsilon, beta, energy_sum, enstrophy_sum


def energy_beta(beta):
    epsilon_N = 1.0 / (2.0 * (beta + lambda_array[0, 0]))
    epsilon = epsilon_N / N  # Energy of the first mode
    tmp = 1 / (beta + lambda_array[:, :])
    energy_sum = 1.0 / (2 * N) * np.sum(tmp)
    enstrophy_sum = 1.0 / (2 * N)  \
        * np.sum(lambda_array[:, :] / (beta + lambda_array[:, :]))
    return epsilon, beta, energy_sum, enstrophy_sum


for i in range(number_beta):
    # epsilon_N = (i+1) * 100.0
    # epsilon, beta, energy_sum, enstrophy_sum = energy_epsilon(epsilon_N)
    beta = - lambda_array[0, 0] + np.exp(-(i+1) * 0.1)
    epsilon, beta, energy_sum, enstrophy_sum = energy_beta(beta)
    epsilon_array[i] = epsilon
    beta_array[i] = beta
    energy_array[i] = energy_sum
    enstrophy_array[i] = enstrophy_sum
    energy_normal_array[i] = energy_sum / enstrophy_sum
    print(i, beta, energy_normal_array[i], epsilon / energy_sum)

print(lambda_array[0, 0], lambda_array[0, 1], 1.0 / lambda_array[0, 0])

fig = plt.figure(figsize=[4, 3])
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9,
                    top=0.9, wspace=0.15, hspace=0.15)
ax = fig.add_subplot(1, 1, 1)

ax.plot(beta_array, energy_normal_array)
ax.set_xlabel("beta")
ax.set_ylabel("energy")
fileeps = figure_dir + 'beta_E.eps'
plt.savefig(fileeps)
plt.show()
plt.close()
