#
# -*- coding: utf-8 -*-
#
from os import remove
from paramiko import SSHClient, AutoAddPolicy

import f90nml
import numpy as np


class ReadClass:
    def __init__(self, local_dir='./data/'):
        self.local_dir = local_dir

    def __del__(self):
        pass

    def set_flag(self, remove_flag):
        self.remove_flag = remove_flag

    def set_ssh(
            self,
            host='ofp.jcahpc.jp',
            user='c24070',
            private_key='/Users/onuki/.ssh/id_rsa.pub',
            port=22,
            remote_dir='/work/gi55/c24070/spectral-model-2D/',
            # remote_data='data/data_default'
            remote_data='data/data2899425/'
    ):
        self.host = host
        self.user = user
        self.private_key = private_key
        self.port = port
        self.remote_dir = remote_dir
        self.remote_data = remote_data
        self._ssh_method()

    def _ssh_method(self):
        self.client = SSHClient()
        self.client.set_missing_host_key_policy(AutoAddPolicy())
        self.client.connect(
            hostname=self.host, port=self.port,
            username=self.user, key_filename=self.private_key)
        self.sftp_client = self.client.open_sftp()

    def _read_data(self, file_name, data_type, shape,
                   ssh_flag=True, remove_flag=True):
        if (ssh_flag):
            self.sftp_client.get(
                self.remote_dir + self.remote_data + file_name,
                self.local_dir + file_name)
        size = shape.prod()
        file_open = open(self.local_dir + file_name, 'r')
        data_tmp = np.fromfile(file_open, dtype=data_type, count=size)
        data = data_tmp.reshape(shape, order='F')
        if (ssh_flag and remove_flag):
            remove(self.local_dir + file_name)
        return data

    def read_setting(self, setting_name, ssh_flag=True):
        if (ssh_flag):
            self.sftp_client.get(
                self.remote_dir + self.remote_data + setting_name + '.out',
                self.local_dir + setting_name + '.out')
        namelist = f90nml.read(self.local_dir + setting_name + '.out')
        tmp = namelist[setting_name]

        self.NK = tmp['NK']
        self.NL = tmp['NL']
        self.NK_truncate = tmp['NK_truncate']
        self.N_process = tmp['N_process']
        # self.N_output_file = tmp['N_output_file']
        self.length_domain_X = tmp['length_domain_X']
        self.length_domain_Y = tmp['length_domain_Y']
        self.time_end = tmp['time_end']
        # self.CFL_factor = tmp['CFL_factor']
        self.time_write = tmp['time_write']
        self.interval_write_budget = tmp['interval_write_budget']
        self.interval_write_variable = tmp['interval_write_variable']
        # self.viscosity = tmp['viscosity']

        self.NL_truncate = self.NK_truncate  #

        self.Del_K = np.pi / self.length_domain_X
        self.Del_L = np.pi / self.length_domain_Y
        self.wavenumber_truncate = self.NK_truncate * self.Del_K
        #

    def read_file_number(self, ssh_flag=True):
        if (ssh_flag):
            self.sftp_client.get(
                self.remote_dir + self.remote_data + 'file_number.out',
                self.local_dir + 'file_number.out')
        file_open = open(self.local_dir + 'file_number.out', 'r')
        self.file_number_list = list([])
        line = file_open.readline()
        while line:
            if (not str.isdecimal(str.strip(line))):
                self.file_number_list.append([])
            else:
                self.file_number_list[-1].append(int(line))
            line = file_open.readline()
        file_open.close()

    def read_budget(self, ssh_flag=True):
        shape = np.array([10, self.file_number_list[-1][-1] + 1])
        data_type = np.dtype('<f4')
        data = np.zeros(shape)
        for ip in range(self.N_process):
            str_ip = '{0:04d}'.format(ip)
            file_name = 'budget' + str_ip + '.out'
            data_tmp = self._read_data(file_name, data_type, shape,
                                       ssh_flag=ssh_flag, remove_flag=False)
            data = data + data_tmp
            print(str_ip)
        data = data / self.N_process
        return data.transpose()

    def read_exp_E(self, beta, ssh_flag=True):
        number_time = self.file_number_list[-1][-1] + 1
        shape = np.array([10, number_time])
        data_type = np.dtype('<f4')
        data = np.zeros(number_time)
        for ip in range(self.N_process):
            str_ip = '{0:04d}'.format(ip)
            file_name = 'budget' + str_ip + '.out'
            data_tmp = self._read_data(file_name, data_type, shape,
                                       ssh_flag=ssh_flag, remove_flag=False)
            data = data + np.exp(-beta*(data_tmp[1, :] - data_tmp[1, 0])
                                 * (self.NK_truncate-1) * (self.NL_truncate-1))
            print(str_ip)
        data = data / self.N_process
        return data

    def read_aspect_ratio(self, ssh_flag=True):
        shape = np.array([self.file_number_list[-1][-1] + 1])
        data_type = np.dtype('<f4')
        data = np.zeros(shape)
        file_name = 'aspect_ratio.out'
        data = self._read_data(file_name, data_type, shape,
                               ssh_flag=ssh_flag, remove_flag=False)
        return data

    def read_real(self, it, ip, var_name, flag=True, ssh_flag=True):
        NL_in = self.NL_truncate
        NK_in = self.NL_truncate

        # Prepare arrays
        shape = np.array([NL_in, NK_in])
        data_type = np.dtype('<f4')

        str_it = '{0:06d}'.format(it)
        str_ip = '{0:04d}'.format(ip)
        print(str_it, str_ip)
        file_name = var_name + str_ip + '_' + str_it + '.out'
        self.R_out = self._read_data(file_name, data_type, shape,
                                     ssh_flag=ssh_flag)
        return self.R_out

    def read_E_ave(self, it, flag=True, ssh_flag=True):
        NL_in = self.NL_truncate
        NK_in = self.NL_truncate

        # Prepare arrays
        shape = np.array([NL_in, NK_in])
        data_type = np.dtype('<f4')

        str_it = '{0:06d}'.format(it)
        print(str_it)
        file_name = 'E' + str_it + '.out'
        self.R_out = self._read_data(file_name, data_type, shape,
                                     ssh_flag=ssh_flag)
        return self.R_out

    def read_Q2_ave(self, it, flag=True, ssh_flag=True):
        NL_in = self.NL_truncate
        NK_in = self.NL_truncate

        # Prepare arrays
        shape = np.array([NL_in, NK_in])
        data_type = np.dtype('<f4')

        str_it = '{0:06d}'.format(it)
        print(str_it)
        file_name = 'Q2' + str_it + '.out'
        self.R_out = self._read_data(file_name, data_type, shape,
                                     ssh_flag=ssh_flag)
        return self.R_out

    def set_axis(self):
        K_axis = np.arange(1, self.NK_truncate+1) * self.Del_K
        L_axis = np.arange(1, self.NL_truncate+1) * self.Del_L
        return K_axis, L_axis

    def set_XY_axis(self):
        LX = self.length_domain_X
        LY = self.length_domain_Y

        X = np.linspace(-LX/2, LX/2, self.NK, dtype='float64', endpoint=False)
        Y = np.linspace(-LY/2, LY/2, self.NL, dtype='float64', endpoint=False)
        # X_mesh, Y_mesh = np.meshgrid(X, Y, indexing='ij')

        return X, Y

    def set_mu(self, aspect):
        m_array = np.linspace(1, self.NK_truncate-1, self.NK_truncate-1)
        n_array = np.linspace(1, self.NL_truncate-1, self.NL_truncate-1)
        m_mesh, n_mesh = np.meshgrid(m_array, n_array)
        mu = np.pi**2 * (m_mesh**2 / aspect + n_mesh**2 * aspect)
        return mu

    def set_free_energy(self, aspect, beta):
        mu = self.set_mu(aspect=aspect)
        free_energy = np.sum(np.log(beta + mu[:, :]) - np.log(mu[:, :]))  \
            / (2 * beta)
        return free_energy

    def set_enstrophy_constant(self, aspect, beta):
        mu = self.set_mu(aspect=aspect)
        enstrophy_constant = np.sum(mu[:, :] / (beta + mu[:, :]))  \
            / (2 * self.NK_truncate * self.NL_truncate)
        return enstrophy_constant

    def set_dissipation_axis(self):
        Del_dis = (self.dissipation_max - self.dissipation_min) / \
            (self.N_dissipation - 2)
        dis_axis = np.linspace(
            self.dissipation_min - Del_dis / 2,
            self.dissipation_max + Del_dis / 2,
            self.N_dissipation, endpoint=True)
        return dis_axis

    def long_name(self, var):
        dictionary = {'T': 'Time',
                      'U': 'Zonal velocity',
                      'V': 'Meridional velocity',
                      'E': 'Energy',
                      'E-diag': 'Energy-diagonal',
                      'E-lower': 'Energy-lower',
                      'E-upper': 'Energy-upper',
                      'Q2': 'Enstrophy',
                      'Q2-diag': 'Enstrophy-diagonal',
                      'Q2-lower': 'Enstrophy-lower',
                      'Q2-upper': 'Enstrophy-upper',
                      'NE': 'Energy transfer rate',
                      'PE': 'Energy production rate',
                      'DE': 'Energy dissipation rate'}
        return dictionary[var]

    def data_PID(self, setting_name):
        dictionary = {'setting': 'data/data_default/',
                      'setting_2': 'data/data_sub_2/'}
        return dictionary[setting_name]
