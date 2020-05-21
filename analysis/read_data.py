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
            remote_dir='/work/gi55/c24070/spectral-model-distortion/',
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
        self.NM = tmp['NM']
        self.N_truncate = tmp['N_truncate']
        self.N_process = tmp['N_process']
        self.NK_local = tmp['NK_local']
        self.N_output_file = tmp['N_output_file']
        self.NK_local_output = tmp['NK_local_output']
        self.length_domain_X = tmp['length_domain_X']
        self.length_domain_Y = tmp['length_domain_Y']
        self.length_domain_Z = tmp['length_domain_Z']
        self.time_end = tmp['time_end']
        self.CFL_factor = tmp['CFL_factor']
        self.time_write = tmp['time_write']
        self.interval_write_budget = tmp['interval_write_budget']
        self.interval_write_variable = tmp['interval_write_variable']
        self.interval_write_B = tmp['interval_write_B']
        self.Coriolis_frequency = tmp['Coriolis_frequency']
        self.buoyancy_frequency = tmp['buoyancy_frequency']
        self.wave_frequency = tmp['wave_frequency']
        self.viscosity = tmp['viscosity']
        self.diffusivity = tmp['diffusivity']
        self.noise_level = tmp['noise_level']
        self.background_shear = tmp['background_shear']
        self.interval_write_spectra = tmp['interval_write_spectra']
        self.N_angle = tmp['N_angle']
        self.interval_write_Richardson = tmp['interval_write_Richardson']
        self.N_Richardson = tmp['N_Richardson']
        self.Richardson_min = tmp['Richardson_min']
        self.Richardson_max = tmp['Richardson_max']
        self.N_dissipation = tmp['N_dissipation']
        self.dissipation_min = tmp['dissipation_min']
        self.dissipation_max = tmp['dissipation_max']
        self.interval_write_profile = tmp['interval_write_profile']
        self.N_profile = tmp['N_profile']

#
        self.NK_half = int(self.NK / 2)
        self.NL_half = int(self.NL / 2)
        self.NM_half = int(self.NM / 2)
#        self.NK_truncate = int((2.0 ** 0.5 / 3 * self.NK))
#        self.NL_truncate = int((2.0 ** 0.5 / 3 * self.NL))
#        self.NM_truncate = int((2.0 ** 0.5 / 3 * self.NM))
        self.NK_truncate = self.N_truncate
        self.NL_truncate = self.N_truncate
        self.NM_truncate = self.N_truncate

        self.NK_effective = self.NK_truncate * 2 + 1
        self.NL_effective = self.NL_truncate * 2 + 1
        self.NM_effective = self.NM_truncate + 1
        #
        self.sine_theta = np.sqrt(
            (self.wave_frequency ** 2 - self.Coriolis_frequency ** 2)
            / (self.buoyancy_frequency ** 2 - self.Coriolis_frequency ** 2))
        self.cosine_theta = np.sqrt(1.0 - self.sine_theta ** 2)
        self.Coriolis_frequency_Z = self.Coriolis_frequency * self.cosine_theta
        self.Del_K = 2 * np.pi / self.length_domain_X
        self.Del_L = 2 * np.pi / self.length_domain_Y
        self.Del_M = 2 * np.pi / self.length_domain_Z
        self.strain_XZ = self.background_shear / self.wave_frequency
        self.wavenumber_truncate = self.N_truncate * self.Del_K
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
        file_name = 'budget.out'
        shape = np.array([10, self.file_number_list[-1][-1] + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape,
                               ssh_flag=ssh_flag, remove_flag=False)
        return data.transpose()

    def read_real(self, IT, var_name,
                  flag=True, NM_out=0, NL_out=0, NK_out=0, ssh_flag=True):
        if (flag):
            NM_out = self.NM_truncate + 1
            NL_out = self.NL_truncate * 2 + 1
            NK_out = self.NK_truncate * 2 + 1
        NM_in = self.NM_truncate + 1
        NL_in = self.NL_truncate * 2 + 1
        NK_in = self.NK_local_output
        NK_in_ALL = NK_in * self.N_output_file
        NK_in_OFF = self.NK_truncate * 2 + 1

        # Prepare arrays
        self.R_out = np.zeros((NM_out, NL_out, NK_out),
                              dtype='float32', order='F')
        R_Eff = np.zeros((NM_in, NL_in, NK_in_ALL),
                         dtype='float32', order='F')
        shape = np.array([NM_in, NL_in, NK_in])
        data_type = np.dtype('<f4')

        str_it = '{0:06d}'.format(IT)
        for ip in range(self.N_output_file):
            str_ip = '{0:04d}'.format(ip)
            print(str_ip)
            file_name = var_name + str_ip + '_' + str_it + '.out'
            R_Div = self._read_data(file_name, data_type, shape,
                                    ssh_flag=ssh_flag)
            ikp = ip * NK_in
            R_Eff[:, :, ikp: ikp + NK_in] = R_Div[:, :, :]

            # Make output array
            NL_in_H = int(NL_in / 2 + 1)
            NK_in_H = int(NK_in_OFF / 2 + 1)
            NL_out_H = NL_out - (NL_in - NL_in_H)
            NK_out_H = NK_out - (NK_in_OFF - NK_in_H)
            self.R_out[0: NM_in, 0: NL_in_H, 0: NK_in_H]  \
                = R_Eff[0: NM_in, 0: NL_in_H, 0: NK_in_H]
            self.R_out[0: NM_in, 0: NL_in_H, NK_out_H: NK_out]  \
                = R_Eff[0: NM_in, 0: NL_in_H, NK_in_H: NK_in_OFF]
            self.R_out[0: NM_in, NL_out_H: NL_out, 0: NK_in_H]  \
                = R_Eff[0: NM_in, NL_in_H: NL_in,  0: NK_in_H]
            self.R_out[0: NM_in, NL_out_H: NL_out, NK_out_H: NK_out]  \
                = R_Eff[0: NM_in, NL_in_H: NL_in,  NK_in_H: NK_in_OFF]

            self.NK_offset_1 = NK_out_H
            self.NL_offset_1 = NL_out_H
            self.NK_offset_2 = NK_out - NK_out_H
            self.NL_offset_2 = NL_out - NL_out_H
        return self.R_out

    def read_complex(self, it, var_name, flag=True,
                     NM_out=0, NL_out=0, NK_out=0, ssh_flag=True):
        if (flag):
            NM_out = self.NM_truncate + 1
            NL_out = self.NL_truncate * 2 + 1
            NK_out = self.NK_truncate * 2 + 1
        NM_in = self.NM_truncate + 1
        NL_in = self.NL_truncate * 2 + 1
        NK_in = self.NK_local_output
        NK_in_ALL = NK_in * self.N_output_file
        NK_in_OFF = self.NK_truncate * 2 + 1

        # Prepare arrays
        self.C_out = np.zeros((NM_out, NL_out, NK_out),
                              dtype='complex64', order='F')
        C_Eff = np.zeros((NM_in, NL_in, NK_in_ALL),
                         dtype='complex64', order='F')
        shape = np.array([NM_in, NL_in, NK_in])
        data_type = np.dtype('<c8')

        str_it = '{0:06d}'.format(it)
        for ip in range(self.N_output_file):
            print('it = ' + str(it) + ', ip = ' + str(ip))
            str_ip = '{0:04d}'.format(ip)
            file_name = var_name + str_ip + '_' + str_it + '.out'
            C_Div = self._read_data(file_name, data_type, shape,
                                    ssh_flag=ssh_flag)
            ikp = ip * NK_in
            C_Eff[:, :, ikp: ikp + NK_in] = C_Div[:, :, :]

        # Make output array
        NL_in_H = int(NL_in / 2 + 1)
        NK_in_H = int(NK_in_OFF / 2 + 1)
        NL_out_H = NL_out - (NL_in - NL_in_H)
        NK_out_H = NK_out - (NK_in_OFF - NK_in_H)
        self.C_out[0: NM_in, 0: NL_in_H, 0: NK_in_H]  \
            = C_Eff[0: NM_in, 0: NL_in_H, 0: NK_in_H]
        self.C_out[0: NM_in, 0: NL_in_H, NK_out_H: NK_out]  \
            = C_Eff[0: NM_in, 0: NL_in_H, NK_in_H: NK_in_OFF]
        self.C_out[0: NM_in, NL_out_H: NL_out, 0: NK_in_H]  \
            = C_Eff[0: NM_in, NL_in_H: NL_in,  0: NK_in_H]
        self.C_out[0: NM_in, NL_out_H: NL_out, NK_out_H: NK_out]  \
            = C_Eff[0: NM_in, NL_in_H: NL_in,  NK_in_H: NK_in_OFF]

        self.NK_offset_1 = NK_out_H
        self.NL_offset_1 = NL_out_H
        self.NK_offset_2 = NK_out - NK_out_H
        self.NL_offset_2 = NL_out - NL_out_H
        return self.C_out

    def read_spec_angle(self, it, var='E', ssh_flag=True):
        num_lat = self.N_angle + 1
        num_long = 2 * self.N_angle
        shape = np.array([num_lat, num_long])
        data_type = np.dtype('<f4')
        str_it = '{0:06d}'.format(it)
        if var == 'E':
            file_name = 'KE_SpecAngle_' + str_it + '.out'
            data_KE = self._read_data(file_name, data_type, shape,
                                      ssh_flag=ssh_flag)
            file_name = 'PE_SpecAngle_' + str_it + '.out'
            data_PE = self._read_data(file_name, data_type, shape,
                                      ssh_flag=ssh_flag)
            data = data_KE + data_PE
        elif var in ('KE', 'PE'):
            file_name = var + '_SpecAngle_' + str_it + '.out'
            data = self._read_data(file_name, data_type, shape,
                                   ssh_flag=ssh_flag)
        else:
            print("Assign var either 'E', 'KE', or 'PE'.")
            return
        return data

    def read_spec_XY_Z(self, it, var='E', ssh_flag=True):
        num_horizontal = self.N_truncate + 1
        num_vertical = self.N_truncate + 1
        shape = np.array([num_vertical, num_horizontal])
        data_type = np.dtype('<f4')
        str_it = '{0:06d}'.format(it)
        if var == 'E':
            file_name = 'KE_SpecXY-Z_' + str_it + '.out'
            data_KE = self._read_data(file_name, data_type, shape,
                                      ssh_flag=ssh_flag)
            file_name = 'PE_SpecXY-Z_' + str_it + '.out'
            data_PE = self._read_data(file_name, data_type, shape,
                                      ssh_flag=ssh_flag)
            data = data_KE + data_PE
        elif var in ('KE', 'PE'):
            file_name = var + '_SpecXY-Z_' + str_it + '.out'
            data = self._read_data(file_name, data_type, shape,
                                   ssh_flag=ssh_flag)
        else:
            print("Assign var either 'E', 'KE', or 'PE'.")
            return
        return data

    def read_profile(self, it, var, ssh_flag=True):
        shape = np.array([self.NM, self.N_profile])
        data_type = np.dtype('<f4')
        str_it = '{0:06d}'.format(it)
        file_name = var + '_profile_' + str_it + '.out'
        data = self._read_data(file_name, data_type, shape,
                               ssh_flag=ssh_flag)
        return data

    def set_axis(self):
        N_tmp = self.NK_truncate * 2 + 1
        K_axis = np.fft.fftfreq(
            N_tmp, self.length_domain_X / N_tmp) * 2 * np.pi
        K_axis = np.fft.fftshift(K_axis)

        N_tmp = self.NL_truncate * 2 + 1
        L_axis = np.fft.fftfreq(
            N_tmp, self.length_domain_Y / N_tmp) * 2 * np.pi
        L_axis = np.fft.fftshift(K_axis)

        N_tmp = self.NM_truncate * 2 + 1
        M_axis = np.fft.fftfreq(
            N_tmp, self.length_domain_Z / N_tmp) * 2 * np.pi
        M_axis = M_axis[0:self.N_truncate+1]
        return K_axis, L_axis, M_axis

    def set_ZG_axis(self):
        length_ZG = self.length_domain_Z / self.cosine_theta
        ZG_axis = np.linspace(0, length_ZG, self.NM, endpoint=False)
        return ZG_axis

    def set_XYZ_G_axis(self, it):
        LX = self.length_domain_X
        LY = self.length_domain_Y
        LZ = self.length_domain_Z

        X = np.linspace(-LX/2, LX/2, self.NK, dtype='float64', endpoint=False)
        Y = np.linspace(-LY/2, LY/2, self.NL, dtype='float64', endpoint=False)
        Z = np.linspace(-LZ/2, LZ/2, self.NM, dtype='float64', endpoint=False)
        X_mesh, Y_mesh, Z_mesh = np.meshgrid(X, Y, Z, indexing='ij')

        strain_XZ = self.cal_strain(it, 0, True)
        strain_YZ = self.cal_strain(it, 0, False)
        X_tmp = X_mesh + strain_XZ * Z_mesh
        YG = Y_mesh + strain_YZ * Z_mesh
        XG = X_tmp * self.cosine_theta + Z_mesh * self.sine_theta
        ZG = - X_tmp * self.sine_theta + Z_mesh * self.cosine_theta
        return XG, YG, ZG

    def cal_strain(self, it, phase_offset=0, flag=True):
        time = it * self.time_write
        if flag:
            strain = self.background_shear / self.wave_frequency * \
                np.sin(self.wave_frequency * time + phase_offset)
        else:
            strain = self.Coriolis_frequency_Z * self.background_shear / \
                self.wave_frequency**2 * \
                np.cos(self.wave_frequency * time + phase_offset)
        return strain

    def read_Richardson(self, it, ssh_flag=True):
        shape = np.array([self.N_Richardson])
        data_type = np.dtype('<f4')
        str_it = '{0:06d}'.format(it)
        file_name = 'Richardson' + str_it + '.out'
        data = self._read_data(file_name, data_type, shape,
                               ssh_flag=ssh_flag)
        return data

    def set_Richardson_axis(self):
        Del_Ri = (self.Richardson_max - self.Richardson_min) / \
            (self.N_Richardson - 2)
        Ri_axis = np.linspace(
            self.Richardson_min - Del_Ri / 2, self.Richardson_max + Del_Ri / 2,
            self.N_Richardson, endpoint=True)
        return Ri_axis

    def read_Richardson_dissipation(self, it, ssh_flag=True):
        shape = np.array([self.N_Richardson, self.N_dissipation])
        data_type = np.dtype('<f4')
        str_it = '{0:06d}'.format(it)
        file_name = 'Richardson_dissipation' + str_it + '.out'
        data = self._read_data(file_name, data_type, shape,
                               ssh_flag=ssh_flag)
        return data

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
                      'W': 'Vertical velocity',
                      'B': 'Buoyancy perturbation',
                      'R': 'Density',
                      'KE': 'Kinetic energy',
                      'PE': 'Potential energy',
                      'TE': 'Total energy',
                      'MC': 'Mixing coefficient',
                      'NKE': 'Kinetic energy transfer rate',
                      'NPE': 'Potential energy transfer rate',
                      'PKE': 'Kinetic energy production rate',
                      'PPE': 'Potential energy production rate',
                      'CKP': 'Energy conversion rate',
                      'DKE': 'Kinetic energy dissipation rate',
                      'DPE': 'Potential energy dissipation rate'}
        return dictionary[var]

    def data_PID(self, setting_name):
        dictionary = {'setting_restart_Fr01': 'data/data2898419/',
                      'setting_restart': 'data/data2892514/',
                      'setting_restart_Fr03': 'data/data2898421/',
                      'setting_restart_Fr04': 'data/data2898422/',
                      'setting_restart_Fr05': 'data/data2898316/',
                      'setting_restart_Fr06': 'data/data2898423/',
                      'setting_restart_Fr07': 'data/data2898424/',
                      'setting_restart_Fr08': 'data/data2898425/',
                      'setting_restart_Fr09': 'data/data2898426/',
                      'setting_restart_Fr10': 'data/data2898317/',
                      'setting_restart_Fr11': 'data/data2899423/',
                      'setting_restart_Fr12': 'data/data2899425/',
                      'setting_Fr01': 'data/data2898419/',
                      'setting': 'data/data2892514/',
                      'setting_Fr03': 'data/data2898421/',
                      'setting_Fr04': 'data/data2898422/',
                      'setting_Fr05': 'data/data2898316/',
                      'setting_Fr06': 'data/data2898423/',
                      'setting_Fr07': 'data/data2898424/',
                      'setting_Fr08': 'data/data2898425/',
                      'setting_Fr09': 'data/data2898426/',
                      'setting_Fr10': 'data/data2898317/',
                      'setting_Fr11': 'data/data2899423/',
                      'setting_Fr12': 'data/data2899425/',
                      'setting_Dossmann_Fr01_Re20000': 'data/data2938881/',
                      'setting_Dossmann_Fr01_Re20000_restart':
                          'data/data2938881/',
                      'setting_Dossmann_Fr01_Re10000': 'data/data2938940/',
                      'setting_Dossmann_Pr1_Re40000': 'data/data2939151/'}
        return dictionary[setting_name]
