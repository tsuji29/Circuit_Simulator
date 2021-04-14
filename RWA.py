#%%
try:
    from global_var import *
except ImportError:
    from .global_var import *

import numpy as np
import copy 

class RWA():

    def __init__(self,H_channel_dict,rwa_freq = 6.0e9):
        # self.H_channel_rwa_dict = {}
        self.H_channel_dict = H_channel_dict
        self.qubit_num = len( H_channel_dict['Channel1']['operator'].split(',') )

        if isinstance(rwa_freq,float):
            self.rwa_freq = rwa_freq * np.ones(self.qubit_num)
        elif len(rwa_freq) == self.qubit_num:
            self.rwa_freq =rwa_freq
        else:
            raise Exception('Incompatible array size')

    def generate_RWA_operators(self):
        for i in range(self.qubit_num):
            operator = ['I']*self.qubit_num
            operator[i] = 'aa'
            self.H_channel_dict[f'Channel_RWA_Q{i+1}'] = {'name':f'Q{i+1}_RWA',
                                                        'operator':','.join(operator),
                                                        'value':-1* self.rwa_freq[i]}

    def generate_rwa_dict(self,omit_freq= 3e9):
        self.get_splitted_H_channel_dict()

        keys = list(self.H_channel_dict.keys())
        for key in keys:
            self.H_channel_dict[key],pulse_omega = self.reform_channel(self.H_channel_dict[key])
            if abs(pulse_omega) >= omit_freq:
                self.H_channel_dict.pop(key)

        self.generate_RWA_operators()

    def reform_channel(self,channel):
        omega = 0 
        for i,oper in enumerate(channel['operator'].split(',')):
            if oper == Destroy:
                omega -= self.rwa_freq[i] 
            elif oper == Create:
                omega += self.rwa_freq[i]

        if omega ==0:
            return channel,omega
        else:
            channel['value'], pulse_omega = self.reform_channel_value(channel['value'],omega)
            return channel,pulse_omega

    def reform_channel_value(self,pulse,omega):
        if isinstance(pulse,(float,int,list,np.ndarray)):
            return {'I':pulse,'Q':0,'lo_freq':omega,'complex':True},omega
        elif isinstance(pulse,dict):
            I = np.asarray(pulse['I']) if isinstance(pulse['I'],list) else pulse['I']
            Q = np.asarray(pulse['Q']) if isinstance(pulse['Q'],list) else pulse['Q']
            if abs(pulse['lo_freq']+omega) < abs(pulse['lo_freq']-omega):
                return {'I':I/2,'Q':Q/2,'lo_freq': omega + pulse['lo_freq'] ,'complex':True},omega + pulse['lo_freq']
            else:
                return {'I':I/2,'Q':-1 * Q / 2,'lo_freq':omega -pulse['lo_freq'] ,'complex':True }, omega -pulse['lo_freq'] 

    def get_splitted_H_channel_dict(self):
        for idx in range(self.qubit_num):
            self.H_channel_dict = self.split_H_channel_dict(self.H_channel_dict,idx)


    def split_H_channel_dict(self,H_channel_dict,index = 0):
        # total_channel = len(H_channel_dict)
        keys = list(H_channel_dict.keys())
        for key in keys:
            # print(key)
            channel = H_channel_dict[key]
            oper_list = channel['operator'].split(',')
            if oper_list[index] in [X,Y]:
                if oper_list[index] == Y:
                    coeff1 = 1j
                    coeff2 = -1j
                else:
                    coeff1 = 1
                    coeff2 = 1

                channel1 = copy.deepcopy(channel)
                oper_list1 = copy.deepcopy(oper_list)
                oper_list1[index] = Create
                channel1['operator'] =','.join(oper_list1)
                channel1['coefficient'] = coeff1 if 'coefficient' not in channel1.keys() else channel1['coefficient']*coeff1
                H_channel_dict[key + f'_s{index}_1'] = channel1

                channel2 = copy.deepcopy(channel)
                oper_list2 = copy.deepcopy(oper_list)
                oper_list2[index] = Destroy
                channel2['operator'] =','.join(oper_list2)
                channel2['coefficient'] = coeff2 if 'coefficient' not in channel2.keys() else channel2['coefficient']*coeff2

                H_channel_dict[key + f'_s{index}_2'] = channel2

                H_channel_dict.pop(key)

        return H_channel_dict




#%%
# address = r'C:\Chuji\Code_and_Data\MyCode\Circuit_Simulator\.json\test_H_dict.json'

# with open(address,'r') as f:
#     dict_ = json.load(f)


# config_rwa = RWA(dict_['Channels'])
# config_rwa.get_rwa_dict()
# print(config_rwa.H_channel_dict.keys())

# config_rwa.get_H_channel_dict()
# print(config_rwa.H_channel_dict.keys())

# print(H_channel_dict.keys())

# %%
