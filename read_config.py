try:
    from global_var import *
    from _util import *
except:
    from .global_var import *
    from ._util import *

import json
import sys


class CircuitParser(object):

    def __init__(self,config_file):
        self.config = self._preprocessing(config_file)
        self.load_general_setting()
        self.Channels={}
        self.channel_count = 0

    def _preprocessing(self,config_file):
        
        if isinstance(config_file,dict):
            self.config_type = 'dict'
            return config_file

        try:
            if config_file.name in ['labber','Labber']:
                self.config_type = 'labber'
        except:
            self.config_type = None

        return config_file

    @property
    def element_num(self):
        return self.qubit_num

    def read_circuit(self):
        self.Channels={}
        self.channel_count = 0
        self.Collapse_Channels={}
        self.collapse_channel_count = 0

    def _get_config_value(self,quantity,datatype=None):
        ## get quantity value, and transfer value to specific datatype
        if self.config_type == 'labber':
            value = self.config.getValue(quantity)
        elif self.config_type =='dict':
            value = self.config.get(quantity)

        if datatype=='str':
            return str(value)
        elif datatype=='float':
            return float(value)
        elif datatype=='int':
            return int(float(value))
        else:
            return value

    def _get_config_trace(self,quantity):
        ## get trace value
        if self.config_type == 'labber':
            try:
                value = self.config.getValueArray(quantity)
            except Exception:
                value = 0
        elif self.config_type =='dict':
            value = self.config.get(quantity)

        if len(value)==0:
            value = 0
        return value

    def _set_config(self,quantity,value):
        if self.config_type=='labber':
            self.config.setValue(quantity,value)

    def _write_channel_config(self,dict_value):
        self.channel_count+=1
        self.Channels.update({f'Channel{self.channel_count}':dict_value})

    def _write_collapse_channel(self,dict_value):
        self.collapse_channel_count+=1
        self.Collapse_Channels.update({f'Collapse Channel{self.collapse_channel_count}':dict_value})

    def write_to_dict(self):
        self.Config_Dict = {'Basic Setting':self.Basic_setting,'Channels':self.Channels,'Collapse Channels':self.Collapse_Channels}

    def write_dict_to_json(self,file_address=r'C:\Users\tsuji\Labber\Drivers\Circuit_Simulator\.json',circuit_name='test'):
        channels = get_checked_dict(self.Channels)
        collapse_channels = get_checked_dict(self.Collapse_Channels)
        basic_settings = get_checked_dict(self.Basic_setting)
        All_Circuit_Parameters = {'Basic Setting':basic_settings,'Channels':channels,'Collapse Channels':collapse_channels}
        # self.config.log(All_Circuit_Parameters)
        with open(file_address+'\\' + circuit_name + '.json', 'w') as f:
            # dict_format=json.dumps(new_channel)
            dict_format=json.dumps(All_Circuit_Parameters,indent=4)
            f.write(dict_format)

    def generate_full_operator_name(self,operator_list,index_list):
        # length of operator-list should be the same as index-list
        if isinstance(index_list,(list,tuple,np.ndarray)):
            ascend_index = np.argsort(index_list)
            index_array=np.asarray(index_list)[ascend_index]
            operator_array=np.asarray(operator_list)[ascend_index]
        else:
            index_array = [index_list]
            operator_array = [operator_list]
        s=''
        k=0
        for index in range(1,self.element_num+1):
            if index in index_array:
                s+=operator_array[k]+','
                k+=1
            else:
                s+=Identity+','
        return s[:-1]

    def _load_channel(self,name='freq',operator=AA,value=None,index=1,isdict=False):
        operator_name=self.generate_full_operator_name(operator,index)
        if isdict:
            I,Q,LO=value
            value= get_IQ_dict(I,Q,LO)
        if isinstance(index,(list,tuple,np.ndarray)):
            for index_value in index:
                name += f'{index_value}' 
            return ValuesToDict( name,operator_name,value)
        else:
            return ValuesToDict(f'Q{index}_'+name,operator_name,value)

    def load_general_setting(self):
        self.Basic_setting={}
        self.circuit_type = self._get_config_value('Circuit Type',datatype='str')
        self.qubit_num = self._get_config_value('Number of Qubits',datatype='int')
        if self._get_config_value('Uniform Truncation',datatype='bool'):
            self.trunc= np.ones(self.qubit_num) *self._get_config_value('Uniform Truncation Degree',datatype='int')
        else:
            self.trunc= np.asarray([ max( 2,self._get_config_value(f'Truncation Q{i}','int') )  for i in range(1,self.qubit_num + 1 )])

        self.Basic_setting.update({'Circuit Type':self.circuit_type})
        self.Basic_setting.update({'Element Number':self.qubit_num})
        self.Basic_setting.update({'Trunc':self.trunc})
        

        self.init_state = self._get_config_value('Initial State',datatype='str')
        if self.init_state in ['','Basis','basis','+Z','ground','Ground']:
            self.init_state = '+Z'*self.qubit_num
        self.Basic_setting.update({'Initial State':self.init_state})


        dict_setting_type = {'Sampling Rate':'float',
                             'Absolute Tolerance':'float',
                             'Relative Tolerance':'float',
                             'Internal Steps':'float',
                             }
        for key,value in dict_setting_type.items():
            self.Basic_setting.update( {key:self._get_config_value(key,datatype=value)} )



class JC_CircuitParser(CircuitParser):

    def __init__(self,config_file):
        super().__init__(config_file)

    def _check_status(self):
        self.modulation_on= self._get_config_value('Modulation On')
        self.decoherence_on = self._get_config_value('Decoherence On')
        self.uniform_LO = self._get_config_value('Uniform XY LO')
        if self.uniform_LO:
            self.LO_Freq=self._get_config_value('XY LO Freq')

    def read_circuit(self,remove_zeros=True):
        self.Channels={}
        self.channel_count = 0
        self.Collapse_Channels={}
        self.collapse_channel_count = 0
        self._check_status()
        for qubit_index in range(1,self.element_num+1):
            self._load_element_freq_and_anhar(qubit_index)
            if self.modulation_on:
                self._load_element_transverse_pulse(qubit_index)
            if self.decoherence_on:
                self._load_element_decoherence(qubit_index)

        for i in range(1,self.element_num+1):
            for j in range(i+1,self.element_num+1):
                self._load_coupling_network(i,j)

        if remove_zeros:
            self.get_cutted_pulses()

    def get_cutted_pulses(self):
        (channel_start,channel_end) = get_start_and_end_index_in_dict(self.Channels)
        # print('start',channel_start,'end',channel_end)
        self.Channels = cut_dict_by_start_and_end_index(self.Channels,channel_start,channel_end)
        (collapse_channel_start,collapse_channel_end) = get_start_and_end_index_in_dict(self.Collapse_Channels)
        self.Collapse_Channels = cut_dict_by_start_and_end_index(self.Collapse_Channels,collapse_channel_start,collapse_channel_end)


    def _load_element_freq_and_anhar(self,qubit_index):
        freq = self._get_config_value(f'Q{qubit_index} f01','float')
        anhar = self._get_config_value(f'Q{qubit_index} Anharmonicity','float')
        if self.modulation_on:
            Z = self._get_config_trace(f'Q{qubit_index} Z')
        else:
            Z = 0
        freq_dict=self._load_channel('freq',operator=AA,value=(freq+Z)*TWO_PI_write,index=qubit_index)
        anhar_dict=self._load_channel('anhar',operator=AAAA,value=anhar/2*TWO_PI_write,index=qubit_index)
        self._write_channel_config(freq_dict)
        self._write_channel_config(anhar_dict)

    def _load_element_decoherence(self,qubit_index):
        T1 = self._get_config_value(f'Q{qubit_index} T1','float')
        Tphi = self._get_config_value(f'Q{qubit_index} Tphi','float')
        T1_dict=self._load_channel('decay',operator=Destroy,value=np.sqrt(1/T1)/TWO_PI_load,index=qubit_index)
        Tphi_dict=self._load_channel('dephase',operator=AA,value=np.sqrt(0.5/Tphi)/TWO_PI_load,index=qubit_index)
        self._write_collapse_channel(T1_dict)
        self._write_collapse_channel(Tphi_dict)

    def _load_element_transverse_pulse(self,qubit_index):
        I= self._get_config_trace(f'Q{qubit_index} XY-I')
        Q= self._get_config_trace(f'Q{qubit_index} XY-Q')

        if self.uniform_LO:
            LO_Freq=self.LO_Freq
        else:
            LO_Freq= self._get_config_value(f'Q{qubit_index} XY LO Freq','float')
        Drive_dict=self._load_channel('drive',operator=X,value=[I*TWO_PI_write,Q*TWO_PI_write,LO_Freq*TWO_PI_write],index=qubit_index,isdict=True)
        self._write_channel_config(Drive_dict)

    def _load_coupling_network(self,i,j):
        g = self._get_config_value(f'g{i}{j}','float')
        coupling_dict = self._load_channel('g',operator=[X,X],value=g*TWO_PI_write,index=[i,j],isdict=False)
        self._write_channel_config(coupling_dict)


class Transmon_CircuitParser(JC_CircuitParser):

    def __init__(self,config_file):
        super().__init__(config_file)
    
    def _check_status(self):
        self.modulation_on= self._get_config_value('Modulation On')
        self.decoherence_on = self._get_config_value('Decoherence On')
        self.uniform_LO = self._get_config_value('Uniform XY LO')
        if self.uniform_LO:
            self.LO_Freq=self._get_config_value('XY LO Freq')
        self.use_capacitance = self._get_config_value('Use Capacitance Network')
        self._get_rmatrix()

    def _load_element_freq_and_anhar(self,qubit_index):
        f01_max = self._get_config_value(f'Q{qubit_index} f01_max')
        f01_min = self._get_config_value(f'Q{qubit_index} f01_min')
        Voltage_Period = self._get_config_value(f'Q{qubit_index} Voltage period')
        Voltage_operating_point = self._get_config_value(f'Q{qubit_index} Voltage operating point')
        if self.modulation_on:
            flux = self._get_config_trace(f'Q{qubit_index} Flux')
        else:
            flux = 0
        
        if self.use_capacitance:
            Ec = electron**2 / h / (2 * self.effective_c_matrix[qubit_index-1][qubit_index-1] ) *1e15
        else:
            Ec = self._get_config_value(f'Q{qubit_index} Ec','float')

        freq,anhar= self.get_transmon_freq(f01_max,f01_min,Ec,Voltage_Period,Voltage_operating_point,flux)
        freq_dict=self._load_channel('freq',operator=AA,value=freq*TWO_PI_write,index=qubit_index)
        anhar_dict=self._load_channel('anhar',operator=AAAA,value=anhar*TWO_PI_write/2,index=qubit_index)
        self._write_channel_config(freq_dict)
        self._write_channel_config(anhar_dict)
        setattr(self,f'Q{qubit_index}_freq',freq)

        self._set_config(f'Q{qubit_index} f01', freq if isinstance(freq,float) else freq[0] )
        self._set_config(f'Q{qubit_index} Anharmonicity', anhar if isinstance(anhar,float) else anhar[0] )
        self._set_config(f'Q{qubit_index} Ec', Ec )

    def _load_coupling_network(self,i,j):
        g_value = self.r_matrix[i-1][j-1] * np.sqrt(getattr(self,f'Q{i}_freq') *getattr(self,f'Q{j}_freq'))
        coupling_dict = self._load_channel('g',operator=[X,X],value=g_value*TWO_PI_write,index=[i,j],isdict=False)
        self._set_config(f'g{i}{j}',g_value if isinstance(g_value,float) else g_value[0] )
        self._write_channel_config(coupling_dict)

    def _get_rmatrix(self):
        if self.use_capacitance:
            C_matrix = np.zeros([self.element_num,self.element_num])
            for i in range(0,self.element_num):
                for j in range(i,self.element_num):
                    C_matrix[i,j] = self._get_config_value(f'C{i+1}{j+1}')
            self.r_matrix,self.effective_c_matrix=self.generate_r_matrix_from_C_matrix(C_matrix)
            for i in range(0,self.element_num):
                for j in range(i+1,self.element_num):
                    self._set_config(f'r{i+1}{j+1}',self.r_matrix[i,j])
        else:
            self.r_matrix = np.zeros([self.element_num,self.element_num])
            for i in range(0,self.element_num):
                for j in range(i+1,self.element_num):
                    self.r_matrix[i,j] = self._get_config_value(f'r{i+1}{j+1}')

    def get_transmon_freq(self,f01_max,f01_min,Ec,Voltage_Period,Voltage_operating_point,flux):
        if Voltage_Period >= 0:
            Ej_max = (f01_max + Ec)**2/(8*Ec)
            d = (f01_min + Ec)**2/(8*Ec*Ej_max)
            phi = (Voltage_operating_point + flux)/ Voltage_Period
            # self.config.log('d',d)
            # self.config.log('Ej_max',phi)
            Ej_flux = Ej_max * np.sqrt( np.cos(np.pi * phi)**2 + np.sin(np.pi*phi)**2 * d**2)
            
            freq= np.sqrt(8*Ej_flux * Ec) -Ec
        elif Voltage_Period == -1:
            freq = (f01_min-f01_max) * flux + f01_max
        elif Voltage_Period == -2:
            freq = (f01_max-f01_min) * flux + f01_min
        anhar = -1 * Ec
        return freq,anhar

    def generate_r_matrix_from_C_matrix(self,C_matrix):
        # get C factor for Kinetic energy terms in Hamiltonian.
        dim = len(C_matrix)
        C_factor = -1*(C_matrix.T+C_matrix)
        for i in range(dim):
            C_factor[i,i]=-1 * ( sum(C_factor[i,:])-C_factor[i,i] ) + C_matrix[i,i]
        # self.config.log('C_matrix',C_matrix)
        # self.config.log('C_factor',C_factor)
        A = np.linalg.inv(C_factor)
        # self.config.log('A_factor',A)
        r_matrix=np.zeros_like(C_matrix)
        for i in range(0,dim):
            for j in range(i+1,dim):
                r_matrix[i,j] = 0.5 * A[i,j] *np.sqrt(C_factor[i,i]*C_factor[j][j] )
        return r_matrix,C_factor












