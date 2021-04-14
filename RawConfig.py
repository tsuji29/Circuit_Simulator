
### this file defines standard dictionary

try:
    from .global_var import *
except:
    from global_var import *


class RawConfig():

    def __init__(self,qubit_num,dimension,circuit_type=1,initial_state='ground',sampling_rate=1e9,atol=1e-9,rtol=1e-7,nsteps=2500):
        self.config = {}
        Circuit_Type_Dict = {0:JC_CIRCUIT,1:TRANSMON_CIRCUIT,2:LCJ_CIRCUIT}
        self.config['Circuit Type'] = Circuit_Type_Dict[circuit_type]
        self.config['Number of Qubits'] = qubit_num
        
        if isinstance(dimension,(float, int)):
            self.config['Uniform Truncation']=True
            self.config['Uniform Truncation Degree']= dimension
        else:
            for i,dims in enumerate(dimension):
                self.config[f'Truncation Q{i}'] = dims

        self.config['Initial State'] = initial_state
        self.config['Sampling Rate'] = sampling_rate
        self.config['Absolute Tolerance'] = atol
        self.config['Relative Tolerance'] = rtol
        self.config['Internal Steps'] = nsteps

    def load_default_value(self,modulation=False, decoherence = False,use_capacitance=False):
        qubit_num = self.config['Number of Qubits']
        self._load_XY_channels(qubit_num,modulation)
        self._load_Decay_channels(qubit_num,decoherence)

        if self.config['Circuit Type']  == JC_CIRCUIT:
            self._load_JC_parameters(qubit_num,modulation)
        elif self.config['Circuit Type'] == TRANSMON_CIRCUIT:
            self._load_transmon_paramters(qubit_num,use_capacitance,modulation)
        else: 
            pass

    def show_items(self):
        for key,value in self.config.items():
            print(key,':',value)

    def setValue(self,quant_name,value):
        if quant_name in self.config.keys():
            if quant_name in ('Circuit Type','Number of Qubits'):
                print(f'setting {quant_name} is not allowed')
            else:
                self.config[quant_name] = value
                if isinstance(value,(int,float)):
                    # print(quant_name,':',value)
                    pass
                else:
                    # print('channel:',quant_name,' is setted')
                    pass
        else:
            raise Exception(f'{quant_name} in not in dictionary keys')

    def get_dict(self):
        return self.config

    def _load_XY_channels(self,qubit_num,modulation):
        self.config['Modulation On'] = modulation
        self.config['Uniform XY LO'] = True
        self.config['XY LO Freq'] = 6e9
        if modulation:
            for q_idx in range(1,qubit_num+1):
                self.config[f'Q{q_idx} XY-I'] = []
                self.config[f'Q{q_idx} XY-Q'] = []
                self.config[f'Q{q_idx} XY LO Freq'] = 6e9


    def _load_Decay_channels(self,qubit_num,decoherence):
        self.config['Decoherence On'] = decoherence
        if decoherence:
            for q_idx in range(1,qubit_num+1):
                self.config[f'Q{q_idx} T1'] = 10e-6
                self.config[f'Q{q_idx} Tphi'] = 10e-6

    def _load_JC_parameters(self,qubit_num,modulation):
        for q_idx in range(1,qubit_num+1):
            self.config[f'Q{q_idx} f01'] = 6e9 + q_idx *150e6
            self.config[f'Q{q_idx} Anharmonicity'] = -250e6
            if modulation:
                self.config[f'Q{q_idx} Z'] = []

        for i in range(1,qubit_num):
            for j in range(i+1,qubit_num+1):
                self.config[f'g{i}{j}'] = 100e6 / ( abs(i-j) )**3

    def _load_transmon_paramters(self,qubit_num,use_capacitance,modulation):
        for q_idx in range(1,qubit_num+1):
            self.config[f'Q{q_idx} f01_max'] = 6.4e9 + q_idx * 100e6 
            self.config[f'Q{q_idx} f01_min'] = 1.8e9
            self.config[f'Q{q_idx} Voltage period'] = 1
            self.config[f'Q{q_idx} Voltage operating point'] = 0
            if modulation:
                self.config[f'Q{q_idx} Flux'] = []

        self.config['Use Capacitance Network'] = use_capacitance
        
        if use_capacitance:
            for i in range(1,qubit_num+1):
                for j in range(i,qubit_num+1):
                    if i==j:
                        self.config[f'C{i}{j}'] = 70 
                    else:
                        self.config[f'C{i}{j}'] = 2.0 / (abs(i-j))**3
        else:
            for i in range(1,qubit_num+1):
                self.config[f'Q{i} Ec'] = 240e6
            for i in range(1,qubit_num):
                for j in range(i+1,qubit_num+1):
                    self.config[f'r{i}{j}'] = 0.016 / (abs(i-j))**3



