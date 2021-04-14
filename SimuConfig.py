
try:
    from .Qobj import *
    from ._util import *
    from .global_var import *
    from .RWA import *
except:
    from Qobj import *
    from _util import *
    from global_var import *
    from RWA import *
import json
import matplotlib.pyplot as plt

class SysH():
    '''
    get circuit information from a dictionary/.json/.mat
    return init_state,tlist,H_idle,H_Mesolve,H_collapse,solver_options,requirements
    '''
    def __init__(self,config,rwa = False,rwa_freq = 6.0e9,omit_freq = 3e9,maximum_points=0):
        self.simu_config = self._preprocess(config)
        
        self.rotation_wave_approx = rwa
        self.rwa_freq = rwa_freq
        self.rwa_omit_freq = omit_freq

        if maximum_points ==0:
            self.max_pulse_len=get_maximum_arraylen_in_dict(self.simu_config['Channels'],max_len=1)
        else:
            self.max_pulse_len = maximum_points

        self._load_basic_setting()
        self.tlist=np.arange(0,(self.max_pulse_len-0.0001)/self.sampling_rate,1/self.sampling_rate)
        if self.init_state_label[0:5] =='RB60_':
            self.init_state = rb_idx2state( int(self.init_state_label[5:]),self.nTrunc )
        elif self.init_state_label[0:8] =='RB_C_60_':
            self.init_state = rb_idx2state_coupler_excited( int(self.init_state_label[8:]),self.nTrunc )
        else:
            self.init_state = str2state( self.init_state_label,self.nTrunc )
        self._load_channel()
        self.channel_plot={}

    def _preprocess(self,config):
        if isinstance(config,str):
            ## a file address is passed
            if config.endswith('.json'):
                return self._load_json(config)
            elif config.endswith('.mat'):
                return 'mat'
        elif isinstance(config,dict):
            return config
        else:
            raise Exception(f'unknown config type')

    def _load_json(self,config_address):
        with open(config_address,'r') as f:
            dict_ = json.load(f)
            return dict_

    def _load_basic_setting(self):
        self.setting_dict = get_dictionary_value(self.simu_config,'Basic Setting')
        self.sampling_rate = get_dictionary_value(self.setting_dict,'Sampling Rate')
        self.nTrunc = list(map(int,get_dictionary_value(self.setting_dict,'Trunc')))
        self.init_state_label = get_dictionary_value(self.setting_dict,'Initial State')

        self.solver_options = {}
        atol = get_dictionary_value(self.setting_dict,'Absolute Tolerance')
        rtol = get_dictionary_value(self.setting_dict,'Relative Tolerance')
        nsteps = get_dictionary_value(self.setting_dict,'Internal Steps')
        if atol:
            self.solver_options.update({'atol':atol})
        if rtol:
            self.solver_options.update({'rtol':rtol})
        if nsteps:
            self.solver_options.update({'nsteps':nsteps})

    def _load_channel(self):
        self.channel_dict = get_dictionary_value(self.simu_config,'Channels')
        if self.rotation_wave_approx:
            RWA_genertor = RWA(self.channel_dict,self.rwa_freq)
            RWA_genertor.generate_rwa_dict(omit_freq=self.rwa_omit_freq)

        self.collapse_dict = get_dictionary_value(self.simu_config,'Collapse Channels')
        self.H_Mesolve = self.get_H_Mesolve(self.channel_dict)
        self.H_collapse = self.get_H_Mesolve(self.collapse_dict)

    def get_H_Mesolve(self,dict_):
        H_list=[]
        H_constant= 0
        for value in dict_.values():
            H = self.get_H(value)
            if not isinstance(H,list):
                H_constant += H
            else:
                H_list.append(H)
        if H_constant:
            H_list.insert(0,H_constant)
        return H_list

    def get_H(self,channel_dict):
        oper_name = channel_dict['operator']
        value = channel_dict['value']
        operator = str2oper(oper_name,self.nTrunc) 
        if 'coefficient' in channel_dict.keys():
            coeff = channel_dict['coefficient'] 
        else:
            coeff = 1

        operator = operator * coeff

        if isinstance(value,(list,np.ndarray,dict)):
            if isinstance(value,dict):
                # if isinstance(value['I'],(float,int)) and isinstance(value['Q'],(float,int)) :
                #     return operator * value['I'] * TWO_PI_load
                try:
                    if value['I']==0 and value['Q']==0:
                        return 0 * operator
                except ValueError:
                    pass
                
                if 'complex' in value.keys():
                    complex_omega = value['complex'] 
                else:
                    complex_omega = False

                IQ_complex,omega = get_format_IQ(value['I'],value['Q'],value['lo_freq'])
                pulse_value=(trunc_array(IQ_complex*TWO_PI_load,self.max_pulse_len),omega*TWO_PI_load+0j if complex_omega else omega*TWO_PI_load)
            else:
                pulse_value = trunc_array(np.asarray(value)*TWO_PI_load,self.max_pulse_len )
            return [operator,pulse_value]
        else:
            return operator*value*TWO_PI_load


    def show_channels(self):
        self.channel_plot = {}
        self.constant_channel = {}
        self.get_array_channel_in_dict(self.channel_dict)
        self.get_array_channel_in_dict(self.collapse_dict)
        print(self.constant_channel)
        if len(self.channel_plot)>0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for key,value in self.channel_plot.items():
                if len(self.tlist) == len(value):
                    ax.plot(self.tlist,value,label=key)
            plt.legend()
        else:
            print('No pulse has been set')

    def get_array_channel_in_dict(self,dict_,prefix=''):
        if ('name' in dict_.keys()) and ('value' in dict_.keys()):
            if isinstance(dict_['value'],(list,np.ndarray)):
                self.channel_plot[prefix+dict_['name']] =trunc_array( dict_['value'],self.max_pulse_len )
                self.constant_channel[prefix+dict_['name']] = str(round(dict_['value'][0]/1e6,4) ) + 'MHz'
            elif isinstance(dict_['value'],float) or isinstance(dict_['value'],int):
                self.constant_channel[prefix+dict_['name']] = str(round(dict_['value']/1e6,4) ) + 'MHz'
            elif isinstance(dict_['value'],dict):
                    if not isinstance(dict_['value']['I'],(int,float)):
                        self.channel_plot[dict_['name']+'_I'] =trunc_array( dict_['value']['I'], self.max_pulse_len )
                    if not isinstance(dict_['value']['Q'],(int,float)):
                        self.channel_plot[dict_['name']+'_Q'] =trunc_array( dict_['value']['Q'], self.max_pulse_len )
        for value in dict_.values():
            if isinstance(value,dict):
                self.get_array_channel_in_dict(value,prefix='')

