#!/usr/bin/env python
#%% 
import InstrumentDriver
import numpy as np
from read_config import *
import time 
import os
import logging
log = logging.getLogger('LabberDriver')


try:
    from Qobj import *
    from Solver import *
    from SimuConfig import *
    from _util import *
    from global_var import *
    from ResultParser import *
except ImportError:
    from .Qobj import *
    from .Solver import *
    from .SimuConfig import *
    from ._util import *
    from .global_var import *
    from .ResultParser import *

current_path=os.path.realpath(__file__)
def mkdir(path):
    folder = os.path.exists(path)
    if not folder:                  
        os.makedirs(path)    


class Driver(InstrumentDriver.InstrumentWorker):
    """ This class implements a Single-qubit simulator"""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        self.name = 'Labber'
        circuit_type = self.getValue('Circuit Type')
        self.log(circuit_type)
        if circuit_type == 'JC Hamiltonian':
            self.circuit = JC_CircuitParser(self)
        elif circuit_type == 'Transmon Circuit': 
            self.circuit = Transmon_CircuitParser(self)
        else:
            self.circuit = CircuitParser(self)


    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation"""
        pass


    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument"""
        # print('name:',quant.name,':',value)
        if quant.name == 'Circuit Type':
            if value == 'JC Hamiltonian':
                self.circuit = JC_CircuitParser(self)
            elif value == 'Transmon Circuit':
                self.circuit = Transmon_CircuitParser(self)
            else:
                self.circuit = CircuitParser(self)
                print('unknown circuit type')
        elif quant.name == 'write json':
            if value:
                self.circuit.load_general_setting()
                self.circuit.read_circuit()
                # print('read circuit times')
                file_address =r'' + os.path.split(current_path)[0] + '\\.json\\' + time.strftime("%Y%m\\%d", time.localtime()) 
                mkdir(file_address)  
                circuit_name = self.getValue('file name') + '_' + str(int(self.getValue('file index')))
                self.circuit.write_dict_to_json(file_address = file_address ,circuit_name=circuit_name )

        quant.setValue(value)
        if self.getValue('auto update values'):
            self.circuit.load_general_setting()
            self.circuit.read_circuit()
        return value


    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation"""
        if quant.name in ('Eigenenergies','Eigenenergies maximum overlap','Ener Gap','Ener Gap Trace'):
            circuit_type = self.getValue('Circuit Type')
            if circuit_type not in ('JC Hamiltonian','Transmon Circuit'):
                return quant.getValue()

            self.circuit.load_general_setting()
            self.circuit.read_circuit()
            self.circuit.write_to_dict()
            self.Simu = SysH(self.circuit.Config_Dict)
            EigenS = EigenSolver(self.Simu.H_Mesolve,self.Simu.H_collapse,self.Simu.tlist,self.Simu.init_state,self.Simu.solver_options,1)
            if quant.name == 'Eigenenergies':
                EigenS.solve(0)
            else:
                EigenS.solve(1)
            EigenS.postprocess()
            EigenResult = EigenStateParser(EigenS.result,EigenS.tlist,list(map(int,self.Simu.nTrunc)))
            # self.log(EigenResult.eigen_ener)
            
            func_str = self.getValue('Level func')
            if quant.name == 'Ener Gap':
                return EigenResult.get_Ener_gap_trace(func_str)[0]
            elif quant.name == 'Ener Gap Trace':
                return EigenResult.get_Ener_gap_trace(func_str)
            else:
                return np.real(EigenResult.eigen_ener)

        return quant.getValue()


if __name__ == '__main__':
    pass
