

from .Qobj import *
from .Solver import *
from .SimuConfig import *
from .Circuit_Simulator import *
from .ResultParser import *
from .RawConfig import *
from .read_config import *


def read_config(config_dict,remove_zeros=True):
    circuit_type = config_dict['Circuit Type']

    if circuit_type == 'JC Hamiltonian':
        circuit = JC_CircuitParser(config_dict)
    elif circuit_type == 'Transmon Circuit':
        circuit = Transmon_CircuitParser(config_dict)
    else:
        circuit = CircuitParser(config_dict)

    circuit.load_general_setting()
    circuit.read_circuit(remove_zeros=remove_zeros)
    circuit.write_to_dict()
    return circuit.Config_Dict

