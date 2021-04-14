TRANSMON_CIRCUIT = 'Transmon Circuit'
JC_CIRCUIT = 'JC Hamiltonian'
LCJ_CIRCUIT = 'LCJ Circuit'
MAX_QUBITS = 9

__version__ = '1.0.0'


if __name__ == "__main__":
    from pathlib import Path
    from driver_config import *
    # General setting
    dir_path = Path(__file__).parent
    f = LDriverDefinition(dir_path/'Labber_API.ini')
    f.add_general_settings(
        name='Circuit_Simulator',
        version=__version__,
        driver_path='Labber_API',
        signal_analyzer=True,
        signal_generator=True,
    )


    qubit_list = list(range(1, MAX_QUBITS+1))

    #region Section: General setting
    f.add_section('General setting')
    #region Group: Basic Circuit Parameter
    f.add_group('Basic Circuit Parameter')

    circuit_type = LCombo(
        'Circuit Type',
        def_value=TRANSMON_CIRCUIT ,
        combo={TRANSMON_CIRCUIT,JC_CIRCUIT,LCJ_CIRCUIT},
    )
    f.add_quantity(circuit_type)

    combo_qubits = LCombo(
        'Number of Qubits',
        combo=qubit_list,
        def_value=qubit_list[1],
    )
    f.add_quantity(combo_qubits)

    f.add_quantity(LDouble(
        'Sampling Rate',
        def_value=1e9,
    )
    )

    f.add_group('Basic Simulation Parameter')

    f.add_quantity(LString(
        'Initial State',
        def_value='basis'
    )
    )

    f.add_quantity(LDouble(
        'Absolute Tolerance',
        def_value=1e-9,
    )
    )

    f.add_quantity(LDouble(
        'Relative Tolerance',
        def_value=1e-7,
    )
    )

    f.add_quantity(LDouble(
        'Internal Steps',
        def_value=2500,
    )
    )

    uniformTruc= LBoolean(
        'Uniform Truncation',
        def_value=True,
    )
    f.add_quantity(uniformTruc)

    f.add_quantity(LDouble(
        'Uniform Truncation Degree',
        def_value=3,
        state_quant=uniformTruc,
        states=True,
    )
    )

    for i in qubit_list:
        f.add_quantity(LDouble(
            f'Truncation Q{i}',
            def_value=3,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )


    f.add_group('General Utility')

    f.add_quantity(LString(
        'file name',
        def_value = 'test',
    )
    )

    f.add_quantity(LDouble(
        'file index',
        def_value = 0,
    )
    )

    f.add_quantity(LButton(
        'write json',
    )
    )

    f.add_quantity(LButton(
        'write mat',
    )
    )

    f.add_quantity(LBoolean(
        'auto update values',
        def_value=False,
    )
    )

    #end Group: Basic Circuit Parameter
    #end Section: General setting

    #region Section: 
    f.add_section('JC Hamiltonian')
    #region Group: 
    f.add_group('Qubits')
    
    for i in qubit_list:
        f.add_quantity(LDouble(
            f'Q{i} f01',
            unit='Hz',
            def_value=6.0e9,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )

    for i in qubit_list:
        f.add_quantity(LDouble(
            f'Q{i} Anharmonicity',
            unit='Hz',
            def_value=-300e6,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )

    f.add_group('Coupling Network')

    for i in qubit_list[:]:
        for j in qubit_list[i:]:        
            f.add_quantity(LDouble(
                f'g{i}{j}',
                unit='Hz',
                def_value=20e6,
                state_quant= combo_qubits,
                states=qubit_list[j-1:]
            )
            )

    f.add_group('Modulation')

    f.add_quantity(LBoolean(
        'Modulation On',
        def_value=True
    )
    )

    uniform_XY_LO=LBoolean('Uniform XY LO',
        def_value=True,
    )
    f.add_quantity(uniform_XY_LO)

    f.add_quantity(LDouble(
        f'XY LO Freq',
        unit='Hz',
        def_value=0,
        state_quant=uniform_XY_LO,
        states=True,
        show_in_measurement_dlg=True,
    )
    )

    
    for i in qubit_list:
        f.add_quantity(LVector(
            f'Q{i} XY-I',
            unit='Hz3',
            x_name='Time',
            x_unit='s',
            permission='WRITE',
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
            show_in_measurement_dlg=True,
        )
        )

        f.add_quantity(LVector(
            f'Q{i} XY-Q',
            unit='Hz',
            x_name='Time',
            x_unit='s',
            permission='WRITE',
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
            show_in_measurement_dlg=True,
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} XY LO Freq',
            unit='Hz',
            def_value=0,
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
            show_in_measurement_dlg=True,
        )
        )

        f.add_quantity(LVector(
            f'Q{i} Z',
            unit='Hz',
            x_name='Time',
            x_unit='s',
            permission='WRITE',
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
            show_in_measurement_dlg=True,
        )
        )


    f.add_group('Decoherence')

    decoherence=LBoolean(
        'Decoherence On',
        def_value=True
    )
    f.add_quantity(decoherence)

    for i in qubit_list:
        f.add_quantity(LDouble(
            f'Q{i} T1',
            unit='s',
            def_value=10e-6,
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} Tphi',
            unit='s',
            def_value=10e-6,
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
        )
        )


    f.add_section('Transmon Circuit')
    #region Group: 
    f.add_group('Qubits')

    for i in qubit_list:
        f.add_quantity(LDouble(
            f'Q{i} f01_max',
            unit='Hz',
            def_value=7.2e9,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} f01_min',
            unit='Hz',
            def_value=4.0e9,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} Ec',
            unit='Hz',
            def_value=200e6,
            state_quant= combo_qubits,
            states=qubit_list[i-1:]
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} Voltage period',
            unit='V',
            def_value=1,
            state_quant=combo_qubits,
            states=qubit_list[i-1:]
        )
        )

        f.add_quantity(LDouble(
            f'Q{i} Voltage operating point',
            unit='V',
            def_value=0,
            state_quant=combo_qubits,
            states=qubit_list[i-1:]
        )
        )


    f.add_group('Coupling_Network')

    for i in qubit_list[:]:
        for j in qubit_list[i:]:        
            f.add_quantity(LDouble(
                f'r{i}{j}',
                def_value=0.01,
                state_quant= combo_qubits,
                states=qubit_list[j-1:]
            )
            )

    f.add_group('Modulation')    

    # uniform_Charge_LO=LBoolean('Uniform Charge LO',
    #     def_value=True,
    # )
    # f.add_quantity(uniform_Charge_LO)

    # f.add_quantity(LDouble(
    #     f'Charge LO Freq',
    #     unit='Hz',
    #     def_value=0,
    #     state_quant=uniform_Charge_LO,
    #     states=True,
    #     show_in_measurement_dlg=True,
    # )
    # )    

    # for i in qubit_list:
    #     f.add_quantity(LVector(
    #         f'Q{i} Charge-I',
    #         unit='V',
    #         x_name='Time',
    #         x_unit='s',
    #         permission='WRITE',
    #         state_quant=combo_qubits,
    #         states=qubit_list[i-1:],
    #         show_in_measurement_dlg=True,
    #     )
    #     )

    #     f.add_quantity(LVector(
    #         f'Q{i} Charge-Q',
    #         unit='V',
    #         x_name='Time',
    #         x_unit='s',
    #         permission='WRITE',
    #         state_quant=combo_qubits,
    #         states=qubit_list[i-1:],
    #         show_in_measurement_dlg=True,
    #     )
    #     )

    #     f.add_quantity(LDouble(
    #         f'Q{i} Charge LO Freq',
    #         unit='Hz',
    #         def_value=0,
    #         state_quant=combo_qubits,
    #         states=qubit_list[i-1:],
    #         show_in_measurement_dlg=True,
    #     )
    #     )

    for i in qubit_list:
        f.add_quantity(LVector(
            f'Q{i} Flux',
            unit='V',
            x_name='Time',
            x_unit='s',
            permission='WRITE',
            state_quant=combo_qubits,
            states=qubit_list[i-1:],
            show_in_measurement_dlg=True,
        )
        )


    f.add_group('Capacitance Network')

    f.add_quantity(LBoolean(
        'Use Capacitance Network',
        def_value=False,
    )
    ) 

    for i in qubit_list[:]:
        for j in qubit_list[i-1:]:        
            f.add_quantity(LDouble(
                f'C{i}{j}',
                unit='fF',
                def_value=70,
                state_quant=combo_qubits,
                states=qubit_list[j-1:]
            )
            )
    
    #region Section: 
    f.add_section('LCJ Circuit')
    #region Group: 

    f.add_quantity(LVector(
        f'L Matrix',
        permission='WRITE',
    )
    )

    f.add_quantity(LVector(
        f'C Matrix',
        permission='WRITE',
    )
    )

    f.add_quantity(LVector(
        f'J Matrix',
        permission='WRITE',
    )
    )

    # Section:
    f.add_section('Output')
    f.add_group('Eigenenergies')

    f.add_quantity(LVector(
        'Eigenenergies',
        unit='Hz',
        x_name='Level',
        permission='READ',
    )
    )

    f.add_quantity(LVector(
        'Eigenenergies maximum overlap',
        unit='Hz',
        x_name='Level',
        permission='READ',
    )
    )

    f.add_quantity(LDouble(
        'Ener Gap',
        unit='Hz',
        permission='READ',
    )
    )

    f.add_quantity(LVector(
        'Ener Gap Trace',
        unit='Hz',
        permission='READ',
    )
    )

    f.add_quantity(LString(
        'Level func',
        permission='WRITE',
    )
    )

# f.add_quantity(LVector(
#     'Raw Hamiltonian',
#     unit='Hz',
#     permission='READ',
# )
# )

# [Diag Hamiltonian]
# label: Diag Hamiltonian
# datatype: VECTOR
# #unit: Hz
# #def_value: 10e6
# permission: READ
# group: Operator
# section: Output

# [Eigenenergies]
# label: Eigenenergies
# datatype: VECTOR
# unit: Hz
# x_name : Level
# x_unit: arb
# permission: READ
# group: Operator
# section: Output

# [Eigenenergies maximum overlap]
# label: Eigenenergies maximum overlap
# datatype: VECTOR
# unit: Hz
# x_name : Level
# x_unit: arb
# permission: READ
# group: Operator
# section: Output

# [Eigen Sorted index]
# datatype: VECTOR
# x_name : Level
# x_unit: arb
# permission: READ
# group: Operator
# section: Output




