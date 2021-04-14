try:
    from ._util import *
    from .global_var import *
    from .Qobj import *
except:
    from _util import *
    from global_var import *
    from Qobj import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import logging
log = logging.getLogger('LabberDriver')

class QuantumParser():

    def __init__(self):
        pass

class EigenStateParser(QuantumParser):

    def __init__(self,eigen,tlist,nTrunc):
        self.eigen_ener,self.eigen_state = eigen
        self.tlist = tlist
        self.nTrunc = nTrunc

    def show_eigenenergy(self,level_min=0,level_max = 3**10):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        total_energy_num = len(self.eigen_ener[0])
        for i in range( max(0,level_min) ,min(level_max,total_energy_num) ):
            ax.plot(self.tlist,self.eigen_ener[:,i])

    def get_Ener_gap_trace(self,func_str):
        # if len(self.eigen_ener) == 1:
        #     return self.get_Ener_gap(self.eigen_ener[0],func_str)
        # else:
        E_arr = np.array([])
        for eig_ener in self.eigen_ener:
            E_arr= np.append(E_arr,self.get_Ener_gap(eig_ener,func_str))
        return E_arr

    def get_Ener_gap(self,eigen_ener,func_str):
        # func_str='-101-100-001+000'
        # split_str_func(func_str) = ['', '-101', '-100', '-001', '+000']
        E=0
        for state in split_str_func(func_str):
            if state == '':
                pass
            elif state.startswith('+'):
                E += eigen_ener[sstr2idx( state[1:] ,self.nTrunc)]
            elif state.startswith('-'):
                E -= eigen_ener[sstr2idx( state[1:] ,self.nTrunc)]
            else:
                E += eigen_ener[sstr2idx( state ,self.nTrunc)]
        return np.real(E)



class UnitaryOperatorParser(QuantumParser):

    def __init__(self,result,nTrunc):
        self.U_trace = result
        self.nTrunc = nTrunc

    def get_U(self,index=-1):
        self.U = self.U_trace[index]
    
    def get_subspace_operator(self,subspace=['000','001','100','111']):
        self.U_sub=qt.Qobj(trunc_to_specific_subspace(self.U.full(),self.nTrunc,subspace))
        
    def set_Target_gate(self,gate_name):
        self.Target_Gate_name=gate_name
        if gate_name=='I':
            A=np.identity(4)
            self.Target_Gate=qt.Qobj(A)
        elif gate_name=='SWAP':
            A=np.identity(4)
            A[1,1]=0;A[2,2]=0;A[1,2]=1;A[2,1]=1
            self.Target_Gate=qt.Qobj(A)
        elif gate_name=='iSWAP':
            A=np.identity(4)*(1+0j)
            A[1,1]=0;A[2,2]=0;A[1,2]=1j;A[2,1]=1j
            self.Target_Gate=qt.Qobj(A)
        elif gate_name=='CZ':
            A=np.identity(4)
            A[3,3]=-1
            self.Target_Gate=qt.Qobj(A)
        elif gate_name=='CZ_I':
            A=np.identity(8)
            A[6,6]=-1
            A[7,7]=-1
            self.Target_Gate=qt.Qobj(A)
        elif gate_name=='I_CZ':
            A=np.identity(8)
            A[3,3]=-1
            A[7,7]=-1
            self.Target_Gate=qt.Qobj(A)

    def get_abs_Gate_Fidelity(self):
        dim = self.Target_Gate.dims[0][0]
        self.Gate_Fidelity_without_phase=qt.Qobj(np.abs((self.Target_Gate.dag()*self.U_sub))**2).tr()/dim

    def get_Gate_Fidelity(self):
        dim = self.Target_Gate.dims[0][0]
        self.Gate_Fidelity=np.abs(1/dim*(self.Target_Gate.dag()*self.U_sub).tr())**2

    def remove_single_qubit_gate(self):
        # self.U_sub = self.eliminate_single_qubit_rot(self.U_sub)
        if self.Target_Gate_name=='CZ':
            self.U_sub = self.eliminate_single_qubit_phase(self.U_sub)
        else:
            self.U_sub = self.remove_local_phases(self.U_sub)

    def remove_local_phases(self,U):
        diag_value=U.diag()
        phaseQ3=np.angle(diag_value[1])-np.angle(diag_value[0])
        phaseQ2=np.angle(diag_value[2])-np.angle(diag_value[0])
        phaseQ1=np.angle(diag_value[4])-np.angle(diag_value[0])
        U001 =qt.tensor(qt.qeye(2),qt.qeye(2),Z_Gate(-phaseQ3).U)
        U010 =qt.tensor(qt.qeye(2),Z_Gate(-phaseQ2).U,qt.qeye(2))
        U100 =qt.tensor(Z_Gate(-phaseQ1).U,qt.qeye(2),qt.qeye(2))
        U001.dims=[[8],[8]]
        U010.dims=[[8],[8]]
        U100.dims=[[8],[8]]
        Operator_after_rot=U001*U010*U100*U
        phase_global=np.angle(Operator_after_rot.diag()[0])
        return np.exp(-1j*phase_global)*Operator_after_rot

    def eliminate_single_qubit_phase(self,U):
        diag_value=U.diag()
        phaseQ2=np.angle(diag_value[1])-np.angle(diag_value[0])
        phaseQ1=np.angle(diag_value[2])-np.angle(diag_value[0])
        # U01=qt.tensor(qt.qeye(2),np.cos(-phaseQ2/2)*qt.qeye(2)-1j*np.sin(-phaseQ2/2)*qt.sigmaz())
        # U10=qt.tensor(np.cos(-phaseQ1/2)*qt.qeye(2)-1j*np.sin(-phaseQ1/2)*qt.sigmaz(),qt.qeye(2))
        U01 =qt.tensor(qt.qeye(2),Z_Gate(-phaseQ2).U)
        U10 =qt.tensor(Z_Gate(-phaseQ1).U,qt.qeye(2))
        U01.dims=[[4],[4]]
        U10.dims=[[4],[4]]
        Operator_after_rot=U10*U01*U
        # remove global phase
        phase_global=np.angle(Operator_after_rot.diag()[0])
        return np.exp(-1j*phase_global)*Operator_after_rot

    def eliminate_single_qubit_rot(self,U):
        Q1_rot_ax_phi,Q1_evolution_angle = self.get_rot_angle(U,1)
        Q2_rot_ax_phi,Q2_evolution_angle = self.get_rot_angle(U,2)
        rotations  = qt.tensor(SingleQubitGate(np.pi/2,Q1_rot_ax_phi,Q1_evolution_angle).U,SingleQubitGate(np.pi/2,Q2_rot_ax_phi,Q2_evolution_angle).U)
        return Qflatten(rotations) * U

    def get_rot_angle(self,U,qubit_index = 1):
        res = minimize(single_rot_cost_func,x0=(0.1,0.5),args=(U,qubit_index),method='Nelder-Mead')
        rot_ax_phi,evolution_angle = res.x
        return rot_ax_phi,evolution_angle


def single_rot_cost_func(params,*args):
    U = args[0]
    qubit_idx =args[1]
    if U.shape[0] == 2:
        return np.abs((SingleQubitGate(np.pi/2,params[0],params[1]).U * U)[0][0][1])
    elif U.shape[0] == 4:
        if qubit_idx ==1:
            rot_oper = qt.tensor(SingleQubitGate(np.pi/2,params[0],params[1]).U,qt.qeye(2))
            return np.abs((Qflatten(rot_oper) * Qflatten(U))[0][0][2])
        else:
            rot_oper = qt.tensor(qt.qeye(2),SingleQubitGate(np.pi/2,params[0],params[1]).U)
            return np.abs((Qflatten(rot_oper) * Qflatten(U))[0][0][1])
    else:
        raise Exception(f'Opearator of shape{U.shape} is not supported now!')


class QutipStateParser(QuantumParser):

    def __init__(self,rho_trace,tlist,nTrunc):
        ## nTrunc = [3,2,4]
        self.tlist = tlist 
        self.rho_trace = rho_trace
        self.state_dim = [ list(nTrunc) , list(nTrunc)]

    def set_measure(self,measure_label):
        rho_sub_trace = self.get_sub_trace(measure_label)
        sName = 'M_'+measure_label
        setattr(self,sName,QubitMeasure(self.tlist,rho_sub_trace))

    def get_sub_trace(self,measure_label):
        trace_seq=[int(i)-1 for i in measure_label.replace('Q','')]
        # sName='rho_trace_'+measure_label
        rho_sub_trace=[]
        for states in self.rho_trace:
            states.dims=self.state_dim
            if len(self.state_dim[0]) == len(trace_seq):
                rho_sub_trace.append(states)
            else:
                rho_sub_trace.append(states.ptrace(trace_seq))
        return rho_sub_trace

def list_to_str(list_):
    a=''
    for value in list_:
        a += str(value) 
    return a

class QubitMeasure():

    def __init__(self,tlist,rho_sub_trace):
        self.rho_trace = rho_sub_trace
        self.tlist = tlist
        self.nTrunc =self.rho_trace[0].dims[0]

        self.calculate_occupation_rate_all()
        self.listPauli = ['I','X','Y','Z']
        self.calculate_pauli()

    def calculate_occupation_rate_all(self):
        self.dict_Population_Trace={}
        for idx in range(multiply_list(self.nTrunc)):
            self.calculate_occupation_rate(idx)

    def calculate_occupation_rate(self,s_idx):
        population_list = []
        sName = 'Occu:' + idx2sstr(s_idx,self.nTrunc)
        for states in self.rho_trace:
            population_list.append( states[s_idx][0][s_idx] )
        self.dict_Population_Trace[sName]=np.real( np.asarray(population_list))

    def calculate_pauli(self):
        List_PauliBasis = self.get_all_PauliBasis( len(self.nTrunc))
        # print(List_PauliBasis)
        self.dict_PauliBasis = {}
        self.dict_Trace_Pauli ={}
        for pauli_basis in List_PauliBasis:
            self.dict_PauliBasis.update({pauli_basis:get_Pauli(pauli_basis)})
            self.dict_Trace_Pauli[pauli_basis] = []

        for rho in self.rho_trace:
            rho_lowest = qt.Qobj(trunc_to_lowest_two_level(rho.full(),self.nTrunc))
            for key,op in self.dict_PauliBasis.items():
                self.dict_Trace_Pauli[key].append( np.real((Qflatten(op)*rho_lowest).tr()) )

    def get_all_PauliBasis(self,nqubit):
        ## 'IXYZ'
        PauliBasis=[]
        for idx in range(4**nqubit):
            sIdx = idx2sstr(idx,4*np.ones(nqubit,dtype=int))
            PauliBasis.append( self.sIdx_2_sPauli(sIdx) )
        return PauliBasis

    def sIdx_2_sPauli(self,sIdx):
        ## '0103' -> 'IXIZ'
        sPauli = ''
        for s in sIdx:
            sPauli += self.listPauli[int(s)]
        return sPauli

    def show_occupation_rate(self):
        self.occu_rate_fig=plt.figure('Occupation')
        ax=self.occu_rate_fig.add_subplot(111)
        for key,value in self.dict_Population_Trace.items():
            ax.plot(self.tlist,value,label=key)
            ax.legend()
        ax.tick_params(labelsize=16)
        plt.show()

    def show_pauli(self):
        self.pauli_fig=plt.figure('Pauli')
        ax=self.pauli_fig.add_subplot(111)
        for key,value in self.dict_Trace_Pauli.items():
            ax.plot(self.tlist,value,label=key)
            ax.legend()
        plt.show()

