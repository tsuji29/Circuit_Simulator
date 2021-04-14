
try:
    from ._util import *
    from .SimuConfig import *
except ImportError:
    from _util import *
    from SimuConfig import *

import qutip as qt
import copy

class Solver():

    def __init__(self,H_Mesolve,H_collapse,tlist,init_state,**kwargs):
        self.H_Mesolve = copy.deepcopy(H_Mesolve)
        self.H_collapse = copy.deepcopy(H_collapse)
        self.tlist = tlist
        self.init_state =copy.deepcopy(init_state)

        if len(tlist) == 1:
            resample_factor = 1
        else:
            try:
                resample_factor = int(kwargs['resample_factor'])
            except Exception:
                print_tips('resample_factor is not defined, set to 1' )
                resample_factor = 1
        
        self._preprocess(resample_factor)

    def _preprocess(self,resample_factor):
        ## try to resampling pulse
        self._resample_pulses(resample_factor)
        ## get Hidle
        self.H_idle = self.H(0,idle = True)
        ## get transformation Unitary Operator and get transformed init state
        self.get_tranformation_U(rot_init_state=True)
        self.init_state = self.U_logic_to_lab * self.init_state

    def solve(self):
        self.raw_result = 0

    def postprocess(self):
        self.result = self.raw_result

    def _resample_pulses(self,resample_factor):
        ## resample pulses by interpolate points
        if resample_factor == 1:
            ## do nothing
            pass
        else:
            dt=(self.tlist[1]-self.tlist[0])/resample_factor
            self.original_tlist = copy.deepcopy(self.tlist)
            self.tlist = np.arange(0,self.tlist[-1],dt)
            for H in self.H_Mesolve:
                if isinstance(H,list):
                    if isinstance(H[1],tuple):
                        real_part = np.real(H[1][0])
                        imag_part = np.imag(H[1][0])
                        H[1] = ( get_Interpolate_pulse(self.original_tlist,real_part,self.tlist) + 1j * get_Interpolate_pulse(self.original_tlist,imag_part,self.tlist),H[1][1] )
                    else:
                        H[1] = get_Interpolate_pulse(self.original_tlist,H[1],self.tlist)

    def H(self,t_index,idle = False):
        H_t = 0
        for H in self.H_Mesolve:
            if isinstance(H,list):
                if isinstance(H[1],tuple):
                    if isinstance(H[1][1],(float,int)):
                        H_t += H[0]*np.real(H[1][0][t_index]*np.exp(1j*H[1][1] * self.tlist[t_index] ))
                    elif isinstance(H[1][1],complex):
                        H_t += H[0]*(H[1][0][t_index]*np.exp(1j*np.real(H[1][1]) * self.tlist[t_index] ))
                else:
                    H_t += H[0]*H[1][t_index]
            else:
                H_t+=H
        return H_t

    def get_tranformation_U(self,rot_init_state = False):
        self.vals_idle, self.vecs_idle = eigensolve_close(self.H_idle.full())
        self.H_idle_diag = qt.Qobj(np.diag(self.vals_idle))
        if rot_init_state:
            self.U_logic_to_lab = qt.Qobj(self.vecs_idle)
        else:
            self.U_logic_to_lab = qt.qeye(len(self.vals_idle))
        self.U_lab_to_logic = self.U_logic_to_lab.dag()


class EigenSolver(Solver):

    def __init__(self,H_Mesolve,H_collapse,tlist,init_state,**kwargs):
        super().__init__(H_Mesolve,H_collapse,tlist,init_state,**kwargs)

        if 'eigen_cloest_to_bare' in kwargs.keys():
            self.eigen_cloest = kwargs['eigen_cloest_to_bare']
        else:
            self.eigen_cloest = False

        if 'sort_by_maximum_overlap' in kwargs.keys():
            self.sort_by_maximum_overlap = kwargs['sort_by_maximum_overlap']
        else:
            self.sort_by_maximum_overlap = False
        
        if 'gap' in kwargs.keys():
            self.gap_threshold = kwargs['gap']
        else:
            self.gap_threshold = 0
        

    def solve(self):
        if self.eigen_cloest:
            self.eigen_eners_cloest,self.eigen_states_cloest = self.get_instantaneous_cloest_eigen()
        else:
            self.eigen_eners_ascend,self.eigen_states_ascend = self.get_instantaneous_eigens()
            if self.gap_threshold:
                for i in range(3):
                    self.rearrangement_eigen_traces_by_ignore_small_gap()

            if self.sort_by_maximum_overlap:
                self.eigen_eners_maximum_overlap = []
                self.eigen_states_maximum_overlap = []
                index_overlap = self.get_maximum_overlap_index()
                for i in range(len(self.eigen_eners_ascend)):
                    self.eigen_eners_maximum_overlap.append(self.eigen_eners_ascend[i][index_overlap])
                    self.eigen_states_maximum_overlap.append(self.eigen_states_ascend[i][index_overlap])

    def postprocess(self):
        if self.eigen_cloest:
            self.result = (self.eigen_eners_cloest,self.eigen_states_cloest)
        elif self.sort_by_maximum_overlap:
            self.result = (np.asarray(self.eigen_eners_maximum_overlap),np.asarray(self.eigen_states_maximum_overlap))
        else:
            self.result = (self.eigen_eners_ascend,self.eigen_states_ascend)

    def get_instantaneous_eigens(self):
        eigenenergies=[]
        eigenstates=[]
        # index_overlap = self.get_maximum_overlap_index()
        for t in range(len(self.tlist)):
            H = self.H(t)
            eigen_ener,eigen_state=eig(H.full())
            eigen_state=eigen_state.T
            index=np.argsort(eigen_ener)
            # if sort_by_maximum_overlap:
            #     eigenenergies.append(eigen_ener[index][index_overlap])
            #     eigenstates.append(eigen_state[index][index_overlap])
            # else:
            eigenenergies.append(eigen_ener[index])
            eigenstates.append(eigen_state[index])
        return np.real(np.asarray(eigenenergies))/2/np.pi,np.asarray(eigenstates)

    def get_instantaneous_cloest_eigen(self):
        eigenenergies=[]
        eigenstates=[]
        # index_overlap = self.get_maximum_overlap_index()
        for t in range(len(self.tlist)):
            H = self.H(t)
            eigen_ener,eigen_state = eigensolve_close( H.full() )
            eigenenergies.append(eigen_ener)
            eigenstates.append(eigen_state)
        return np.real(np.asarray(eigenenergies))/2/np.pi,np.asarray(eigenstates)

    def get_maximum_overlap_index(self):
        ## be careful using this function, it may fail in degenerate case !!!!
        H = self.H(0)
        eigenvalues = eigensolve_close(H.full())[0]
        position_index = np.argsort(eigenvalues)
        return np.argsort(position_index)

    def rearrangement_eigen_traces_by_ignore_small_gap(self):
        for i in range(len(self.eigen_eners_ascend[0])-6):
            for k in range(1,5):
                self.swap_two_eigen_trace(self.eigen_eners_ascend[:,i],self.eigen_eners_ascend[:,i+k],self.eigen_states_ascend[:,i],self.eigen_states_ascend[:,i+k],self.gap_threshold )

    def swap_two_eigen_trace(self,eigen_ener1,eigen_ener2,eigen_state1,eigen_state2,gap):
        ener_diff = eigen_ener2 - eigen_ener1
        anticross_idx = np.where( np.abs(ener_diff) < gap )[0]
        if len(anticross_idx) == 0 or isinstance(ener_diff,float):
            # return eigen_ener1,eigen_ener2,eigen_state1,eigen_state2
            pass
        else:
            extreme_points  = self.get_extreme_points(ener_diff,anticross_idx)
            for point in extreme_points:
                eigen_ener1_temp = copy.deepcopy(eigen_ener1)
                eigen_state1_temp = copy.deepcopy(eigen_state1)
                eigen_ener1[point:] = eigen_ener2[point:]
                eigen_ener2[point:] = eigen_ener1_temp[point:]
                eigen_state1[point:] = eigen_state2[point:]
                eigen_state2[point:] = eigen_state1_temp[point:]
            # return eigen_ener1,eigen_ener2,eigen_state1,eigen_state2
            pass

    def get_extreme_points(self,ener_diff,anticross_idx):
        # print(anticross_idx)
        start_idxs = [anticross_idx[0]]
        end_idxs = []
        for idx_count,idx in enumerate(anticross_idx):
            if idx+1 in anticross_idx:
                continue
            else:
                end_idxs.append(idx)
                if idx_count != len(anticross_idx)-1:
                    start_idxs.append(anticross_idx[idx_count+1])
        # print(start_idxs,end_idxs)

        extreme_points = []
        for i in range(len(start_idxs)):
            if ener_diff[start_idxs[i]]*ener_diff[end_idxs[i]] > 0:
                if start_idxs[i] == end_idxs[i]:
                    extreme_points.append(start_idxs[i])
                else:
                    extreme_points.append( np.argmin(ener_diff[start_idxs[i]:end_idxs[i]])+start_idxs[i] )    
        return extreme_points



class UnitaryOperator(Solver):

    def __init__(self,H_Mesolve,H_collapse,tlist,init_state,**kwargs):
        super().__init__(H_Mesolve,H_collapse,tlist,init_state,**kwargs)
        if 'save_trace' in kwargs.keys():
            self.save_trace = kwargs['save_trace']
        else:
            self.save_trace = False

    def solve(self):
        dt = self.tlist[1] -self.tlist[0]
        self.H_trace = []
        self.U_trace_raw = []
        for index in range(len(self.tlist)):
            H = self.H(index)
            if index == 0:
                U = (-1j * H * 0).expm()
            else:
                U = (-1j * H * dt).expm() * U
            if self.save_trace:
                self.H_trace.append(H)
                self.U_trace_raw.append(U)
        if not self.save_trace:
            self.H_trace.append(H)
            self.U_trace_raw.append(U)

    def postprocess(self):
        self.U_trace=[]
        for U in self.U_trace_raw:

            self.U_trace.append( self.U_lab_to_logic*U*self.U_lab_to_logic.dag() )
        self.result = self.U_trace

    # def solve_eigen(self):
    #     self.eigen_ener,self.eigen_state = self.get_instantaneous_eigens()

    # def get_instantaneous_eigens(self):
    #     eigenenergies=[]
    #     eigenstates=[]
    #     for H in range(len(self.H_trace)):
    #         # H = self.H(t)
    #         eigen_ener,eigen_state=eig(H.full())
    #         eigen_state=eigen_state.T
    #         index=np.argsort(eigen_ener)
    #         eigenenergies.append(eigen_ener[index])
    #         eigenstates.append(eigen_state[index])
    #     return np.asarray(eigenenergies)/2/np.pi , np.asarray(eigenstates)



def qutip_fun0(t,args):
    index=0
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))

def qutip_fun1(t,args):
    index=1
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))
def qutip_fun2(t,args):
    index=2
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))
def qutip_fun3(t,args):
    index=3
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))
def qutip_fun4(t,args):
    index=4
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))
def qutip_fun5(t,args):
    index=5
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))

def qutip_fun6(t,args):
    index=6
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))

def qutip_fun7(t,args):
    index=7
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))

def qutip_fun8(t,args):
    index=8
    if isinstance(args[f'omega_{index}'],(complex)):
        return (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*np.real(args[f'omega_{index}'])*t)
    else:
        return np.real( (interpolate.splev(t,args[f'I_{index}']) +1j * interpolate.splev(t,args[f'Q_{index}']) )* np.exp(1j*args[f'omega_{index}']*t))



qutip_func_list = [qutip_fun0,qutip_fun1,qutip_fun2,qutip_fun3,qutip_fun4,qutip_fun5,qutip_fun6,qutip_fun7,qutip_fun8]

class QutipMeSolver(Solver):

    def __init__(self,H_Mesolve,H_collapse,tlist,init_state,**kwargs):
        super().__init__(H_Mesolve,H_collapse,tlist,init_state,**kwargs)
        self.func_args={}
        self.Hamiltonian_reform()

        qutip_atol = 1e-8 if 'atol' not in kwargs.keys() else kwargs['atol'] 
        qutip_rtol = 1e-6 if 'rtol' not in kwargs.keys() else kwargs['rtol'] 
        qutip_nsteps = 1000 if 'nsteps' not in kwargs.keys() else kwargs['nsteps'] 
        self.qutip_solver_options = qt.Options(atol = qutip_atol,rtol= qutip_rtol,nsteps=qutip_nsteps)

    def solve(self):
        self.raw_result = qt.mesolve(H=self.H_Mesolve,rho0 = self.init_state, tlist =self.tlist, c_ops = self.H_collapse, args=self.func_args,options=self.qutip_solver_options)

    def Hamiltonian_reform(self):
        i=0
        for H in self.H_Mesolve:
            if isinstance(H,list):
                if isinstance(H[1],tuple):
                    I_spline = get_Interpolate_pulse(self.tlist,H[1][0].real,new_tlist=None)
                    Q_spline = get_Interpolate_pulse(self.tlist,H[1][0].imag,new_tlist=None)
                    self.func_args.update({f'I_{i}':I_spline,f'Q_{i}':Q_spline, f'omega_{i}':H[1][1]})
                    H[1] = qutip_func_list[i]
                    i += 1

    def postprocess(self):
        self.result=[]
        for k,t in enumerate(self.tlist):
            state_lab_rot = self.raw_result.states[k]
            if not state_lab_rot.isoper:
                rho_lab_rot = state_lab_rot * state_lab_rot.dag()
            else: 
                rho_lab_rot = state_lab_rot
            rho_logic_rot = T(rho_lab_rot, self.U_lab_to_logic)
            rho_logic = T(rho_logic_rot, U(self.H_idle_diag,t).dag() )
            self.result.append(rho_logic)



