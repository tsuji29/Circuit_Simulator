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

import copy

class Simulator():

    def __init__(self,config,rwa = False,rwa_freq = 6.0e9,omit_freq = 3e9,maximum_points=0):
        self.Simu = SysH(config,rwa,rwa_freq,omit_freq,maximum_points)
        self.solvers = {0:'qutip mesolver',1:'unitary operator',2:'eigen solver'}
        # print_tips('supported solvers:',self.solvers)
                
        self.solver_options = copy.deepcopy(self.Simu.solver_options)

    def show_pulse(self):
        self.Simu.show_channels()

    def performsimulation(self,solver_type=0,resample_factor=1,save_trace = False,sort_by_maximum_overlap= True,**other_kwargs):
        self.solver_type = solver_type
        self.solver_options.update(other_kwargs)
        self.solver_options.update({'resample_factor':resample_factor})
        if solver_type in [0,'qutip mesolver']:
            # print_tips(solver_options)
            self.StateSolver = QutipMeSolver( self.Simu.H_Mesolve,self.Simu.H_collapse,self.Simu.tlist,self.Simu.init_state,**self.solver_options )
            self.StateSolver.solve()
            self.StateSolver.postprocess()
            self.StateResult = QutipStateParser(self.StateSolver.result,self.StateSolver.tlist,list(map(int,self.Simu.nTrunc)))
        elif solver_type in [1,'unitary operator']:
            self.solver_options.update({'save_trace':save_trace})
            # print_tips(solver_options)
            self.UnitarySolver = UnitaryOperator(self.Simu.H_Mesolve,self.Simu.H_collapse,self.Simu.tlist,self.Simu.init_state,**self.solver_options)
            self.UnitarySolver.solve()
            self.UnitarySolver.postprocess()
            self.UnitaryResult = UnitaryOperatorParser(self.UnitarySolver.result,list(map(int,self.Simu.nTrunc)))
        elif solver_type in [2,'eigen solver']:
            self.solver_options.update({'sort_by_maximum_overlap':sort_by_maximum_overlap})
            # print_tips(solver_options)
            self.EigenS= EigenSolver(self.Simu.H_Mesolve,self.Simu.H_collapse,self.Simu.tlist,self.Simu.init_state,**self.solver_options)
            self.EigenS.solve()
            self.EigenS.postprocess()
            self.EigenResult = EigenStateParser(self.EigenS.result,self.EigenS.tlist,list(map(int,self.Simu.nTrunc)))


