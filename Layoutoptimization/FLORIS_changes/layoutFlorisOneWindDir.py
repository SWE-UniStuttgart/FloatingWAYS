# Copyright 2022 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See https://floris.readthedocs.io for documentation


import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from scipy.spatial.distance import cdist
import pdb

from .layout_optimization_base import LayoutOptimization

class LayoutOptimizationPyOptSparseOneDir(LayoutOptimization):
    def __init__(
        self,
        fi,
        boundaries,
        min_dist=None,
        freq=None,
        solver=None,
        optOptions=None,
        timeLimit=None,
        storeHistory='hist.hist',
        hotStart=None
    ):
        super().__init__(fi, boundaries, min_dist=min_dist, freq=freq)
        
        self.ymin = self.ymin-120
        self.ymax = self.ymax+120

        self.y0 = self._norm(self.fi.layout_y, self.ymin, self.ymax)
        self.x0 = self.fi.layout_x

        self.storeHistory = storeHistory
        self.timeLimit = timeLimit
        self.hotStart = hotStart

        try:
            import pyoptsparse
        except ImportError:
            err_msg = (
                "It appears you do not have pyOptSparse installed. "
                + "Please refer to https://pyoptsparse.readthedocs.io/ for "
                + "guidance on how to properly install the module."
            )
            self.logger.error(err_msg, stack_info=True)
            raise ImportError(err_msg)

        # Insantiate ptOptSparse optimization object with name and objective function
        self.optProb = pyoptsparse.Optimization('layout', self._obj_func)

        self.optProb = self.add_var_group(self.optProb)
        self.optProb = self.add_con_group(self.optProb)
        self.optProb.addObj("obj")

        if solver is not None:
            self.solver = solver
            print("Setting up optimization with user's choice of solver: ", self.solver)
        else:
            self.solver = "SLSQP"
            print("Setting up optimization with default solver: SLSQP.")
        if optOptions is not None:
            self.optOptions = optOptions
        else:
            if self.solver == "SNOPT":
                self.optOptions = {"Major optimality tolerance": 1e-7}
            else:
                self.optOptions = {}

        exec("self.opt = pyoptsparse." + self.solver + "(options=self.optOptions)")

    def _optimize(self):
        if hasattr(self, "_sens"):
            self.sol = self.opt(self.optProb, sens=self._sens)
        else:
            if self.timeLimit is not None:
                self.sol = self.opt(self.optProb, sens="CDR", storeHistory=self.storeHistory, timeLimit=self.timeLimit, hotStart=self.hotStart)
            else:
                self.sol = self.opt(self.optProb, sens="CDR", storeHistory=self.storeHistory, hotStart=self.hotStart)
        return self.sol

    def _obj_func(self, varDict):
        # Parse the variable dictionary
        self.parse_opt_vars(varDict)

        # Update turbine map with turbince locations
        self.fi.reinitialize(layout_x = self.x0, layout_y = self.y)

        # Compute the objective function
        funcs = {}

        funcs["obj"] = (
            -1 * (self.fi.get_farm_AEP(self.freq,freq_warn=0) - self.initial_AEP)/ self.initial_AEP
        )
        

        funcs = self.compute_cons(funcs)

        fail = False
        return funcs, fail


    def parse_opt_vars(self, varDict):

        self.y = self._unnorm(varDict["y"], self.ymin, self.ymax)


    def parse_sol_vars(self, sol):
        self.y = list(self._unnorm(sol.getDVs()["y"], self.ymin, self.ymax))[1]

    def add_var_group(self, optProb):

        optProb.addVarGroup(
            "y", self.nturbs, varType="c", lower=0.0, upper=1.0, value=self.y0
        )

        return optProb

    def add_con_group(self, optProb):
        optProb.addConGroup("perp_con", self.nturbs, lower=-1.0, upper=1.0)

        return optProb

    def compute_cons(self, funcs):
        funcs["perp_con"] = self.perp_constraint()

        return funcs


    
    def perp_constraint(self):
        y_initial = self._unnorm(self.y0, self.ymin, self.ymax)

        disty = [
            (self.y[i] - y_initial[i])/120
            for i in range(self.nturbs)
        ]

        return disty



    def _get_initial_and_final_locs(self):
        x_initial = self.x0
        y_initial = self._unnorm(self.y0, self.ymin, self.ymax)

        x_opt, y_opt = self.get_optimized_locs()
        return x_initial, y_initial, x_opt, y_opt

    def get_optimized_locs(self):
        x_opt = self.x0
        y_opt = self._unnorm(self.sol.getDVs()["y"], self.ymin, self.ymax)

        return x_opt,y_opt