
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.io import load
from seisflows.seistools.shared import getpar, setpar

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class specfem3d_cubit(loadclass('solver', 'specfem3d')):
    """ For SPECFEM3D simulations with CUBIT models
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(specfem3d_cubit, self).check()


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)
        unix.cd(self.getpath)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.data_wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type=''):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(exists(model_path))

        self.initialize_solver_directories()
        unix.cd(self.getpath)
        unix.cp(glob(model_path +'/'+ '*'), self.getpath +'/'+ 'OUTPUT_FILES')
        self.export_model(PATH.OUTPUT +'/'+ model_name)


    def forward(self):
        """ Calls SPECFEM3D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        #self.mpirun('bin/xgenerate_databases')
        self.mpirun('bin/xspecfem3D')


    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        self.mpirun('bin/xspecfem3D')


    def initialize_io_machinery(self):
        """ Writes mesh files expected by input/output methods
        """
        if system.getnode() == 0:
            parts = self.load(PATH.MODEL_INIT +'/'+ 'DATABASES_MPI')
            path = PATH.GLOBAL +'/'+ 'mesh'

            if not exists(path):
                for key in self.model_parameters:
                    if key not in self.inversion_parameters:
                        unix.mkdir(path +'/'+ key)
                        for proc in range(PAR.NPROC):
                            with open(path +'/'+ key +'/'+ '%06d' % proc, 'w') as file:
                                np.save(file, parts[key][proc])

            if 'OPTIMIZE' in PATH:
                if not exists(PATH.OPTIMIZE +'/'+ 'm_new'):
                    savenpy(PATH.OPTIMIZE +'/'+ 'm_new', self.merge(parts))


