
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.shared import getpar, setpar

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterError, ParameterObj

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

        if PAR.FORMAT != 'ascii_specfem3d':
            raise Exception ('Currently, CUBIT models must be used in combination with ASCII data ouput format.')

        if PAR.WORKFLOW == 'inversion':
            raise Exception ('The machinery in SPECFEM3D for iteratively updating CUBIT models is in disrepair. As a result, it is only possible to run migration, not inversions.')


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

        print 'DEBUGGING INFO'
        print 'model_name:', model_name
        print 'model_path:', model_path
        print ''

        self.initialize_solver_directories()
        unix.touch('OUTPUT_FILES/values_from_mesh.txt')
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


    def load(self):
        raise Exception('Not implemented for CUBIT models.')


    def save(self):
        raise Exception('Not implemented for CUBIT models.')


    def import_model(self, path):
        src = join(path, 'model')
        dst = self.getpath +'/'+ 'OUTPUT_FILES'
        unix.cp(src, dst)

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/*.sem*')

