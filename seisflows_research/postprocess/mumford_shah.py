
from glob import glob

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import call, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import
from seisflows.tools.math import nabla

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize


class mumford_shah(custom_import('postprocess', 'regularize')):

    def check(self):
        """ Checks parameters and paths
        """
        super(mumford_shah, self).check()

        # check parameters
        if 'FIXRADIUS' not in PAR:
            setattr(PAR, 'FIXRADIUS', 7.5)

        if 'GAMMA' not in PAR:
            setattr(PAR, 'GAMMA', 0.)

        if 'ETA' not in PAR:
            setattr(PAR, 'ETA', 0.)

        # check paths
        if 'MUMFORD_SHAH_INPUT' not in PATH:
            raise ParameterError

        if 'MUMFORD_SHAH_OUTPUT' not in PATH:
            raise ParameterError

        if 'MUMFORD_SHAH_CONFIG' not in PATH:
            raise ParameterError

        if 'MUMFORD_SHAH_BIN' not in PATH:
            raise ParameterError

        assert PAR.GAMMA >= 0
        assert PAR.ETA >= 0 


    def write_gradient(self, path):
        super(mumford_shah, self).write_gradient(path)

        g = solver.load(path +'/'+ 'gradient', suffix='_kernel')
        m = solver.load(path +'/'+ 'model')
        mesh = self.getmesh()

        for parameter in solver.parameters:
            for iproc in range(PAR.NPROC):
                g[parameter][iproc] += PAR.GAMMA *\
                    self.get_damping_term(parameter, iproc)

        # save gradient
        self.save(path, solver.merge(g), backup='noregularize')

        # save edges
        src = PATH.MUMFORD_SHAH_OUTPUT
        dst = PATH.OUTPUT+'/'+'mumford_shah'+('_%04d' % optimize.iter)
        unix.mv(src, dst)


    def detect_edges(self):
        from seisflows.seistools.io.sem import _read_bin

        path_input = PATH.MUMFORD_SHAH_INPUT
        path_output = PATH.MUMFORD_SHAH_OUTPUT
        path_run = PATH.SUBMIT

        unix.cp(glob(PATH.GRAD+'/'+'model/*'), path_input)
        unix.mkdir(path_output)
        unix.cd(path_run)

        # writes damping term to disk
        with open('mumford_shah.log', 'w') as fileobj:
            call(PATH.MUMFORD_SHAH_BIN+'/'+'psemimage ' +
                 PATH.MUMFORD_SHAH_CONFIG +
                 ' -ksp_type fgmres ' +
                 ' -pc_type asm ' +
                 ' -ksp_gmres_restart 300 ',
                 stdout=fileobj)


    def get_damping_term(self, parameter, iproc):
        from seisflows.seistools.io.sem import _read_bin

        # reads damping term from disk
        filename = '%s/proc%06d_%s.bin.dm' % \
            (PATH.MUMFORD_SHAH_OUTPUT, iproc, parameter)
        return _read_bin(filename)


    def process_kernels(self, path, parameters):
        """ Processes kernels in accordance with parameter settings
        """
        fullpath = path +'/'+ 'kernels'
        assert exists(path)

        if exists(fullpath +'/'+ 'sum'):
            unix.mv(fullpath +'/'+ 'sum', fullpath +'/'+ 'sum_nofix')

        system.run('postprocess', 'fix_near_field', 
                   hosts='all', 
                   path=fullpath)

        system.run('solver', 'combine',
                   hosts='head',
                   path=fullpath,
                   parameters=parameters)

        system.run('postprocess', 'detect_edges',
                   hosts='head')


