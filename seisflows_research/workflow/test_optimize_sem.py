
from os.path import abspath, exists
import os
import sys

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.array import loadnpy, savenpy
from seisflows.config import , \
    ParameterError

from seisflows.plugins.io import sem

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

import optimize


class test_optimize_sem(object):

    parameters = []
    parameters += ['vp']
    parameters += ['vs']

    def check(self):
        # check parameters
        if 'OPTIMIZE' not in PAR: 
            setattr(PAR, 'OPTIMIZE', 'base')

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH,'SCRATCH',abspath('.'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',abspath('.'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH,'OPTIMIZE',abspath('./scratch'))


    def main(self):
        wd = PATH.OPTIMIZE
        unix.mkdir(wd)

        print 'Checking mesh properties...'
        self.check_mesh_properties(PATH.MODEL)

        print 'Reading model...'
        model = self.load(PATH.MODEL)
        savenpy(wd+'/'+'m_new', self.merge(model))

        print 'Reading gradient...'
        model = self.load(PATH.GRADIENT)
        savenpy(wd+'/'+'g_new', self.merge(model))

        print 'Testing optimization machinery...'
        optimize.setup()
        optimize.compute_direction()

        print 'Finished...\n'


    def load(self, path):
        model = {}
        for key in self.parameters:
            model[key] = []
            for iproc in range(self.mesh_properties.nproc):
                model[key] += [sem.read(path, key, iproc)]
        return model


    def merge(self, model):
        """ Converts model from dictionary to vector representation
        """
        v = np.array([])
        for key in self.parameters:
            for iproc in range(self.mesh_properties.nproc):
                v = np.append(v, model[key][iproc])
        return v


    def check_mesh_properties(self, path):
            nproc = 0
            ngll = []
            while True:
                dummy = sem.read(path, self.parameters[0], nproc)
                ngll += [len(dummy)]
                nproc += 1
                if not exists('%s/proc%06d_%s.bin' % (path, nproc, self.parameters[0])):
                    break

            self.mesh_properties = Struct([
                ['nproc', nproc],
                ['ngll', ngll]])


