
from os.path import join
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.optimize.lib.LBFGS import LBFGS
from seisflows.optimize.lib.NLCG import NLCG
from seisflows.optimize.lib.io import Writer, StepWriter


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import solver


def fix(A):
    nrow = A.shape[0]
    ncol = A.shape[1]
    for i in range(ncol):
       if A[i,i] < 0.:
           A[:,i] *= -1
    return A


def mat2str(A):
    nrow = A.shape[0]
    ncol = A.shape[1]

    lines = []
    for ir in range(nrow):
        line = ''
        for ic in range(ncol):
            line += '%10.3e ' % A[ir,ic]
        lines += [line + '\n']
    return lines


class pca(loadclass('optimize', 'base')):

    def check(self):
        super(pca, self).check()

        if 'PCADAMP' not in PAR:
            setattr(PAR, 'PCADAMP', 0.)


    def setup(self):
        """ Sets up nonlinear optimization machinery
        """
        unix.mkdir(PATH.OPTIMIZE)

        # prepare output writers
        self.writer = Writer(
                path=PATH.OUTPUT)

        self.stepwriter = StepWriter(
                path=PATH.SUBMIT)

        # prepare algorithm machinery
        if PAR.SCHEME in ['NLCG']:
            self.NLCG = NLCG(
                path=PATH.OPTIMIZE,
                maxiter=PAR.NLCGMAX,
                thresh=PAR.NLCGTHRESH,
                precond=self.precond,
                load=self.load,
                save=self.save)

        elif PAR.SCHEME in ['LBFGS']:
            self.LBFGS = LBFGS(
                path=PATH.OPTIMIZE,
                memory=PAR.LBFGSMEM,
                maxiter=PAR.LBFGSMAX,
                thresh=PAR.LBFGSTHRESH,
                precond=self.precond)


        # write initial model
        if exists(PATH.MODEL_INIT):
            src = PATH.MODEL_INIT
            dst = join(PATH.OPTIMIZE, 'm_new')
            savenpy(dst, solver.merge(solver.load(src)))


    def load(self, filename):
        W = self.change_of_variables()

        v = loadnpy(filename)

        with open(join(PATH.SUBMIT, 'output.pca'), 'a') as f:
            f.write(filename)
            f.write('\n')
            f.writelines(mat2str(W))
            f.write('\n')
            f.write('\n')

        old = solver.split(v)
        new = {}

        # convert from original to pca parameters
        for ii,ikey in enumerate(solver.parameters):
            new[ikey] = []

            # initialize with zeros
            for iproc in range(solver.mesh.nproc):
                ngll = solver.mesh.ngll[iproc]
                new[ikey] += [np.zeros(ngll)]

            # change of variables
            for iproc in range(solver.mesh.nproc):
                for jj,jkey in enumerate(solver.parameters):
                        new[ikey][iproc] += W[ii,jj]*old[jkey][iproc]

        return solver.merge(new)


    def save(self, filename, v):
        W = np.linalg.inv(self.change_of_variables())

        # convert from pca to original parameters
        new = solver.split(v)
        old = {}

        for ii,ikey in enumerate(solver.parameters):
            old[ikey] = []

            # initialize with zeros
            for iproc in range(solver.mesh.nproc):
                ngll = solver.mesh.ngll[iproc]
                old[ikey] += [np.zeros(ngll)]

            # change of variables
            for iproc in range(solver.mesh.nproc):
                for jj,jkey in enumerate(solver.parameters):
                            old[ikey][iproc] += W[ii,jj]*new[jkey][iproc]

        savenpy(filename, solver.merge(old))



    def covariance(self):
        path = join(PATH.GRAD, 'gradient')
        g = solver.load(path, suffix='_kernel')

        # compute covariance
        n = len(solver.parameters)
        cov = np.zeros((n,n))
        for ii,ikey in enumerate(solver.parameters):
            for jj,jkey in enumerate(solver.parameters):
                for iproc in range(solver.mesh.nproc):
                    cov[ii,jj] += np.dot(g[ikey][iproc], g[jkey][iproc])

        return cov


    def change_of_variables(self):
        n = len(solver.parameters)
        I = np.eye(n,n)

        C = self.covariance()
        d,E = np.linalg.eig(C)
        E = fix(E)
        D = np.diag(d)
        W = np.dot(D, E)
        W += PAR.PCADAMP * np.diag(W).max() * I
        W /= sum(np.diag(W))/n

        return W

