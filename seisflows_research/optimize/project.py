
from os.path import join
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.tools.math import angle, polyfit2, backtrack2
from seisflows.optimize.lib.LBFGS import LBFGS
from seisflows.optimize.lib.NLCG import NLCG
from seisflows.optimize.lib.io import Writer, StepWriter


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import solver


class project(loadclass('optimize', 'base')):

    def initialize_search(self):
        """ Determines initial step length for line search
        """
        unix.cd(PATH.OPTIMIZE)

        m = self.load('m_new')
        p = self.load('p_new')
        f = loadtxt('f_new')

        norm_m = max(abs(solver.merge(solver.load(PATH.MODEL_INIT))))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        # reset search history
        self.search_history = [[0., f]]
        self.step_count = 0
        self.isdone = 0
        self.isbest = 0
        self.isbrak = 0

        # determine initial step length
        if self.iter == 1:
            alpha = p_ratio*PAR.STEPINIT
        elif self.restarted:
            alpha = p_ratio*PAR.STEPINIT
        elif PAR.SCHEME in ['LBFGS']:
            alpha = 1.
        else:
            alpha = self.initial_step()

        # optional ad hoc scaling
        if PAR.STEPOVERSHOOT:
            alpha *= PAR.STEPOVERSHOOT

        # optional maximum step length safegaurd
        if PAR.STEPTHRESH:
            if alpha > p_ratio * PAR.STEPTHRESH and \
                self.iter > 1:
                alpha = p_ratio * PAR.STEPTHRESH

        # write trial model corresponding to chosen step length
        savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)

        # upate log
        self.stepwriter(steplen=0., funcval=f)


    def update_status(self):
        """ Updates line search status

            Maintains line search history by keeping track of step length and
            function value from each trial model evaluation. From line search
            history, determines whether stopping criteria have been satisfied.
        """
        unix.cd(PATH.OPTIMIZE)

        x_ = loadtxt('alpha')
        f_ = loadtxt('f_try')
        if np.isnan(f_):
            raise ValueError

        # update search history
        self.search_history += [[x_, f_]]
        self.step_count += 1
        x = self.step_lens()
        f = self.func_vals()

        fmin = f.min()
        imin = f.argmin()

        # is current step length the best so far?
        vals = self.func_vals(sort=False)
        if np.all(vals[-1] < vals[:-1]):
            self.isbest = 1

        # are stopping criteria satisfied?
        if PAR.LINESEARCH == 'Fixed':
            if (fmin < f[0]) and any(fmin < f[imin:]):
                self.isdone = 1

        elif PAR.LINESEARCH == 'Bracket' or \
            self.iter == 1 or self.restarted:
            if self.isbrak:
                self.isdone = 1
            elif (fmin < f[0]) and any(fmin < f[imin:]):
                self.isbrak = 1

        elif PAR.LINESEARCH == 'Backtrack':
            if fmin < f[0]:
                self.isdone = 1

        # update log
        self.stepwriter(steplen=x_, funcval=f_)

        return self.isdone


    def compute_step(self):
        """ Computes next trial step length
        """
        unix.cd(PATH.OPTIMIZE)

        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        s = loadtxt('s_new')

        norm_m = max(abs(solver.merge(solver.load(PATH.MODEL_INIT))))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        x = self.step_lens()
        f = self.func_vals()

        # compute trial step length
        if PAR.LINESEARCH == 'Fixed':
            alpha = p_ratio*(self.step_count + 1)*PAR.STEPINIT

        elif PAR.LINESEARCH == 'Bracket' or \
            self.iter==1 or self.restarted:
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                alpha = polyfit2(x, f)

            elif any(f[1:] <= f[0]):
                alpha = loadtxt('alpha')*PAR.STEPFACTOR**-1
            else:
                alpha = loadtxt('alpha')*PAR.STEPFACTOR

        elif PAR.LINESEARCH == 'Backtrack':
            # calculate slope along 1D profile
            slope = s/self.dot(g,g)**0.5
            if PAR.ADHOCFACTOR:
                slope *= PAR.ADHOCFACTOR            

            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)

        # write trial model corresponding to chosen step length
        savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)


