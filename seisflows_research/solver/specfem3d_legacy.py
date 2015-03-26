import subprocess
import glob

from seisflows.tools import unix
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class specfem3d_legacy(loadclass('solver', 'specfem3d')):

    def mpirun(self, runfile, args='', outfile='/dev/null'):
        """ Wrapper for mpirun
        """
        unix.cd('bin')

        with open(outfile) as f:
            subprocess.call(
                system.mpiargs() +
                unix.basename(runfile) +
                args,
                shell=True,
                stdout=f)
        unix.cd('..')

    def combine(self, path=''):
        """ combines SPECFEM3D kernels
        """
        unix.cd(self.getpath)

        # create temporary files and directories needed by xsum_kernels
        dirs = unix.ls(path)
        with open('../kernels_list.txt', 'w') as file:
            file.write('\n'.join(dirs) + '\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path + '/' + dir
            dst = unix.pwd() + '/' + 'INPUT_KERNELS' + '/' + dir
            unix.ln(src, dst)

        # sum kernels
        self.mpirun(PATH.SPECFEM_BIN + '/' + 'xsum_kernels')
        unix.mv('OUTPUT_SUM', path + '/' + 'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('../kernels_list.txt')

        unix.cd(path)

