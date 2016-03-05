
from glob import glob
from os.path import join
import subprocess

from seisflows.tools import unix
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import preprocess


class specfem3d_legacy(custom_import('solver', 'specfem3d')):

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


    ### workaround legacy code issues

    def export_kernels(self, path):
        # export kernels
        unix.mkdir_gpfs(join(path, 'kernels'))
        unix.mkdir(join(path, 'kernels', basename(self.path)))
        src = join(glob(self.model_databases +'/'+ '*kernel.bin'))
        dst = join(path, 'kernels', basename(self.path))
        unix.mv(src, dst)


    def combine(self, path=''):
        """ Sums individual source contributions. Wrapper over xcombine_sem
            utility.
        """
        super(specfem3d_legacy, self).combine(path)

        try:
            files = glob(self.model_databases +'/'+ '*alpha*_kernel.bin')
            unix.rename('alpha', 'vp', files)

            files = glob(self.model_databases +'/'+ '*beta*_kernel.bin')
            unix.rename('beta', 'vs', files)
        except:
            pass

