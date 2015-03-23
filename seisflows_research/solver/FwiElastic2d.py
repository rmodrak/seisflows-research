
from seisflows.tools.config import loadclass

specfem2d = loadclass('solver', 'specfem2d_binary')


class FwiElastic2d(FwiElastic, specfem2d):
    """ Adds elastic inversion machinery to SPECFEM2D
    """
    pass

