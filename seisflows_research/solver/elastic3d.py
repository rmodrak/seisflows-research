
from seisflows.tools.config import loadclass

class elastic3d(loadclass('solver', 'elastic'), loadclass('solver', 'specfem3d')):
    """ Adds elastic inversion machinery to SPECFEM2D
    """
    pass

