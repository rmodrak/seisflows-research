
from seisflows.tools.config import loadclass

class elastic2d(loadclass('solver', 'elastic'), loadclass('solver', 'specfem2d')):
    """ Adds elastic inversion machinery to SPECFEM2D
    """
    pass

