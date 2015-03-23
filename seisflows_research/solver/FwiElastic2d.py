
from seisflows.tools.config import loadclass

class FwiElastic2d(loadclass('solver', 'FwiElastic'), loadclass('solver', 'specfem2d_binary')):
    """ Adds elastic inversion machinery to SPECFEM2D
    """
    pass

