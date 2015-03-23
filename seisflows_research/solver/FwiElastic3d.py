
from seisflows.tools.config import loadclass

class FwiElastic3d(loadclass('solver', 'FwiElastic'), loadclass('solver', 'specfem3d')):
    """ Adds elastic inversion machinery to SPECFEM2D
    """
    pass

