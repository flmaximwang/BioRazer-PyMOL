from .calculator import *
from .fitter import *
from .selectors import *
from .transformers import *
from .visualizers import *

ADVANCE_MODULE_WARNING = "WARNING: Advanced modules are not available in the default PyMOL installation. Some packages used in advanced modules are missing in the default PyMOL, please install pymol and related packages with conda. This doesn't affect normal usage of normal modules."

try:
    # from .analyzers_advanced import *
    from .visualizers_advanced import *
except ImportError as e:
    print(f"ADVANCE MODULE IMPORT ERROR: {e}")
    print(ADVANCE_MODULE_WARNING)