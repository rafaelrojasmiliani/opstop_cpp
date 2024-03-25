try:
    import inspect
    import sys
    from .pyopstop import *  # noqa
    # ---
    from .plot import *  # noqa
except ImportError:
    import pyopstop
    from pyopstop import *  # noqa
    from .plot import *  # noqa


submodules = inspect.getmembers(pyopstop, inspect.ismodule)
for module_info in submodules:
    sys.modules['opstop.' + module_info[0]] = module_info[1]
