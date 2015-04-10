from os import listdir

__all__ = [pkg for pkg in listdir(__path__[0]) if "." not in pkg]

try:
    for pkg in __all__:
        __import__(".".join([__name__, pkg]))
    del pkg, listdir

except ImportError:
    # this is repeated in each __init__ to avoid circular calls to functions
    from sys import platform
    from os import path

    missing = set()
    for pkg in __all__:
        slash = "\\" if 'win' in platform else "/"
        full_pkg = slash.join([__path__[0], pkg])
        if not path.exists(slash.join([full_pkg, "__init__.py"])):
            missing.add(full_pkg)
    exit('Missing __init__.py in:\n  {:s}'.format("\n  ".join(missing)))
