from pkg_resources import get_distribution

try:
    __version__ = get_distribution("simutator").version
except:
    __version__ = "local"


__all__ = ["utils"]

from simutator import *
