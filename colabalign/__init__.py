'''ColabAlign: all-vs-all US-align structural comparison, UPGMA dendrograms,
tree clustering and cluster-representative selection.'''

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('colabalign')
except PackageNotFoundError:  # running from a source checkout without installing
    __version__ = '0+unknown'

__all__ = ['ColabAlign', 'main', '__version__']


def __getattr__(name):
    # Lazy import keeps `import colabalign` (and `colabalign --version`) light.
    if name in ('ColabAlign', 'main'):
        from . import pipeline  # pylint: disable=import-outside-toplevel
        return getattr(pipeline, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
