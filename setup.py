from   setuptools import find_packages, setup

from   fuchsia import __author__, __author_email__, __version__

setup(
    author = __author__,
    author_email = __author_email__,
    name = "fuchsia",
    packages = find_packages(),
    version = __version__,
)
