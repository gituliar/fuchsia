from   setuptools import setup

from   fuchsia import __author__, __author_email__, __version__

with open('README') as f:
    long_description = f.read()

setup(
    author = __author__,
    author_email = __author_email__,
    description = "A tool for reducing differential equations for Feynman master integrals",
    download_url = "http://gituliar.org/fuchsia/fuchsia-"+__version__+".tar.gz",
    entry_points = {'console_scripts': ['fuchsia = fuchsia:main']},
    install_requires = ['sage>=7.0'],
    license = "ISC License",
    long_description = long_description,
    name = "fuchsia",
    py_modules = ['fuchsia'],
    test_suite = 'test.full_test_suite',
    url = "http://www.gituliar.org/fuchsia/",
    version = __version__,
)
