from   setuptools import setup

from   fuchsia import __author__, __author_email__, __version__

with open('README.md') as f:
    long_description = f.read()

setup(
    author = __author__,
    author_email = __author_email__,
    description = "A tool for reducing differential equations for multiloop master integrals",
    entry_points = {'console_scripts': ['fuchsia = fuchsia:entry_point']},
    install_requires = ['docopt', 'sage==7.0'],
    license = "ISC License",
    long_description = long_description,
    name = "fuchsia",
    py_modules = ['fuchsia'],
    test_suite = 'test',
    url = "http://www.gituliar.org/fuchsia/",
    version = __version__,
)
