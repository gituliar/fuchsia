from   setuptools import setup

from   fuchsia import __author__, __author_email__, __version__

with open('README') as f:
    long_description = f.read()

setup(
    author = __author__,
    author_email = __author_email__,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
    ],
    description = "Fuchsia reduces differential equations for Feynman master integrals to canonical form",
    download_url = "http://gituliar.net/fuchsia/fuchsia-"+__version__+".tar.gz",
    entry_points = {'console_scripts': ['fuchsia = fuchsia:main']},
    install_requires = ['sage>=7.0'],
    keywords = "differential equations, Feynman integrals",
    license = "ISC License",
    long_description = long_description,
    name = "fuchsia",
    platforms = ["any"],
    py_modules = ['fuchsia'],
    test_suite = 'test.test_suite_maxima',
    url = "http://www.gituliar.net/fuchsia/",
    version = __version__,
)
