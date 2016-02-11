import logging
import os.path as path

import sage.all

logging.basicConfig(
    format='\033[32m%(levelname)s [%(asctime)s]\033[0m %(message)s',
#    format='\033[32m%(levelname)s [%(asctime)s]\033[0m %(pathname)s:%(funcName)s, line %(lineno)s %(message)s\n',
#    format='%(levelname)s [%(asctime)s] %(message)s %(pathname)s:%(funcName)s, line %(lineno)s',
    datefmt='%Y-%m-%d %I:%M:%S',
    level=logging.INFO,
)

__author__ = "Oleksandr Gituliar"
__author_email__ = "oleksandr.gituliar@ifj.edu.pl"
__version__ = open(path.join(path.dirname(__file__),"VERSION")).readline().rstrip('\n')
