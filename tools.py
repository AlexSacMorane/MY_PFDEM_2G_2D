# -*- encoding=utf-8 -*-

from pathlib import Path
import shutil, os

#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_folder(name):
    '''
    Create a new folder. If it already exists, it is erased.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)
