#!/usr/bin/env python


from setuptools import setup, Extension
import numpy

PACKAGE_NAME = "embrace_eep"

setup (
    ext_modules = [
        Extension ( name = 'embrace_eep.raweep',
          sources = [ 'src/embrace_eep/src/raweep_module.c' ] ) ],

    package_data = { "": [ "*.h"] } )
