#!/usr/bin/env python
###################################################################
#                                                                 #
#   Setup file for compiling and installing the python wrapper    #
#   for dasslc                                                    #
#                                                                 #
#   Author: Ataide Neto                                           #
#   email: ataide@peq.coppe.ufrj.br                               #
#   Universidade Federal do Rio de Janeiro                        #
#   Version: 0.1-6                                                #
#                                                                 #
###################################################################

from setuptools import setup, Extension
import numpy
import os

if os.path.exists(os.path.join("sparse","lib","sparse.a")):
    extraObjects = [os.path.join("sparse","lib","sparse.a")]
    defineMacros = [('SPARSE', None)]
else:
    extraObjects = []
    defineMacros = []

setup(
    name = 'Dasslc2py',
    version = '0.1-6',
    license = 'MIT',
    author = 'ataide@peq.coppe.ufrj.br',
    ext_modules=[
        Extension("dasslc",
            sources = ["dasslcmodule.c",os.path.join("dasslc", "dasslc.c")],
            include_dirs = [numpy.get_include()],
            define_macros = defineMacros,
            extra_objects = extraObjects,
            extra_compile_args = ["-O3"]
        )
    ]
)
