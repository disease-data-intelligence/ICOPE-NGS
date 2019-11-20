#! /usr/bin/env python

import types

def imports(global_modules):
    for name, val in global_modules.items():
        if isinstance(val, types.ModuleType):
            yield val


def print_modules(modules):
    print("# Loaded modules:\nModule \tVersion")
    for x in modules: 
        try:
          print(x.__name__, "\t", x.__version__)
        except: 
           print(x.__name__)

