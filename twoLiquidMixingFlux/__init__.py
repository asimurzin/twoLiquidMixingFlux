#!/usr/bin/env python

#----------------------------------------------------------------------------
## Copyright (C) 2011 Dhasthagheer
##
## Author : Andrey Simurzin


#----------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "010701" ):
    from twoLiquidMixingFlux.r1_7_1 import *
    pass


#----------------------------------------------------------------------------
def entry_point():
    try:
       engine = main_standalone
       pass
    except NameError:
       print
       print "There is no implementation of the current OpenFOAM version"
       print
       import os; os._exit( os.EX_OK )
       pass
    
    import sys; argv = sys.argv
    return engine( len( argv ), argv )

#----------------------------------------------------------------------------
if __name__ == "__main__" :
    entry_point()
    pass
    
    
#----------------------------------------------------------------------------
