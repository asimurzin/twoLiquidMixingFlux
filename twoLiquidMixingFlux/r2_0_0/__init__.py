#!/usr/bin/env python

#----------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey Petrov, Andrey Simurzin
##


#----------------------------------------------------------------------------
from Foam import ref, man


#----------------------------------------------------------------------------
def _createFields( runTime, mesh, g ):
    ref.ext_Info() << "Reading field p_rgh\n" << ref.nl
    p_rgh = man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.MUST_READ,
                                              ref.IOobject.AUTO_WRITE ),
                                mesh )

    ref.ext_Info() << "Reading field alpha1\n" << ref.nl
    alpha1 = man.volScalarField( man.IOobject( ref.word( "alpha1" ),
                                               ref.fileName( runTime.timeName() ),
                                               mesh,
                                               ref.IOobject.MUST_READ,
                                               ref.IOobject.AUTO_WRITE ),
                                 mesh )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    phi = man.createPhi( runTime, mesh, U )

    ref.ext_Info() << "Reading transportProperties\n" << ref.nl
    twoPhaseProperties = man.twoPhaseMixture( U, phi )

    rho1 = twoPhaseProperties.rho1()
    rho2 = twoPhaseProperties.rho2()

    Dab = ref.dimensionedScalar( twoPhaseProperties.lookup( ref.word( "Dab" ) ) )

    # Read the reciprocal of the turbulent Schmidt number
    alphatab = ref.dimensionedScalar( twoPhaseProperties.lookup( ref.word( "alphatab" ) ) )

    # Need to store rho for ddt(rho, U)
    
    rho = man.volScalarField ( ref.word( "rho" ), alpha1 * rho1 + ( ref.scalar(1) - alpha1 ) * rho2 )
    rho.oldTime()

    # Mass flux
    # Initialisation does not matter because rhoPhi is reset after the
    # alpha1 solution before it is used in the U equation.
    rhoPhi = man.surfaceScalarField( man.IOobject( ref.word( "rho*phi" ),
                                                   ref.fileName( runTime.timeName() ),
                                                   mesh,
                                                   ref.IOobject.NO_READ,
                                                   ref.IOobject.NO_WRITE ),
                                     rho1 * phi )

    # Construct incompressible turbulence model
    turbulence = man.incompressible.turbulenceModel.New( U, phi, twoPhaseProperties )

    ref.ext_Info() << "Calculating field g.h\n" << ref.nl
    gh = man.volScalarField ( ref.word( "gh" ), g & man.volVectorField( mesh.C(), man.Deps( mesh ) ) )
    ghf = man.surfaceScalarField ( ref.word( "ghf" ), g & man.surfaceVectorField( mesh.Cf(), man.Deps( mesh ) ) )

    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.NO_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            p_rgh + rho * gh )

    pRefCell = 0
    pRefValue = 0.0
    pRefCell, pRefValue = ref.setRefCell( p, p_rgh, mesh.solutionDict().subDict( ref.word( "PIMPLE" ) ), pRefCell, pRefValue )

    if p_rgh.needReference():
        p += ref.dimensionedScalar( ref.word( "p" ),
                                    p.dimensions(), 
                                    pRefValue - ref.getRefCellValue( p, pRefCell ) )
        p_rgh << p - rho * gh
        pass
     

    return p_rgh, alpha1, U, phi, twoPhaseProperties, rho1, rho2, Dab, alphatab, rho, rhoPhi, turbulence, gh, ghf, p, pRefCell, pRefValue
    

#--------------------------------------------------------------------------------------
def alphaEqn( mesh, phi, alpha1, alphatab, Dab, rhoPhi, rho, rho1, rho2, turbulence ):
    alpha1Eqn = ref.fvm.ddt( alpha1 ) + ref.fvm.div( phi, alpha1 ) - ref.fvm.laplacian( Dab + alphatab * turbulence.ext_nut(), alpha1, ref.word( "laplacian(Dab,alpha1)" ) )

    alpha1Eqn.solve()

    rhoPhi << alpha1Eqn.flux() * ( rho1 - rho2 ) + phi * rho2

    rho << alpha1 * rho1 + ( ref.scalar( 1 ) - alpha1 ) * rho2
    
    ref.ext_Info() << "Phase 1 volume fraction = " << alpha1.weightedAverage( mesh.V() ).value()  \
               << "  Min(alpha1) = " << alpha1.ext_min().value() << "  Max(alpha1) = " << alpha1.ext_max().value() << ref.nl
    
    pass

#----------------------------------------------------------------------------------------
def fun_UEqn( mesh, U, p_rgh, ghf, rho, rhoPhi, turbulence, twoPhaseProperties, pimple ):
    muEff = man.surfaceScalarField( ref.word( "muEff" ),
                                    man.surfaceScalarField( twoPhaseProperties.muf(), man.Deps( twoPhaseProperties ) )
                                    + man.fvc.interpolate( rho * man.volScalarField( turbulence.ext_nut(), man.Deps( turbulence ) ) ) )

    UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( rhoPhi, U ) - man.fvm.laplacian( muEff, U ) - ( man.fvc.grad( U ) & man.fvc.grad( muEff ) )

    UEqn.relax()
        
    if pimple.momentumPredictor():
        ref.solve( UEqn  == ref.fvc.reconstruct( ( - ghf * ref.fvc.snGrad( rho ) - ref.fvc.snGrad( p_rgh ) ) * mesh.magSf() ) )
        pass

    return UEqn


#--------------------------------------------------------------------------------------
def fun_pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, rho, pimple, corr, pRefCell, pRefValue, cumulativeContErr ):
    rAU = 1.0/UEqn.A()
     
    rAUf = ref.fvc.interpolate( rAU )
    
    U << rAU * UEqn.H()

    phiU = ref.surfaceScalarField( ref.word( "phiU" ),
                                   ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) )
                               
    ref.adjustPhi( phiU, U, p )
    
    phi << phiU - ghf * ref.fvc.snGrad( rho ) * rAUf * mesh.magSf()

    for nonOrth in range( pimple.nNonOrthCorr() + 1 ):
        p_rghEqn = ref.fvm.laplacian( rAUf, p_rgh ) == ref.fvc.div( phi ) 
        p_rghEqn.setReference(pRefCell, ref.getRefCellValue( p_rgh, pRefCell ) )

        p_rghEqn.solve( mesh.solver( p_rgh.select( pimple.finalInnerIter( corr, nonOrth ) ) ) )
        
        if nonOrth == pimple.nNonOrthCorr():
           phi -= p_rghEqn.flux()
           pass
        pass
    
    U += rAU * ref.fvc.reconstruct( ( phi() - phiU ) / rAUf ) # mixed calculations
    U.correctBoundaryConditions()

    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )
    
    p == p_rgh + rho * gh

    if p_rgh.needReference():
       p += ref.dimensionedScalar( ref.word( "p" ),
                                   p.dimensions(),
                                   pRefValue - ref.getRefCellValue( p, pRefCell ) )
       p_rgh << p - rho * gh
       pass
    
    return cumulativeContErr
    
    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    g = man.readGravitationalAcceleration( runTime, mesh )

    cumulativeContErr = ref.initContinuityErrs()
    
    p_rgh, alpha1, U, phi, twoPhaseProperties, rho1, rho2, Dab, \
    alphatab, rho, rhoPhi, turbulence, gh, ghf, p, pRefCell, pRefValue = _createFields( runTime, mesh, g )

    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
    
    CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
    
    runTime = ref.setInitialDeltaT( runTime, adjustTimeStep, maxCo, CoNum )
    
    pimple = man.pimpleControl( mesh )
    
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl

    while runTime.run():
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
        CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
        
        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )

        runTime.increment()
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        # --- Pressure-velocity PIMPLE corrector loop
        pimple.start()
        while pimple.loop():

            twoPhaseProperties.correct()

            alphaEqn( mesh, phi, alpha1, alphatab, Dab, rhoPhi, rho, rho1, rho2, turbulence )
            
            UEqn = fun_UEqn( mesh, U, p_rgh, ghf, rho, rhoPhi, turbulence, twoPhaseProperties, pimple )
            
            # --- PISO loop
            for corr in range( pimple.nCorr() ):
                cumulativeContErr = fun_pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, rho, pimple, \
                                              corr, pRefCell, pRefValue, cumulativeContErr )
                pass

            if pimple.turbCorr():
                turbulence.correct()
                pass
            
            pimple.increment()    
            pass

        runTime.write()
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl

        pass
    
    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher \n "     
   pass


#--------------------------------------------------------------------------------------
