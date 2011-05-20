#!/usr/bin/env python

#----------------------------------------------------------------------------
## Copyright (C) 2011 Dhasthagheer
##
## Author : Andrey Simurzin


#----------------------------------------------------------------------------
def _createFields( runTime, mesh, g ):
    from Foam.OpenFOAM import ext_Info, nl
    from Foam.OpenFOAM import IOdictionary, IOobject, word, fileName
    from Foam.finiteVolume import volScalarField
    ext_Info() << "Reading field p_rgh\n" << nl
    p_rgh = volScalarField( IOobject( word( "p_rgh" ),
                                      fileName( runTime.timeName() ),
                                      mesh,
                                      IOobject.MUST_READ,
                                      IOobject.AUTO_WRITE ),
                            mesh )

    ext_Info() << "Reading field alpha1\n" << nl
    alpha1 = volScalarField( IOobject( word( "alpha1" ),
                                       fileName( runTime.timeName() ),
                                       mesh,
                                       IOobject.MUST_READ,
                                       IOobject.AUTO_WRITE ),
                             mesh )

    ext_Info() << "Reading field U\n" << nl
    from Foam.finiteVolume import volVectorField
    U = volVectorField( IOobject( word( "U" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    from Foam.finiteVolume.cfdTools.incompressible import createPhi
    phi = createPhi( runTime, mesh, U )

    ext_Info() << "Reading transportProperties\n" << nl
    from Foam.transportModels import twoPhaseMixture
    twoPhaseProperties = twoPhaseMixture( U, phi )

    rho1 = twoPhaseProperties.rho1()
    rho2 = twoPhaseProperties.rho2()

    from Foam.OpenFOAM import dimensionedScalar
    Dab = dimensionedScalar( twoPhaseProperties.lookup( word( "Dab" ) ) )

    # Read the reciprocal of the turbulent Schmidt number
    alphatab = dimensionedScalar( twoPhaseProperties.lookup( word( "alphatab" ) ) )

    # Need to store rho for ddt(rho, U)
    
    from Foam.OpenFOAM import scalar
    rho = volScalarField ( word( "rho" ), alpha1 * rho1 + ( scalar(1) - alpha1 ) * rho2 )
    rho.oldTime()

    # Mass flux
    # Initialisation does not matter because rhoPhi is reset after the
    # alpha1 solution before it is used in the U equation.
    from Foam.finiteVolume import surfaceScalarField
    rhoPhi = surfaceScalarField( IOobject( word( "rho*phi" ),
                                           fileName( runTime.timeName() ),
                                           mesh,
                                           IOobject.NO_READ,
                                           IOobject.NO_WRITE ),
                                 rho1 * phi )

    # Construct incompressible turbulence model
    from Foam import incompressible
    turbulence = incompressible.turbulenceModel.New( U, phi, twoPhaseProperties )

    ext_Info() << "Calculating field g.h\n" << nl
    gh = volScalarField ( word( "gh" ), g & mesh.C() )
    ghf = surfaceScalarField ( word( "ghf" ), g & mesh.Cf() )

    p = volScalarField( IOobject( word( "p" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.NO_READ,
                                  IOobject.AUTO_WRITE ),
                        p_rgh + rho * gh )

    pRefCell = 0
    pRefValue = 0.0
    from Foam.finiteVolume import setRefCell, getRefCellValue
    pRefCell, pRefValue = setRefCell( p, p_rgh, mesh.solutionDict().subDict( word( "PIMPLE" ) ), pRefCell, pRefValue )

    if p_rgh.needReference():
        p.ext_assign( p + dimensionedScalar( word( "p" ),
                                             p.dimensions(), 
                                             pRefValue - getRefCellValue( p, pRefCell ) ) )
        p_rgh.ext_assign( p - rho * gh )
        pass
     

    return p_rgh, alpha1, U, phi, twoPhaseProperties, rho1, rho2, Dab, alphatab, rho, rhoPhi, turbulence, gh, ghf, p, pRefCell, pRefValue
    

#--------------------------------------------------------------------------------------
def alphaEqn( mesh, phi, alpha1, alphatab, Dab, rhoPhi, rho, rho1, rho2, turbulence ):
    from Foam import fvm
    from Foam.OpenFOAM import word
    alpha1Eqn = fvm.ddt( alpha1 ) + fvm.div( phi, alpha1 ) - fvm.laplacian( Dab + alphatab * turbulence.ext_nut(), alpha1, word( "laplacian(Dab,alpha1)" ) )

    alpha1Eqn.solve()

    rhoPhi.ext_assign( alpha1Eqn.flux() * ( rho1 - rho2 ) + phi * rho2 )
    from Foam.OpenFOAM import scalar
    rho.ext_assign( alpha1 * rho1 + ( scalar( 1 ) - alpha1 ) * rho2 )
    
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "Phase 1 volume fraction = " << alpha1.weightedAverage( mesh.V() ).value()  \
               << "  Min(alpha1) = " << alpha1.ext_min().value() << "  Max(alpha1) = " << alpha1.ext_max().value() << nl
    
    pass

#----------------------------------------------------------------------------------------
def fun_UEqn( mesh, U, p_rgh, ghf, rho, rhoPhi, turbulence, twoPhaseProperties, momentumPredictor, finalIter ):
    from Foam.OpenFOAM import word
    from Foam.finiteVolume import surfaceScalarField
    from Foam import fvc
    muEff = surfaceScalarField( word( "muEff" ),
                                twoPhaseProperties.muf() + fvc.interpolate( rho * turbulence.ext_nut() ) )

    from Foam import fvm, fvc

    UEqn = fvm.ddt( rho, U ) + fvm.div( rhoPhi, U ) - fvm.laplacian( muEff, U ) - ( fvc.grad( U ) & fvc.grad( muEff ) )

    if finalIter:
        UEqn.relax( 1.0 )
        pass
    else:
        UEqn.relax()
        pass
        
    if momentumPredictor:
        from Foam.finiteVolume import solve
        solve( UEqn  == fvc.reconstruct( ( - ghf * fvc.snGrad( rho ) - fvc.snGrad( p_rgh ) ) * mesh.magSf() ),
               mesh.solver( U.select( finalIter ) ) )
        pass

    return UEqn


#--------------------------------------------------------------------------------------
def fun_pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, rho, finalIter, corr, nCorr, nNonOrthCorr, pRefCell, pRefValue, cumulativeContErr ):
    rAU = 1.0/UEqn.A()
     
    from Foam import fvc
    rAUf = fvc.interpolate( rAU )
    
    U.ext_assign( rAU * UEqn.H() )
    
    from Foam.finiteVolume import surfaceScalarField
    from Foam.OpenFOAM import word
    phiU = surfaceScalarField( word( "phiU" ),
                               ( fvc.interpolate( U ) & mesh.Sf() ) + fvc.ddtPhiCorr( rAU, rho, U, phi ) )
                               
    from Foam.finiteVolume import adjustPhi
    adjustPhi(phiU, U, p)
    
    phi.ext_assign( phiU - ghf * fvc.snGrad( rho ) * rAUf * mesh.magSf() )

    from Foam import fvm
    for nonOrth in range( nNonOrthCorr + 1 ):
        p_rghEqn = fvm.laplacian( rAUf, p_rgh ) == fvc.div( phi ) 
        p_rghEqn.setReference( pRefCell, pRefValue )

        p_rghEqn.solve( mesh.solver( p_rgh.select( finalIter and corr == nCorr-1 and nonOrth == nNonOrthCorr) ) )
        
        if nonOrth == nNonOrthCorr:
           phi.ext_assign( phi - p_rghEqn.flux() )
           pass
        pass
    
    U.ext_assign( U + rAU * fvc.reconstruct( ( phi - phiU ) / rAUf ) )
    U.correctBoundaryConditions()

    from Foam.finiteVolume.cfdTools.incompressible import continuityErrs
    cumulativeContErr = continuityErrs( mesh, phi, runTime, cumulativeContErr )
    
    p == p_rgh + rho * gh

    if p_rgh.needReference():
       from Foam.OpenFOAM import pRefValue
       p.ext_assign( p + dimensionedScalar( word( "p" ),
                                            p.dimensions(),
                                            pRefValue - getRefCellValue(p, pRefCell) ) )
       p_rgh.ext_assign( p - rho * gh )
       pass
    
    return cumulativeContErr
    
    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    from Foam.OpenFOAM.include import setRootCase
    args = setRootCase( argc, argv )

    from Foam.OpenFOAM.include import createTime
    runTime = createTime( args )

    from Foam.OpenFOAM.include import createMesh
    mesh = createMesh( runTime )
    
    from Foam.finiteVolume.cfdTools.general.include import readGravitationalAcceleration
    g = readGravitationalAcceleration( runTime, mesh )

    from Foam.finiteVolume.cfdTools.general.include import readPIMPLEControls
    pimple, nOuterCorr, nCorr, nNonOrthCorr, momentumPredictor, transonic = readPIMPLEControls( mesh )
    
    from Foam.finiteVolume.cfdTools.general.include import initContinuityErrs
    cumulativeContErr = initContinuityErrs()
    
    p_rgh, alpha1, U, phi, twoPhaseProperties, rho1, rho2, Dab, \
    alphatab, rho, rhoPhi, turbulence, gh, ghf, p, pRefCell, pRefValue = _createFields( runTime, mesh, g )

    from Foam.finiteVolume.cfdTools.general.include import readTimeControls
    adjustTimeStep, maxCo, maxDeltaT = readTimeControls( runTime )
    
    from Foam.finiteVolume.cfdTools.incompressible import CourantNo
    CoNum, meanCoNum = CourantNo( mesh, phi, runTime )
    
    from Foam.finiteVolume.cfdTools.general.include import setInitialDeltaT
    runTime = setInitialDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
    
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "\nStarting time loop\n" << nl

    while runTime.run():
        pimple, nOuterCorr, nCorr, nNonOrthCorr, momentumPredictor, transonic = readPIMPLEControls( mesh )
        adjustTimeStep, maxCo, maxDeltaT = readTimeControls( runTime )
        CoNum, meanCoNum = CourantNo( mesh, phi, runTime )
        
        from Foam.finiteVolume.cfdTools.general.include import setDeltaT
        runTime = setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )

        runTime.increment()
        ext_Info() << "Time = " << runTime.timeName() << nl << nl
        
        # --- Pressure-velocity PIMPLE corrector loop
        for oCorr in range( nOuterCorr ):
            finalIter = oCorr == nOuterCorr-1

            twoPhaseProperties.correct()

            alphaEqn( mesh, phi, alpha1, alphatab, Dab, rhoPhi, rho, rho1, rho2, turbulence )
            
            UEqn = fun_UEqn( mesh, U, p_rgh, ghf, rho, rhoPhi, turbulence, twoPhaseProperties, momentumPredictor, finalIter )
            
            # --- PISO loop
            for corr in range( nCorr ):
                cumulativeContErr = fun_pEqn( runTime, mesh, UEqn, U, p, p_rgh, gh, ghf, phi, rho, \
                                              finalIter, corr, nCorr, nNonOrthCorr, pRefCell, pRefValue, cumulativeContErr )
                pass

            turbulence.correct()
            pass

        runTime.write()
        ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << nl

        pass
    
    ext_Info() << "End\n" << nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "010701" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam1.7.1 or higher \n "     
   pass


#--------------------------------------------------------------------------------------
