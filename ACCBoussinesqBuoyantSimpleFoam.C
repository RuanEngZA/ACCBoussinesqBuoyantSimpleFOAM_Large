/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantSimpleFoam

Description
    Steady-state solver for buoyant, turbulent flow of incompressible fluids

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "ActuatorDisk.H"
#include "HeatExchanger.H"
#include "singlePhaseTransportModel.H"
#include "powerLawXVelocityFvPatchVectorField.H"
#include "powerLawYVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar FanSpeed, Power31, Power33, Power35, TotTorque31, TotTorque33, TotTorque35, mFanTop31, mFanBot31, mFan31, mFanTot31, PowerTot31, PowerTot33, PowerTot35, mFanSum31, mFanCounter31, vForceSum, energySourceSum, mFanTop33, mFanBot33, mFan33, mFanTot33, mFanSum33, mFanCounter33, mFanTop35, mFanBot35, mFan35, mFanTot35, mFanSum35, mFanCounter35;

scalar mFanTop1, mFanBot1, mFan1, mFanTot1, mFanSum1, mFanCounter1, TotTorque1, Power1, PowerTot1;
scalar mFanTop2, mFanBot2, mFan2, mFanTot2, mFanSum2, mFanCounter2, TotTorque2, Power2, PowerTot2;
scalar mFanTop3, mFanBot3, mFan3, mFanTot3, mFanSum3, mFanCounter3, TotTorque3, Power3, PowerTot3;
scalar mFanTop4, mFanBot4, mFan4, mFanTot4, mFanSum4, mFanCounter4, TotTorque4, Power4, PowerTot4;
scalar mFanTop5, mFanBot5, mFan5, mFanTot5, mFanSum5, mFanCounter5, TotTorque5, Power5, PowerTot5;

scalar volumeForceDeltaX, volumeForceDeltaY, volumeForceDeltaZ, volumeForceDeltaCounter, Power_previous;
vector U_up, U_down;
vector VolumeForce_previous;
scalar HeatSource_previous;
scalar percentChange, heatOnOff;

scalar rho = 1.0857;
scalar mu = 1.7948e-5;
scalar cpa = 1006.609;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createDiskFields.H"
    #include "createHEFields.H"
    #include "createTransFields.H"

    simpleControl simple(mesh);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //========================================================================
    //				Preliminary Field Loops
    //========================================================================
    
    PreCalcs Pre;
    Pre.readDataB2(mesh);
    
    ActuatorDiskModel AcDisk;
    AcDisk.readDataB2(mesh);
    AcDisk.readDataA(mesh);
    
    forAll(mesh.C(), cellI) {
      if(diskID[cellI] == 2) {
// 	DiskRad[cellI] = AcDisk.getRadius(mesh.C()[cellI]);
	DiskTheta[cellI] = AcDisk.getTheta(mesh.C()[cellI], DiskRad[cellI]);
	DiskGama[cellI] = AcDisk.getBladeAngleB(DiskRad[cellI]);
	DiskChord[cellI] = AcDisk.getChordLengthB(DiskRad[cellI]);
      }
      if(diskID[cellI] == 4) { 
// 	DiskRad[cellI] = AcDisk.getRadiusA(mesh.C()[cellI]);
	DiskTheta[cellI] = AcDisk.getTheta(mesh.C()[cellI], DiskRad[cellI]);
	DiskGama[cellI] = AcDisk.getBladeAngleA(DiskRad[cellI]);
	DiskChord[cellI] = 0.722;
      }
    }
    
    //Read in Airfoil profile data for NASA 0413-LS profile
    AcDisk.readCSV_NASA0413LS();
    double timeHolder = 0;
    
    HeatExchanger HE;
    HE.readDataHE(mesh);
     
    tensor heAxisTransform = HE.axesRotation();
     
    heLeftDExpl =  heAxisTransform.T() & (heLeftDExpl & heAxisTransform);
    heLeftFExpl =  heAxisTransform.T() & (heLeftFExpl & heAxisTransform);
     
    heRightDExpl =  heAxisTransform.T() & (heRightDExpl & heAxisTransform);
    heRightFExpl =  heAxisTransform.T() & (heRightFExpl & heAxisTransform);
    
    scalar relax;
    
    Istream& is1 = mesh.solutionDict().subDict("General").lookup("VolForceRelax");
    is1.format(IOstream::ASCII); 
    is1 >> relax;
    
    Istream& is2 = mesh.solutionDict().subDict("General").lookup("HeatOnOff");
    is2.format(IOstream::ASCII); 
    is2 >> heatOnOff;
           
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
	    #include "TransportEqn.H"
            #include "UEqn.H"
            #include "TEqn.H"
            #include "pEqn.H" 
	    #include "Post.H"
        }

        turbulence->correct();

	runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	    
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
