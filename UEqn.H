    // Solve the Momentum equation
    Istream& is3 = mesh.solutionDict().subDict("General").lookup("VolForceRelax");
    is3.format(IOstream::ASCII); 
    is3 >> relax;
        
    timeHolder = timeHolder + 1;
    vector U_up, U_down;
    VolumeForce_previous.x() = 0;
    VolumeForce_previous.y() = 0;
    VolumeForce_previous.z() = 0;
    volumeForceDeltaX = 0;
    volumeForceDeltaY = 0;
    volumeForceDeltaZ = 0;
    volumeForceDeltaCounter = 1;
    
//     if(timeHolder >= 480) { 
//       relax = 0.001;
//     } else {
//       relax = relax;
//     }
            
    forAll(mesh.C(), cellI) { 
      if(diskID[cellI] == 2 || diskID[cellI] == 4) { 
	//Store previous iteration value
	VolumeForce_previous.x() = VolumeForce[cellI].x();
	VolumeForce_previous.y() = VolumeForce[cellI].y();
	VolumeForce_previous.z() = VolumeForce[cellI].z();
	
	U_up.x() = UForwardX[cellI];
	U_up.y() = UForwardY[cellI];
	U_up.z() = UForwardZ[cellI];
	
	U_down.x() = UarbBackward[cellI].x();
	U_down.y() = UarbBackward[cellI].y();
 	U_down.z() = UarbBackward[cellI].z();
		
	scalar volume = mesh.V()[cellI];
	
	VolumeForce[cellI] = AcDisk.calcVolumeForce(mesh.C()[cellI], DiskRad[cellI], U_up, U_down, DiskTheta[cellI], DiskGama[cellI], DiskChord[cellI], relativeVelocity[cellI], betaAC[cellI], alpha[cellI], Cl[cellI], Cd[cellI], Torque[cellI], volume, AverageVelocity[cellI], diskID[cellI], rhok[cellI], mu);
		
	//Explicitly relax source terms
	VolumeForce[cellI].x() = VolumeForce_previous.x() + relax*(VolumeForce[cellI].x() - VolumeForce_previous.x());
	VolumeForce[cellI].y() = VolumeForce_previous.y() + relax*(VolumeForce[cellI].y() - VolumeForce_previous.y());
	VolumeForce[cellI].z() = VolumeForce_previous.z() + relax*(VolumeForce[cellI].z() - VolumeForce_previous.z());
	
	volumeForceDeltaX = VolumeForce[cellI].x() - VolumeForce_previous.x() + volumeForceDeltaX;
	volumeForceDeltaY = VolumeForce[cellI].y() - VolumeForce_previous.y() + volumeForceDeltaY;
	volumeForceDeltaZ = VolumeForce[cellI].z() - VolumeForce_previous.z() + volumeForceDeltaZ;
	volumeForceDeltaCounter = volumeForceDeltaCounter + 1;
	
	VolumeForce_previous.x() = 0;
	VolumeForce_previous.y() = 0;
	VolumeForce_previous.z() = 0;
      }
    }
    
    if(Pstream::master() == true) { 
	  Sout << "\nVolumeForce Calculated " << nl << endl;
	}
        
    HEExForce = HEExLeft*mu*(heLeftDExpl & U) + HEExLeft*0.5*mag(U)*rho*(heLeftFExpl & U) + HEExRight*mu*(heRightDExpl & U) + HEExRight*0.5*mag(U)*rho*(heRightFExpl & U);
         
    HE.PressureDrop(mesh, HEExLeft, HEExRight, rho, U, HEImForce, rhoDum, muDum, UDum, SpDum, KthetaTField, SuDum);
        
    scalar ma_Tot = HE.massFlowRate(mesh, ma, projectedArea, U, rhok, HEEnergyLeft, HEEnergyRight);
    
        
    if(ma_Tot < 300) { 
      ma_Tot = 300;
    }
    
    forAll(mesh.C(), cellI) { 
      if (HEExLeft[cellI] == 1) {
	scalar pDum = p[cellI];
	scalar Tup = HeLeftUpstreamTrans[cellI];
	scalar Tdown = HeLeftDownstreamTrans[cellI];
	
	HE.TheoreticalPressureDrop(mesh, Tup, Tdown, U[cellI], HEImForce[cellI], pDum, rho, KthetaTField[cellI], mu);

      } else {
	if(HEExRight[cellI] == 1) { 
	scalar pDum= p[cellI];
	scalar Tup = HeRightUpstreamTrans[cellI];
	scalar Tdown = HeRightDownstreamTrans[cellI];
	
	HE.TheoreticalPressureDrop(mesh, Tup, Tdown, U[cellI], HEImForce[cellI], pDum, rho, KthetaTField[cellI], mu); 

	}
      }
           

    if(HEExLeft[cellI] == 1) { 
//     HeatSource_previous = energySource[cellI];
    
    scalar Tup = HeLeftUpstreamTrans[cellI];
    scalar Tdown = HeLeftDownstreamTrans[cellI];
    scalar pDum = p[cellI];
    			
    energySource[cellI] = HE.EnergySource(mesh, rho, Tup, Tdown, mesh.V()[cellI], pDum, phi[cellI], U[cellI], HEExLeft[cellI], HEExRight[cellI], T[cellI], ma_Tot);
    
//     energySource[cellI] = HeatSource_previous + heatRelax*(energySource[cellI] - HeatSource_previous);
//     HeatSource_previous = 0;

  } else { 
    if(HEExRight[cellI] == 1) {
//       HeatSource_previous = energySource[cellI];
      
      scalar Tup = HeRightUpstreamTrans[cellI];
      scalar Tdown = HeRightDownstreamTrans[cellI];
      scalar pDum = p[cellI];
		      
      energySource[cellI] = HE.EnergySource(mesh, rho, Tup, Tdown, mesh.V()[cellI], pDum, phi[cellI], U[cellI], HEExLeft[cellI], HEExRight[cellI], T[cellI], ma_Tot);
      
//       energySource[cellI] = HeatSource_previous + heatRelax*(energySource[cellI] - HeatSource_previous);
//       HeatSource_previous = 0;

    }
	  
    }
    
  }
  
    if(Pstream::master() == true) { 
	Sout << "\nHeat Exchanger Calculated (maTot = " << ma_Tot << ")" <<  nl << endl;
      
    }
    
    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U) - VolumeForce + fvm::SuSp(HEImForce,U) + HEExForce
     ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn()
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
    
    
//    //Create list of size = number of processors
    label thisProcNb = Pstream::myProcNo();
//     labelList VX(Pstream::nProcs(), 0.0);
//     labelList VY(Pstream::nProcs(), 0.0);
//     labelList VZ(Pstream::nProcs(), 0.0);
//     
//     //Populate with power output of each processor
//     VX[thisProcNb] = volumeForceDeltaX;
//     VY[thisProcNb] = volumeForceDeltaY;
//     VZ[thisProcNb] = volumeForceDeltaZ;
//     
//     //Reduce to make accessible to all processors
//     reduce(VX, sumOp<labelList>());
//     reduce(VY, sumOp<labelList>());
//     reduce(VZ, sumOp<labelList>());
//     
//     
//     scalar VXSum = 0;
//     scalar VXCounter = 1;
//     scalar VYSum = 0;
//     scalar VYCounter = 1;
//     scalar VZSum = 0;
//     scalar VZCounter = 1;
//     
//     for(int i = 0; i < Pstream::nProcs(); i++) {
// 	    if(VX[i] != 0) { 
// 	      VXSum = VXSum + VX[i];
// 	      VXCounter = VXCounter + 1;
// 	    }
// 	    if(VY[i] != 0) { 
// 	      VYSum = VYSum + VY[i];
// 	      VYCounter = VYCounter + 1;
// 	    }
// 	    if(VZ[i] != 0) { 
// 	      VZSum = VZSum + VZ[i];
// 	      VZCounter = VZCounter + 1;
// 	    }
//       }
//           
//     //Output total power if on master processor
//     if(Pstream::master() == true) { 
//       for(int i = 0; i < Pstream::nProcs(); i++) {
// 	    if(VX[i] != 0) { 
// 	      VXSum = VXSum + VX[i];
// 	      VXCounter = VXCounter + 1;
// 	    }
// 	    if(VY[i] != 0) { 
// 	      VYSum = VYSum + VY[i];
// 	      VYCounter = VYCounter + 1;
// 	    }
// 	    if(VZ[i] != 0) { 
// 	      VZSum = VZSum + VZ[i];
// 	      VZCounter = VZCounter + 1;
// 	    }
//       }
//       
//       Sout << "\nVolumeForce DeltaX: " << VXSum/(VXCounter) << endl;
//       Sout << "VolumeForce DeltaY: " << VYSum/(VYCounter) << endl;
//       Sout << "VolumeForce DeltaZ: " << VZSum/(VZCounter) << endl;
// //       Sout << "energy Source: " << gSum(energySource) << endl;
//       Sout << "Relax: " << relax << nl << endl;
//       
//  }
    
    
    
    
    
