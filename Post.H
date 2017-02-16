   /*---------------------------------------------------------------------------------------------
   * 					Post Results
   * Returns the power output as P = M*W where M = SUM(Ftheta*radius) and W is omega
   * --------------------------------------------------------------------------------------------*/
   const static scalar Pi = 3.141592654;

   Istream& is1 = mesh.solutionDict().subDict("B2Fan").lookup("FanSpeed");
   is1.format(IOstream::ASCII); 
   is1 >> FanSpeed;
   
   Power31 = 0;
   TotTorque31 = 0;
   Power33 = 0;
   TotTorque33 = 0;
   Power35 = 0;
   TotTorque35 = 0;
   
   Power1 = 0;
   TotTorque1 = 0;
   Power2 = 0;
   TotTorque2 = 0;
   Power3 = 0;
   TotTorque3 = 0; 
   Power4 = 0;
   TotTorque4 = 0;
   Power5 = 0;
   TotTorque5 = 0;
   
  const label cellZoneIDB31 = mesh.cellZones().findZoneID("Bfan15");
  const cellZone& zoneB31 = mesh.cellZones()[cellZoneIDB31];
  const cellZoneMesh& zoneMeshB31 = zoneB31.zoneMesh();
  const labelList& cellsZoneB31 = zoneMeshB31[cellZoneIDB31];
  
  const label cellZoneIDB33 = mesh.cellZones().findZoneID("Bfan35");
  const cellZone& zoneB33 = mesh.cellZones()[cellZoneIDB33];
  const cellZoneMesh& zoneMeshB33 = zoneB33.zoneMesh();
  const labelList& cellsZoneB33 = zoneMeshB33[cellZoneIDB33];
  
  const label cellZoneIDB35 = mesh.cellZones().findZoneID("Bfan45");
  const cellZone& zoneB35 = mesh.cellZones()[cellZoneIDB35];
  const cellZoneMesh& zoneMeshB35 = zoneB35.zoneMesh();
  const labelList& cellsZoneB35 = zoneMeshB35[cellZoneIDB35];
  
  const label cellZoneIDB1 = mesh.cellZones().findZoneID("Bfan65");
  const cellZone& zoneB1 = mesh.cellZones()[cellZoneIDB1];
  const cellZoneMesh& zoneMeshB1 = zoneB1.zoneMesh();
  const labelList& cellsZoneB1 = zoneMeshB1[cellZoneIDB1];
  
  const label cellZoneIDB2 = mesh.cellZones().findZoneID("Bfan63");
  const cellZone& zoneB2 = mesh.cellZones()[cellZoneIDB2];
  const cellZoneMesh& zoneMeshB2 = zoneB2.zoneMesh();
  const labelList& cellsZoneB2 = zoneMeshB2[cellZoneIDB2];
  
  const label cellZoneIDB3 = mesh.cellZones().findZoneID("Bfan11");
  const cellZone& zoneB3 = mesh.cellZones()[cellZoneIDB3];
  const cellZoneMesh& zoneMeshB3 = zoneB3.zoneMesh();
  const labelList& cellsZoneB3 = zoneMeshB3[cellZoneIDB3];
  
  const label cellZoneIDB4 = mesh.cellZones().findZoneID("Bfan41");
  const cellZone& zoneB4 = mesh.cellZones()[cellZoneIDB4];
  const cellZoneMesh& zoneMeshB4 = zoneB4.zoneMesh();
  const labelList& cellsZoneB4 = zoneMeshB4[cellZoneIDB4];
  
  const label cellZoneIDB5 = mesh.cellZones().findZoneID("Bfan61");
  const cellZone& zoneB5 = mesh.cellZones()[cellZoneIDB5];
  const cellZoneMesh& zoneMeshB5 = zoneB5.zoneMesh();
  const labelList& cellsZoneB5 = zoneMeshB5[cellZoneIDB5];
   
   forAll(cellsZoneB31, cIB) {
      if(diskID[cellsZoneB31[cIB]] == 2 || diskID[cellsZoneB31[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB31[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB31[cIB]]) + VolumeForce[cellsZoneB31[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB31[cIB]]);
	TotTorque31 = TotTorque31 + FTAN*mesh.V()[cellsZoneB31[cIB]]*mag(DiskRad[cellsZoneB31[cIB]]);
	 }
   }
   
   forAll(cellsZoneB33, cIB) {
      if(diskID[cellsZoneB33[cIB]] == 2 || diskID[cellsZoneB33[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB33[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB33[cIB]]) + VolumeForce[cellsZoneB33[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB33[cIB]]);
	TotTorque33 = TotTorque33 + FTAN*mesh.V()[cellsZoneB33[cIB]]*mag(DiskRad[cellsZoneB33[cIB]]);
	 }
   }
   
   forAll(cellsZoneB35, cIB) {
      if(diskID[cellsZoneB35[cIB]] == 2 || diskID[cellsZoneB35[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB35[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB35[cIB]]) + VolumeForce[cellsZoneB35[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB35[cIB]]);
	TotTorque35 = TotTorque35 + FTAN*mesh.V()[cellsZoneB35[cIB]]*mag(DiskRad[cellsZoneB35[cIB]]);
	 }
   }
   
   forAll(cellsZoneB1, cIB) {
      if(diskID[cellsZoneB1[cIB]] == 2 || diskID[cellsZoneB1[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB1[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB1[cIB]]) + VolumeForce[cellsZoneB1[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB1[cIB]]);
	TotTorque1 = TotTorque1 + FTAN*mesh.V()[cellsZoneB1[cIB]]*mag(DiskRad[cellsZoneB1[cIB]]);
	 }
   }
   
   forAll(cellsZoneB2, cIB) {
      if(diskID[cellsZoneB2[cIB]] == 2 || diskID[cellsZoneB2[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB2[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB2[cIB]]) + VolumeForce[cellsZoneB2[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB2[cIB]]);
	TotTorque2 = TotTorque2 + FTAN*mesh.V()[cellsZoneB2[cIB]]*mag(DiskRad[cellsZoneB2[cIB]]);
	 }
   }
   
   forAll(cellsZoneB3, cIB) {
      if(diskID[cellsZoneB3[cIB]] == 2 || diskID[cellsZoneB3[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB3[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB3[cIB]]) + VolumeForce[cellsZoneB3[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB3[cIB]]);
	TotTorque3 = TotTorque3 + FTAN*mesh.V()[cellsZoneB3[cIB]]*mag(DiskRad[cellsZoneB3[cIB]]);
	 }
   }
   
   forAll(cellsZoneB4, cIB) {
      if(diskID[cellsZoneB4[cIB]] == 2 || diskID[cellsZoneB4[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB4[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB4[cIB]]) + VolumeForce[cellsZoneB4[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB4[cIB]]);
	TotTorque4 = TotTorque4 + FTAN*mesh.V()[cellsZoneB4[cIB]]*mag(DiskRad[cellsZoneB4[cIB]]);
	 }
   }
   
   forAll(cellsZoneB5, cIB) {
      if(diskID[cellsZoneB5[cIB]] == 2 || diskID[cellsZoneB5[cIB]] == 4) { 
	scalar FTAN = VolumeForce[cellsZoneB5[cIB]].x()*-Foam::sin(DiskTheta[cellsZoneB5[cIB]]) + VolumeForce[cellsZoneB5[cIB]].y()*Foam::cos(DiskTheta[cellsZoneB5[cIB]]);
	TotTorque5 = TotTorque5 + FTAN*mesh.V()[cellsZoneB5[cIB]]*mag(DiskRad[cellsZoneB5[cIB]]);
	 }
   }
   
  Power31 = TotTorque31*(FanSpeed*2*Pi/60);
  Power33 = TotTorque33*(FanSpeed*2*Pi/60);
  Power35 = TotTorque35*(FanSpeed*2*Pi/60);
  Power1 = TotTorque1*(FanSpeed*2*Pi/60);
  Power2 = TotTorque2*(FanSpeed*2*Pi/60);
  Power3 = TotTorque3*(FanSpeed*2*Pi/60);
  Power4 = TotTorque4*(FanSpeed*2*Pi/60);
  Power5 = TotTorque5*(FanSpeed*2*Pi/60);
  
  label patchIDTop31 = mesh.boundary().findPatchID("fanTop15_shell");
  label patchIDBot31 = mesh.boundary().findPatchID("fanBottom15_shell");
  
  label patchIDTop33 = mesh.boundary().findPatchID("fanTop35_shell");
  label patchIDBot33 = mesh.boundary().findPatchID("fanBottom35_shell");
     
  label patchIDTop35 = mesh.boundary().findPatchID("fanTop45_shell");
  label patchIDBot35 = mesh.boundary().findPatchID("fanBottom45_shell");
  
  label patchIDTop1 = mesh.boundary().findPatchID("fanTop65_shell");
  label patchIDBot1 = mesh.boundary().findPatchID("fanBottom65_shell");
  
  label patchIDTop2 = mesh.boundary().findPatchID("fanTop63_shell");
  label patchIDBot2 = mesh.boundary().findPatchID("fanBottom63_shell");
     
  label patchIDTop3 = mesh.boundary().findPatchID("fanTop11_shell");
  label patchIDBot3 = mesh.boundary().findPatchID("fanBottom11_shell");
  
  label patchIDTop4 = mesh.boundary().findPatchID("fanTop41_shell");
  label patchIDBot4 = mesh.boundary().findPatchID("fanBottom41_shell");
  
  label patchIDTop5 = mesh.boundary().findPatchID("fanTop61_shell");
  label patchIDBot5 = mesh.boundary().findPatchID("fanBottom61_shell");
  
  mFanTop31 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop31]*gSum(phi.boundaryField()[patchIDTop31])));
  
  mFanTop33 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop33]*gSum(phi.boundaryField()[patchIDTop33])));
  
  mFanTop35 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop35]*gSum(phi.boundaryField()[patchIDTop35])));
  
  mFanBot31 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot31]*gSum(phi.boundaryField()[patchIDBot31])));
  
  mFanBot33 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot33]*gSum(phi.boundaryField()[patchIDBot33])));
  
  mFanBot35 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot35]*gSum(phi.boundaryField()[patchIDBot35])));
  
  mFanTop1 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop1]*gSum(phi.boundaryField()[patchIDTop1])));
  
  mFanTop2 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop2]*gSum(phi.boundaryField()[patchIDTop2])));
  
  mFanTop3 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop3]*gSum(phi.boundaryField()[patchIDTop3])));
  
  mFanTop4 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop4]*gSum(phi.boundaryField()[patchIDTop4])));
  
  mFanTop5 = Foam::sqrt(gSum(phi.boundaryField()[patchIDTop5]*gSum(phi.boundaryField()[patchIDTop5])));
  
  mFanBot1 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot1]*gSum(phi.boundaryField()[patchIDBot1])));
  
  mFanBot2 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot2]*gSum(phi.boundaryField()[patchIDBot2])));
  
  mFanBot3 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot3]*gSum(phi.boundaryField()[patchIDBot3])));
  
  mFanBot4 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot4]*gSum(phi.boundaryField()[patchIDBot4])));
  
  mFanBot5 = Foam::sqrt(gSum(phi.boundaryField()[patchIDBot5]*gSum(phi.boundaryField()[patchIDBot5])));
    
  mFan31 = (mFanTop31 + mFanBot31)/2; 
  mFan33 = (mFanTop33 + mFanBot33)/2;   
  mFan35 = (mFanTop35 + mFanBot35)/2;  
  mFan1 = (mFanTop1 + mFanBot1)/2;
  mFan2 = (mFanTop2 + mFanBot2)/2; 
  mFan3 = (mFanTop3 + mFanBot3)/2;
  mFan4 = (mFanTop4 + mFanBot4)/2;
  mFan5 = (mFanTop5 + mFanBot5)/2;
        
  //Create list of size = number of processors
  thisProcNb = Pstream::myProcNo();
  labelList mFanPerProc31(Pstream::nProcs(), 0.0);
  labelList mFanPerProc33(Pstream::nProcs(), 0.0);
  labelList mFanPerProc35(Pstream::nProcs(), 0.0);
  labelList mFanPerProc1(Pstream::nProcs(), 0.0);
  labelList mFanPerProc2(Pstream::nProcs(), 0.0);
  labelList mFanPerProc3(Pstream::nProcs(), 0.0);
  labelList mFanPerProc4(Pstream::nProcs(), 0.0);
  labelList mFanPerProc5(Pstream::nProcs(), 0.0);
  
  labelList pntFieldPowerPerProc31(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc33(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc35(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc1(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc2(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc3(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc4(Pstream::nProcs(), 0.0);
  labelList pntFieldPowerPerProc5(Pstream::nProcs(), 0.0);
  
  //Populate with power output of each processor
  mFanPerProc31[thisProcNb] = mFan31*1e3;
  mFanPerProc33[thisProcNb] = mFan33*1e3;
  mFanPerProc35[thisProcNb] = mFan35*1e3;
  mFanPerProc1[thisProcNb] = mFan1*1e3;
  mFanPerProc2[thisProcNb] = mFan2*1e3;
  mFanPerProc3[thisProcNb] = mFan3*1e3;
  mFanPerProc4[thisProcNb] = mFan4*1e3;
  mFanPerProc5[thisProcNb] = mFan5*1e3;
  
  pntFieldPowerPerProc31[thisProcNb] = Power31;
  pntFieldPowerPerProc33[thisProcNb] = Power33;
  pntFieldPowerPerProc35[thisProcNb] = Power35;
  pntFieldPowerPerProc1[thisProcNb] = Power1;
  pntFieldPowerPerProc2[thisProcNb] = Power2;
  pntFieldPowerPerProc3[thisProcNb] = Power3;
  pntFieldPowerPerProc4[thisProcNb] = Power4;
  pntFieldPowerPerProc5[thisProcNb] = Power5;
  
  //Reduce to make accessible to all processors
  reduce(mFanPerProc31, sumOp<labelList>());
  reduce(mFanPerProc33, sumOp<labelList>());
  reduce(mFanPerProc35, sumOp<labelList>());
  reduce(mFanPerProc1, sumOp<labelList>());
  reduce(mFanPerProc2, sumOp<labelList>());
  reduce(mFanPerProc3, sumOp<labelList>());
  reduce(mFanPerProc4, sumOp<labelList>());
  reduce(mFanPerProc5, sumOp<labelList>());
  
  reduce(pntFieldPowerPerProc31, sumOp<labelList>());
  reduce(pntFieldPowerPerProc33, sumOp<labelList>());
  reduce(pntFieldPowerPerProc35, sumOp<labelList>());
  reduce(pntFieldPowerPerProc1, sumOp<labelList>());
  reduce(pntFieldPowerPerProc2, sumOp<labelList>());
  reduce(pntFieldPowerPerProc3, sumOp<labelList>());
  reduce(pntFieldPowerPerProc4, sumOp<labelList>());
  reduce(pntFieldPowerPerProc5, sumOp<labelList>());
    
  mFanSum31 = 0;
  mFanCounter31 = 0;
  mFanTot31 = 0;
  mFanSum33 = 0;
  mFanCounter33 = 0;
  mFanTot33 = 0;
  mFanSum35 = 0;
  mFanCounter35 = 0;
  mFanTot35 = 0;
  
  mFanSum1 = 0;
  mFanCounter1 = 0;
  mFanTot1 = 0;
  mFanSum2 = 0;
  mFanCounter2 = 0;
  mFanTot2 = 0;
  mFanSum3 = 0;
  mFanCounter3 = 0;
  mFanTot3 = 0;
  mFanSum4 = 0;
  mFanCounter4 = 0;
  mFanTot4 = 0;
  mFanSum5 = 0;
  mFanCounter5 = 0;
  mFanTot5 = 0;
  
  PowerTot31 = 0;
  PowerTot33 = 0;
  PowerTot35 = 0;
  PowerTot1 = 0;
  PowerTot2 = 0;
  PowerTot3 = 0;
  PowerTot4 = 0;
  PowerTot5 = 0;
  
  vForceSum = 0;
  energySourceSum = 0;
  
  vForceSum = mag(gSum(VolumeForce));
  energySourceSum = gSum(energySource);
  
  //Output total power if on master processor
  if(Pstream::master() == true) { 
    for(int i = 0; i < Pstream::nProcs(); i++) {
	  if(mFanPerProc31[i] > 0) { 
	    mFanSum31 = mFanSum31 + mFanPerProc31[i];
	    mFanCounter31 = mFanCounter31 + 1;
	  }
	  if(mFanPerProc33[i] > 0) { 
	    mFanSum33 = mFanSum33 + mFanPerProc33[i];
	    mFanCounter33 = mFanCounter33 + 1;
	  }
	  if(mFanPerProc35[i] > 0) { 
	    mFanSum35 = mFanSum35 + mFanPerProc35[i];
	    mFanCounter35 = mFanCounter35 + 1;
	  }
	  if(mFanPerProc1[i] > 0) { 
	    mFanSum1 = mFanSum1 + mFanPerProc1[i];
	    mFanCounter1 = mFanCounter1 + 1;
	  }
	  if(mFanPerProc2[i] > 0) { 
	    mFanSum2 = mFanSum2 + mFanPerProc2[i];
	    mFanCounter2 = mFanCounter2 + 1;
	  }
	  if(mFanPerProc3[i] > 0) { 
	    mFanSum3 = mFanSum3 + mFanPerProc3[i];
	    mFanCounter3 = mFanCounter3 + 1;
	  }
	  if(mFanPerProc4[i] > 0) { 
	    mFanSum4 = mFanSum4 + mFanPerProc4[i];
	    mFanCounter4 = mFanCounter4 + 1;
	  }
	  if(mFanPerProc5[i] > 0) { 
	    mFanSum5 = mFanSum5 + mFanPerProc5[i];
	    mFanCounter5 = mFanCounter5 + 1;
	  }
    }
    
    if(mFanCounter31 == 0) { 
      mFanTot31 = 0; 
    } else {
      mFanTot31 = mFanSum31/(1e3*mFanCounter31);
    }
    if(mFanCounter33 == 0) { 
      mFanTot33 = 0; 
    } else {
      mFanTot33 = mFanSum33/(1e3*mFanCounter33);
    }
    if(mFanCounter35 == 0) { 
      mFanTot35 = 0; 
    } else {
      mFanTot35 = mFanSum35/(1e3*mFanCounter35);
    }
    
    if(mFanCounter1 == 0) { 
      mFanTot1 = 0; 
    } else {
      mFanTot1 = mFanSum1/(1e3*mFanCounter1);
    }
    if(mFanCounter2 == 0) { 
      mFanTot2 = 0; 
    } else {
      mFanTot2 = mFanSum2/(1e3*mFanCounter2);
    }
    if(mFanCounter3 == 0) { 
      mFanTot3 = 0; 
    } else {
      mFanTot3 = mFanSum3/(1e3*mFanCounter3);
    }
    if(mFanCounter4 == 0) { 
      mFanTot4 = 0; 
    } else {
      mFanTot4 = mFanSum4/(1e3*mFanCounter4);
    }
    if(mFanCounter5 == 0) { 
      mFanTot5 = 0; 
    } else {
      mFanTot5 = mFanSum5/(1e3*mFanCounter5);
    }
    
    PowerTot31 = sum(pntFieldPowerPerProc31);
    PowerTot33 = sum(pntFieldPowerPerProc33);
    PowerTot35 = sum(pntFieldPowerPerProc35);
    PowerTot1 = sum(pntFieldPowerPerProc1);
    PowerTot2 = sum(pntFieldPowerPerProc2);
    PowerTot3 = sum(pntFieldPowerPerProc3);
    PowerTot4 = sum(pntFieldPowerPerProc4);
    PowerTot5 = sum(pntFieldPowerPerProc5);
              
    Sout << " " << nl << endl;
    Sout << "Volume Flow Rate through Fan15: " << mFanTot31 << endl;
    Sout << "Volume Flow Rate through Fan35: " << mFanTot33 << endl;
    Sout << "Volume Flow Rate through Fan45: " << mFanTot35 << endl;
    Sout << "Volume Flow Rate through Fan65: " << mFanTot1 << endl;
    Sout << "Volume Flow Rate through Fan63: " << mFanTot2 << endl;
    Sout << "Volume Flow Rate through Fan11: " << mFanTot3 << endl;
    Sout << "Volume Flow Rate through Fan41: " << mFanTot4 << endl;
    Sout << "Volume Flow Rate through Fan61: " << mFanTot5 << nl << endl;
    
    Sout << "Power of Fan15: " << PowerTot31 << endl;
    Sout << "Power of Fan35: " << PowerTot33 << endl;
    Sout << "Power of Fan45: " << PowerTot35 << endl;
    Sout << "Power of Fan65: " << PowerTot1 << endl;
    Sout << "Power of Fan63: " << PowerTot2 << endl;
    Sout << "Power of Fan11: " << PowerTot3 << endl;
    Sout << "Power of Fan41: " << PowerTot4 << endl;
    Sout << "Power of Fan61: " << PowerTot5 << nl << endl;
    Sout << "Relax: " << relax << nl << endl;
    
  }

  scalar TOutlet = 0;
  scalar counter_Inlet = 1;
  scalar counter_Outlet = 1;
  scalar pUp = 0;
  scalar pDown = 0;
  
  forAll(mesh.C(), cellI) { 
      //Upstream heat exchanger zone/Inlet
      if(HEID[cellI] == 1) { 
	pUp = pUp + p_rgh[cellI];
	counter_Inlet = counter_Inlet + 1;  
      }
      
      //Downstream heat exchanger zone/outlet
      if(HEID[cellI] == 3) { 
	TOutlet = TOutlet + T[cellI];
	pDown = pDown + p_rgh[cellI];
	counter_Outlet = counter_Outlet + 1;
      }
    }
    
    if(counter_Inlet >= 2) { 
      pUp = pUp/(counter_Inlet - 1);
    } else { 
      pUp = 0;
    }
    
    if(counter_Outlet >= 2) { 
      pDown = pDown/(counter_Outlet - 1);
      TOutlet = TOutlet/(counter_Outlet-1);
    } else {
      pDown = 0;
      TOutlet = 0;
    }
    
      
  //Create list of size = number of processors
  labelList TOutletPerProc(Pstream::nProcs(), 0.0);
  labelList PUpPerProc(Pstream::nProcs(), 0.0);
  labelList PDownPerProc(Pstream::nProcs(), 0.0);
  
  //Populate with power output of each processor
  TOutletPerProc[thisProcNb] = TOutlet;
  PUpPerProc[thisProcNb] = pUp;
  PDownPerProc[thisProcNb] = pDown;
  
  //Reduce to make accessible to all processors
  reduce(TOutletPerProc, sumOp<labelList>());
  reduce(PUpPerProc, sumOp<labelList>());
  reduce(PDownPerProc, sumOp<labelList>());
  
  scalar TOutletSum = 0;
  scalar TOutletSumCounter = 0;
  scalar PUpSum = 0;
  scalar PUpSumCounter = 0;
  scalar PDownSum = 0;
  scalar PDownSumCounter = 0;
  
  if(Pstream::master() == true) {
      for(int i = 0; i < Pstream::nProcs(); i++) {
	if(TOutletPerProc[i] > 0) { 
	  TOutletSum = TOutletSum + TOutletPerProc[i];
	  TOutletSumCounter = TOutletSumCounter + 1;
	}
	if(PUpPerProc[i] > 0) { 
	  PUpSum = PUpSum + PUpPerProc[i];
	  PUpSumCounter = PUpSumCounter + 1;
	}
	if(PDownPerProc[i] > 0) { 
	  PDownSum = PDownSum + PDownPerProc[i];
	  PDownSumCounter = PDownSumCounter + 1;
	}
      }
            
    scalar PDiff = PUpSum/PUpSumCounter - PDownSum/PDownSumCounter;
    Sout << "Heat Exchanger Outlet Temperature: " << TOutletSum/TOutletSumCounter << endl;
    Sout << "Heat Exchanger Mass Flow Rate: " << ma_Tot << endl;
    Sout << "Heat Exchanger Pressure Difference: " << PDiff << nl << endl;
    
  }
    