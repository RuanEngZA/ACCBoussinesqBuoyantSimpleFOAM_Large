{ 
  forAll(mesh.C(), cellI) {
        
    if(diskID[cellI] == 1 || diskID[cellI] == 3) { 
    UTransport[cellI].x() = U[cellI].x();
    UTransport[cellI].y() = U[cellI].y();
    UTransport[cellI].z() = U[cellI].z();
    
    UX[cellI] = U[cellI].x();
    UY[cellI] = U[cellI].y();
    UZ[cellI] = U[cellI].z();
    }
    
  }

// 	  tmp <fvVectorMatrix> UarbForEqn
// 	  (
// 	    fvm::div(phiFor, UarbForward) 
// 	    ==
// 	    fvm::laplacian(DSFor, UarbForward) + 
// 	    fvm::Sp(Large_for*1e30, UarbForward) -
// 	    Large_for*U*1e30 - 
// 	    fvm::Sp(DatStop*1e30, UarbForward)
// 	  );
// 	  
// 	  solve(UarbForEqn());
	  
	  tmp <fvVectorMatrix> UarbBackEqn
	  (
	    fvm::div(phiBack, UarbBackward) 
	    ==
	    fvm::laplacian(DSBack, UarbBackward) + 
	    fvm::Sp(Large_back*1e3, UarbBackward) -
	    Large_back*UTransport*1e3- 
	    fvm::Sp(DatStop*1e3, UarbBackward)
	  );
	  
	  solve(UarbBackEqn());
	  
	  tmp <fvScalarMatrix> UForwardEqnX
	  (
	    fvm::div(phiFor, UForwardX)
	    ==
	    fvm::laplacian(DSX, UForwardX) + 
	    fvm::Sp(Large_for*1e10, UForwardX) - 
	    Large_for*1e10*UX - 
	    fvm::Sp(DatStop*1e10, UForwardX)
	  );
	  
	  solve(UForwardEqnX());
	  
	  tmp <fvScalarMatrix> UForwardEqnY
	  (
	    fvm::div(phiFor, UForwardY)
	    ==
	    fvm::laplacian(DSY, UForwardY) + 
	    fvm::Sp(Large_for*1e10, UForwardY) - 
	    Large_for*1e10*UY - 
	    fvm::Sp(DatStop*1e10, UForwardY)
	  );
	  
	  solve(UForwardEqnY());
	  
	  tmp <fvScalarMatrix> UForwardEqnZ
	  (
	    fvm::div(phiFor, UForwardZ)
	    ==
	    fvm::laplacian(DSZ, UForwardZ) + 
	    fvm::Sp(Large_for*1e10, UForwardZ) - 
	    Large_for*1e10*UZ - 
	    fvm::Sp(DatStop*1e10, UForwardZ)
	  );
	  
	  solve(UForwardEqnZ());
	  
	  
	  forAll(mesh.C(), cellI) {
	    
	    HeRightDownstreamTrans[cellI] = 0;
	    HeRightUpstreamTrans[cellI] = 0;
	    HeLeftDownstreamTrans[cellI] = 0;
	    HeLeftUpstreamTrans[cellI] = 0;
	    
	    if(Left_for[cellI] == 1 || Left_back[cellI] == 1) { 
	    UTransportLeft[cellI] = T[cellI];
	    } else {
	      if(Right_for[cellI] == 1 || Right_back[cellI] == 1) { 
		UTransportRight[cellI] = T[cellI];
	      }
	    }
	  }

	  tmp <fvScalarMatrix> UTForEqnLeft
	  (
	    fvm::div(phiForT_Left, HeLeftUpstreamTrans) 
	    ==
	    fvm::laplacian(DSHE, HeLeftUpstreamTrans) + 
	    fvm::Sp(Left_for*1e6, HeLeftUpstreamTrans) -
	    Left_for*UTransportLeft*1e6 - 
	    fvm::Sp(DatStopHE*1e6, HeLeftUpstreamTrans)
	  );
	  
	  solve(UTForEqnLeft());
	  
	  tmp <fvScalarMatrix> UTBackEqnLeft
	  (
	    fvm::div(phiBackT_Left, HeLeftDownstreamTrans) 
	    ==
	    fvm::laplacian(DSHE, HeLeftDownstreamTrans) + 
	    fvm::Sp(Left_back*1e6, HeLeftDownstreamTrans) -
	    Left_back*UTransportLeft*1e6 - 
	    fvm::Sp(DatStopHE*1e6, HeLeftDownstreamTrans)
	  );
	 	  
	  solve(UTBackEqnLeft());
	  
	  tmp <fvScalarMatrix> UTForEqnRight
	  (
	    fvm::div(phiForT_Right, HeRightUpstreamTrans) 
	    ==
	    fvm::laplacian(DSHE, HeRightUpstreamTrans) + 
	    fvm::Sp(Right_for*1e6, HeRightUpstreamTrans) -
	    Right_for*UTransportRight*1e6 - 
	    fvm::Sp(DatStopHE*1e6, HeRightUpstreamTrans)
	  );
	  
	  solve(UTForEqnRight());
	  
	  tmp <fvScalarMatrix> UTBackEqnRight
	  (
	    fvm::div(phiBackT_Right, HeRightDownstreamTrans) 
	    ==
	    fvm::laplacian(DSHE, HeRightDownstreamTrans) + 
	    fvm::Sp(Right_back*1e6, HeRightDownstreamTrans) -
	    Right_back*UTransportRight*1e6 - 
	    fvm::Sp(DatStopHE*1e6, HeRightDownstreamTrans)
	  );
	 	  
	  solve(UTBackEqnRight());

}