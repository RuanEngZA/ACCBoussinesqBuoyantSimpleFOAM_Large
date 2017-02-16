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
    HeatExchanger

Description
    Calculates the corresponding pressure drop across the heat exchanger bundle
\*---------------------------------------------------------------------------*/

#include "HeatExchanger.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include <fstream>
#include <stdlib.h> 

namespace Foam {
  
  defineTypeNameAndDebug(HeatExchanger, 0);
  HeatExchanger::HeatExchanger(){}	//Default constructor
  HeatExchanger::~HeatExchanger(){}	//Default destructor
    
  void HeatExchanger::readDataHE(const fvMesh &iMesh) { 
    
    Istream& is1 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HEApexAngle");
    is1.format(IOstream::ASCII); 
    is1 >> HEApexAngle;
    
    Istream& is2 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HESigma_21");
    is2.format(IOstream::ASCII); 
    is2 >> HESigma_21;
    
    Istream& is3 = iMesh.solutionDict().subDict("HeatExchanger").lookup("a");
    is3.format(IOstream::ASCII); 
    is3 >> a;
    
    Istream& is4 = iMesh.solutionDict().subDict("HeatExchanger").lookup("b");
    is4.format(IOstream::ASCII); 
    is4 >> b;
    
    Istream& is5 = iMesh.solutionDict().subDict("HeatExchanger").lookup("nb");
    is5.format(IOstream::ASCII); 
    is5 >> nb;
    
    Istream& is6 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HEWidth");
    is6.format(IOstream::ASCII); 
    is6 >> HEWidth;
    
    Istream& is7 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HESigma");
    is7.format(IOstream::ASCII); 
    is7 >> HESigma;
    
    Istream& is8 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HECellWidthNormal");
    is8.format(IOstream::ASCII); 
    is8 >> HECellWidth;
    
    Istream& is9 = iMesh.solutionDict().subDict("HeatExchanger").lookup("Ny1_a");
    is9.format(IOstream::ASCII); 
    is9 >> Ny1_a;
    
    Istream& is10 = iMesh.solutionDict().subDict("HeatExchanger").lookup("Ny1_b");
    is10.format(IOstream::ASCII); 
    is10 >> Ny1_b;
    
    Istream& is11 = iMesh.solutionDict().subDict("HeatExchanger").lookup("steamhA");
    is11.format(IOstream::ASCII); 
    is11 >> steamhA;
    
    Istream& is12 = iMesh.solutionDict().subDict("HeatExchanger").lookup("steamFrontalAreaRatio");
    is12.format(IOstream::ASCII); 
    is12 >> steamFrontalAreaRatio;
    
    Istream& is13 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HEZoneLongCellNo");
    is13.format(IOstream::ASCII); 
    is13 >> HECellNo;
    
    Istream& is14 = iMesh.solutionDict().subDict("HeatExchanger").lookup("T_steam");
    is14.format(IOstream::ASCII); 
    is14 >> T_steam;
    
    Istream& is15 = iMesh.solutionDict().subDict("HeatExchanger").lookup("Ny2_a");
    is15.format(IOstream::ASCII); 
    is15 >> Ny2_a;
    
    Istream& is16 = iMesh.solutionDict().subDict("HeatExchanger").lookup("Ny2_b");
    is16.format(IOstream::ASCII); 
    is16 >> Ny2_b;
    
    Istream& is17 = iMesh.solutionDict().subDict("HeatExchanger").lookup("ntb1");
    is17.format(IOstream::ASCII); 
    is17 >> ntb1;
    
    Istream& is18 = iMesh.solutionDict().subDict("HeatExchanger").lookup("ntb2");
    is18.format(IOstream::ASCII); 
    is18 >> ntb2;
    
    Istream& is19 = iMesh.solutionDict().subDict("HeatExchanger").lookup("Afr");
    is19.format(IOstream::ASCII); 
    is19 >> Afr;
    
    Istream& is20 = iMesh.solutionDict().subDict("HeatExchanger").lookup("HEWidth");
    is20.format(IOstream::ASCII); 
    is20 >> HEWidth;
    
  }
    
  
  void HeatExchanger::TheoreticalPressureDrop(const fvMesh &iMesh, scalar &Tup, scalar &Tdown, const vector &U, scalar &HeatImForce, scalar &p, const scalar &rho, scalar &KthetaTField, const scalar &mu) { 
    /*--------------------------------------------------------------------------
			Theoretical Pressure Drop Calculations
			    *used for debuging purposes*
     Calculates the pressure loss dependant on flow through the heat exchanger
    --------------------------------------------------------------------------*/     
    vector UDum(U.x(), U.y(), U.z());
    
    scalar T_average = (Tup + Tdown)/2;
        
    //Properties evaluated at heat exchanger average temperature
    scalar rho_air = Properties_DArho(T_average, p);
    scalar mu_air = Properties_DAmu(T_average);
           
    scalar Ry = rho_air*mag(UDum)/mu_air;
             
    if(mag(UDum) == 0) { 
      Khe = 0;
    } else { 
      Khe = a*pow(Ry,b);
    }
    
    //Loss coefficient calculations
    Kdj = 1.664;			//Calculated from spreadsheet
    Ko = 7.117;				//Calculated from spreadsheet
    Kturn = 4.393396;
        
    Kup = 0.28;				//Upstream loss coefficient
    Kdown = 0.35;			//Downstream loss coefficient
    Kts = 1.6;				//ACC platform structure loss coefficient
    
    KthetaTField = Khe + Kts + Ko /*+ Kdj*/ + Kturn;        
                
    if(debug >= 2) { 
      Sout << "Heat exchanger loss coefficient: " << Khe << endl;
      Sout << "Heat exchanger loss coefficient (thetaT): " << KthetaTField << endl;
    }
    
  }
  
  void HeatExchanger::PressureDrop(const fvMesh &iMesh, const volScalarField &HEExLeft, const volScalarField &HEExRight, scalar &rho, const volVectorField &U, volScalarField &HEImForce, const dimensionedScalar rhoDum, const dimensionedScalar muDum, const dimensionedScalar UDum, const dimensionedScalar SpDum, const volScalarField &KthetaTField, const dimensionedScalar SuDum) {
    /*--------------------------------------------------------------------------
			Implicit Momentum Sink Calculation
     Calculates a volScalarField based on the draft equation for use as an 
     implicit source term in the momentum equation defined as Sp = f(U)
    --------------------------------------------------------------------------*/
    scalar Lz = HEWidth;
    
    HEImForce = HEExLeft*SpDum/Lz*((KthetaTField/2)*rho*mag(U)/UDum + ((Kup + Kdown)/2)*rho*mag(U)/UDum) + 
    HEExRight*SpDum/Lz*((KthetaTField/2)*rho*mag(U)/UDum + ((Kup + Kdown)/2)*rho*mag(U)/UDum);
        
}

scalar HeatExchanger::massFlowRate(const fvMesh &iMesh, volScalarField &ma, const volScalarField &projectedArea, const volVectorField &U, const volScalarField &rhok, const volScalarField &HEEnergyLeft, const volScalarField &HEEnergyRight) { 
  
  forAll(iMesh.C(), cellI) { 
    if(HEEnergyLeft[cellI] == 1) {  
      vector UDum(U[cellI].x(), U[cellI].y(), U[cellI].z());
      ma[cellI] = projectedArea[cellI]*rhok[cellI]*mag(UDum);
       } else {  
      if(HEEnergyRight[cellI] == 1) { 
	vector UDum(U[cellI].x(), U[cellI].y(), U[cellI].z()); 
	ma[cellI] = projectedArea[cellI]*rhok[cellI]*mag(UDum); 
      }
    }
   }
  
  scalar ma_Tot = gSum(ma)/30;
  return ma_Tot;
}

scalar HeatExchanger::EnergySource(const fvMesh &iMesh, scalar &rho, scalar &Tup, scalar &Tdown, const double &volume, scalar &P, scalar &phi, const vector &U, const scalar iHEExLeft, const scalar iHEExRight, scalar &T, scalar &ma_Tot) {
  
  scalar Tai1 = 0;
  scalar Tao1 = 0;
  scalar Tao2 = 0;
  scalar Ta1 = 0;
  scalar Ta2 = 0;
  scalar UHE = 0;
        
  vector UDum(U.x(), U.y(), U.z());
  
  if(mag(UDum) == 0) { 
    UDum.x() = 1;
    UDum.y() = 1;
    UDum.z() = 1;
  }
  
  if(ma_Tot < 300) { 
    ma_Tot = 300;
  }
  
  if(iHEExRight > 0 && iHEExLeft == 0) { 
    UHE = -U.z()*Foam::cos(HEApexAngle*deg2Rad) + U.x()*Foam::sin(HEApexAngle*deg2Rad);
  } else { 
    if(iHEExLeft > 0 && iHEExRight == 0) { 
      UHE = -U.z()*Foam::cos(HEApexAngle*deg2Rad) - U.x()*Foam::sin(HEApexAngle*deg2Rad);
    }
  }

  if(UHE <= 0) { 
    Tai1 = Tdown;
    Tao2 = Tup;	
  } else {
    Tai1 = Tup;
    Tao2 = Tdown;
  }
          
  for (int i = 0; i <= 5; i++) { 
    
    Ta1 = (Tao1 + Tai1)/2;
    
    //Properties evaluated at inlet temperature
    scalar mu_air1 = Properties_DAmu(Ta1);
    scalar cp_air1 = Properties_DAcp(Ta1);
    scalar k_air1 = Properties_DAk(Ta1);
    scalar rho_air1 = Properties_DArho(Ta1,P);
      
    scalar Pr_air1 = cp_air1*mu_air1/k_air1;
    scalar Ry1 = rho_air1*mag(UDum)/mu_air1;
    scalar Ny1 = Ny1_a*pow(Ry1, Ny1_b);
	    
    scalar hA1 = k_air1*pow(Pr_air1, 0.333)*Afr*ntb1/ntb2*Ny1*nb;
     
    scalar UA1 = 1/(1/hA1);
    
    scalar e1 = 1-Foam::exp(-UA1/(ma_Tot*cp_air1));
    
    Tao1 = e1*T_steam + (1-e1)*Tai1;
    
    Ta2 = (Tao2 + Tao1)/2;

    //Properties evaluated at outlet temperature
    scalar mu_air2 = Properties_DAmu(Ta2);
    scalar cp_air2 = Properties_DAcp(Ta2);
    scalar k_air2 = Properties_DAk(Ta2);
    scalar rho_air2 = Properties_DArho(Ta2, P);
    
    scalar Pr_air2 = cp_air2*mu_air2/k_air2;
    scalar Ry2 = rho_air2*mag(UDum)/mu_air2;
    scalar Ny2 = Ny2_a*pow(Ry2, Ny2_b);
		    
    scalar hA2 = k_air2*pow(Pr_air2, 0.333)*Afr*Ny2*nb;
    scalar UA2 = 1/(1/hA2 + 1/(steamhA*Afr*steamFrontalAreaRatio));
    scalar e2 = 1-Foam::exp(-UA2/(ma_Tot*cp_air2));

    Tao2 = e2*T_steam + (1-e2)*Tao1;
  }
  
  scalar T_average = (Tai1 + Tao2)/2;
        
   scalar energySource = Properties_DArho(T_average,P) * mag(UDum) * Properties_DAcp(T_average) * (Tao2 - Tai1)/HEWidth;
 
  return energySource;
  
}
    
  scalar HeatExchanger::Properties_DAmu(scalar &T) { 
      T = max(T, 220);
      T = min(T, 380);    
      scalar da_mu = 2.287973e-6+(6.259793e-8*T)-(3.131956e-11*pow(T,2))+(8.15038e-15*pow(T,3));
      
      return da_mu;
  }

  scalar HeatExchanger::Properties_DArho(scalar &T, scalar &P) { 
    T = max(T, 220);
    T = min(T, 380);    
    scalar rho = P/(287.08*T);
    
    return rho;
  }
  
  scalar HeatExchanger::Properties_DAcp(scalar &T) { 
    T = max(T, 220);
    T = min(T, 380);   
    scalar da_cp =1045.356-(0.3161783*T)+(0.0007083814*pow(T,2))-(0.0000002705209*pow(T,3));
    
    return da_cp;

  }
  
  scalar HeatExchanger::Properties_DAk(scalar &T) { 
    T = max(T, 220);
    T = min(T, 380);
    scalar da_k = -0.0004937787+(0.0001018087*T)-(0.00000004627937*pow(T,2))+(0.00000000001250603*pow(T,3));
    
    return da_k;
  }
  
  tensor HeatExchanger::axesRotation() {
    
    vector axis1(Foam::cos(HEApexAngle*deg2Rad), 0, Foam::sin(HEApexAngle*deg2Rad));
    vector axis2(0, 1, 0);

    vector a = axis1/mag(axis1);
    vector b = axis2;

    b = b - (b & a)*a;

    b = b/mag(b);
    vector c = a^b;

    Rglob_loc = tensor(a, b, c);

    // the global->local transformation
    Rglob_loc = Rglob_loc;
    
    // the local->global transformation
    Rloc_glob = Rglob_loc.T();
    
    return Rglob_loc;
  }
}
    
   
    
    
    
    
    
    
  
 