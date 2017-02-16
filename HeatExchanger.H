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
    ActuatorDiskModel_

Description
    Actuator Disk Model implemented on a simplified mesh

\*---------------------------------------------------------------------------*/

#ifndef HEATEX
#define HEATEX

#include "fvCFD.H"
#include "polyMesh.H"
#include <vector>

namespace Foam {
    
  class HeatExchanger {
  public: 
    
    TypeName("HeatExchanger");
    HeatExchanger();
    virtual ~HeatExchanger();
    
    void readDataHE(const fvMesh &iMesh);
    
    void TheoreticalPressureDrop(const fvMesh &iMesh, scalar &Tup, scalar &Tdown, const vector &U, scalar &HeatImForce, scalar &p, const scalar &rho, scalar &KthetaTField, const scalar &mu);
    
    void PressureDrop(const fvMesh &iMesh, const volScalarField &HEExLeft, const volScalarField &HEExRight, scalar &rho, const volVectorField &U, volScalarField &HEImForce, const dimensionedScalar rhoDum, const dimensionedScalar muDum, const dimensionedScalar UDum, const dimensionedScalar SpDum, const volScalarField &KthetaTField, const dimensionedScalar SuDum);
    
    scalar massFlowRate(const fvMesh &iMesh, volScalarField &ma, const volScalarField &projectedArea, const volVectorField &U, const volScalarField &rhok, const volScalarField &HEEnergyLeft, const volScalarField &HEEnergyRight);
    
    scalar EnergySource(const fvMesh &iMesh, scalar &rho, scalar &Tup, scalar &Tdown, const double &volume, scalar &P, scalar &phi, const vector &U, const scalar iHEExLeft, const scalar iHEExRight, scalar &T, scalar &ma_Tot);
    
    scalar Properties_DAmu(scalar &T);
    scalar Properties_DArho(scalar &T, scalar &P);
    scalar Properties_DAcp(scalar &T);
    scalar Properties_DAk(scalar &T);
    
    tensor axesRotation();
    
    scalar Kci, Kdj, Ko, Kturn, Kup, Kdown, Khe, Kts, Ny1_a, Ny1_b, steamhA, steamFrontalAreaRatio, Ny2_a, Ny2_b, ntb1, ntb2;
    
    tensor Rloc_glob, Rglob_loc;
   
    //Circle Values
    const static scalar Pi = 3.141592654;
    const static scalar deg2Rad = 3.141592654/180; 
    const static scalar rad2Deg = 180/3.141592654;
    
    scalar HEApexAngle, meanHEIncidence, HEContractionRatio, HESigma, HESigma_21, a, b, Afr, nb, HECellWidth, HECellNo, T_steam;
    scalar HEWidth;
  private:
    
    
  };
  
}//end namespace Foam
#endif
