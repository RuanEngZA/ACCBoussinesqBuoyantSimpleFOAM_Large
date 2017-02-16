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
    ActuatorDiskModel

Description
    Calculates the corresponding forces imparted by the blade onto the fluid
\*---------------------------------------------------------------------------*/

#include "ActuatorDisk.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include <fstream>
#include <stdlib.h> 

namespace Foam {
  
  defineTypeNameAndDebug(ActuatorDiskModel, 0);
  ActuatorDiskModel::ActuatorDiskModel(){}	//Default constructor
  ActuatorDiskModel::~ActuatorDiskModel(){}	//Default destructor
    
    void ActuatorDiskModel::readDataB2(const fvMesh &iMesh) { 
  
  Istream& is1 = iMesh.solutionDict().subDict("B2Fan").lookup("AcDiskLocation");
  is1.format(IOstream::ASCII); 
  is1 >> GenesisAcDiskLocation;
  
  Istream& is2 = iMesh.solutionDict().subDict("B2Fan").lookup("EPS");
  is2.format(IOstream::ASCII); 
  is2 >> EPS;
  
  Istream& is3 = iMesh.solutionDict().subDict("B2Fan").lookup("HubSize");
  is3.format(IOstream::ASCII); 
  is3 >> HubSize;
  
  Istream& is4 = iMesh.solutionDict().subDict("B2Fan").lookup("AcDiskWidth");
  is4.format(IOstream::ASCII); 
  is4 >> AcDiskWidth;
  
  Istream& is5 = iMesh.solutionDict().subDict("B2Fan").lookup("FlowDirection");
  is5.format(IOstream::ASCII); 
  is5 >> FlowDir;
  
  Istream& is6 = iMesh.solutionDict().subDict("B2Fan").lookup("FanDiameter");
  is6.format(IOstream::ASCII); 
  is6 >> FanDiameter;
  
  Istream& is7 = iMesh.solutionDict().subDict("B2Fan").lookup("ChordLength");
  is7.format(IOstream::ASCII); 
  is7 >> ChordLength;
  
  Istream& is8 = iMesh.solutionDict().subDict("B2Fan").lookup("DiskSize");
  is8.format(IOstream::ASCII); 
  is8 >> DiskSize;
  
  Istream& is9 = iMesh.solutionDict().subDict("B2Fan").lookup("UpDiskLocation");
  is9.format(IOstream::ASCII); 
  is9 >> GenesisUpDiskLocation;
  
  Istream& is10 = iMesh.solutionDict().subDict("B2Fan").lookup("DownDiskLocation");
  is10.format(IOstream::ASCII); 
  is10 >> GenesisDownDiskLocation;
  
  Istream& is11 = iMesh.solutionDict().subDict("B2Fan").lookup("GamaHub");
  is10.format(IOstream::ASCII); 
  is11 >> GamaHubB;
  
  Istream& is12 = iMesh.solutionDict().subDict("B2Fan").lookup("FanSpeed");
  is10.format(IOstream::ASCII); 
  is12 >> FanSpeed;
  
  Istream& is13 = iMesh.solutionDict().subDict("B2Fan").lookup("BladeNo");
  is10.format(IOstream::ASCII); 
  is13 >> BladeNo;
  
      Sout << "**B2Fan Properties** " << endl;
      Sout << "Fan Diameter: " << FanDiameter << endl;
      Sout << "Hub Size: " << HubSize << endl;
      Sout << "Chord Length: " << ChordLength  << endl;
      Sout << "Fan Speed: " << FanSpeed << endl;
      Sout << "Numbers of Blade: " << BladeNo << endl;
      Sout << "Blade Angle at Hub: " << GamaHubB << "\n" << endl;
      
      Sout << "**Actuator Disk Properties** " << endl;
      Sout << "First Actuator Disk Location: " << GenesisAcDiskLocation << endl;
      Sout << "First Upstream Disk Location: " << GenesisUpDiskLocation << endl;
      Sout << "First Downstream Disk Location: " << GenesisDownDiskLocation << endl;
      Sout << "Actuator Disk Width: " << AcDiskWidth << endl;
      Sout << "Flow Direction: " << FlowDir << endl;
      Sout << "Tolerance: " << EPS << "\n" << endl;
        }
  
  void ActuatorDiskModel::readDataA(const fvMesh &iMesh) { 
  
  Istream& is1 = iMesh.solutionDict().subDict("AFan").lookup("AcDiskLocation");
  is1.format(IOstream::ASCII); 
  is1 >> GenesisAcDiskLocationA;
    
  Istream& is3 = iMesh.solutionDict().subDict("AFan").lookup("HubSize");
  is3.format(IOstream::ASCII); 
  is3 >> HubSizeA;
      
  Istream& is6 = iMesh.solutionDict().subDict("AFan").lookup("FanDiameter");
  is6.format(IOstream::ASCII); 
  is6 >> FanDiameterA;
  
  Istream& is7 = iMesh.solutionDict().subDict("AFan").lookup("ChordLength");
  is7.format(IOstream::ASCII); 
  is7 >> ChordLengthA;
    
  Istream& is9 = iMesh.solutionDict().subDict("AFan").lookup("UpDiskLocation");
  is9.format(IOstream::ASCII); 
  is9 >> GenesisUpDiskLocationA;
  
  Istream& is10 = iMesh.solutionDict().subDict("AFan").lookup("DownDiskLocation");
  is10.format(IOstream::ASCII); 
  is10 >> GenesisDownDiskLocationA;
  
  Istream& is11 = iMesh.solutionDict().subDict("AFan").lookup("GamaHub");
  is10.format(IOstream::ASCII); 
  is11 >> GamaHubA;
  
  Istream& is12 = iMesh.solutionDict().subDict("AFan").lookup("FanSpeed");
  is10.format(IOstream::ASCII); 
  is12 >> FanSpeedA;
  
  Istream& is13 = iMesh.solutionDict().subDict("AFan").lookup("BladeNo");
  is10.format(IOstream::ASCII); 
  is13 >> BladeNoA;
  
      Sout << "**AFan Properties** " << endl;
      Sout << "Fan Diameter: " << FanDiameterA << endl;
      Sout << "Hub Size: " << HubSizeA << endl;
      Sout << "Chord Length: " << ChordLengthA  << endl;
      Sout << "Fan Speed: " << FanSpeedA << endl;
      Sout << "Numbers of Blade: " << BladeNoA << endl;
      Sout << "Blade Angle at Hub: " << GamaHubA << "\n" << endl;
      
      Sout << "**Actuator Disk Properties** " << endl;
      Sout << "First Actuator Disk Location: " << GenesisAcDiskLocationA << endl;
      Sout << "First Upstream Disk Location: " << GenesisUpDiskLocationA << endl;
      Sout << "First Downstream Disk Location: " << GenesisDownDiskLocationA << endl;
  }
  
  vector ActuatorDiskModel::getRadius(const vector &iMeshC) { 
    /*--------------------------------------------------------------------------
			      Get Radius
     Calculate the radius of the cells and return as a scalar
    --------------------------------------------------------------------------*/  
  
   vector fanCentre(GenesisAcDiskLocation);
   
    //Find Radius of element 
    scalar radiusX = iMeshC.x() - fanCentre.x();
    scalar radiusY = iMeshC.y() - fanCentre.y();
    scalar radiusZ = iMeshC.z() - fanCentre.z();
    
    vector radius(radiusX, radiusY, radiusZ);
    return (radius);
  }
  
    vector ActuatorDiskModel::getRadiusA(const vector &iMeshC) { 
    /*--------------------------------------------------------------------------
			      Get Radius
     Calculate the radius of the cells and return as a scalar
    --------------------------------------------------------------------------*/  
  
   vector fanCentre(GenesisAcDiskLocationA);
   
    //Find Radius of element 
    scalar radiusX = iMeshC.x() - fanCentre.x();
    scalar radiusY = iMeshC.y() - fanCentre.y();
    scalar radiusZ = iMeshC.z() - fanCentre.z();
    
    vector radius(radiusX, radiusY, radiusZ);
    return (radius);
  }
  
  scalar ActuatorDiskModel::getTheta(const vector &iMeshC, vector &DiskRad) {
    /*--------------------------------------------------------------------------
			      Get Theta
     Calculates theta using ArcTan2 function. Moving from +X in anti-clockwise 
     direction.
    --------------------------------------------------------------------------*/ 
  
    //Find angular distance from 0 rads of element
    scalar theta;
    
    theta = Foam::atan2(DiskRad.y(), DiskRad.x());
    //theta = theta + Pi;
    if(theta < 0.0) { 
      theta = theta + 2*Pi;
    }
        
    return(theta);
  
  }
  
  scalar ActuatorDiskModel::getBladeAngleB(const vector &radiusVec) { 
    /*--------------------------------------------------------------------------
			      Get Blade Angle
     Calculates blade angle at radial station using empirical data. 
     Value is added to initial blade angle setting at hub.
    --------------------------------------------------------------------------*/   
  scalar radius = mag(radiusVec);
  scalar rR = 2*radius/FanDiameter;
  scalar GamaHubRad = GamaHubB*deg2Rad;
  
  if(rR > 1) { 
    Sout << "BFan Error: Radius ratio is larger than 1\n" << endl;
  }
  
  scalar gamaTemp = (- 5.9120646e1*rR*rR*rR   + 1.7184723e2*rR*rR  - 1.8066539e2*rR + 7.8602940e1 - 30)*deg2Rad;
  scalar gamaAct = gamaTemp + GamaHubRad;
     
  return(gamaAct);
  
}

  scalar ActuatorDiskModel::getChordLengthB (const vector &radiusVec) { 
    /*--------------------------------------------------------------------------
			      Get Chord Length
     Calculates Chord length at each radial station using empirical data.
    --------------------------------------------------------------------------*/     
    scalar radius = mag(radiusVec);
    scalar rR = 2*radius/FanDiameter;
    
    if(rR > 1) { 
      Sout << "BFan Error: Radius ratio is larger than 1\n" << endl;
    }
    
    scalar chordLengthTemp = (-5.16667E+01*rR + 2.04667E+02)/1000*(FanDiameter/(2*0.771)); 
    return (chordLengthTemp);
    
  }
  
  scalar ActuatorDiskModel::getBladeAngleA(const vector &radiusVec) { 
    /*--------------------------------------------------------------------------
			      Get Blade Angle
     Calculates blade angle at radial station using empirical data. 
     Value is added to initial blade angle setting at hub.
    --------------------------------------------------------------------------*/   
  scalar radius = mag(radiusVec);
  scalar rR = 2*radius/FanDiameter;
  scalar GamaHubRad = GamaHubA*deg2Rad;
  
  if(rR > 1) { 
    Sout << "AFan Error: Radius ratio is larger than 1\n" << endl;
  }
  
  scalar gamaTemp = (79.6662498129*rR*rR*rR*rR - 222.5210124696*rR*rR*rR + 234.3688530132*rR*rR - 133.8034253129*rR + 42.1651302853)*deg2Rad;
//   scalar gamaTemp = (90.8936000237*rR*rR*rR*rR - 251.3832092838*rR*rR*rR + 260.8724014201*rR*rR - 144.0854890221*rR + 46.1165586222)*deg2Rad;
  scalar gamaAct = gamaTemp + GamaHubRad;
     
  return(gamaAct);
  }
  
  scalar ActuatorDiskModel::getClA_fan(const scalar alpha, scalar &ClField) { 
    
    scalar Cl = 0;
    if(alpha > -180 && alpha <= -90) { 
      Cl = 5.46805129531068e-8*alpha*alpha*alpha*alpha + 2.95403178508562e-5*alpha*alpha*alpha + 0.005385767*alpha*alpha + 0.3771675547*alpha + 8.2696835957;
    }
    if(alpha > -90 && alpha <= -20) { 
      Cl = 3.76837906764838e-8*alpha*alpha*alpha*alpha + 1.00486461185167e-5*alpha*alpha*alpha + 0.0011700772*alpha*alpha + 0.0531441458*alpha + 0.0232915317;
    } 
    if(alpha > -20 && alpha <= 16) { 
      Cl = 2.72237063959306e-6*alpha*alpha*alpha*alpha - 0.0001577441*alpha*alpha*alpha - 0.0009522191*alpha*alpha + 0.1201661795*alpha + 0.516735376;
    }  
    if(alpha > 16 && alpha <= 90) {
      Cl = 1.28790818398648e-7*alpha*alpha*alpha*alpha - 3.20497440300299e-5*alpha*alpha*alpha + 0.0027330317*alpha*alpha - 0.1118481625*alpha + 2.9207100635;
    }
    if(alpha > 90 && alpha < 180) { 
      Cl = - 5.46805129531068e-8*alpha*alpha*alpha*alpha + 2.95403178508562e-5*alpha*alpha*alpha - 0.005385767*alpha*alpha + 0.3771675547*alpha - 8.2696835957;
    }
    ClField = Cl;
    return Cl;
  }
  
  scalar ActuatorDiskModel::getCdA_fan(const scalar alpha, scalar &CdField) { 
    
    scalar Cd = 0;
    if(alpha > -180 && alpha <= -90) { 
      Cd = -4.78392499918153e-9*alpha*alpha*alpha*alpha - 8.73439832177165e-6*alpha*alpha*alpha - 0.0030080502*alpha*alpha - 0.3474974626*alpha - 11.0022320603;
    }
    if(alpha > -90 && alpha <= -14) { 
      Cd = -0.000000081*alpha*alpha*alpha*alpha  - 1.60647353960824e-5*alpha*alpha*alpha - 0.0012103674*alpha*alpha - 0.0595391282*alpha - 0.5632715352;
    } 
    if(alpha > -14 && alpha <= 14) { 
      Cd = 1.064922920829e-6*alpha*alpha*alpha*alpha - 1.17938166788587e-6*alpha*alpha*alpha - 1.33401909287259e-5*alpha*alpha - 0.0002556615*alpha + 0.0216712551;
    }  
    if(alpha > 14 && alpha <= 90) {
      Cd = -2.49009180487414e-8*alpha*alpha*alpha*alpha + 4.38241449538984e-6*alpha*alpha*alpha - 0.0004623645*alpha*alpha + 0.0483843997*alpha - 0.6125908767;
    }
    if(alpha > 90 && alpha < 180) { 
      Cd = -4.78392499918153e-9*alpha*alpha*alpha*alpha + 8.73439832177165e-6*alpha*alpha*alpha - 0.0030080502*alpha*alpha + 0.3474974626*alpha - 11.0022320603;
    }
    CdField = Cd;
    return Cd;
  }
  
  
  void ActuatorDiskModel::readCSV_NASA0413LS() { 
    /*--------------------------------------------------------------------------
			      Read NASA 0413-LS2 CSV Data
     Reads CSV files containing lift and drag coefficient data for Reynolds numbers ranging from 2e6 to 9e6.
     IMPORTANT: Ensure correct file path for each data set
    --------------------------------------------------------------------------*/   
  std::ifstream file, file1, file2, file3;
  Sout << "\nReading CSV files" << endl;
  file.open("/home/ruan/Dropbox/Masters/Model_Development/Solvers/ACC/ACCBoussinesqBuoyantSimpleFoam_Large/NASA0413LS/Re2e6.csv");
  file1.open("/home/ruan/Dropbox/Masters/Model_Development/Solvers/ACC/ACCBoussinesqBuoyantSimpleFoam_Large/NASA0413LS/Re4e6.csv");
  file2.open("/home/ruan/Dropbox/Masters/Model_Development/Solvers/ACC/ACCBoussinesqBuoyantSimpleFoam_Large/NASA0413LS/Re6e6.csv");
  file3.open("/home/ruan/Dropbox/Masters/Model_Development/Solvers/ACC/ACCBoussinesqBuoyantSimpleFoam_Large/NASA0413LS/Re9e6.csv");
  string alpha, cl, cd, alpha1, cl1, cd1, alpha2, cl2, cd2, alpha3, cl3, cd3;
  char* sz;
  int counter = 0;
  
  Sout << "NASA 0413-LS2 2e6" << endl;
  while (file.good()) {
    getline(file, alpha , ',');
    getline(file, cl , ',');
    getline(file, cd , '\n');
    
    if(!alpha.empty() && !cl.empty() && !cd.empty()) {
      alpha_2e6_Arr_0413.push_back(strtod(alpha.c_str(), &sz));
      cl_2e6_Arr_0413.push_back(strtod(cl.c_str(), &sz));
      cd_2e6_Arr_0413.push_back(strtod(cd.c_str(), &sz));
   
      if(debug >= 2) { 
	scalar alphaDisp = alpha_2e6_Arr_0413[counter];
	scalar clDisp = cl_2e6_Arr_0413[counter];
	scalar cdDisp = cd_2e6_Arr_0413[counter];
	Sout << "Alpha at 2e6: " << alphaDisp << " Cl: " << clDisp << " Cd: " << cdDisp << endl;
	counter++;
      }
    } 
  }
  file.close();
  
  counter = 0;
  Sout << "\nNASA 0413-LS2 4e6" << endl;
  while(file1.good()) { 
    getline(file1, alpha1 , ',');
    getline(file1, cl1 , ',');
    getline(file1, cd1 , '\n');
    
    if(!alpha1.empty() && !cl1.empty() && !cd1.empty()) {
      alpha_4e6_Arr_0413.push_back(strtod(alpha1.c_str(), &sz));
      cl_4e6_Arr_0413.push_back(strtod(cl1.c_str(), &sz));
      cd_4e6_Arr_0413.push_back(strtod(cd1.c_str(), &sz));
   
      if(debug >= 2) { 
	scalar alphaDisp = alpha_4e6_Arr_0413[counter];
	scalar clDisp = cl_4e6_Arr_0413[counter];
	scalar cdDisp = cd_4e6_Arr_0413[counter];
	Sout << "Alpha at 4e6: " << alphaDisp << " Cl1: " << clDisp << " Cd1: " << cdDisp << endl;
	counter++;
      }
    } 
  }
  file1.close();
  
  counter = 0;
  Sout << "\nNASA 0413-LS2 6e6" << endl;
  while(file2.good()) { 
    getline(file2, alpha2 , ',');
    getline(file2, cl2 , ',');
    getline(file2, cd2 , '\n');
    
    if(!alpha2.empty() && !cl2.empty() && !cd2.empty()) {
      alpha_6e6_Arr_0413.push_back(strtod(alpha2.c_str(), &sz));
      cl_6e6_Arr_0413.push_back(strtod(cl2.c_str(), &sz));
      cd_6e6_Arr_0413.push_back(strtod(cd2.c_str(), &sz));
   
      if(debug >= 2) { 
	scalar alphaDisp = alpha_6e6_Arr_0413[counter];
	scalar clDisp = cl_6e6_Arr_0413[counter];
	scalar cdDisp = cd_6e6_Arr_0413[counter];
	Sout << "Alpha at 6e6: " << alphaDisp << " Cl1: " << clDisp << " Cd1: " << cdDisp << endl;
	counter++;
      }
    } 
  }
  file2.close();
  
  counter = 0;
  Sout << "\nNASA 0413-LS2 9e6" << endl;
  while(file3.good()) { 
    getline(file3, alpha3 , ',');
    getline(file3, cl3 , ',');
    getline(file3, cd3 , '\n');
    
    if(!alpha3.empty() && !cl3.empty() && !cd3.empty()) {
      alpha_9e6_Arr_0413.push_back(strtod(alpha3.c_str(), &sz));
      cl_9e6_Arr_0413.push_back(strtod(cl3.c_str(), &sz));
      cd_9e6_Arr_0413.push_back(strtod(cd3.c_str(), &sz));
   
      if(debug >= 2) { 
	scalar alphaDisp = alpha_9e6_Arr_0413[counter];
	scalar clDisp = cl_9e6_Arr_0413[counter];
	scalar cdDisp = cd_9e6_Arr_0413[counter];
	Sout << "Alpha at 9e6: " << alphaDisp << " Cl1: " << clDisp << " Cd1: " << cdDisp << endl;
	counter++;
      }
    } 
  }
  file3.close();
  
}

scalar ActuatorDiskModel::getLiftNASA0413(const scalar &alpha, const scalar &Re, scalar &ClField) { 
    /*--------------------------------------------------------------------------
			      Get Lift Coefficient
     Calculates lift coefficient using linear interpolation and extrapolation
     for values at different Reynolds values.
    --------------------------------------------------------------------------*/   
  scalar Cl;
  
  if(Re < 2e6) { 
    scalar Cl_ReHigh = getCl(alpha, alpha_4e6_Arr_0413, cl_4e6_Arr_0413);
    scalar Cl_ReLow = getCl(alpha, alpha_2e6_Arr_0413, cl_2e6_Arr_0413);
    
    //Cl and Cd lower extrapolation scheme
    scalar m = (Cl_ReHigh - Cl_ReLow)/(4e6 - 2e6);
    scalar b = Cl_ReHigh - m*4e6;
    Cl = m*Re + b;
    
    if(debug >= 2) {
      Sout << "Cl_ReHigh at 4e6: " << Cl_ReHigh << endl;
      Sout << "Cl_ReLow at 2e6: " << Cl_ReLow << endl;
      Sout << "Cl (Actual): " << Cl << endl;
    }
    ClField = Cl;
    return(Cl);
  }
  
  if(Re >= 2e6 && Re <= 4e6) {
    scalar Cl_ReHigh = getCl(alpha, alpha_4e6_Arr_0413, cl_4e6_Arr_0413);
    scalar Cl_ReLow = getCl(alpha, alpha_2e6_Arr_0413, cl_2e6_Arr_0413);
        
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cl = (Re - 2e6)/(4e6 - 2e6)*(Cl_ReHigh - Cl_ReLow) + Cl_ReLow;
  
    if(debug >= 2) {
      Sout << "Cl_ReHigh at 4e6: " << Cl_ReHigh << endl;
      Sout << "Cl_ReLow at 2e6: " << Cl_ReLow << endl;
      Sout << "Cl (Actual): " << Cl << endl;
    }
    ClField = Cl;
    return(Cl);
  }
  
  if(Re > 4e6 && Re <= 6e6) {
    scalar Cl_ReHigh = getCl(alpha, alpha_6e6_Arr_0413, cl_6e6_Arr_0413);
    scalar Cl_ReLow = getCl(alpha, alpha_4e6_Arr_0413, cl_4e6_Arr_0413);
        
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cl = (Re - 4e6)/(6e6 - 4e6)*(Cl_ReHigh - Cl_ReLow) + Cl_ReLow;
  
    if(debug >= 2) {
      Sout << "Cl_ReHigh at 6e6: " << Cl_ReHigh << endl;
      Sout << "Cl_ReLow at 4e6: " << Cl_ReLow << endl;
      Sout << "Cl (Actual): " << Cl << endl;
    }
    ClField = Cl;
    return(Cl);
  }
  
  if(Re > 6e6 && Re <= 9e6) {
    scalar Cl_ReHigh = getCl(alpha, alpha_9e6_Arr_0413, cl_9e6_Arr_0413);
    scalar Cl_ReLow = getCl(alpha, alpha_6e6_Arr_0413, cl_6e6_Arr_0413);
        
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cl = (Re - 6e6)/(9e6 - 6e6)*(Cl_ReHigh - Cl_ReLow) + Cl_ReLow;
  
    if(debug >= 2) {
      Sout << "Cl_ReHigh at 9e6: " << Cl_ReHigh << endl;
      Sout << "Cl_ReLow at 6e6: " << Cl_ReLow << endl;
      Sout << "Cl (Actual): " << Cl << endl;
    }
    ClField = Cl;
    return(Cl);
  }
  
  if(Re > 9e6) { 
    scalar Cl_ReHigh = getCl(alpha, alpha_4e6_Arr_0413, cl_4e6_Arr_0413);
    scalar Cl_ReLow = getCl(alpha, alpha_2e6_Arr_0413, cl_2e6_Arr_0413);
    
    //Cl and Cd upper extrapolation scheme
    scalar m = (Cl_ReHigh - Cl_ReLow)/(9e6 - 6e6);
    scalar b = Cl_ReHigh - m*9e6;
    Cl = m*Re + b;
    
    if(debug >= 2) {
      Sout << "Cl_ReHigh at 9e6: " << Cl_ReHigh << endl;
      Sout << "Cl_ReLow at 6e6: " << Cl_ReLow << endl;
      Sout << "Cl (Actual): " << Cl << endl;
    }
    ClField = Cl;
    return(Cl);
  } else {
    Sout << "Reynolds Error" << endl;
    Cl = 1;
    ClField = Cl;
    return(Cl);
  }
}

scalar ActuatorDiskModel::getDragNASA0413(const scalar &alpha, const scalar &Re, scalar &CdField) { 
    /*--------------------------------------------------------------------------
			      Get Drag Coefficient
     Calculates drag coefficient using linear interpolation and extrapolation
     for values at different Reynolds values.
    --------------------------------------------------------------------------*/   
  scalar Cd;
  
  if(Re < 2e6) { 
    scalar Cd_ReHigh = getCd(alpha, alpha_4e6_Arr_0413, cd_4e6_Arr_0413);
    scalar Cd_ReLow = getCd(alpha, alpha_2e6_Arr_0413, cd_2e6_Arr_0413);
    
    //Cl and Cd lower extrapolation scheme
    scalar m = (Cd_ReHigh - Cd_ReLow)/(4e6 - 2e6);
    scalar b = Cd_ReHigh - m*4e6;
    Cd = m*Re + b;
    
    if(debug >= 2) {
      Sout << "Cd_ReHigh at 4e6: " << Cd_ReHigh << endl;
      Sout << "Cd_ReLow at 2e6: " << Cd_ReLow << endl;
      Sout << "Cd (Actual): " << Cd << endl;
    }
    CdField = Cd;
    return(Cd);
  }
  if(Re >= 2e6 && Re <= 4e6) {
       
    scalar Cd_ReHigh = getCd(alpha, alpha_4e6_Arr_0413, cd_4e6_Arr_0413);
    scalar Cd_ReLow = getCd(alpha, alpha_2e6_Arr_0413, cd_2e6_Arr_0413);
    
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cd = (Re - 2e6)/(4e6 - 2e6)*(Cd_ReHigh - Cd_ReLow) + Cd_ReLow; 
    
    if(debug >= 2) {
      Sout << "Cd_ReHigh at 4e6: " << Cd_ReHigh << endl;
      Sout << "Cd_ReLow at 2e6: " << Cd_ReLow << endl;
      Sout << "Cd (Actual): " << Cd << endl;
     }
     CdField = Cd;
     return(Cd);
  }
  
  if(Re > 4e6 && Re <= 6e6) {
    scalar Cd_ReHigh = getCd(alpha, alpha_6e6_Arr_0413, cd_6e6_Arr_0413);
    scalar Cd_ReLow = getCd(alpha, alpha_4e6_Arr_0413, cd_4e6_Arr_0413);
        
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cd = (Re - 4e6)/(6e6 - 4e6)*(Cd_ReHigh - Cd_ReLow) + Cd_ReLow;
  
    if(debug >= 2) {
      Sout << "Cd_ReHigh at 6e6: " << Cd_ReHigh << endl;
      Sout << "Cd_ReLow at 4e6: " << Cd_ReLow << endl;
      Sout << "Cd (Actual): " << Cd << endl;
    }
    CdField = Cd;
    return(Cd);
  }
  
  if(Re > 6e6 && Re <= 9e6) {
    scalar Cd_ReHigh = getCd(alpha, alpha_9e6_Arr_0413, cd_9e6_Arr_0413);
    scalar Cd_ReLow = getCd(alpha, alpha_6e6_Arr_0413, cd_6e6_Arr_0413);
        
    //Cl and Cd Interpolation Scheme between Reynolds values
    Cd = (Re - 6e6)/(9e6 - 6e6)*(Cd_ReHigh - Cd_ReLow) + Cd_ReLow;
  
    if(debug >= 2) {
      Sout << "Cd_ReHigh at 9e6: " << Cd_ReHigh << endl;
      Sout << "Cd_ReLow at 6e6: " << Cd_ReLow << endl;
      Sout << "Cd (Actual): " << Cd << endl;
    }
    CdField = Cd;
    return(Cd);
  }
  
  if(Re > 9e6) { 
    scalar Cd_ReHigh = getCd(alpha, alpha_4e6_Arr_0413, cd_4e6_Arr_0413);
    scalar Cd_ReLow = getCd(alpha, alpha_2e6_Arr_0413, cd_2e6_Arr_0413);
    
    //Cl and Cd upper extrapolation scheme
    scalar m = (Cd_ReHigh - Cd_ReLow)/(9e6 - 6e6);
    scalar b = Cd_ReHigh - m*9e6;
    Cd = m*Re + b;
    
    if(debug >= 2) {
      Sout << "Cd_ReHigh at 9e6: " << Cd_ReHigh << endl;
      Sout << "Cd_ReLow at 6e6: " << Cd_ReLow << endl;
      Sout << "Cd (Actual): " << Cd << endl;
    }
    CdField = Cd;
    return(Cd);
  } else {
    Sout << "Reynolds Error" << endl;
    Cd = 1;
    CdField = Cd;
    return(Cd);
  }
}

scalar ActuatorDiskModel::getCl(const scalar &alpha, const std::vector<long double> &alpha_Arr, const std::vector<long double> &cl_Arr) {
    /*--------------------------------------------------------------------------
			      Get Lift at Alpha
     Interpolates between two alpha values to calculate correct lift value
    --------------------------------------------------------------------------*/   
    int alphaArrSize =  alpha_Arr.size();
    int alphaCounter_Low = 0;
    int alphaCounter_High = 0;
    
    for(int i = 1; i < alphaArrSize; i++) { 
      //Determine two points between which alpha lies
      if(alpha >= alpha_Arr[i-1] && alpha <= alpha_Arr[i]) {
	alphaCounter_Low = i-1;
	alphaCounter_High = i;
      }
    }
    if(alphaCounter_Low == 0 && alphaCounter_High == 0) { 
      Sout << "Alpha is out of range" << endl;
      scalar Cl = 1;
      return Cl;
    }
    
    //Interpolation variables for Cl and Cd values between alpha values
    scalar cl_AlphaHigh = cl_Arr[alphaCounter_High];
    scalar cl_AlphaLow = cl_Arr[alphaCounter_Low];
        
    scalar alpha_High = alpha_Arr[alphaCounter_High];
    scalar alpha_Low = alpha_Arr[alphaCounter_Low];
    
    //Cl and Cd Interpolation Scheme implemented between alpha values
    scalar Cl = (alpha - alpha_Low)/(alpha_High - alpha_Low) * (cl_AlphaHigh - cl_AlphaLow) + cl_AlphaLow; 
    
      if(debug >= 2) {
	Sout << "\nalphaHigh: " << alpha_High << " alphaLow: " << alpha_Low << endl;
	Sout << "Cl_AlphaHigh: " << cl_AlphaHigh << " Cl_AlphaLow: " << cl_AlphaLow << endl;
	Sout << "Cl (Alpha Level): " << Cl << endl;  
      }
    return Cl;
}

scalar ActuatorDiskModel::getCd(const scalar &alpha, const std::vector<long double> &alpha_Arr, const std::vector<long double> &cd_Arr) {
      /*--------------------------------------------------------------------------
			      Get Drag at Alpha
     Interpolates between two alpha values to calculate correct drag value
    --------------------------------------------------------------------------*/ 
    int alphaArrSize =  alpha_Arr.size();
    int alphaCounter_Low = 0;
    int alphaCounter_High = 0;
    
    for(int i = 1; i < alphaArrSize; i++) { 
      //Determine two points between which alpha lies
      if(alpha >= alpha_Arr[i-1] && alpha <= alpha_Arr[i]) {
	alphaCounter_Low = i-1;
	alphaCounter_High = i;
      }
    }
    if(alphaCounter_Low == 0 && alphaCounter_High == 0) { 
      Sout << "Alpha is out of range" << endl;
      scalar Cd = 1;
      return Cd;
    }
    
    //Interpolation variables for Cl and Cd values between alpha values
    scalar cd_AlphaHigh = cd_Arr[alphaCounter_High];
    scalar cd_AlphaLow = cd_Arr[alphaCounter_Low];
        
    scalar alpha_High = alpha_Arr[alphaCounter_High];
    scalar alpha_Low = alpha_Arr[alphaCounter_Low];
    
    //Cl and Cd Interpolation Scheme implemented between alpha values
    scalar Cd = (alpha - alpha_Low)/(alpha_High - alpha_Low) * (cd_AlphaHigh - cd_AlphaLow) + cd_AlphaLow; 
    
      if(debug >= 2) {
	Sout << "\nalphaHigh: " << alpha_High << " alphaLow: " << alpha_Low << endl;
	Sout << "Cd_AlphaHigh: " << cd_AlphaHigh << " Cd_AlphaLow: " << cd_AlphaLow << endl;
	Sout << "Cd (Alpha Level): " << Cd << endl;  
      }
    return Cd;
}
  
vector ActuatorDiskModel::calcVolumeForce(const vector &iMeshC, const vector &radiusVec, const vector &U_up, const vector &U_down, const scalar &theta, const scalar &gamma, const scalar &chordLength, vector &relativeVelocity, scalar &betaField, scalar &alphaField, scalar &ClField, scalar &CdField, scalar &Torque, scalar &volume, vector &AverageVelocity, const scalar &diskID, const scalar &rho, const scalar &mu) {
  
  /*---------------------------------------------------------------------------------------------
   * 					Calculate Volume Force
   * Calculates the volume forces to be implemented into the momentum equation for use as source terms
   * Computed in the Axial, Radial and Tangential directions first and transformed back to local coordinate systems as Fx, Fy and Fz
   * --------------------------------------------------------------------------------------------*/
  
  /*===============================================================================================
   *					Calculate relative velocity angle Beta
   * =============================================================================================*/
  scalar beta = 0;
  
  //Calculate average velocity at midspan of the airfoil using the values on the upstream and downstream disks
  scalar Ux = 0.5*U_up.x() + 0.5*U_down.x();
  scalar Uy = 0.5*U_up.y() + 0.5*U_down.y();
  scalar Uz = 0.5*U_up.z() + 0.5*U_down.z();
  
  AverageVelocity.x() = Ux;
  AverageVelocity.y() = Uy;
  AverageVelocity.z() = Uz;
  
  //Compute the position of rotor element relative to rotor centre
  scalar radiusX = radiusVec.x();
  scalar radiusY = radiusVec.y();
  scalar radiusZ = radiusVec.z();
  
  scalar radius = mag(radiusVec);
  
  scalar URADIAL, C1x, Ca;
  
  //Compute X and Y velocity components in cylindrical coordinates of element
  URADIAL = Ux*cos(theta) + sin(theta)*Uy;
  C1x = -sin(theta)*Ux + cos(theta)*Uy;
  Ca = Uz;
  
  //Blade speed
  scalar UFan = (FanSpeed*2*Pi/60)*radius;
  
  //Relative velocity magnitude in 3 dimensions
  scalar W1x = UFan - C1x;
  scalar W1 = sqrt(Ca*Ca + W1x*W1x + URADIAL*URADIAL);
  
  relativeVelocity.x() = W1x;
  relativeVelocity.y() = URADIAL;
  relativeVelocity.z() = Ca;
  
  //Compute the beta angle 
  if(W1x >= 0) { 
    beta = atan(-Ca/W1x);
  }
  if(W1x < 0) {
    beta = Pi + atan(-Ca/W1x);
  }
    
    betaField = beta * rad2Deg;
        
  /*===============================================================================================
   *					    Calculate angle of attack alpha
   * =============================================================================================*/
  scalar alpha = gamma - beta;
  
  if(alpha > Pi) { 
    alpha = alpha - 2*Pi;
  }
  
  if(alpha <= -Pi) { 
    alpha = alpha + 2*Pi;
  }
  
  scalar alpha_deg = alpha*rad2Deg; 	//Convert to degrees 
  alphaField = alpha*rad2Deg;
  
 /*===============================================================================================
  *					      Calculate blade solidity
  * =============================================================================================*/
  scalar sigma = BladeNo*chordLength/(2*Pi*radius);
   
  if(debug >= 1) { 
    Sout << "\nMesh Coordinate: " << iMeshC << endl;
    Sout << "Radius X: " << radiusX << " Y: " << radiusY  << " Z: " << radiusZ << endl;
    Sout << "Theta: " << theta << endl;
    Sout << "Gama: " << gamma << endl;
    Sout << "Chord Length: " << chordLength  << endl;
    Sout << "Solidity: " << sigma << endl;
    Sout << "\nUpstream Velocity X: " << U_up.x() << " Y: " << U_up.y() << " Z: " << U_up.z() << endl;
    Sout << "Downstream Velocity X: " << U_down.x() << " Y: " << U_down.y() << " Z: " << U_down.z() << endl;
    Sout << "\nAverage Velocity in X: " << Ux << " in Y " << Uy << " in Z " << Uz << endl;
    Sout << "Radius relative to centre in X: " << radiusX << " in Y: " << radiusY << " in Z: " << radiusZ << endl;
    Sout << "Radial Velocity: " << URADIAL << endl;
    Sout << "Tangential Velocity: " << C1x << endl;
    Sout << "Axial Velocity: " << Ca << endl;
    Sout << "Circumferential Speed: " << UFan << endl;
    Sout << "Relative Velocity: " << W1 << endl;
    Sout << "Beta: " << betaField << endl;
    Sout << "Alpha: " << alpha_deg << endl;
  }
  
  /*===============================================================================================
   *						Calculate Cl and Cd
   * =============================================================================================*/
  scalar Reynolds = W1*chordLength*rho/mu;
  scalar Cl = 0;
  scalar Cd = 0;
  
  if(diskID == 2) { 
    Cl = getLiftNASA0413(alpha_deg, Reynolds, ClField);
    Cd = getDragNASA0413(alpha_deg, Reynolds, CdField);
  }
  if(diskID == 4) { 
    Cl = getClA_fan(alpha_deg, ClField);
    Cd = getCdA_fan(alpha_deg, CdField);
  }

  /*===============================================================================================
   *					Calculate force source terms
   * =============================================================================================*/
  
  //Local cylindrical coordinate system
  scalar FAXIAL = 0.5*rho*W1*W1*sigma/AcDiskWidth*(Cl*cos(beta) - Cd*sin(beta));
  scalar FTAN = 0.5*rho*W1*W1*sigma/AcDiskWidth*(Cl*sin(beta) + Cd*cos(beta));
  scalar FRAD = 0;
  
  //Transformed to local cartesian coordinate system
  scalar Fx = -sin(theta)*FTAN; 
  scalar Fy = cos(theta)*FTAN;
  scalar Fz = -FAXIAL;
  
  //Calculate Torque 
  scalar FTAN_Vol = FTAN*volume;
  Torque = FTAN_Vol * radius;
    
  if(debug >= 1) {   
    Sout << "\nRe: " << Reynolds << endl;
    Sout << "Cl: " << Cl << endl;
    Sout << "Cd: " << Cd << endl;
    
    Sout << "\nAxial Force: " << FAXIAL << endl;
    Sout << "Tangential Force: " << FTAN << endl;
    Sout << "Radial Force: " << FRAD << endl;
    Sout << "Force in X: " << Fx << endl;
    Sout << "Force in Y: " << Fy << endl;
    Sout << "Force in Z: " << Fz  << endl;
    
    Sout <<"\nCell Volume: " << volume << endl;
    Sout << "Torque: " << Torque << endl;
    Sout <<"====================" << endl;
  }
  
  vector VolumeForce(Fx, Fy, Fz);
  
  return (VolumeForce);
  } 
}  