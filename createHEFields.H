Info << "Reading Fields Related to Heat Exchanger" << endl;
    volScalarField HEID
	(
	  IOobject
	  (
	    "HEID",
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::AUTO_WRITE
	  ),
	mesh
	);
	
    volScalarField HEExLeft
    (
      IOobject
      (
	"HEExLeft",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
      ),
    mesh
    );
    
    volScalarField HEExRight
    (
      IOobject
      (
	"HEExRight",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
      ),
    mesh
    );
	
    Info << "Reading HeatExchangerForce" << endl;
    volVectorField HEExForce
      (
	IOobject
	(
	  "HEExForce",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedVector("HEExForce", dimensionSet(0,1,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
      );
      
      volVectorField HEExPres
      (
	IOobject
	(
	  "HEExPres",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedVector("HEExPres", dimensionSet(1,-2,-2,0,0,0,0), Foam::vector(0.0,0.0,0.0)), 
      BoundaryDummy.boundaryField().types()
      );
      
      volScalarField HEImForce
      (
	IOobject
	(
	  "HEImForce",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("HEImForce", dimensionSet(0,0,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      volScalarField HEEnergyLeft
      (
	IOobject
	(
	  "HEEnergyLeft",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("HEEnergyLeft", dimensionSet(1,-3,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      volScalarField HEEnergyRight
      (
	IOobject
	(
	  "HEEnergyRight",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("HEEnergyRight", dimensionSet(1,-3,-1,0,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      
      volScalarField energySource
      (
	IOobject
	(
	  "energySource",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("energySource", dimensionSet(0,0,-1,1,0,0,0), 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      //Meyer's Implementation for flow conditioning 
      dimensionedTensor heLeftDExpl
      (
	"heLeftDExpl",
       dimensionSet (0, 0, -1, 0, 0, 0, 0),
       tensor (0.0, 0.0, 0.0, 0.0, 1e7, 0.0, 0.0, 0.0, 1e7)
      );
      
      dimensionedTensor heLeftFExpl
      (
	"heLeftFExpl",
       dimensionSet (0, -1, 0, 0, 0, 0, 0),
       tensor (0.0, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 100)
      );
      
      dimensionedTensor heRightDExpl
      (
	"heRightDExpl",
       dimensionSet (0, 0, -1, 0, 0, 0, 0),
       tensor (1e7, 0.0, 0.0, 0.0, 1e7, 0.0, 0.0, 0.0, 0.0)
      );
      
      dimensionedTensor heRightFExpl
      (
	"heRightFExpl",
       dimensionSet (0, -1, 0, 0, 0, 0, 0),
       tensor (100, 0.0, 0.0, 0.0, 100, 0.0, 0.0, 0.0, 0)
      );
                      
      volScalarField KthetaTField
      (
	IOobject
	(
	  "KthetaTField",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("KthetaTField", dimless, 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      volScalarField ma
      (
	IOobject
	(
	  "ma",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
      mesh,
      dimensionedScalar("ma", dimless, 0), 
      BoundaryDummy.boundaryField().types()
      );
      
      volScalarField projectedArea
      (
	IOobject
	(
	  "projectedArea",
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::AUTO_WRITE
	),
	mesh
      );
                  
      dimensionedScalar muDum
      (
	"muDum",
	dimensionSet (1, -1, -1, 0, 0, 0, 0),
	scalar (1.0)
      );
      
      dimensionedScalar rhoDum
      (
	"rhoDum",
	dimensionSet (1, -3, 0, 0, 0, 0, 0),
	scalar (1.0)
      );
      
      dimensionedScalar UDum
      (
	"UDum",
	dimensionSet (0, 1, -1, 0, 0, 0, 0),
	scalar (1.0)
      );
      
      dimensionedScalar SpDum
      (
	"SpDum",
	dimensionSet (0, 0, -1, 0, 0, 0, 0),
	scalar (1.0)
      );
      
      dimensionedScalar SuDum
      (
	"SuDum",
	dimensionSet (1, -2, -2, 0, 0, 0, 0),
	scalar (1.0)
      );