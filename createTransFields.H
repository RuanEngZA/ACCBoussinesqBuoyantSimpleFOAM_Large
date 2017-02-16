Info << "Reading fields related to Transport Equations" << endl;

volVectorField UTForward_Left
    (
        IOobject
        (
            "UTForward_Left",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UTForward_Left", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(-0.766044,0.0,-0.766044)), 
	BoundaryDummy.boundaryField().types()
    );
    
    volScalarField HeLeftUpstreamTrans
    (
        IOobject
        (
            "HeLeftUpstreamTrans",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("HeLeftUpstreamTrans", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
	BoundaryDummy.boundaryField().types()
    );
        
    volVectorField UTBackward_Left
    (
        IOobject
        (
            "UTBackward_Left",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UTBackward_Left", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.766044,0.0,0.766044)), 
	BoundaryDummy.boundaryField().types()
    );
    
    volScalarField HeLeftDownstreamTrans
    (
        IOobject
        (
            "HeLeftDownstreamTrans",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("HeLeftDownstreamTrans", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
	BoundaryDummy.boundaryField().types()
    );
    
    volVectorField UTForward_Right
    (
        IOobject
        (
            "UTForward_Right",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UTForward_Right", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0.766044,0.0,-0.766044)), 
	BoundaryDummy.boundaryField().types()
    );
    
    volScalarField HeRightUpstreamTrans
    (
        IOobject
        (
            "HeRightUpstreamTrans",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("HeRightUpstreamTrans", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
	BoundaryDummy.boundaryField().types()
    );
    
    volVectorField UTBackward_Right
    (
        IOobject
        (
            "UTBackward_Right",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UTBackward_Right", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(-0.766044,0.0,-0.766044)), 
	BoundaryDummy.boundaryField().types()
    );
    
    volScalarField HeRightDownstreamTrans
    (
        IOobject
        (
            "HeRightDownstreamTrans",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("HeRightDownstreamTrans", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
	BoundaryDummy.boundaryField().types()
    );
    
    volScalarField Left_for
    (
      IOobject
      (
	"Left_for",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    
    volScalarField Left_back
    (
      IOobject
      (
	"Left_back",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    
    volScalarField Right_for
    (
      IOobject
      (
	"Right_for",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    
    volScalarField Right_back
    (
      IOobject
      (
	"Right_back",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    
    Info << "Reading field DatStop\n" << endl;
    volScalarField DatStopHE
    (
      IOobject
      (
	"DatStopHE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    
    dimensionedScalar DSHE
    (
      "DSHE", 
      dimensionSet(0, 2, -1, 0, 0, 0, 0), 
      1e-8
    ); 
    
    Info << "Reading surfaceScalarField phi\n"<<endl;
    surfaceScalarField phiForT_Left
    (
	IOobject
	(
	    "phiForT_Left",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UTForward_Left) & mesh.Sf()
    );
    
    surfaceScalarField phiBackT_Left
    (
	IOobject
	(
	    "phiBackT_Left",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UTBackward_Left) & mesh.Sf()
    );
    
    surfaceScalarField phiForT_Right
    (
	IOobject
	(
	    "phiForT_Right",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UTForward_Right) & mesh.Sf()
    );
    
    surfaceScalarField phiBackT_Right
    (
	IOobject
	(
	    "phiBackT_Right",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(UTBackward_Right) & mesh.Sf()
    );
    
    volScalarField UTransportLeft
    (
      IOobject
      (
	"UTransportLeft",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("UTransportLeft", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
      BoundaryDummy.boundaryField().types()
    );
    
    volScalarField UTransportRight
    (
      IOobject
      (
	"UTransportRight",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("UTransportRight", dimensionSet(0,1,-1,0,0,0,0), 0.0), 
      BoundaryDummy.boundaryField().types()
    );