//  -*- mode:c++; tab-width:4;  -*-
#include "WCSimDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4UserLimits.hh"
#include "G4ReflectionFactory.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "WCSimTuningParameters.hh" //jl145

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


/***********************************************************
 *
 * This file containts the functions which construct 
 * the nuPRISM detector.  It is called in the Construct()
 * method in WCSimDetectorConstruction.cc.
 *
 * This has been taken almost entirely from WCConstructCylinder
 * but we can modify it at will once the nuPRISM detector design
 * becomes more finalised - same as the HyperK class
 *
 * I have removed the Cap construction, as this should be unchanged from the cylinder constructor, so is called from there
 *
 ***********************************************************/


G4LogicalVolume* WCSimDetectorConstruction::ConstructNuPrism()
{
    G4cout << "**** Building NuPrism Detector ****" << G4endl;

    //-----------------------------------------------------
    // Positions
    //-----------------------------------------------------

    debugMode = false;

    WCIDRadius = WCIDDiameter/2.;
    // the number of regular cell in phi direction:
    WCBarrelRingNPhi     = (G4int)(WCBarrelNumPMTHorizontal/WCPMTperCellHorizontal); 
    // the part of one ring, that is covered by the regular cells: 
    totalAngle  = 2.0*pi*rad*(WCBarrelRingNPhi*WCPMTperCellHorizontal/WCBarrelNumPMTHorizontal) ;
    // angle per regular cell:
    dPhi        =  totalAngle/ WCBarrelRingNPhi;
    // it's height:
    barrelCellHeight  = (WCIDHeight-2.*WCBarrelPMTOffset)/WCBarrelNRings;
    // the height of all regular cells together:
    mainAnnulusHeight = WCIDHeight -2.*WCBarrelPMTOffset -2.*barrelCellHeight;

    innerAnnulusRadius = WCIDRadius - WCPMTExposeHeight-1.*mm;
    outerAnnulusRadius = WCIDRadius + WCBlackSheetThickness + 1.*mm;//+ Steel structure etc.
    // the radii are measured to the center of the surfaces
    // (tangent distance). Thus distances between the corner and the center are bigger.
    WCLength    = WCShaftHeight;	//jl145 - reflects top veto blueprint, cf. Farshid Feyzi
    // Add WCTotalDiameter/2. to radius - nuPRISM shaft will be ~10m in diameter, but might change
    WCRadius    = (WCShaftDiameter/2. + WCBlackSheetThickness)/cos(dPhi/2.) ; // TODO: OD 

    // now we know the extend of the detector and are able to tune the tolerance
    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(WCLength > WCRadius ? WCLength : WCRadius);
    G4cout << "Computed tolerance = "
        << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
        << " mm" << G4endl;

    //Decide if adding Gd
    water = "Water";
    if (WCAddGd)
    {water = "Doped Water";}

    //-----------------------------------------------------
    // Volumes
    //-----------------------------------------------------

    // The water barrel is placed in an tubs of air

    G4Tubs* solidWC = new G4Tubs("WC",
            0.0*m,
            WCRadius+1.*m, 
            .5*WCLength+1*m,	//jl145 - per blueprint
            0.*deg,
            360.*deg);

    G4LogicalVolume* logicWC = 
        new G4LogicalVolume(solidWC,
                G4Material::GetMaterial("Air"),
                "WC",
                0,0,0);


    G4VisAttributes* showColor = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
    logicWC->SetVisAttributes(showColor);

    logicWC->SetVisAttributes(G4VisAttributes::Invisible); //amb79

    //-----------------------------------------------------
    // Create water tub for nuPRISM shaft
    //-----------------------------------------------------
    G4Tubs* solidWCShaft = new G4Tubs("WCShaft",
            0.0*m,
            WCRadius,
            .5*WCLength,  //jl145 - per blueprint
            0.*deg,
            360.*deg);

    G4LogicalVolume* logicWCShaft = 
        new G4LogicalVolume(solidWCShaft,
                G4Material::GetMaterial(water),
                "WCShaft",
                0,0,0);

    G4VPhysicalVolume* physiWCShaft = 
        new G4PVPlacement(0,
                G4ThreeVector(0.,0.,0.),
                logicWCShaft,
                "WCShaft",
                logicWC,
                false,
                0); 


    //-----------------------------------------------------
    // Create instrumented region of detector 
    //-----------------------------------------------------
    G4Tubs* solidWCBarrel = new G4Tubs("WCBarrel",
            0.0*m,
            WCIDRadius, 
            .5*WCIDHeight,
            0.*deg,
            360.*deg);

    G4LogicalVolume* logicWCBarrel =
        new G4LogicalVolume(solidWCBarrel,
                G4Material::GetMaterial(water),
                "WCBarrel",
                0,0,0);

    //Places instrumented region at desired position in shaft 
    G4VPhysicalVolume* physiWCBarrel;
    physiWCBarrel = new G4PVPlacement(0,
            G4ThreeVector(0.,0.,WCIDOffset),
            logicWCBarrel,
            "WCBarrel",
            logicWCShaft,
            false,
            0);

    // if(!debugMode)
    //    logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible); 

    //-----------------------------------------------------
    // Form annular section of barrel to hold PMTs 
    //----------------------------------------------------


    G4double mainAnnulusZ[2] = {-mainAnnulusHeight/2., mainAnnulusHeight/2};
    G4double mainAnnulusRmin[2] = {innerAnnulusRadius, innerAnnulusRadius};
    G4double mainAnnulusRmax[2] = {outerAnnulusRadius, outerAnnulusRadius};

    G4Polyhedra* solidWCBarrelAnnulus = new G4Polyhedra("WCBarrelAnnulus",
            0.*deg, // phi start
            totalAngle, 
            (G4int)WCBarrelRingNPhi, //NPhi-gon
            2,
            mainAnnulusZ,
            mainAnnulusRmin,
            mainAnnulusRmax);

    G4LogicalVolume* logicWCBarrelAnnulus = 
        new G4LogicalVolume(solidWCBarrelAnnulus,
                G4Material::GetMaterial(water),
                "WCBarrelAnnulus",
                0,0,0);
    // G4cout << *solidWCBarrelAnnulus << G4endl; 
    G4VPhysicalVolume* physiWCBarrelAnnulus = 
        new G4PVPlacement(0,
                G4ThreeVector(0.,0.,0.),
                logicWCBarrelAnnulus,
                "WCBarrelAnnulus",
                logicWCBarrel,
                false,
                0,true);
    if(!debugMode)
        logicWCBarrelAnnulus->SetVisAttributes(G4VisAttributes::Invisible); //amb79
    //-----------------------------------------------------
    // Subdivide the BarrelAnnulus into rings
    //-----------------------------------------------------
    G4double RingZ[2] = {-barrelCellHeight/2.,
        barrelCellHeight/2.};

    G4Polyhedra* solidWCBarrelRing = new G4Polyhedra("WCBarrelRing",
            0.*deg,//+dPhi/2., // phi start
            totalAngle, //phi end
            (G4int)WCBarrelRingNPhi, //NPhi-gon
            2,
            RingZ,
            mainAnnulusRmin,
            mainAnnulusRmax);

    G4LogicalVolume* logicWCBarrelRing = 
        new G4LogicalVolume(solidWCBarrelRing,
                G4Material::GetMaterial(water),
                "WCBarrelRing",
                0,0,0);

    G4VPhysicalVolume* physiWCBarrelRing = 
        new G4PVReplica("WCBarrelRing",
                logicWCBarrelRing,
                logicWCBarrelAnnulus,
                kZAxis,
                (G4int)WCBarrelNRings-2,
                barrelCellHeight);

    if(!debugMode)
        logicWCBarrelRing->SetVisAttributes(G4VisAttributes::Invisible);
    else {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
        tmpVisAtt->SetForceWireframe(true);
        logicWCBarrelRing->SetVisAttributes(tmpVisAtt);
    }

    //-----------------------------------------------------
    // Subdivisions of the BarrelRings are cells
    //------------------------------------------------------


    G4Polyhedra* solidWCBarrelCell = new G4Polyhedra("WCBarrelCell",
            -dPhi/2.+0.*deg, // phi start
            dPhi, //total Phi
            1, //NPhi-gon
            2,
            RingZ,
            mainAnnulusRmin,
            mainAnnulusRmax); 
    //G4cout << *solidWCBarrelCell << G4endl; 
    G4LogicalVolume* logicWCBarrelCell = 
        new G4LogicalVolume(solidWCBarrelCell,
                G4Material::GetMaterial(water),
                "WCBarrelCell",
                0,0,0);

    G4VPhysicalVolume* physiWCBarrelCell = 
        new G4PVReplica("WCBarrelCell",
                logicWCBarrelCell,
                logicWCBarrelRing,
                kPhi,
                (G4int)WCBarrelRingNPhi,
                dPhi,
                0.); 

    if(!debugMode)
        logicWCBarrelCell->SetVisAttributes(G4VisAttributes::Invisible);
    else {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
        tmpVisAtt->SetForceWireframe(true);
        logicWCBarrelCell->SetVisAttributes(tmpVisAtt);
    }

    //-----------------------------------------------------------
    // The Blacksheet, a daughter of the cells containing PMTs,
    // and also some other volumes to make the edges light tight
    //-----------------------------------------------------------

    //-------------------------------------------------------------
    // add barrel blacksheet to the normal barrel cells 
    // ------------------------------------------------------------
    G4double annulusBlackSheetRmax[2] = {(WCIDRadius+WCBlackSheetThickness),
        WCIDRadius+WCBlackSheetThickness};
    G4double annulusBlackSheetRmin[2] = {(WCIDRadius),
        WCIDRadius};

    G4Polyhedra* solidWCBarrelCellBlackSheet = new G4Polyhedra("WCBarrelCellBlackSheet",
            -dPhi/2., // phi start
            dPhi, //total phi
            1, //NPhi-gon
            2,
            RingZ,
            annulusBlackSheetRmin,
            annulusBlackSheetRmax);

    logicWCBarrelCellBlackSheet =
        new G4LogicalVolume(solidWCBarrelCellBlackSheet,
                G4Material::GetMaterial("Blacksheet"),
                "WCBarrelCellBlackSheet",
                0,0,0);

    G4VPhysicalVolume* physiWCBarrelCellBlackSheet =
        new G4PVPlacement(0,
                G4ThreeVector(0.,0.,0.),
                logicWCBarrelCellBlackSheet,
                "WCBarrelCellBlackSheet",
                logicWCBarrelCell,
                false,
                0,true);

    G4LogicalBorderSurface * WaterBSBarrelCellSurface 
        = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
                physiWCBarrelCell,
                physiWCBarrelCellBlackSheet, 
                OpWaterBSSurface);


    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
        = new G4VisAttributes(G4Colour(0.2,0.9,0.2));
    if(debugMode)
        logicWCBarrelCellBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    else
        logicWCBarrelCellBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);


    //-----------------------------------------------------------
    // add extra tower if nessecary
    // ---------------------------------------------------------

    // we have to declare the logical Volumes 
    // outside of the if block to access it later on 
    G4LogicalVolume* logicWCExtraTowerCell;
    G4LogicalVolume* logicWCExtraBorderCell;
    if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

        // as the angles between the corners of the main annulus 
        // and the corners extra tower are different, we need to adjust the 
        // tangent distance the surfaces of the extra tower. Otherwise
        // the corners of the main annulus and the extra tower would 
        // not match. 
        G4double extraTowerRmin[2];
        G4double extraTowerRmax[2];
        for(int i = 0; i < 2 ; i++){
            extraTowerRmin[i] = mainAnnulusRmin[i] != 0 ? mainAnnulusRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
            extraTowerRmax[i] = mainAnnulusRmax[i] != 0 ? mainAnnulusRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
        }
        G4Polyhedra* solidWCExtraTower = new G4Polyhedra("WCextraTower",
                totalAngle-2.*pi,//+dPhi/2., // phi start
                2.*pi -  totalAngle // total angle.
                -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m),
                // we need this little Gap between the extra tower and the main annulus
                // to avoid a shared surface. Without the gap the photons stuck
                // at this place for mare than 25 steps and the howl simulation
                // crashes.
                1, //NPhi-gon
                2,
                mainAnnulusZ,
                extraTowerRmin,
                extraTowerRmax);

        G4LogicalVolume* logicWCExtraTower = 
            new G4LogicalVolume(solidWCExtraTower,
                    G4Material::GetMaterial(water),
                    "WCExtraTower",
                    0,0,0);
        G4VPhysicalVolume* physiWCExtraTower = 
            new G4PVPlacement(0,
                    G4ThreeVector(0.,0.,0.),
                    logicWCExtraTower,
                    "WCExtraTower",
                    logicWCBarrel,
                    false,
                    0,true);


        logicWCExtraTower->SetVisAttributes(G4VisAttributes::Invisible);
        //-------------------------------------------
        // subdivide the extra tower into cells  
        //------------------------------------------

        G4Polyhedra* solidWCExtraTowerCell = new G4Polyhedra("WCExtraTowerCell",
                totalAngle-2.*pi,//+dPhi/2., // phi start
                2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
                1, //NPhi-gon
                2,
                RingZ,
                extraTowerRmin,
                extraTowerRmax); 
        //G4cout << * solidWCExtraTowerCell << G4endl;
        logicWCExtraTowerCell = 
            new G4LogicalVolume(solidWCExtraTowerCell,
                    G4Material::GetMaterial(water),
                    "WCExtraTowerCell",
                    0,0,0);
        G4VPhysicalVolume* physiWCTowerCell = 
            new G4PVReplica("extraTowerCell",
                    logicWCExtraTowerCell,
                    logicWCExtraTower,
                    kZAxis,
                    (G4int)WCBarrelNRings-2,
                    barrelCellHeight);
        logicWCExtraTowerCell->SetVisAttributes(G4VisAttributes::Invisible);

        //---------------------------------------------
        // add blacksheet to this cells
        //--------------------------------------------

        G4double towerBSRmin[2];
        G4double towerBSRmax[2];
        for(int i = 0; i < 2; i++){
            towerBSRmin[i] = annulusBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
            towerBSRmax[i] = annulusBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
        }
        G4Polyhedra* solidWCTowerBlackSheet = new G4Polyhedra("WCExtraTowerBlackSheet",
                totalAngle-2.*pi,//+dPhi/2., // phi start
                2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
                1, //NPhi-gon
                2,
                RingZ,
                towerBSRmin,
                towerBSRmax);
        //G4cout << * solidWCTowerBlackSheet << G4endl;
        logicWCTowerBlackSheet =
            new G4LogicalVolume(solidWCTowerBlackSheet,
                    G4Material::GetMaterial("Blacksheet"),
                    "WCExtraTowerBlackSheet",
                    0,0,0);

        G4VPhysicalVolume* physiWCTowerBlackSheet =
            new G4PVPlacement(0,
                    G4ThreeVector(0.,0.,0.),
                    logicWCTowerBlackSheet,
                    "WCExtraTowerBlackSheet",
                    logicWCExtraTowerCell,
                    false,
                    0,true);

        G4LogicalBorderSurface * WaterBSTowerCellSurface 
            = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
                    physiWCTowerCell,
                    physiWCTowerBlackSheet, 
                    OpWaterBSSurface);


        if(debugMode)
            logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
        else
            logicWCTowerBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    }


    //jl145------------------------------------------------
    // Add top veto volume
    //-----------------------------------------------------

    G4bool WCTopVeto = (WCSimTuningParams->GetTopVeto());

    G4LogicalVolume* logicWCTopVeto;

    if(WCTopVeto){

        G4double WCTyvekThickness = 1.0*mm; //completely made up

        G4VSolid* solidWCTopVeto;
        solidWCTopVeto =
            new G4Tubs(			"WCTopVeto",
                    0.0*m,
                    WCIDRadius + WCTyvekThickness,
                    0.5*m + WCTyvekThickness,
                    0.*deg,
                    360.*deg);

        logicWCTopVeto = 
            new G4LogicalVolume(solidWCTopVeto,
                    G4Material::GetMaterial(water),
                    "WCTopVeto",
                    0,0,0);

        G4VPhysicalVolume* physiWCTopVeto =
            new G4PVPlacement(	0,
                    G4ThreeVector(0.,0.,WCIDHeight/2
                        +1.0*m),
                    logicWCTopVeto,
                    "WCTopVeto",
                    logicWCBarrel,
                    false,0,true);

        //Add the top veto Tyvek
        //-----------------------------------------------------

        G4VSolid* solidWCTVTyvek;
        solidWCTVTyvek =
            new G4Tubs(			"WCTVTyvek",
                    0.0*m,
                    WCIDRadius,
                    WCTyvekThickness/2,
                    0.*deg,
                    360.*deg);


        G4LogicalVolume* logicWCTVTyvek =
            new G4LogicalVolume(solidWCTVTyvek,
                    G4Material::GetMaterial("Tyvek"),
                    "WCTVTyvek",
                    0,0,0);

        //Bottom
        G4VPhysicalVolume* physiWCTVTyvekBot =
            new G4PVPlacement(	0,
                    G4ThreeVector(0.,0.,-0.5*m
                        -WCTyvekThickness/2),
                    logicWCTVTyvek,
                    "WCTVTyvekBot",
                    logicWCTopVeto,
                    false,0,true);

        G4LogicalBorderSurface * WaterTyTVSurfaceBot =
            new G4LogicalBorderSurface(	"WaterTyTVSurfaceBot",
                    physiWCTopVeto,
                    physiWCTVTyvekBot,
                    OpWaterTySurface);

        //Top
        G4VPhysicalVolume* physiWCTVTyvekTop =
            new G4PVPlacement(	0,
                    G4ThreeVector(0.,0.,0.5*m
                        +WCTyvekThickness/2),
                    logicWCTVTyvek,
                    "WCTVTyvekTop",
                    logicWCTopVeto,
                    false,0,true);

        G4LogicalBorderSurface * WaterTyTVSurfaceTop =
            new G4LogicalBorderSurface(	"WaterTyTVSurfaceTop",
                    physiWCTopVeto,
                    physiWCTVTyvekTop,
                    OpWaterTySurface);

        //Side
        G4VSolid* solidWCTVTyvekSide;
        solidWCTVTyvekSide =
            new G4Tubs(			"WCTVTyvekSide",
                    WCIDRadius,
                    WCIDRadius + WCTyvekThickness,
                    0.5*m + WCTyvekThickness,
                    0.*deg,
                    360.*deg);


        G4LogicalVolume* logicWCTVTyvekSide =
            new G4LogicalVolume(solidWCTVTyvekSide,
                    G4Material::GetMaterial("Tyvek"),
                    "WCTVTyvekSide",
                    0,0,0);

        G4VPhysicalVolume* physiWCTVTyvekSide =
            new G4PVPlacement(	0,
                    G4ThreeVector(0.,0.,0.),
                    logicWCTVTyvekSide,
                    "WCTVTyvekSide",
                    logicWCTopVeto,
                    false,0,true);

        G4LogicalBorderSurface * WaterTyTVSurfaceSide =
            new G4LogicalBorderSurface(	"WaterTyTVSurfaceSide",
                    physiWCTopVeto,
                    physiWCTVTyvekSide,
                    OpWaterTySurface);

    }

    //
    //
    //jl145------------------------------------------------


    //////////// M Fechner : I need to  declare the PMT  here in order to
    // place the PMTs in the truncated cells
    //-----------------------------------------------------
    // The PMT
    //-----------------------------------------------------

    ////////////J Felde: The PMT logical volume is now constructed separately 
    // from any specific detector geometry so that any geometry can use the same definition. 
    // K.Zbiri: The PMT volume and the PMT glass are now put in parallel. 
    // The PMT glass is the sensitive volume in this new configuration.

    G4LogicalVolume* logicWCPMT = ConstructPMT(WCPMTRadius, WCPMTExposeHeight);



    //jl145------------------------------------------------
    // Add top veto PMTs
    //-----------------------------------------------------

    if(WCTopVeto){

        G4double WCTVPMTSpacing = (WCSimTuningParams->GetTVSpacing())*cm;
        G4double WCTVEdgeLimit = WCCapEdgeLimit;
        G4int TVNCell = WCTVEdgeLimit/WCTVPMTSpacing + 2;

        int icopy = 0;

        for ( int i = -TVNCell ; i <  TVNCell; i++) {
            for (int j = -TVNCell ; j <  TVNCell; j++)   {

                G4double xoffset = i*WCTVPMTSpacing + WCTVPMTSpacing*0.5;
                G4double yoffset = j*WCTVPMTSpacing + WCTVPMTSpacing*0.5;

                G4ThreeVector cellpos =
                    G4ThreeVector(	xoffset, yoffset, -0.5*m);

                if ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCTVEdgeLimit) {

                    G4VPhysicalVolume* physiCapPMT =
                        new G4PVPlacement(	0,						// no rotation
                                cellpos,				// its position
                                logicWCPMT,				// its logical volume
                                "WCPMT",				// its name 
                                logicWCTopVeto,			// its mother volume
                                false,					// no boolean os
                                icopy);					// every PMT need a unique id.

                    icopy++;
                }
            }
        }

        G4double WCTVEfficiency = icopy*WCPMTRadius*WCPMTRadius/((WCIDRadius)*(WCIDRadius));
        G4cout << "Total on top veto: " << icopy << "\n";
        G4cout << "Coverage was calculated to be: " << WCTVEfficiency << "\n";

    }

    //
    //
    //jl145------------------------------------------------


    ///////////////   Barrel PMT placement
    G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
    WCPMTRotation->rotateY(90.*deg);

    G4double barrelCellWidth = 2.*WCIDRadius*tan(dPhi/2.);
    G4double horizontalSpacing   = barrelCellWidth/WCPMTperCellHorizontal;
    G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

    for(G4double i = 0; i < WCPMTperCellHorizontal; i++){
        for(G4double j = 0; j < WCPMTperCellVertical; j++){
            G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius,
                    -barrelCellWidth/2.+(i+0.5)*horizontalSpacing,
                    -barrelCellHeight/2.+(j+0.5)*verticalSpacing);


            G4VPhysicalVolume* physiWCBarrelPMT =
                new G4PVPlacement(WCPMTRotation,              // its rotation
                        PMTPosition, 
                        logicWCPMT,                // its logical volume
                        "WCPMT",             // its name
                        logicWCBarrelCell,         // its mother volume
                        false,                     // no boolean operations
                        (int)(i*WCPMTperCellVertical+j),
                        true);                       

            // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
            // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
            // this is still the case.
        }
    }
    //-------------------------------------------------------------
    // Add PMTs in extra Tower if necessary
    //------------------------------------------------------------
    if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

        G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
        WCPMTRotation->rotateY(90.*deg);
        WCPMTRotation->rotateX((2*pi-totalAngle)/2.);//align the PMT with the Cell

        G4double towerWidth = WCIDRadius*tan(2*pi-totalAngle);

        G4double horizontalSpacing   = towerWidth/(WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal);
        G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

        for(G4double i = 0; i < (WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal); i++){
            for(G4double j = 0; j < WCPMTperCellVertical; j++){
                G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.),
                        towerWidth/2.-(i+0.5)*horizontalSpacing,
                        -barrelCellHeight/2.+(j+0.5)*verticalSpacing);
                PMTPosition.rotateY(-(2*pi-totalAngle)/2.); // align with the symmetry 
                //axes of the cell 

                G4VPhysicalVolume* physiWCBarrelPMT =
                    new G4PVPlacement(WCPMTRotation,             // its rotation
                            PMTPosition, 
                            logicWCPMT,                // its logical volume
                            "WCPMT",             // its name
                            logicWCExtraTowerCell,         // its mother volume
                            false,                     // no boolean operations
                            (int)(i*WCPMTperCellVertical+j),
                            true);                       

                // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
                // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
                // this is still the case.
            }
        }

    }


    G4LogicalVolume* logicTopCapAssembly = ConstructCaps(-1);
    G4LogicalVolume* logicBottomCapAssembly = ConstructCaps(1);

    G4VPhysicalVolume* physiTopCapAssembly =
        new G4PVPlacement(0,
                G4ThreeVector(0.,0.,(mainAnnulusHeight/2.+ capAssemblyHeight/2.)),
                logicTopCapAssembly,
                "TopCapAssembly",
                logicWCBarrel,
                false, 0,true);

    G4VPhysicalVolume* physiBottomCapAssembly =
        new G4PVPlacement(0,
                G4ThreeVector(0.,0.,(-mainAnnulusHeight/2.- capAssemblyHeight/2.)),
                logicBottomCapAssembly,
                "BottomCapAssembly",
                logicWCBarrel,
                false, 0,true);

    return logicWC;
}

