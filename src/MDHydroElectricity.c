/******************************************************************************
NIMA
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2013, UNH - CCNY/CUNY

MDHydroElectricity.c

nehsani00@ccny.cuny.edu

 *******************************************************************************/


#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input
static int _MDInResCapacityID  = MFUnset;
static int _MDInResStorageID   = MFUnset;
static int _MDInResMaxHeightID = MFUnset;
static int _MDInMaxHydroCapID  = MFUnset;
static int _MDInResReleaseID   = MFUnset;
static int _MDInAvgNStepsID    = MFUnset;
// Output
static int _MDOutMegaWattID = MFUnset;

static void _MDHydroPower (int itemID) {
    float resstorage; // Reservoir storage [m3]
    float resrelease; // Reservoir release [m3/s]
    float resCapacity; // Reservoir capacity [m3]
    float resmaxH;
    float resH;
    float maxhydropcap;
    float hydrogen;
//    int nSteps;
    
//    if ((nSteps     = MFVarGetInt   (_MDInAvgNStepsID,       itemID,   0)) <= 2) {
//        return;
//    }

    if ((maxhydropcap = MFVarGetFloat(_MDInMaxHydroCapID,  itemID, 0.0)) > 0.0) {
        resCapacity   = MFVarGetFloat(_MDInResCapacityID,  itemID, 0.0);
        resrelease    = MFVarGetFloat(_MDInResReleaseID,   itemID, 0.0);
        resstorage    = MFVarGetFloat(_MDInResStorageID,   itemID, 0.0);
        resmaxH       = 1.25*MFVarGetFloat(_MDInResMaxHeightID, itemID, 0.0); // Consider 25% overhead

//////////////////////////////////////////////////////////////////////
//////////Reservoir as a paraboloid////////////////////////////    
///////////////////////////////////////////////////////////////////
  /*     
        float a;
        a = 2 * resCapacity / (3.14 * (pow(resmaxH,2)));
        resH = sqrt(2 * resstorage / (3.14 * a));
  */      
/////////////////////////////////////////////
//////Reservoir As a Cone///////////////////
////////////////////////////////////////////  
       
        float R; 
        float tga;
        resstorage=resstorage+0.3*resCapacity; //Adding back dead capacity for head calculation
        R=sqrt(3*resCapacity/(3.14*resmaxH));
        tga=R/resmaxH;
        resH=pow((3*resstorage/(3.14*pow(tga,2))),1/3);
///////////////////////////////////////////////
//////////////////////////////////////////////////
/*        float delta;
        delta=(resCapacity-resstorage)/resCapacity;
        resH=0.2*delta*resmaxH;
      
*///////////////////////////////////////////////
////////////////////////////////////////////////
            
            hydrogen =  9810 * resH * resrelease / 1000000; // Power Generation in MEga Watt
        if (hydrogen > maxhydropcap) {
            hydrogen = maxhydropcap;
        }
        MFVarSetFloat(_MDOutMegaWattID, itemID, hydrogen);
    }
    else {
        MFVarSetFloat(_MDOutMegaWattID, itemID, 0.0);
        return;
    }
}

enum { MDnone, MDcalculate };

int MDHydroPowerDef() {
	int  optID = MFUnset;
	const char *optStr, *optName = MDOptHydroPower;
	const char *options [] = { MDNoneStr, MDCalculateStr, (char *) NULL };

	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);
        if ((optID == MDnone) || (_MDOutMegaWattID != MFUnset)) return (_MDOutMegaWattID);
        
    MFDefEntering("HydroElectricity");
    
    switch (optID) {
        case MDcalculate:
            if (    //(                         MDWaterBalanceDef() == CMfailed) ||
                    ((_MDInAvgNStepsID      = MDAvgNStepsDef() )  == CMfailed) ||
                    ((_MDInResStorageID     = MDReservoirDef() )  == CMfailed) ||
                    ((_MDInResReleaseID     = MDReservoirDef() )  == CMfailed) ||
                    ((_MDInMaxHydroCapID    = MFVarGetID(MDVarMaxHydroCap,                  "MW",   MFInput,  MFState,  MFBoundary)) == CMfailed) ||
                    ((_MDInResMaxHeightID   = MFVarGetID(MDVarResMaxHeight,                 "m",    MFInput,  MFState,  MFBoundary)) == CMfailed) ||
                    ((_MDInResCapacityID    = MFVarGetID(MDVarReservoirCapacity,            "m3",   MFInput,  MFState,  MFBoundary)) == CMfailed) ||
                    ((_MDOutMegaWattID      = MFVarGetID(MDVarMegaWatt,                     "MW",   MFOutput, MFState,  MFInitial))  == CMfailed) ||
                    ((MFModelAddFunction(_MDHydroPower) == CMfailed))
                    ) return (CMfailed);
			break;
		default: MFOptionMessage (optName, optStr, options); return (CMfailed);
    }
    MFDefLeaving("HydroElectricity");
    return (_MDOutMegaWattID);
}

