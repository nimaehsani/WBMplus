/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDRunoff.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>


// Input
static int _MDInRunoffPoolReleaseID = MFUnset;	// RJS 042612
static int _MDInStormRunoffTotalID  = MFUnset;	// RJS 082812
static int _MDInBaseFlowID  	    = MFUnset;
static int _MDInRunoffCorrID	    = MFUnset;
static int _MDInTotalSurfRunoffID   = MFUnset;	// RJS 082812
// Output
static int _MDOutRunoffID   	      = MFUnset;
static int _MDOutTotalSurfRunoffID    = MFUnset;  // RJS 082812
static int _MDOutPropROStormWaterID   = MFUnset;  // RJS 100313
static int _MDOutPropROSurfaceWaterID = MFUnset;  // RJS 100313
static int _MDOutPropROGroundWaterID  = MFUnset;  // RJS 100313

static void _MDRunoff (int itemID) {
// Input

	float baseFlow;
	float stormRunoffTotal;		// RJS 082812
	float runoffPoolRelease;	// RJS 082812
	float surfaceRO;
	float runoffCorr;
        float propStW;                  // RJS 100313
        float propSuW;                  // RJS 100313
        float propGrW;                  // RJS 100313
        float totalRO;                  // RJS 100313

	baseFlow 	   = MFVarGetFloat (_MDInBaseFlowID,  itemID, 0.0);
	runoffPoolRelease  = MFVarGetFloat (_MDInRunoffPoolReleaseID, itemID, 0.0);		// RJS 042712
	stormRunoffTotal   = MFVarGetFloat (_MDInStormRunoffTotalID,  itemID, 0.0);		// RJS 082812
	surfaceRO	   = runoffPoolRelease + stormRunoffTotal;							// RJS 082812
	runoffCorr	   = _MDInRunoffCorrID == MFUnset ? 1.0 : MFVarGetFloat (_MDInRunoffCorrID, itemID, 1.0);

        totalRO = baseFlow + runoffPoolRelease + stormRunoffTotal;
        propStW = totalRO > 0.0 ? stormRunoffTotal / totalRO : 0.33333;
        propSuW = totalRO > 0.0 ? runoffPoolRelease / totalRO : 0.33333;
        propGrW = totalRO > 0.0 ? baseFlow / totalRO : 0.33333;
        
//	if ((itemID == 293)) {
//			printf("m = %d, d = %d --yes-- id = %d, area = %f,  baseFlow = %f, runoffPoolRelease = %f, stormRunoff = %f, surfaceRO = %f\n", MFDateGetCurrentMonth(), MFDateGetCurrentDay(), itemID, MFModelGetArea (itemID), baseFlow * MFModelGetArea (itemID) / (1000 * 86400), runoffPoolRelease * MFModelGetArea (itemID) / (1000 * 86400), stormRunoffTotal * MFModelGetArea (itemID) / (1000 * 86400),surfaceRO * MFModelGetArea (itemID) / (1000 * 86400));
//			printf("m = %d, d = %d --yes-- id = %d, area = %f,  baseFlow = %f, runoffPoolRelease = %f, stormRunoff = %f, surfaceRO = %f\n", MFDateGetCurrentMonth(), MFDateGetCurrentDay(), itemID, MFModelGetArea (itemID), baseFlow, runoffPoolRelease, stormRunoffTotal,surfaceRO);
//	}

	MFVarSetFloat (_MDOutRunoffID,          itemID, (baseFlow + surfaceRO) * runoffCorr);
	MFVarSetFloat (_MDOutTotalSurfRunoffID, itemID, surfaceRO);                                     // RJS 082812
        MFVarSetFloat (_MDOutPropROStormWaterID,  itemID, propStW);                                       // RJS 100313
        MFVarSetFloat (_MDOutPropROSurfaceWaterID,   itemID, propSuW);                                       // RJS 100313
        MFVarSetFloat (_MDOutPropROGroundWaterID, itemID, propGrW);                                       // RJS 100313
}


static void _MDRunoffInput (int itemID) {														// RJS 061312  ADDED THIS WHOLE FUNCTION
// Input
	float baseFlow;																				// "
	float surfaceRO;
	float prop = 0.33333;																			// "

	baseFlow  = MFVarGetFloat (_MDInBaseFlowID,          itemID, 0.0);							// "
//	surfaceRO = MFVarGetFloat (_MDInRunoffPoolReleaseID, itemID, 0.0);							// "
	surfaceRO = MFVarGetFloat (_MDInTotalSurfRunoffID, itemID, 0.0);							// RJS 082812, replaces line above
	

	MFVarSetFloat (_MDOutRunoffID, itemID, (baseFlow + surfaceRO));								// "
	MFVarSetFloat (_MDOutPropROStormWaterID, itemID, prop);
	MFVarSetFloat (_MDOutPropROSurfaceWaterID, itemID, prop);
	MFVarSetFloat (_MDOutPropROGroundWaterID, itemID, prop);
}																								// "

 
enum { MDinput, MDcalculate, MDcorrected };

int MDRunoffDef () {
	int  optID = MFUnset;
	const char *optStr, *optName = MDVarRunoff;
	const char *options [] = { MDInputStr, MDCalculateStr, "corrected", (char *) NULL };

	if (_MDOutRunoffID != MFUnset) return (_MDOutRunoffID);

	MFDefEntering ("Runoff");
	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);
	switch (optID) {
//		case MDinput: _MDOutRunoffID = MFVarGetID (MDVarRunoff,         "mm",     MFInput,  MFFlux, MFBoundary); break;						// RJS commented out 061312
		case MDinput: 																										       			// RJS 061312
			if (((_MDInBaseFlowID          = MFVarGetID (MDVarBaseFlow,          "mm",    MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	// RJS 061312
//				((_MDInRunoffPoolReleaseID = MFVarGetID (MDVarRunoffPoolRelease, "mm",    MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	// RJS 061312, commented out 082812
				((_MDInTotalSurfRunoffID   = MFVarGetID (MDVarTotalSurfRunoff,   "mm",    MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	// RJS 082812		
                                ((_MDOutPropROStormWaterID   = MFVarGetID (MDVarPropROStormWater,    "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
                                ((_MDOutPropROSurfaceWaterID = MFVarGetID (MDVarPropROSurfaceWater,  "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313 
                                ((_MDOutPropROGroundWaterID  = MFVarGetID (MDVarPropROGroundWater,   "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||  
				
				((_MDOutRunoffID           = MFVarGetID (MDVarRunoff,            "mm",    MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	// RJS 061312
				(MFModelAddFunction (_MDRunoffInput) == CMfailed)) return (CMfailed);														// RJS 061312
			break;																															// RJS 061312
		case MDcorrected:
			if ((_MDInRunoffCorrID  = MFVarGetID (MDVarRunoffCorretion, MFNoUnit, MFInput,  MFState, MFBoundary)) == CMfailed)
				return (CMfailed);
			break;	// RJS 082812
		case MDcalculate:		
			if (((_MDInBaseFlowID   = MDBaseFlowDef   ()) == CMfailed) ||
//			    ((_MDInSurfRunoffID = MDSurfRunoffDef ()) == CMfailed) ||				//commented out RJS 042612
				((_MDInRunoffPoolReleaseID = MDSurfRunoffPoolDef ()) == CMfailed) ||		// RJS 042612
//				((_MDInStormRunoffTotalID  = MDStormRunoffDef ()) == CMfailed) ||		// RJS 082812
				((_MDInStormRunoffTotalID  = MFVarGetID (MDVarStormRunoffTotal, "mm",    MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
                                ((_MDOutPropROStormWaterID   = MFVarGetID (MDVarPropROStormWater,    "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
                                ((_MDOutPropROSurfaceWaterID = MFVarGetID (MDVarPropROSurfaceWater,  "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313 
                                ((_MDOutPropROGroundWaterID  = MFVarGetID (MDVarPropROGroundWater,   "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||  
				((_MDOutRunoffID           = MFVarGetID (MDVarRunoff,          "mm",     MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
				((_MDOutTotalSurfRunoffID  = MFVarGetID (MDVarTotalSurfRunoff, "mm",     MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
				(MFModelAddFunction (_MDRunoff) == CMfailed)) return (CMfailed);
			break;
		default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}
	MFDefLeaving  ("Runoff");
	return (_MDOutRunoffID);
}
