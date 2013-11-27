/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

WBMMain.c

balazs.fekete@unh.edu

*******************************************************************************/
#include <wbm.h>

enum { MDpet, MDsurplus, MDinfiltration, MDrunoff, MDdischarge,  MDbalance, MDwatertemp, MDthermal, MDgeometry, MDbgc, MDbgc_DOC, MDbgc_DIN, MDbgc_DINPLUSBIOMASS, MDHydroElectricity/*, MDfecal*/};

int main (int argc,char *argv []) {
	int argNum;
	int  optID = MDbalance;
	const char *optStr, *optName = MDOptModel;
	const char *options [] = { "pet", "surplus", "infiltration", "runoff", "discharge",  "balance", "watertemp", "thermal",  "geometry", "bgc", "bgc_DOC", "bgc_DIN", "bgc_DINPLUSBIOMASS", "HydroElectricity", /*"fecal",*/ (char *) NULL };

	argNum = MFOptionParse (argc,argv);

	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);

	switch (optID) {
		case MDpet:          return (MFModelRun (argc,argv,argNum,MDRainPotETDef));
		case MDsurplus:      return (MFModelRun (argc,argv,argNum,MDRainWaterSurplusDef));
		case MDinfiltration: return (MFModelRun (argc,argv,argNum,MDRainInfiltrationDef));
		case MDrunoff:       return (MFModelRun (argc,argv,argNum,MDRunoffDef));
		case MDdischarge:    return (MFModelRun (argc,argv,argNum,MDDischargeDef));
		case MDbalance:      return (MFModelRun (argc,argv,argNum,MDWaterBalanceDef));
		case MDwatertemp:    return (MFModelRun (argc,argv,argNum,MDWTempRiverRouteDef));
		case MDthermal:	     return (MFModelRun (argc,argv,argNum,MDThermalInputsDef));		// RJS 013112
		case MDgeometry:     return (MFModelRun (argc,argv,argNum,MDRiverWidthDef));
		case MDbgc:          return (MFModelRun (argc,argv,argNum,MDBgcRoutingDef));
		case MDbgc_DOC:      return (MFModelRun (argc,argv,argNum,MDBgcDOCRoutingDef));
		case MDbgc_DIN:      return (MFModelRun (argc,argv,argNum,MDBgcDINRoutingDef));
		case MDbgc_DINPLUSBIOMASS:    return (MFModelRun (argc,argv,argNum,MDBgcDINPlusBiomassRoutingDef));
                case MDHydroElectricity:      return (MFModelRun (argc,argv,argNum,MDHydroPowerDef));
		default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}
	return (CMfailed);
}
