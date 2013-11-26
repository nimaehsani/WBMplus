/******************************************************************************
NIMA
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2013, UNH - CCNY/CUNY

MDHydroPower.c

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
static int _MDInMAxHydroCapID  = MFUnset;
static int _MDInResReleaseID   = MFUnset;

// Output

static int _MDOutHydroPowerGenerationID = MFUnset;

static void _MDHydroPowerGen(int itemID) {


    float resstorage; // Reservoir storage [m3]
    float resrelease; // Reservoir release [m3/s]
    float resCapacity; // Reservoir capacity [m3]
    float resmaxH;
    float resH;
    float maxhydropcap;
    float hydrogen;
    float a; // Y = a X^2

    if ((maxhydropcap = MFVarGetFloat(_MDInMAxHydroCapID,  itemID, 0.0)) <= 0.0) {
        resCapacity   = MFVarGetFloat(_MDInResCapacityID,  itemID, 0.0);
        resrelease    = MFVarGetFloat(_MDInResReleaseID,   itemID, 0.0);
        resstorage    = MFVarGetFloat(_MDInResStorageID,   itemID, 0.0);
        resmaxH       = MFVarGetFloat(_MDInResMaxHeightID, itemID, 0.0);

        a = sqrt(2 * resCapacity / (3.14 * (pow(resmaxH,2))));
        resH = sqrt(2 * resstorage / (3.14 * (pow(a,2))));
        hydrogen = 0.9 * 9810 * resH * resrelease / 1000000; // Power Generation in MEga Watt
        if (hydrogen > maxhydropcap) {
            hydrogen = maxhydropcap;
        }
        MFVarSetFloat(_MDOutHydroPowerGenerationID, itemID, hydrogen);

    }
}


enum { MDnone, MDinput, MDcalculate };


int MDHydroPowerGenDef() {
    int optID = MFUnset;
    const char *optStr, *optName = MDOptHydroPower;
    const char *options [] = {MDNoneStr, MDCalculateStr, (char *) NULL};

    if ((optStr = MFOptionGet(optName)) != (char *) NULL) optID = CMoptLookup(options, optStr, true);

    if ((optID == MDnone) || (_MDOutHydroPowerGenerationID != MFUnset)) return (_MDOutHydroPowerGenerationID);

    MFDefEntering("HydroPowerGen");
    switch (optID) {
        case MDcalculate:
            if (    ((MFModelAddFunction(_MDHydroPowerGen) == CMfailed)) || 
                    ((_MDInMAxHydroCapID    = MFVarGetID(MDVarMAxHydroCap,                  "MW",   MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                    ((_MDInResMaxHeightID   = MFVarGetID(MDVarResMaxHeight,                 "m",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInResCapacityID    = MFVarGetID(MDVarReservoirCapacity,            "m3",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInResStorageID     = MFVarGetID(MDVarReservoirStorage,             "m3",   MFInput, MFState, MFInitial))  == CMfailed) ||
                    ((_MDInResReleaseID     = MFVarGetID(MDVarReservoirRelease,             "m3/s", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ((_MDOutHydroPowerGenerationID     = MFVarGetID(MDVarHydroPowerGeneration,             "MW", MFOutput, MFFlux,  MFInitial)) == CMfailed) 
                    ) return (CMfailed);
            break;
        default: MFOptionMessage(optName, optStr, options);
            return (CMfailed);
    }
    MFDefLeaving("HydroPowerGen");
    return (_MDOutHydroPowerGenerationID);
}

