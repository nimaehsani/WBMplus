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

static void _MDHydroGeneration(int itemID) {


    float resstorage; // Reservoir storage [m3]
    float resrelease; // Reservoir release [m3/s]
    float resCapacity; // Reservoir capacity [m3]
    float resmaxH;
    float resH;
    float maxhydropcap;
    float hydrogen;
    float a; // Y = a X^2

    if ((maxhydropcap = MFVarGetFloat(_MDInMAxHydroCapID,  itemID, 0.0)) <= 0.0) {
        resCapacity   = MFVarGetFloat(_MDInResCapacityID,  itemID, 0.0)
        resrelease    = MFVarGetFloat(_MDInResReleaseID,   itemID, 0.0);
        resstorage    = MFVarGetFloat(_MDInResStorageID,   itemID, 0.0);
        resmaxH       = MFVarGetFloat(_MDInResMaxHeightID, itemID, 0.0);

        a = (2 * resCapacity / (3.14 * resmaxH^2))^0.5;
        resH = (2 * resstorage / (3.14 * a^2))^0.5;
        hydrogen = 0.9 * 9810 * resH * resrelease / 1000000; // Power Generation in MEga Watt
        if (hydrogen > maxhydropcap) {
            hydrogen = maxhydropcap;
        }
        MFVarSetFloat(_MDOutHydroPowerGenerationID, itemID, hydrogen);

    }
}


enum {
    MDnone, MDcalculate
};


int MDHydroPowerDef() {
    int optID = MFUnset;
    const char *optStr, *optName = MDOptHydroPower;
    const char *options [] = {MDNoneStr, MDCalculateStr, (char *) NULL};

    if ((optStr = MFOptionGet(optName)) != (char *) NULL) optID = CMoptLookup(options, optStr, true);

    if ((optID == MDnone) || (_MDOutHydroPowerGenerationID != MFUnset)) return (_MDOutHydroPowerGenerationID);

    MFDefEntering("HydroPower");
    switch (optID) {
        case MDcalculate:
            if (    ((MFModelAddFunction(_MDHydroGeneration) == CMfailed)) || 
                    ((_MDInMAxHydroCapID    = MFVarGetID(MDVarMAxHydroCap,                  "MW",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInResMaxHeightID   = MFVarGetID(MDVarResMaxHeight,                 "m",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInResCapacityID    = MFVarGetID(MDVarReservoirCapacity,            "m3",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInResStorageID     = MFVarGetID(MDVarReservoirStorage,             "m3",   MFOutput, MFState, MFInitial))  == CMfailed) ||
                    ((_MDInResReleaseID     = MFVarGetID(MDVarReservoirRelease,             "m3/s", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ) return (CMfailed);
            break;
        default: MFOptionMessage(optName, optStr, options);
            return (CMfailed);
    }
    MFDefLeaving("HydroPower");
    return (_MDOutHydroPowerGenerationID);
}

