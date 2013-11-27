/******************************************************************************
NIMA&BALAZS
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDReservoirs.c

nehsani00@ccny.cuny.edu
dominik.wisser@unh.edu

Updated with a Neural Network function for 
Reservoir Operation.
2013 NE

 *******************************************************************************/


#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input
static int _MDInDischargeID          = MFUnset;
static int _MDInDischMeanID          = MFUnset;
static int _MDInResCapacityID        = MFUnset;
static int _MDInMegaWattID = MFUnset;

// Output
static int _MDOutResStorageID        = MFUnset;
static int _MDOutResStorageChgID     = MFUnset;
static int _MDOutResReleaseID        = MFUnset;
static int _MDOutResRelease_t_1_ID   = MFUnset;
static int _MDOutResRelease_t_2_ID   = MFUnset;
static int _MDOutResRelease_t_3_ID   = MFUnset;
static int _MDOutDisch_t_1_ID        = MFUnset;
static int _MDOutDisch_t_2_ID        = MFUnset;
static int _MDOutDisch_t_3_ID        = MFUnset;
static int _MDOutDischMinID          = MFUnset; // -MDIn or MDOut? Value Updates in Each State and should be Used in Next Step
static int _MDOutDischMaxID          = MFUnset; // ...
static int _MDOutReleaseMinID        = MFUnset; // ...
static int _MDOutReleaseMaxID        = MFUnset; // ...
static int _MDOutLastMonthID         = MFUnset;
static int _MDOutMonthToDayInFlowID  = MFUnset;
static int _MDOutMonthToDayReleaseID = MFUnset;
//static int _MDOutReservoirReleaseID  = MFUnset;
static int _MDInAvgNStepsID          = MFUnset;

static float ANNOUTPUT(float I1[3][1], float I2[2][1], float I3) {

    float FirstLayerBias[6][1] = {
        {-1.8349638599010885},
        {0.31681519142609904},
        {0.90820782843797254},
        {0.24315763055718279},
        {0.22750233763000977},
        {-5.9638615037974105},
    };

    float FirstLayerWeight [6][3] = {
        {0.334525270188955, 2.32729460498139, -1.08218040335814},
        {-0.132942683641029, -2.07435066955086, -7.71723103154886},
        {-0.961228914800186, -1.59451723505862, 7.64275651580055},
        {-0.138149619471609, -0.369268995299719, -1.05878021476956},
        {-0.830754365521219, -1.14297810944283, -0.591263860604682},
        {-2.35234764862494, -0.170592638429000, 6.45623554696750},
    };

    float SecondLayerWeight[4][2] = {
        {-0.461376663064363, -2.29680500565283},
        {1.14189579706809, 3.60122216280990},
        {5.13965726884100, -8.54024357801218},
        {0.797602695519510, 9.04132111622455},
    };

    float SecondLayerBias [4][1] = {
        {1.7667185407719592},
        {-0.44451918724448725},
        {0.27782065349273022},
        {0.95132682701557358},
    };

    float ThirdLayerWeight[2][1] = {
        {-3.80846357210734},
        {6.09704986254896},
    };

    float ThirdLayerBias [2][1] = {
        {1.9804896408537742},
        {-1.4736095994793379},
    };

    float FourthLayerWeight1[6][6] = {
        {-0.682425987793126, 2.20351531749950, 1.39097166906628, -2.39124203637190, 2.27986610311803, -0.633039777101503},
        {-0.478344057060439, 2.22542158997640, 1.32425417280671, -2.91086448864336, 2.55366825518159, 0.915264601641264},
        {-0.330288477684471, 0.636650781140816, -0.267485765559433, -0.0845957816699171, -0.0311750961112513, 0.519857878801029},
        {0.871647922418761, -0.0363653064870689, -1.23489938556299, 3.48735527686833, -4.09204564273278, -1.29124333772506},
        {-1.14578627733860, 0.362029723412182, -1.49653052173892, -2.14316546336101, 2.53719571843746, 1.03545166271711},
        {-3.93424835569323, -0.140516448398705, -0.287183959118834, 2.81056586729428, -1.13796052928068, -1.61032679050170},
    };

    float FourthLayerWeight2[6][4] = {
        {1.73940117579596, 4.31257163704923, 2.00908995797761, 0.542330032098119},
        {1.82716477758809, 3.90502202733331, 1.71260564791546, 0.673923447124240},
        {-0.297770956962901, 0.353599592055875, -0.0335628624334481, 1.90359407919045},
        {0.211761429085206, -1.37612121135069, -0.303662471752540, 1.13329378872693},
        {-0.299078423419769, 0.649961868699184, 0.173313681427546, 0.573538752543163},
        {-3.38710042074078, 1.16448748228937, -0.0607349889977066, -2.01182147710900},
    };

    float FourthLayerWeight3[6][2] = {
        {-1.84202629698469, -0.877541434013040},
        {-1.96090394140442, -0.946074755836152},
        {0.0173356328650234, -0.00969615110828418},
        {-0.0839568134797032, -0.129727081135696},
        {0.0356547772242401, 0.0691169740242024},
        {0.330775400271462, 0.0723975371228814},
    };

    float FourthLayerBias[6][1] = {
        {0.0051805882973190165},
        {1.7022493850334097},
        {-1.493864901771131},
        {-0.49352259943031906},
        {0.78319807803209185},
        {-0.86188699136698088},
    };

    float FifthLayerWeight [1][6] = {
        {-4.60687439047288, 4.68810658158703, 7.48854158670334, -2.67349645890795, -5.84593641782589, -1.13313370598710},
    };

    float FifthLayerBias = 0.22258525600927653;



    // LAYER ONE 

    float FirstLayerOut [6][1];
    int i;
    float n;

    for (i = 0; i < 6; i++) {
        n = I1[0][0] * FirstLayerWeight[i][0] + I1[1][0] * FirstLayerWeight[i][1] + I1[2][0] * FirstLayerWeight[i][2] + FirstLayerBias [i][0];
        FirstLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }

    // LAYER TWO 
    float SecondLayerOut[4][1];

    for (i = 0; i < 4; i++) {
        n = I2[0][0] * SecondLayerWeight[i][0] + I2[1][0] * SecondLayerWeight[i][1] + SecondLayerBias [i][0];
        SecondLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }

    // LAYER THREE 
    float ThirdLayerOut[2][1];

    for (i = 0; i < 2; i++) {
        n = I3 * ThirdLayerWeight[i][0] + ThirdLayerBias [i][0];
        ThirdLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }

    // LAYER FOUR //

    float FourthLayerOut[6][1];

    for (i = 0; i < 6; i++) {
        n = FirstLayerOut[0][0] * FourthLayerWeight1[i][0] + FirstLayerOut[1][0] * FourthLayerWeight1[i][1] + FirstLayerOut[2][0] * FourthLayerWeight1[i][2] + FirstLayerOut[3][0] * FourthLayerWeight1[i][3] + FirstLayerOut[4][0] * FourthLayerWeight1[i][4] + FirstLayerOut[5][0] * FourthLayerWeight1[i][5]
                + SecondLayerOut[0][0] * FourthLayerWeight2[i][0] + SecondLayerOut[1][0] * FourthLayerWeight2[i][1] + SecondLayerOut[2][0] * FourthLayerWeight2[i][2] + SecondLayerOut[3][0] * FourthLayerWeight2[i][3]
                + ThirdLayerOut[0][0] * FourthLayerWeight3[i][0] + ThirdLayerOut[1][0] * FourthLayerWeight3[i][1]
                + FourthLayerBias [i][0];
        FourthLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }

    // LAYER FIVE 

    float FifthLayerOut;

    n = FourthLayerOut[0][0] * FifthLayerWeight[0][0] + FourthLayerOut[1][0] * FifthLayerWeight[0][1] + FourthLayerOut[2][0] * FifthLayerWeight[0][2] + FourthLayerOut[3][0] * FifthLayerWeight[0][3] + FourthLayerOut[4][0] * FifthLayerWeight[0][4] + FourthLayerOut[5][0] * FifthLayerWeight[0][5]
            + FifthLayerBias;
    FifthLayerOut = 1 / (1 + exp(-n));

    float ANNOUTPUT;

    ANNOUTPUT = FifthLayerOut;

    return (ANNOUTPUT);

}
 

static void _MDReservoirNeuralNet(int itemID) {
    
    int   m             = MFDateGetCurrentMonth();
    int   d             = MFDateGetCurrentDay();
    int   lastmonth;
    float discharge;      // Current discharge [m3/s]
    float resCapacity;    // Reservoir capacity [m3]
    float minresStorage;
    float discharge_t_1;
    float discharge_t_2;
    float discharge_t_3;
    float discharge_min;
    float discharge_max;
    float release_min;
    float release_max;
    float mtdInflow;       // Total Monthly (Up to that Date) InFlow
    float avmtdInflow;     // Average Monthly (Up to that Date) InFlow
    float mtdRelease;      // Total Monthly (Up to that Date) Release
    float avmtdRelease;    // Average Monthly (Up to that Date) Release
    float I1[3][1];        // Input to ANN (ANNOUTPUT.c) 
    float I2[2][1];        // Input to ANN (ANNOUTPUT.c) 
    float I3;              // Input to ANN (ANNOUTPUT.c) 
    float ANN;
    float resStorage;      // Reservoir storage [m3]
    float resStorageChg;   // Reservoir storage change [m3/dt]
    float resRelease;      // Reservoir release [m3/s]
    float res_release_t_1; // STANDARD Reservoir release (resRelease) [m3/s]
    float res_release_t_2;
    float res_release_t_3;
    float SIMOUT;
    float prevResStorage;
    float SD_t_3;
    float SD_t_2;
    float SD_t_1;
    float SR_t_3;
    float SR_t_2;
    int   nSteps;

    discharge = MFVarGetFloat(_MDInDischargeID, itemID, 0.0);
   
    if ((resCapacity = MFVarGetFloat(_MDInResCapacityID,    itemID, 0.0)) <= 0.0) {
                       MFVarSetFloat(_MDOutResStorageID,    itemID, 0.0);
                       MFVarSetFloat(_MDOutResStorageChgID, itemID, 0.0);
                       MFVarSetFloat(_MDOutResReleaseID,    itemID, discharge);
                       //printf("discharge NO Res= %f \n" , discharge);

        //		if (itemID == 25014) printf("@@@ m= %d, d= %d, balance = %f, resCapacity = %f, Q = %f, meanQ = %f, resRelease = %f, resStorage = %f, prevResStorage = %f\n", MFDateGetCurrentMonth(), MFDateGetCurrentDay(), balance, resCapacity, discharge, meanDischarge, resRelease, resStorage*1000000000, prevResStorage*1000000000);
        return;
    } //else {
        resCapacity   = MFVarGetFloat(_MDInResCapacityID, itemID, 0.0);
        discharge_max = MFVarGetFloat(_MDOutDischMaxID,   itemID, 0.0);
        discharge_min = MFVarGetFloat(_MDOutDischMinID,   itemID, 0.0);
        //if (/*discharge > 0 &&*/ discharge < discharge_min) { // Chekking MAx-Min Flow Into Reservoir
         //   discharge_min = discharge;
        //} else {
        nSteps     = MFVarGetInt   (_MDInAvgNStepsID,       itemID,   0);
        
        if (nSteps <=2){
            discharge_max = 0.001;
            discharge_min = 999;
        } else {
            if (discharge > discharge_max) {
                discharge_max = discharge;
            } else {
            if (discharge < discharge_min) {
                discharge_min = discharge;
            }
            }
        }

 
        release_max = MFVarGetFloat(_MDOutReleaseMaxID,       itemID, 0);
        release_min = MFVarGetFloat(_MDOutReleaseMinID,       itemID, 0.0);
        
        if (nSteps <=2){
            release_max = 0.001;
            discharge_min = 999;
        }


        lastmonth   = MFVarGetFloat(_MDOutLastMonthID,        itemID, 1.0);
        mtdInflow   = MFVarGetFloat(_MDOutMonthToDayInFlowID, itemID, 2.0);

        /*
        printf("discharge WT Res= %f \n" , discharge);
        printf("discharge_max= %f \n" , discharge_max);
        printf("discharge_min= %f \n" , discharge_min);
        printf("release_max= %f \n" , release_max);
        printf("release_min= %f \n" , release_min);
        */

        if(m==0){
            m=2;
        }
        if (d==0){
            d=1;
        }
        if (lastmonth == m) {
            mtdInflow   = mtdInflow + discharge;
            avmtdInflow = mtdInflow / d;
        } else {
            mtdInflow =discharge;
            avmtdInflow = discharge;
        } 
        MFVarSetFloat(_MDOutMonthToDayInFlowID, itemID, mtdInflow);

        // Discharge and Release Are Standardized [0, 1]  

        //discharge_t_1   = (avmtdInflow - discharge_min) / (discharge_max - discharge_min); // This Month
        //discharge_t_1   = (discharge - discharge_min) / (discharge_max - discharge_min); // This Month Use Daily Discharge Unstead of Average Monthly
        discharge_t_2   = MFVarGetFloat(_MDOutDisch_t_2_ID,      itemID, 0.0001); // Last Month
        discharge_t_3   = MFVarGetFloat(_MDOutDisch_t_3_ID,      itemID, 0.0001); // Two Month Ago
        res_release_t_2 = MFVarGetFloat(_MDOutResRelease_t_2_ID, itemID, 0.0001); // Last Month
        res_release_t_3 = MFVarGetFloat(_MDOutResRelease_t_3_ID, itemID, 0.0001); // Two Month Ago
/*
        I1[0][0] = discharge_t_3;
        I1[1][0] = discharge_t_2;
        I1[2][0] = discharge_t_1;

        I2[0][0] = res_release_t_3;
        I2[1][0] = res_release_t_2;
*/
        discharge_t_1   = discharge;
        SD_t_1 = (avmtdInflow   - discharge_min) / (discharge_max - discharge_min);
        SD_t_2 = (discharge_t_2 - discharge_min) / (discharge_max - discharge_min);
        SD_t_3 = (discharge_t_3 - discharge_min) / (discharge_max - discharge_min);
        
        SR_t_3 = (res_release_t_3 - release_min) / (release_max - release_min);
        SR_t_2 = (res_release_t_2 - release_min) / (release_max - release_min);
        
        
        I1[0][0] = SD_t_3;
        I1[1][0] = SD_t_2;
        I1[2][0] = SD_t_1;

        I2[0][0] = SR_t_3;
        I2[1][0] = SR_t_2;
        
        I3 = m;

       /* printf ("First Input: I1 \n");
        printf ("%f \n",I1[0][0]);
        printf ("%f \n",I1[1][0]);
        printf ("%f \n",I1[2][0]);
        printf ("Second Input: I2 \n");
        printf ("%f \n",I2[0][0]);
        printf ("%f \n",I2[1][0]);
        printf ("Third Input: I3 \n");
        printf ("%f \n",I3); */
      
        ANN = ANNOUTPUT (I1, I2, I3)* (release_max - release_min) + release_min;
        if (nSteps <=2){
            if (ANN<0) {
                ANN= -ANN;
            }
        }        

       // printf("ANN: %f \n", ANNOUTPUT);
       // printf("EXP: %f \n", exp (2));

        prevResStorage = MFVarGetFloat(_MDOutResStorageID, itemID, 0);
        if (prevResStorage == 0){
            prevResStorage= 0.5*resCapacity;
        }

        resStorageChg = (discharge - ANN)*3600 * 24;
        minresStorage = resCapacity * 0.05;

        if (prevResStorage + resStorageChg < resCapacity && prevResStorage + resStorageChg > minresStorage) {
            SIMOUT = ANN;
            if (SIMOUT < 0) {
                printf("Error: Negative release (1)! \n");
                printf("%f %f %f %f %f %f %f\n", SIMOUT, release_max, release_min, resCapacity, resStorageChg, minresStorage, discharge);
            } 
            resStorage = prevResStorage + resStorageChg;
        } else {
            if (prevResStorage + resStorageChg > resCapacity) {
                SIMOUT = ((discharge * 3600 * 24)-(resCapacity - prevResStorage)) / (3600 * 24);
                if (SIMOUT < 0) {
                    printf("Error: Negative release (2)! \n");
                    printf("%f %f %f %f %f %f %f\n", SIMOUT, release_max, release_min, resCapacity, resStorageChg, minresStorage, discharge);
                }
                resStorage = resCapacity;
            } else {
                SIMOUT = (prevResStorage - minresStorage + (discharge * 3600 * 24)) / (3600 * 24);
                if (SIMOUT < 0) {
                    printf("Error: Negative release (3)! \n");
                    printf("%f %f %f %f %f %f %f %f %f\n", ANN, SIMOUT, release_max, release_min, prevResStorage, resCapacity, resStorageChg, minresStorage, discharge);
                }//                                         0.00  -0.000031  0.000000    0.000000       0.000000      270132.5     0.000000       2.701325      0.000000
                resStorage = minresStorage;
            }
        }
        
        if (nSteps <=2){
            if (SIMOUT<0) {
                SIMOUT= -SIMOUT;
            }
        }   

        mtdRelease = MFVarGetFloat(_MDOutMonthToDayReleaseID, itemID, 0.0);
        if (lastmonth == m) {
            mtdRelease   = mtdRelease + SIMOUT;
            avmtdRelease = mtdRelease / d;
        } else {
            mtdRelease   = SIMOUT;
            avmtdRelease = SIMOUT;
        }
        MFVarSetFloat(_MDOutMonthToDayReleaseID, itemID, mtdRelease);

            resRelease = SIMOUT;

        if (nSteps <=2){
            release_max = 0.001;
            release_min = 999;
        } else {
            if (resRelease > release_max) {
                release_max = resRelease;
            } else {
            if (resRelease < release_min) {
                release_min = resRelease;
            }
            }
        }


        res_release_t_1 = resRelease;
        MFVarSetFloat(_MDOutDischMaxID, itemID, discharge_max);
        MFVarSetFloat(_MDOutDischMinID, itemID, discharge_min);
        
        MFVarSetFloat(_MDOutReleaseMaxID,       itemID, release_max);
        MFVarSetFloat(_MDOutReleaseMinID,       itemID, release_min);
        MFVarSetFloat(_MDOutResReleaseID,       itemID, resRelease);
        MFVarSetFloat(_MDOutLastMonthID,        itemID, m);
        //MFVarSetFloat(_MDOutReservoirReleaseID, itemID, res_release_t_1);
        MFVarSetFloat(_MDOutDisch_t_3_ID,       itemID, discharge_t_2);   //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutDisch_t_2_ID,       itemID, discharge_t_1);   //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_3_ID,  itemID, res_release_t_2); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_2_ID,  itemID, res_release_t_1); //You Should Set t-2 before t-1
        
        MFVarSetFloat(_MDOutResStorageChgID,    itemID, resStorageChg);
        MFVarSetFloat(_MDOutResStorageID,       itemID, resStorage);

    //}
    

}

static void _MDReservoirDW(int itemID) {

    // Input
    float discharge; // Current discharge [m3/s]
    float meanDischarge; // Long-term mean annual discharge [m3/s]
    float resCapacity; // Reservoir capacity [km3]

    // Output
    float resStorage; // Reservoir storage [km3]
    float resStorageChg; // Reservoir storage change [km3/dt]
    float resRelease; // Reservoir release [m3/s] 

    // local
    float prevResStorage; // Reservoir storage from the previous time step [km3]
    float dt; // Time step length [s]
    float balance; // water balance [m3/s]

    // Parameters
    float drySeasonPct = 0.60; // RJS 071511
    float wetSeasonPct = 0.16; // RJS 071511
    float year = 0; // RJS 082311

    discharge = MFVarGetFloat(_MDInDischargeID, itemID, 0.0);
    //printf("discharge= %f \n", discharge);
    meanDischarge = MFVarGetFloat(_MDInDischMeanID, itemID, discharge);
    year = MFDateGetCurrentYear();

    if ((resCapacity = MFVarGetFloat(_MDInResCapacityID, itemID, 0.0)) <= 0.0) {
        MFVarSetFloat(_MDOutResStorageID, itemID, 0.0);
        MFVarSetFloat(_MDOutResStorageChgID, itemID, 0.0);
        MFVarSetFloat(_MDOutResReleaseID, itemID, discharge);
        printf("discharge NO Res= %f \n" , discharge);
        return;
    }
    printf("discharge WT Res= %f \n" , discharge);
    dt = MFModelGet_dt();
    prevResStorage = MFVarGetFloat(_MDOutResStorageID, itemID, 0.0);

    resRelease = discharge > meanDischarge ? wetSeasonPct * discharge : drySeasonPct * discharge + (meanDischarge - discharge);

    resStorage = prevResStorage + (discharge - resRelease) * 86400.0 / 1e9;

    if (resStorage > resCapacity) {
        resRelease = discharge * dt / 1e9 + prevResStorage - resCapacity;
        resRelease = resRelease * 1e9 / dt;
        resStorage = resCapacity;
    } else if (resStorage < 0.0) {
        resRelease = prevResStorage + discharge * dt / 1e9;
        resRelease = resRelease * 1e9 / dt;
        resStorage = 0;
    }

    resStorageChg = resStorage - prevResStorage;

    balance = discharge - resRelease - (resStorageChg / 86400 * 1e9); // water balance

    MFVarSetFloat(_MDOutResStorageID, itemID, resStorage);
    MFVarSetFloat(_MDOutResStorageChgID, itemID, resStorageChg);
    MFVarSetFloat(_MDOutResReleaseID, itemID, resRelease);
}

enum {
    MDnone, MDcalculate, MDneuralnet
};

int MDReservoirDef() {
    int optID = MFUnset;
    const char *optStr, *optName = MDOptReservoirs;
    const char *options [] = {MDNoneStr, MDCalculateStr, "neuralnet", (char *) NULL};

    if ((optStr = MFOptionGet(optName)) != (char *) NULL) optID = CMoptLookup(options, optStr, true);

    if ((optID == MDnone) || (_MDOutResReleaseID != MFUnset)) return (_MDOutResReleaseID);

    MFDefEntering("Reservoirs");
    switch (optID) {
        case MDcalculate:
            if (    ((_MDInDischMeanID      = MDDischMeanDef())   == CMfailed) ||
                    ((_MDInDischargeID      = MDDischLevel2Def()) == CMfailed) ||
                    ((_MDInResCapacityID    = MFVarGetID(MDVarReservoirCapacity,            "km3",  MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutResStorageID    = MFVarGetID(MDVarReservoirStorage,             "km3",  MFOutput, MFState, MFInitial))  == CMfailed) ||
                    ((_MDOutResStorageChgID = MFVarGetID(MDVarReservoirStorageChange,       "km3",  MFOutput, MFState, MFBoundary)) == CMfailed) || //RJS, changed MFBoundary o MFIniial
                    ((_MDOutResReleaseID    = MFVarGetID(MDVarReservoirRelease,             "m3/s", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirDW) == CMfailed)) return (CMfailed);
            break;
        case MDneuralnet:

            if (    ((_MDInAvgNStepsID          = MDAvgNStepsDef ())  == CMfailed) ||
                    ((_MDInDischargeID          = MDDischLevel2Def()) == CMfailed) ||
                    ((_MDInResCapacityID        = MFVarGetID(MDVarReservoirCapacity,      "m3",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutDisch_t_1_ID        = MFVarGetID(MDVarDisch_t_1_,             "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_2_ID        = MFVarGetID(MDVarDisch_t_2_,             "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_3_ID        = MFVarGetID(MDVarDisch_t_3_,             "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutDischMinID          = MFVarGetID(MDVarDischMin,               "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutDischMaxID          = MFVarGetID(MDVarDischMax,               "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutReleaseMinID        = MFVarGetID(MDVarReleaseMin,             "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutReleaseMaxID        = MFVarGetID(MDVarReleaseMax,             "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutLastMonthID         = MFVarGetID(MDVarLastMonth,              "m"   , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutMonthToDayInFlowID  = MFVarGetID(MDVarMonthToDayInFlow,       "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutMonthToDayReleaseID = MFVarGetID(MDVarMonthToDayRelease,      "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    //((_MDOutReservoirReleaseID  = MFVarGetID(MDVarReservoirRelease,       "m3/s", MFOutput, MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutResStorageID        = MFVarGetID(MDVarReservoirStorage,       "m3"  , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageChgID     = MFVarGetID(MDVarReservoirStorageChange, "m3"  , MFOutput, MFState, MFInitial)) == CMfailed) || //RJS, changed MFBoundary o MFIniial
                    ((_MDOutResReleaseID        = MFVarGetID(MDVarReservoirRelease,       "m3/s", MFOutput, MFFlux,  MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_1_ID   = MFVarGetID(MDVarResRelease_t_1_,        "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_2_ID   = MFVarGetID(MDVarResRelease_t_2_,        "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_3_ID   = MFVarGetID(MDVarResRelease_t_3_,        "m3/s", MFOutput, MFFlux, MFInitial)) == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirNeuralNet) == CMfailed)) return (CMfailed);
            
            if (((optStr = MFOptionGet (MDOptHydroPower)) != (char *) NULL) && (CMoptLookup (options,optStr,true) == CMfailed)) {
		if (((_MDInMegaWattID       = MFVarGetID(MDVarMegaWatt,                     "MW",   MFOutput, MFState,  MFInitial))  == CMfailed))
	    	return (CMfailed);
            
            
  /*           if (((optStr = MFOptionGet (MDOptHydroPower)) != (char *) NULL) && (CMoptLookup (options,optStr,true) == CMfailed)) {
		if ((_MDInMegaWattID =MDHydroPowerDef()) == CMfailed) return (CMfailed);
            }
  */        break;
        default: MFOptionMessage(optName, optStr, options);
            return (CMfailed);
    }
    MFDefLeaving("Reservoirs");
    return (_MDOutResReleaseID);
}
