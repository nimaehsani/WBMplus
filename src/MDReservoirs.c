/******************************************************************************
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDReservoirs.c

nehsani@ccny.cuny.edu
dominik.wisser@unh.edu

Updated with a Neural Network function for
Reservoir Operation.
2013 NE
*******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input
static int _MDInDischargeID   = MFUnset;
static int _MDInDischMeanID   = MFUnset;
static int _MDInResCapacityID = MFUnset;

// Output
static int _MDOutResStorageID      = MFUnset;
static int _MDOutPreResStorageID   = MFUnset;
static int _MDOutResStorageChgID   = MFUnset;
static int _MDOutResReleaseID      = MFUnset;
static int _MDOutResRelease_t_1_ID = MFUnset;
static int _MDOutResRelease_t_2_ID = MFUnset;
static int _MDOutResRelease_t_3_ID = MFUnset;
static int _MDOutDisch_t_1_ID      = MFUnset;
static int _MDOutDisch_t_2_ID      = MFUnset;
static int _MDOutDisch_t_3_ID      = MFUnset;


static float ANNOUTPUT(float I1[3][1], float I2[2][1], float I3) { //Do Not MAke Any Changes Unless You Have a New ANN

  	//load('ThirdNet.mat')//

    float FirstLayerBias[6][1] = {
		{-4.8895325263214628},
		{-4.2422210723319198},
		{2.0905793315974699},
		{0.062997739670340558},
		{3.4464068164930848},
		{-4.1917667542064212},
	 };

    float FirstLayerWeight [6][3] = {
        {2.6794945718292671, -1.8860215868651271, 3.8920537842478833},
        {5.0518600975427157, 0.58440858895092118, -0.1809806525515216},
        {-1.9107162268023963, 4.4304306759814951, 3.245958891312041},
        {5.2568153836149412, 5.640102302513899, 0.54757652946826185},
        {1.7919171951657034, -4.6431706287207151, -1.0662755569236482},
        {-2.7452104159420112, 3.2105699689452489, 2.8364581448264503},
    };

    float SecondLayerWeight[4][2] = {
        {3.5942082756296103, 4.2943610222677089},
        {4.0705919037941269, 3.8446565271485422},
        {3.8085937875413247, 4.3569429583894914},
        {-4.8095967403969526, 2.804018069823321},
    };

    float SecondLayerBias [4][1] = {
		{-6.7443900577024927},
		{-4.8989411723704768},
		{-1.2629061907389167},
		{-2.0349323312021075},
    };

    float ThirdLayerWeight[2][1] = {
		{5.7306281576798348},
		{-5.6983330355644579},
    };

    float ThirdLayerBias [2][1] = {
		{-5.674311196800347},
		{-1.414839271870252},
    };

    float FourthLayerWeight1[6][6] = {
		{0.59567976588053761, -0.42919386892406064, -3.1046255138463863, -2.031376203619728, -0.003169071524200952, 0.46422920611620572},    
		{-0.38860966694126992, 0.71546117494102701, 0.63256474302958321, 2.0631335661157193, -0.7436704516276611, -0.29807803523773141},
		{0.44644882591245294, 0.42693726550109096, 0.37924395004334754, 0.32838150484835671, -0.75012892309493773, 0.52495411888003063},
		{-0.90755022140724995, -0.42085199334330003, 0.39145131689240348, 0.31729932166605329, 0.10605839509496631, -0.52660332238125251},
		{-1.2593056682703245, -0.27880615471243519, -1.1357584470045672, -1.1255215407785255, 1.6019643123649754, -0.92553464260621088},
		{-0.83784993451087131, -0.85867478295740773, 0.63309931052669299, 4.8385639509044935, -0.23696517402363465, 0.29606380533166787},
	};

    float FourthLayerWeight2[6][4] = {
		{1.026847186296679, 0.93084659386990265, -0.46559859650904062, 1.0087503984655652},
		{0.94986023154975252, -0.24515885817439695, 0.15965582247548105, -0.084507368370036529},
		{0.17064610716527154, 0.63444226840100992, 0.29585412006070971, 0.80420007394364079},
		{-0.72531745825869021, -0.54636510122067761, -0.15646125993642854, -0.46040383494843767},
		{-1.4452237986760499, -0.30964228054349008, -1.1455804967879082, -0.50303387975699898},
		{-0.47269390814290246, -0.30727184799941598, 1.3920403947490985, -0.17152102598998564},
    };

    float FourthLayerWeight3[6][2] = {
		{1.4439172198184893, 0.45388239499909921},
		{0.2866904468317571, -0.91626401711486505},
		{-0.60786043266814516, 0.20994888632583841},
		{-0.7060221589545238, -0.16849240104849023},
		{0.6342608930164797, -0.54649519225093968},
		{-0.0066559337628386561, -0.23313098033012913},
    };

    float FourthLayerBias[6][1] = {
		{-1.9399653673860706},
		{0.80567524388644829},
		{-0.4137462842858215},
		{-0.21639239964056947},
		{-0.04059493473381641},
		{-1.572847505807266},
	};

    float FifthLayerWeight [1][6] = {
		{-2.050037504773798, 3.0105771035841253, 2.2787786536727679, -2.8542908258666677, -3.918050145213408, 3.9155779322676665},
    };

    float FifthLayerBias = -3.186701263519236;
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


static float ANNOUTPUT2(float I1[3][1], float I2[2][1], float I3) { //Do Not MAke Any Changes Unless You Have a New ANN
	//load('NET2.mat')//

    float FirstLayerBias[6][1] = {
		{-5.4809843376405594},
		{-5.3356092603822143},
		{3.9943069103529747},
		{-0.41793111778343439},
		{0.52164853846445347},
		{-2.1258056711126687},
	 };

    float FirstLayerWeight [6][3] = {
		{3.5701652993129285, 3.4421731676111036, -1.1372085688206515},
		{3.2662789867946485, 3.8713942775562606, 0.47468008881130808},
		{-2.9641176055273681, -4.1339224270062536, 0.11388697197324653},
		{1.0882340809083537, 5.1160331710319573, -1.0196907015138945},
		{1.3329392275433336, 5.0054780204444134, 3.5827525670996208},
		{-6.3291251837696016, -2.3864787278331745, -1.4964509294577444},
    };

    float SecondLayerWeight[4][2] = {
		{2.9931599106323565, 4.7329174543236894},
		{-5.5151222107223044, -0.94710233332376947},
		{-1.1530846954498826, -13.59986436567403},
		{-4.021429027144638, 3.9193131082822963},
    };

    float SecondLayerBias [4][1] = {
		{-6.6634218043389994},
		{4.2100854406186343},
		{-0.051557487525120342},
		{-2.5910494853166908},
    };

    float ThirdLayerWeight[2][1] = {
		{5.3749171820088932},
		{-5.6016850634657329},
    };

    float ThirdLayerBias [2][1] = {
		{-5.6826819780009536},
		{-0.80279829403308589},
    };

    float FourthLayerWeight1[6][6] = {
		{-0.64107839974325942, -0.91291101194364632, 1.1106113365389383, 0.77423919372195538, -1.0815297220582916, 3.8672288925820628},
		{-0.6705578811857601, 0.31756410975285798, 0.43365183416967029, -0.39514061389827898, 0.63864150643258821, 1.4554404306158881},
		{-0.17499129214751333, -0.97525172843551489, 0.86715279984608118, -1.0361623211741953, 0.031543319482775542, 0.072873960341362579},
		{-0.38357133556058426, 0.26195993302273723, 0.71647486388234172, -1.0959860319811008, -2.2306541344421387, 0.88907968671398852},
		{-0.47675687715152387, -0.3047210727563715, 0.46924659749425307, 0.43068700512880737, -0.13399867000639926, -0.66017545403691869},
		{0.32804211761772728, -0.46820037178215129, -0.30486713298573137, 1.1339205377051609, -0.35445240528960298, 0.30941389388713286},
	};

    float FourthLayerWeight2[6][4] = {
		{0.12638269026648186, 0.34670904705605676, 12.332749891806046, 0.18034786075562786},
		{0.59527544903200835, 0.085197202949219356, -2.1124124593738705, -0.39934434343486314},
		{-1.5822316212431684, 0.24685500786798323, 0.29075763191452569, -0.93491035346606588},
		{0.33570393969884593, 0.84281873387360395, -0.57533719925859261, -0.97659102839170753},
		{-0.62025928743149672, 0.075046426729973664, -0.52650894431957684, -0.27834146727155462},
		{0.52482280355612221, 0.52414512037275385, 0.40257177869249311, -0.27688524274849602},
	};

    float FourthLayerWeight3[6][2] = {
		{0.043932815991353125, 0.24911272559419509},
		{0.37028729736592131, 0.17238131148621758},
		{-1.2561529007116925, -0.95726503529769447},
		{-0.59807596692332354, -0.35805181279549236},
		{0.12227542506276243, -0.84849888933061446},
		{0.1945436450183344, 0.12667948525302025},
	};

    float FourthLayerBias[6][1] = {
		{2.082569098546124},
		{0.92746776505107764},
		{0.47749248410557088},
		{-0.12759069416577937},
		{-0.78119797669207691},
		{1.3131580011884232},
	};

    float FifthLayerWeight [1][6] = {
		{-5.9855754805412191, 2.0568290993389731, -2.0452036306293744, -3.2894942805034328, -1.5141235197668719, -4.6282617655087588},
    };

    float FifthLayerBias = -4.37912338405821;
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

    float ANNOUTPUT2;
    ANNOUTPUT2 = FifthLayerOut;

    return (ANNOUTPUT2);
}

 
 
static void _MDReservoirNeuralNet(int itemID) {

    float discharge;      // Current discharge [m3/s]
    float resCapacity;    // Reservoir capacity [m3]
    float minresStorage;
    float discharge_t_1;
    float discharge_t_2;
    float discharge_t_3;

    float I1[3][1];        // Input to ANN (ANNOUTPUT.c) 
    float I2[2][1];        // Input to ANN (ANNOUTPUT.c) 
    float I3;              // Input to ANN (ANNOUTPUT.c) 
    float ANN;
    float ANN1;
    float ANN2;
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
    int   y             = MFDateGetCurrentYear();

    discharge = MFVarGetFloat(_MDInDischargeID,   itemID, 0.0);
   
    if (((resCapacity = MFVarGetFloat(_MDInResCapacityID,    itemID, 0.0)) <= 0.0) ||  y<=1900 ){
                       MFVarSetFloat(_MDOutResStorageID,    itemID, 0.0);
                       MFVarSetFloat(_MDOutResStorageChgID, itemID, 0.0);
                       MFVarSetFloat(_MDOutResReleaseID,    itemID, discharge);
        return;
    } 
        resCapacity   = MFVarGetFloat(_MDInResCapacityID, itemID, 0.0);  
        prevResStorage = MFVarGetFloat(_MDOutPreResStorageID, itemID, 0*resCapacity);

	 if (prevResStorage<=0*resCapacity){
           prevResStorage=0*resCapacity;
        }

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * !!!!!! 0.2*resCapacity should be added to Reservoir storage result file !!!!!!!!
 * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


        discharge_t_1   = discharge;       
        discharge_t_2   = MFVarGetFloat(_MDOutDisch_t_2_ID,      itemID, 0); // Last Month
        discharge_t_3   = MFVarGetFloat(_MDOutDisch_t_3_ID,      itemID, 0); // Two Month Ago
        res_release_t_2 = MFVarGetFloat(_MDOutResRelease_t_2_ID, itemID, 0); // Last Month
        res_release_t_3 = MFVarGetFloat(_MDOutResRelease_t_3_ID, itemID, 0); // Two Month Ago

        SD_t_1 = discharge_t_1 *24*3600/resCapacity;
        SD_t_2 = discharge_t_2 *24*3600/resCapacity;
        SD_t_3 = discharge_t_3 *24*3600/resCapacity;
        
        SR_t_3 = res_release_t_3 *24*3600/resCapacity;
        SR_t_2 = res_release_t_2 *24*3600/resCapacity;
        
        I1[0][0] = SD_t_3;
        I1[1][0] = SD_t_2;
        I1[2][0] = SD_t_1;

        I2[0][0] = SR_t_3;
        I2[1][0] = SR_t_2;
        

        I3 = prevResStorage/resCapacity;  

 
        ANN1 = ANNOUTPUT (I1, I2, I3)*resCapacity/(24*3600) ;
        ANN2 = ANNOUTPUT2 (I1, I2, I3)*resCapacity/(24*3600) ;
        ANN=(ANN1+ANN2)/2;
        


        

        resStorageChg = (discharge - ANN)*3600 * 24;
////////
        minresStorage =0*resCapacity;
////////
        

        if ((prevResStorage + resStorageChg <= resCapacity) && (prevResStorage + resStorageChg >= minresStorage)) {
            SIMOUT = ANN;
             
            if (SIMOUT < 0) {
                printf("Error: Negative release (1)! \n");
                printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge);
            } 
            resStorage = prevResStorage + resStorageChg;
        } else {
            if ((prevResStorage + resStorageChg) > resCapacity) {
                SIMOUT = ((discharge * 3600 * 24)-resCapacity + prevResStorage) / (3600 * 24);
                
                if (SIMOUT < 0) {
                    printf("Error: Negative release (2)! \n");
                        printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge);                }
                resStorage = resCapacity;
            } else {
                SIMOUT = (prevResStorage - minresStorage + (discharge * 3600 * 24)) / (3600 * 24);
               
                if (SIMOUT < 0) {
                    printf("Error: Negative release (3)! \n");
                        printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge);                }
                resStorage = minresStorage;
            }
        }
                //  printf("%f %f %f %f %f  %f %f\n",SIMOUT, ANN, discharge, resCapacity, resStorageChg, minresStorage, resStorage); 
     
///////////////////////////////////////////////////////////////////
/////////////////// Maximum Discharge ///////////////////////
/*        if (SIMOUT > (resCapacity/(3.21*24*3600))){
            resRelease=(resCapacity/(3.21*24*3600));
            resStorage = prevResStorage + (discharge * 3600 * 24)-(resRelease * 3600 * 24);
            if (resStorage > resCapacity) {
                resRelease = ((discharge * 3600 * 24)-(resCapacity - prevResStorage)) / (3600 * 24);
                resStorage = resCapacity;
            }
                
        } else {
            resRelease = SIMOUT;
        }
 */
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////// 

                      resRelease = SIMOUT;

        res_release_t_1 = resRelease;

        MFVarSetFloat(_MDOutResReleaseID,       itemID, resRelease);
        MFVarSetFloat(_MDOutDisch_t_3_ID,       itemID, discharge_t_2);   //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutDisch_t_2_ID,       itemID, discharge_t_1);   //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_3_ID,  itemID, res_release_t_2); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_2_ID,  itemID, res_release_t_1); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResStorageChgID,    itemID, resStorageChg);
        MFVarSetFloat(_MDOutResStorageID,       itemID, resStorage);
        MFVarSetFloat(_MDOutPreResStorageID,    itemID, resStorage);
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
                    ((_MDOutResStorageChgID = MFVarGetID(MDVarReservoirStorageChange,       "km3",  MFOutput, MFState, MFBoundary)) == CMfailed) || 
                    ((_MDOutResReleaseID    = MFVarGetID(MDVarReservoirRelease,             "m3/s", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     (                        MDHydroPowerDef()   == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirDW) == CMfailed)) return (CMfailed);
            break;
        case MDneuralnet:

            if (    ((_MDInDischargeID          = MDDischLevel2Def()) == CMfailed) ||
                    ((_MDInResCapacityID        = MFVarGetID(MDVarReservoirCapacity,      "m3",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutDisch_t_1_ID        = MFVarGetID(MDVarDisch_t_1_,             "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_2_ID        = MFVarGetID(MDVarDisch_t_2_,             "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_3_ID        = MFVarGetID(MDVarDisch_t_3_,             "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageID        = MFVarGetID(MDVarReservoirStorage,       "m3"  , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutPreResStorageID     = MFVarGetID(MDVarPreResStorage,          "m3"  , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageChgID     = MFVarGetID(MDVarReservoirStorageChange, "m3"  , MFOutput, MFState, MFInitial)) == CMfailed) || 
                    ((_MDOutResReleaseID        = MFVarGetID(MDVarReservoirRelease,       "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_1_ID   = MFVarGetID(MDVarResRelease_t_1_,        "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_2_ID   = MFVarGetID(MDVarResRelease_t_2_,        "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_3_ID   = MFVarGetID(MDVarResRelease_t_3_,        "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ( MDHydroPowerDef()   == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirNeuralNet) == CMfailed)) return (CMfailed);
              break;
        default: MFOptionMessage(optName, optStr, options);
            return (CMfailed);
    }
    MFDefLeaving("Reservoirs");
    return (_MDOutResReleaseID);
}


