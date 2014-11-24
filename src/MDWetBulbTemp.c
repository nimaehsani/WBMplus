/******************************************************************************
NIMA EHSANI
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2014, UNH - CCNY/CUNY

MDWetBulbTemp.c

nehsani@ccny.cuny.edu
brosenzweig@ccny.cuny.edu
amiara@ccny.cuny.edu

 *******************************************************************************/
/********************************************************************************
 * Calculates wet-bulb temperature from temperature, pressure and specific
humidity using the approximation described in Chappell, et al 1973. 
This approximation is used as a first guess and the tradition Clausius- 
Clapeyron solution using Newton-Raphson Iteration Method is applied for 
better results.  If this solution does not converge,the Chappell 
approximation is returned. 
 * ******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input
static int _MDInAirTemperatureID = MFUnset;
static int _MDInSpecificHumidityID = MFUnset;
static int _MDInAirPressureID = MFUnset;

// Output
static int _MDOutWetBulbTempID = MFUnset;

static void _MDWetBulbTemp(int itemID) {
    float airtemp;
    float specifichumidity;
    float airpressure;
    float wetbulbtemp;
    float eps = 0.1; //threshold value for iteration 
    float cp = 1005; //specific heat of air (J/(kg*K)
    float e; //e (HPa)
    float es; //Sat. vapor press. (hPa) 
    float w; //Mixing ratio from specific humidity (kg/kg)
    float ws; //Mixing ratio at saturation (kg/kg)
    float td; //dewpoint (Deg. C)
    float a; //Dewpoint depression
    float b; //Mixing ratio depression
    float c;
    float twn; //Chappell Approximation for Wet Bulb Temperature
    float twn2;
    float lv; //latent heat of vaporization (J/kg)
    float firstsol;
    float delta;
    int x;
    float esw; //Sat. vapor press. (hPa) at Twn


    airtemp          = MFVarGetFloat(_MDInAirTemperatureID,   itemID, 0.0);
    airpressure      = MFVarGetFloat(_MDInAirPressureID,      itemID, 0.0) / 100; //pressure (HPa)
    specifichumidity = MFVarGetFloat(_MDInSpecificHumidityID, itemID, 0.0) * 1000; //specific humidity (g/kg)

    e  = (specifichumidity * airpressure) / 622; //e (HPa)
    es = 6.112 * exp(17.27 * (airtemp / (237.3 + airtemp))); //Sat. vapor press. (hPa) 
    
    if (e >= es) {//at or above saturation, wetbulbtemp=airtemp
        wetbulbtemp = airtemp;
    } else {
        if (airtemp == 32800 || airpressure == 32800 || specifichumidity == 32800) { //missing data
            wetbulbtemp = 32800;
        } else {
            w   = specifichumidity / (1000) / (1 - (specifichumidity / 1000)); //Mixing ratio from specific humidity (kg/kg)
            ws  = 0.62197 * (es / (airpressure - es)); //Mixing ratio at saturation (kg/kg)
            td  = (243.5 * log(exp(1) / 6.11)) / (17.67 - log(exp(1) / 6.11)); //dewpoint (Deg. C)
            a   = (airtemp - td); //Dewpoint depression
            b   = (ws - w); //Mixing ratio depression
            c   = (597.3 - 0.566 * (airtemp - 273.16)) / 0.24;
            twn = airtemp - ((a * b * c) / (a + b * c)); //Chappell Approximation for Wet Bulb Temperature

            //Iterate for a better solution
            w = specifichumidity / (1000) / (1 - (specifichumidity / 1000)); //Mixing ratio from specific humidity (kg/kg)
            lv = 1918460*pow(((airtemp + 273) / ((airtemp + 273) - 33.91)),2); //latent heat of vaporization (J/kg)
            firstsol = twn;
            delta = eps + 1;

            for (x = 0; x < 6; x++) {
                esw = 611 * exp(17.27 * (twn / (237.3 + twn))) / 100; //Sat. vapor press. (hPa) at Twn
                ws = 0.62197 * (esw / (airpressure - esw)); //dimensionless mixing ratio at saturation
                twn2 = (twn + airtemp + lv / cp * (w - ws)) / 2; //average with prev. value for stability
                delta = abs(twn2 - twn);
                if (delta <= eps) {
                    break;
                } else {
                    if (x == 6) { //%alternate solution for nonconvergence
                        twn2 = firstsol; //Return the Chappell approximation
                    } else {
                        twn = twn2;
                    }
                }
            }
            wetbulbtemp = twn2; //Return Final Wet Bulb Temperature
        }
    }

    MFVarSetFloat(_MDOutWetBulbTempID, itemID, wetbulbtemp);
}

    enum {
        MDnone, MDcalculate
    };

    int MDWetBulbTempDef() {
        int optID = MFUnset;
        const char *optStr, *optName = MDOptWetBulbTemp;
        const char *options [] = {MDNoneStr, MDCalculateStr, (char *) NULL};

        if ((optStr = MFOptionGet(optName)) != (char *) NULL) optID = CMoptLookup(options, optStr, true);
        if ((optID == MDnone) || (_MDOutWetBulbTempID != MFUnset)) return (_MDOutWetBulbTempID);

        MFDefEntering("WetBulbTemp");

        switch (optID) {
            case MDcalculate: //_MDInAirPressureID
                if (    ((_MDInAirTemperatureID =   MFVarGetID(MDVarAirTemperature, "DegC", MFInput, MFState, MFBoundary)) == CMfailed) || //((_MDInResReleaseID        = MFVarGetID(MDVarReservoirRelease,       "m3/s", MFInput, MFState, MFInitial)) == CMfailed) ||
                        ((_MDInSpecificHumidityID = MFVarGetID(MDVarSpecificHumidity, "kg/kg", MFInput, MFState, MFBoundary)) == CMfailed) ||
                        ((_MDInAirPressureID = MFVarGetID(MDVarAirPressure, "Pa", MFInput, MFState, MFBoundary)) == CMfailed) ||
                        ((_MDOutWetBulbTempID = MFVarGetID(MDVarWetBulbTemp, "DegC", MFOutput, MFState, MFInitial)) == CMfailed) ||
                        ((MFModelAddFunction(_MDWetBulbTemp) == CMfailed))
                        ) return (CMfailed);
                break;
            default: MFOptionMessage(optName, optStr, options);
                return (CMfailed);
        }
        MFDefLeaving("WetBulbTemp");
        return (_MDOutWetBulbTempID);
    }
