#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "essentials.h"
#include "prms.h"
#include "assert.h"



// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double* yic;  
extern prms* ppc;
extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;
extern FILE* OutFile;
extern double G_CLO_TF;
extern int G_CLO_STEPS_PER_DAY;
extern string G_CLO_LOCATION;
extern double G_CLO_P_HOSP_TO_ICU;
extern bool G_B_DIAGNOSTIC_MODE;
extern bool G_B_CHECKPOP_MODE;
extern bool G_B_USESOCIALCONTACTMATRIX;
extern bool G_B_REMOVE_99COLUMNS;
extern bool G_B_BINARY_OUTPUT;
extern bool G_B_DEATHRATE_OUTPUT;

extern double G_CLO_SYMP_FRAC;          //SHOULD BE DEPRECATED SOON
extern double G_CLO_SYMP_FRAC_EQUAL;
extern double G_B_SYMP_FRAC_DAVIES_20200616;
extern double G_B_SYMP_FRAC_SIMPLEAVERAGE;



extern double G_CLO_HOSPFRAC_YOUNG_DEV; //DEPRECATED
extern double G_CLO_HOSPFRAC_OLD_DEV;   //DEPRECATED
extern double G_CLO_HOSPFRAC_MID_DEV;   //DEPRECATED

extern double G_CLO_HOSPFRAC_10;
extern double G_CLO_HOSPFRAC_20;
extern double G_CLO_HOSPFRAC_30;
extern double G_CLO_HOSPFRAC_40;
extern double G_CLO_HOSPFRAC_50;
extern double G_CLO_HOSPFRAC_60;
extern double G_CLO_HOSPFRAC_70;
extern double G_CLO_HOSPFRAC_80;

// extern double G_CLO_MIXINGLEVEL[];
// extern double G_CLO_MIXINGLEVEL_POSTLD[];
extern int G_CLO_FIRSTLOCKDOWN_ENDDAY;
extern bool G_CLO_POSTLD_MIXING_SET;


extern double G_CLO_DEATHPROB_HOME_60;
extern double G_CLO_DEATHPROB_HOME_70;
extern double G_CLO_DEATHPROB_HOME_80;

extern double G_CLO_DEATHPROB_POSTVENT;

extern double G_CLO_PROB_NONICU_DEATH_80;

extern double G_CLO_MIN_MIXINGLEVEL_00;
extern double G_CLO_MIN_MIXINGLEVEL_10;
extern double G_CLO_MIN_MIXINGLEVEL_20;
extern double G_CLO_MIN_MIXINGLEVEL_30;
extern double G_CLO_MIN_MIXINGLEVEL_40;
extern double G_CLO_MIN_MIXINGLEVEL_50;
extern double G_CLO_MIN_MIXINGLEVEL_60;
extern double G_CLO_MIN_MIXINGLEVEL_70;
extern double G_CLO_MIN_MIXINGLEVEL_80;

extern double G_CLO_RELSUSC_0;
extern double G_CLO_RELSUSC_10;
extern double G_CLO_RELSUSC_20;
extern double G_CLO_RELSUSC_30;
extern double G_CLO_RELSUSC_40;
extern double G_CLO_RELSUSC_50;
extern double G_CLO_RELSUSC_60;
extern double G_CLO_RELSUSC_70;
extern double G_CLO_RELSUSC_80;

extern double G_CLO_TIME_SYMPTOHOSP;
extern double G_CLO_SELFISOLATION_FACTOR;

extern double G_CLO_EARLYMARCH_HOSPRATE;     // scaling factor showing how much more likely hospitalization was in early March than later in the epidemic
extern double G_CLO_EARLYMARCH_ENDDAY;       // the day that the high-hosp rate ended (probably because the patient load became too high)


extern double G_CLO_ICUFRAC_DEV;
extern double G_CLO_ICUFRAC_DEV_SECONDPHASE;
extern int G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY;
extern double G_CLO_PROB_ICU_TO_VENT;
extern double G_CLO_VENTDEATH_MID_DEV;
extern double G_CLO_VENTDEATH_70_DEV;
extern double G_CLO_VENTDEATH_80_DEV;
extern double G_CLO_MEANTIME_ON_VENT_SURV;
extern double G_CLO_RELATIVE_BETA_HOSP;
extern double G_CLO_DEV_LEN_HOSPSTAY;
extern double G_CLO_LEN_NAT_IMMU;
extern double G_CLO_REPORTING_RATE;


extern int G_CLO_NORMALCY_BEGINDAY;

//
// vaccination-related
//
extern double G_CLO_VAC_1_FOI;
extern bool G_B_TEST_BEFORE_VACCINATE;
extern bool G_B_USE_VAC_FRAC;
extern bool G_B_USE_VAC_RATIO;

extern double G_CLO_VAC_1_PROTECT_DURATION;
extern double G_CLO_VAC_1_EFF_HALFLIFE;
extern double G_CLO_VAC_1_EFF_SLOPE;
// extern int G_CLO_VAC_1_PHASE1_BEGINDAY;
// extern int G_CLO_VAC_1_PHASE1_ENDDAY;
// extern int G_CLO_VAC_1_PHASE1_DPD;
extern int G_CLO_VAC_1_PHASE2_BEGINDAY;
extern int G_CLO_VAC_1_PHASE2_ENDDAY;
// extern double G_CLO_VAC_1_FRAC_00;
// extern double G_CLO_VAC_1_FRAC_10;
// extern double G_CLO_VAC_1_FRAC_20;
// extern double G_CLO_VAC_1_FRAC_30;
// extern double G_CLO_VAC_1_FRAC_40;
// extern double G_CLO_VAC_1_FRAC_50;
// extern double G_CLO_VAC_1_FRAC_60;
// extern double G_CLO_VAC_1_FRAC_70;
// extern double G_CLO_VAC_1_FRAC_80;
extern double G_CLO_VAC_1_REL_EFF_00;
extern double G_CLO_VAC_1_REL_EFF_10;
extern double G_CLO_VAC_1_REL_EFF_20;
extern double G_CLO_VAC_1_REL_EFF_30;
extern double G_CLO_VAC_1_REL_EFF_40;
extern double G_CLO_VAC_1_REL_EFF_50;
extern double G_CLO_VAC_1_REL_EFF_60;
extern double G_CLO_VAC_1_REL_EFF_70;
extern double G_CLO_VAC_1_REL_EFF_80;

/* extern double G_CLO_VAC_2_FOI;
extern double G_CLO_VAC_2_PROTECT_DURATION;
extern double G_CLO_VAC_2_EFF_HALFLIFE;
extern double G_CLO_VAC_2_EFF_SLOPE;
extern int G_CLO_VAC_2_PHASE1_BEGINDAY;
extern int G_CLO_VAC_2_PHASE1_ENDDAY;
extern int G_CLO_VAC_2_PHASE1_DPD;
extern int G_CLO_VAC_2_PHASE2_BEGINDAY;
extern int G_CLO_VAC_2_PHASE2_ENDDAY;
// extern double G_CLO_VAC_2_FRAC_00;
extern double G_CLO_VAC_2_FRAC_10;
extern double G_CLO_VAC_2_FRAC_20;
extern double G_CLO_VAC_2_FRAC_30;
extern double G_CLO_VAC_2_FRAC_40;
extern double G_CLO_VAC_2_FRAC_50;
extern double G_CLO_VAC_2_FRAC_60;
extern double G_CLO_VAC_2_FRAC_70;
extern double G_CLO_VAC_2_FRAC_80;
extern double G_CLO_VAC_2_REL_EFF_00;
extern double G_CLO_VAC_2_REL_EFF_10;
extern double G_CLO_VAC_2_REL_EFF_20;
extern double G_CLO_VAC_2_REL_EFF_30;
extern double G_CLO_VAC_2_REL_EFF_40;
extern double G_CLO_VAC_2_REL_EFF_50;
extern double G_CLO_VAC_2_REL_EFF_60;
extern double G_CLO_VAC_2_REL_EFF_70;
extern double G_CLO_VAC_2_REL_EFF_80; */

extern double G_C_COMM[NUMAC][NUMAC];

// coefficients for lockdown contact matrix
extern double G_CLO_CONTACT_COEFF_00;
extern double G_CLO_CONTACT_COEFF_10;
extern double G_CLO_CONTACT_COEFF_20;
extern double G_CLO_CONTACT_COEFF_30;
extern double G_CLO_CONTACT_COEFF_40;
extern double G_CLO_CONTACT_COEFF_50;
extern double G_CLO_CONTACT_COEFF_60;
extern double G_CLO_CONTACT_COEFF_70;
extern double G_CLO_CONTACT_COEFF_80;

// coefficients for post-lockdown contact matrix
extern double G_CLO_CONTACT_COEFF_POSTLD_00;
extern double G_CLO_CONTACT_COEFF_POSTLD_10;
extern double G_CLO_CONTACT_COEFF_POSTLD_20;
extern double G_CLO_CONTACT_COEFF_POSTLD_30;
extern double G_CLO_CONTACT_COEFF_POSTLD_40;
extern double G_CLO_CONTACT_COEFF_POSTLD_50;
extern double G_CLO_CONTACT_COEFF_POSTLD_60;
extern double G_CLO_CONTACT_COEFF_POSTLD_70;
extern double G_CLO_CONTACT_COEFF_POSTLD_80;

// coefficients for POLYMOD contact matrix
extern double G_CLO_CONTACT_COEFF_NORM_00;
extern double G_CLO_CONTACT_COEFF_NORM_10;
extern double G_CLO_CONTACT_COEFF_NORM_20;
extern double G_CLO_CONTACT_COEFF_NORM_30;
extern double G_CLO_CONTACT_COEFF_NORM_40;
extern double G_CLO_CONTACT_COEFF_NORM_50;
extern double G_CLO_CONTACT_COEFF_NORM_60;
extern double G_CLO_CONTACT_COEFF_NORM_70;
extern double G_CLO_CONTACT_COEFF_NORM_80;

//  END  ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void PrintUsageModes();
void SetLocationData( string loc );

bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

void PrintUsageModes()
{
    printf("\n\tUSAGE: ./odesim  outfilename   [-options] \n\n");
}


//
// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=2;
    
    // this counter is here to check if we have set any of the G_CLO_MIXINGLEVEL_POSTLD parameters
    // this number has to be either 0 or 8; anything else leads to an error
    int num_postld_contact_rates_set = 0;

    if( argc<2 )
    { 
        PrintUsageModes(); 
	exit(-1);
    }
	
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to outfile\n\n");
    }
    else 
    {
        if( argv[1][0]=='-'){
            fprintf(stderr, "You probably meant to use '%s' as an option, not the output file.\n", argv[1]);
            exit(-1);
        }
        OutFile = fopen( argv[1], "w" );
    }

    string str;

    i=start;
    
    // read in options
    while(i<argc)
    {
	str = argv[i];
        
        // ### 1 ### IF BLOCK FOR BETA
        if( str == "-beta" )
        {
            ppc->v_betas.clear();
            ppc->v_betatimes.clear();
            i++;
            bool bFirstBetaPushedBack = false;
            
            //BEGIN LOOPING THROUGH THE LIST OF BETAS
            while(i<argc)
            {
                string s( argv[i] );
                if( isFloat(s) ) // if the current string is a floating point number, write it into the v_betas array
                {
                    // if the command line argument is <0, just set it back to zero
                    double d = atof( argv[i] );
                    if( d < 0.0 ) d = 0.0;
                    
                    ppc->v_betas.push_back( d );
                    
                    // if it's the first beta value being read in, then push it back twice, as it will be used as 
                    // the Jan 1 to March 1 beta (when there were almost no cases) as well as the beta for the 
                    // first time period we are invesigating
                    if(!bFirstBetaPushedBack) {ppc->v_betas.push_back( d );bFirstBetaPushedBack=true;}
                    
                    // increment and move on in this sub-loop
                    i++;
                }
                else
                {
                    // if the current string is NOT a float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }
                 
            } 
            //END OF LOOPING THROUGH THE LIST OF BETAS
             
            // MAKE SURE AT LEAST ONE BETA VALUE WAS READ IN 
            if( ppc->v_betas.size() == 0 )
            {
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-beta' option must be followed by at least one floating point number.\n\n");
            }
            else
            {
                double total_time_period = G_CLO_TF - 60.0;
                double time_step = total_time_period / ( (double) ppc->v_betas.size() - 1 );
                
                ppc->v_betatimes.push_back(0.0);
                
                for(int jj=1; jj<ppc->v_betas.size(); jj++)
                {
                    ppc->v_betatimes.push_back( 60.0 + time_step*((double)(jj-1)) );
                }
                
            }
                         
        }
        //
        // ### - ### END IF BLOCK FOR BETA
        //
        else if( str == "-binary-output" )  {           G_B_BINARY_OUTPUT = true;                   }   
        else if( str == "-checkpop" )
        {
            G_B_CHECKPOP_MODE = true;
        }
        /*else if( str == "-dev-hosp-young" )
        {
            G_CLO_HOSPFRAC_YOUNG_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_YOUNG_DEV > 3.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 3.0;
            if( G_CLO_HOSPFRAC_YOUNG_DEV < 0.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 0.0;
        }
        else if( str == "-dev-hosp-mid" )
        {
            G_CLO_HOSPFRAC_MID_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_MID_DEV > 3.0 ) G_CLO_HOSPFRAC_MID_DEV = 3.0;
            if( G_CLO_HOSPFRAC_MID_DEV < 0.0 ) G_CLO_HOSPFRAC_MID_DEV = 0.0;
        }
        else if( str == "-dev-hosp-old" )
        {
            G_CLO_HOSPFRAC_OLD_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_OLD_DEV > 1.3 ) G_CLO_HOSPFRAC_OLD_DEV = 1.3;
            if( G_CLO_HOSPFRAC_OLD_DEV < 0.0 ) G_CLO_HOSPFRAC_OLD_DEV = 0.0;
        }*/
        /* 
        else if( str == "-contact-rate-10" )        {     G_CLO_MIXINGLEVEL[1] = atof( argv[++i] );        }
        else if( str == "-contact-rate-20" )        {     G_CLO_MIXINGLEVEL[2] = atof( argv[++i] );        }
        else if( str == "-contact-rate-30" )        {     G_CLO_MIXINGLEVEL[3] = atof( argv[++i] );        }
        else if( str == "-contact-rate-40" )        {     G_CLO_MIXINGLEVEL[4] = atof( argv[++i] );        }
        else if( str == "-contact-rate-50" )        {     G_CLO_MIXINGLEVEL[5] = atof( argv[++i] );        }
        else if( str == "-contact-rate-60" )        {     G_CLO_MIXINGLEVEL[6] = atof( argv[++i] );        }
        else if( str == "-contact-rate-70" )        {     G_CLO_MIXINGLEVEL[7] = atof( argv[++i] );        }
        else if( str == "-contact-rate-80" )        {     G_CLO_MIXINGLEVEL[8] = atof( argv[++i] );        }
        else if( str == "-contact-rate-postld-10" ) {     G_CLO_MIXINGLEVEL_POSTLD[1] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-20" ) {     G_CLO_MIXINGLEVEL_POSTLD[2] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-30" ) {     G_CLO_MIXINGLEVEL_POSTLD[3] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-40" ) {     G_CLO_MIXINGLEVEL_POSTLD[4] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-50" ) {     G_CLO_MIXINGLEVEL_POSTLD[5] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-60" ) {     G_CLO_MIXINGLEVEL_POSTLD[6] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-70" ) {     G_CLO_MIXINGLEVEL_POSTLD[7] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-80" ) {     G_CLO_MIXINGLEVEL_POSTLD[8] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
         */
        else if( str == "-dev-icu-frac" )
        {
            G_CLO_ICUFRAC_DEV = atof( argv[++i] );
            if( G_CLO_ICUFRAC_DEV > 2.0 ) G_CLO_ICUFRAC_DEV = 2.0;
            if( G_CLO_ICUFRAC_DEV < 0.0 ) G_CLO_ICUFRAC_DEV = 0.0;
        }
        else if( str == "-dev-icu-frac-phase2" )
        {
            G_CLO_ICUFRAC_DEV_SECONDPHASE = atof( argv[++i] );
            if( G_CLO_ICUFRAC_DEV_SECONDPHASE > 2.0 ) G_CLO_ICUFRAC_DEV_SECONDPHASE = 2.0;
            if( G_CLO_ICUFRAC_DEV_SECONDPHASE < 0.0 ) G_CLO_ICUFRAC_DEV_SECONDPHASE = 0.0;
        }
        else if( str == "-dev-icu-frac-phase2beginday" )   {     G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY = atoi( argv[++i] );       }
        else if( str == "-dev-len-hospstay" )
        {
            G_CLO_DEV_LEN_HOSPSTAY = atof( argv[++i] );
            if( G_CLO_DEV_LEN_HOSPSTAY > 2.5 ) G_CLO_DEV_LEN_HOSPSTAY = 2.5;
            if( G_CLO_DEV_LEN_HOSPSTAY < 0.0 ) G_CLO_DEV_LEN_HOSPSTAY = 0.0;
        }
        else if( str == "-dev-ventdeath-mid" )
        {
            G_CLO_VENTDEATH_MID_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_MID_DEV > 1.7 ) G_CLO_VENTDEATH_MID_DEV = 1.7;
            if( G_CLO_VENTDEATH_MID_DEV < 0.0 ) G_CLO_VENTDEATH_MID_DEV = 0.0;
        }
        else if( str == "-dev-ventdeath-70" )
        {
            G_CLO_VENTDEATH_70_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_70_DEV > 1.4 ) G_CLO_VENTDEATH_70_DEV = 1.4;
            if( G_CLO_VENTDEATH_70_DEV < 0.8 ) G_CLO_VENTDEATH_70_DEV = 0.8;
        }
        else if( str == "-dev-ventdeath-80" )
        {
            G_CLO_VENTDEATH_80_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_80_DEV > 1.1 ) G_CLO_VENTDEATH_80_DEV = 1.1;
            if( G_CLO_VENTDEATH_80_DEV < 0.7 ) G_CLO_VENTDEATH_80_DEV = 0.7;
        }
        else if( str == "-death-prob-home-60" )     {     G_CLO_DEATHPROB_HOME_60 = atof( argv[++i] );          }
        else if( str == "-death-prob-home-70" )     {     G_CLO_DEATHPROB_HOME_70 = atof( argv[++i] );          }
        else if( str == "-death-prob-home-80" )     {     G_CLO_DEATHPROB_HOME_80 = atof( argv[++i] );          } 
        else if( str == "-death-prob-nonicu-80" )   {     G_CLO_PROB_NONICU_DEATH_80 = atof( argv[++i] );       }
        else if( str == "-death-prob-postvent" )    {     G_CLO_DEATHPROB_POSTVENT   = atof( argv[++i] );       }
        else if( str == "-death-rate-output" )      {     G_B_DEATHRATE_OUTPUT = true;                          }
        else if( str == "-diag" )                   {     G_B_DIAGNOSTIC_MODE = true;                           }
        else if( str == "-earlymarch-hosp-factor" ) {     G_CLO_EARLYMARCH_HOSPRATE = atof( argv[++i] );        }
        else if( str == "-firstlockdown-endday" )   {     G_CLO_FIRSTLOCKDOWN_ENDDAY = atoi( argv[++i] );       }
        else if( str == "-introday" )       {           G_CLO_INTRODUCTION_TIME = atof( argv[++i] );        }
        else if( str == "-loc" )            {           G_CLO_LOCATION = argv[++i];                         }
        else if( str == "-hosp-frac-10" )   {           G_CLO_HOSPFRAC_10 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-20" )   {           G_CLO_HOSPFRAC_20 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-30" )   {           G_CLO_HOSPFRAC_30 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-40" )   {           G_CLO_HOSPFRAC_40 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-50" )   {           G_CLO_HOSPFRAC_50 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-60" )   {           G_CLO_HOSPFRAC_60 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-70" )   {           G_CLO_HOSPFRAC_70 = atof( argv[++i] );              }
        else if( str == "-hosp-frac-80" )   {           G_CLO_HOSPFRAC_80 = atof( argv[++i] );              }
        else if( str == "-mean-time-vent" ) {           G_CLO_MEANTIME_ON_VENT_SURV = atof( argv[++i] );    }        
        else if( str == "-min-mixinglevel-00" ) {       G_CLO_MIN_MIXINGLEVEL_00 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-10" ) {       G_CLO_MIN_MIXINGLEVEL_10 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-20" ) {       G_CLO_MIN_MIXINGLEVEL_20 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-30" ) {       G_CLO_MIN_MIXINGLEVEL_30 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-40" ) {       G_CLO_MIN_MIXINGLEVEL_40 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-50" ) {       G_CLO_MIN_MIXINGLEVEL_50 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-60" ) {       G_CLO_MIN_MIXINGLEVEL_60 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-70" ) {       G_CLO_MIN_MIXINGLEVEL_70 = atof( argv[++i] );       }        
        else if( str == "-min-mixinglevel-80" ) {       G_CLO_MIN_MIXINGLEVEL_80 = atof( argv[++i] );       }
        else if( str == "-printIndices" )
        {
            printf("STARTS %d\n", 0 );          // start index of S
            printf("NUMAC %d\n", NUMAC );       // number age groups
            printf("NUME %d\n", NUME );         // number of Exposed-class
            printf("STARTE %d\n", STARTE );     // start index of E
            printf("NUMA %d\n", NUMA );         // Asymptomatic-class
            printf("STARTA %d\n", STARTA );
            printf("NUMI %d\n", NUMI );         // I (symptomatic)-class
            printf("STARTI %d\n", STARTI );
            printf("NUMHA %d\n", NUMHA );       // Hospotalized Acute-class
            printf("STARTHA %d\n", STARTHA );
            printf("STARTCA %d\n", STARTCA );   // Critical care (ICU) Acute-class
            printf("NUMV %d\n", NUMV );         // mechanical Ventilator-class
            printf("STARTV %d\n", STARTV );
            printf("STARTCR %d\n", STARTCR );   // Critical care (ICU) Recovering-class
            printf("STARTHR %d\n", STARTHR );   // Hospitalized Recovering-class
            printf("STARTD %d\n", STARTD );             // Died-at-home class (age-stratified)
            printf("STARTDHOSP %d\n", STARTDHOSP );     // Died-in-hospital class (age-stratified)
            printf("STARTR %d\n", STARTR );             // Recovered-at-home-class (age-stratified)
            printf("STARTRHOSP %d\n", STARTRHOSP );         // Recovered-in-hospital-and-discharged class (age-stratified)
            printf("STARTJ %d\n", STARTJ );   // Cumulative I-class (age-stratified)
            printf("STARTK %d\n", STARTK );   // Cumulative hospitalization incidence-class (age-stratified)
            
            printf("STARTJTEST %d\n", STARTJTEST );
            printf("STARTZ_1 %d\n", STARTZ_1 );
            printf("NUMZ_1 %d\n", NUMZ_1 );
            printf("STARTZ_2 %d\n", STARTZ_2 );
            printf("NUMZ_2 %d\n", NUMZ_2 );

            printf("STARTRR %d\n", STARTRR );
            printf("NUMRR %d\n", NUMRR );
            printf("STARTRRHOSP %d\n", STARTRRHOSP );
            printf("NUMRRHOSP %d\n", NUMRRHOSP );

            printf("STARTJZ_1 %d\n", STARTJZ_1 );
            printf("STARTJZ_2 %d\n", STARTJZ_2 );

            printf("DIMENSION %d\n", DIMENSION );

            exit(0);
        }
        else if( str == "-prob-icu-vent" )          {           G_CLO_PROB_ICU_TO_VENT     = atof( argv[++i] ); }
        
        else if( str == "-remove99" )               {           G_B_REMOVE_99COLUMNS = true;                    }   
        else if( str == "-version" ) 
        {
            printf("%d\n", 6);
            exit(0);
        }
        else if( str == "-rel-beta-hosp" )
        {
            G_CLO_RELATIVE_BETA_HOSP = atof( argv[++i] );
            if( G_CLO_RELATIVE_BETA_HOSP > 1.0 ) G_CLO_RELATIVE_BETA_HOSP = 1.0;
            if( G_CLO_RELATIVE_BETA_HOSP < 0.0 ) G_CLO_RELATIVE_BETA_HOSP = 0.0;
        }
        else if( str == "-scm" )
        {
            G_B_USESOCIALCONTACTMATRIX = true;
            //G_C_COMM[0][3] = 0.5;
        }
        else if( str == "-self-isolation-factor" )  {            G_CLO_SELFISOLATION_FACTOR = atof( argv[++i] );    }
        else if( str == "-steps-per-day" )          {            G_CLO_STEPS_PER_DAY = atof( argv[++i] );           }
        else if( str == "-susc-0-20" )
        {
            G_CLO_RELSUSC_0  = atof( argv[++i] );
            G_CLO_RELSUSC_10 = G_CLO_RELSUSC_0;
        }
        else if( str == "-susc-60-100" )
        {
            G_CLO_RELSUSC_60 = atof( argv[++i] );
            G_CLO_RELSUSC_70 = G_CLO_RELSUSC_60;
            G_CLO_RELSUSC_80 = G_CLO_RELSUSC_60;            
        }
        else if( str == "-symp-frac" ) // what you're really setting here is the symp fraction for 30-39 year-olds .. should be DEPRECATED
        {
            //G_CLO_SYMP_FRAC = atof( argv[++i] );
            //if( G_CLO_SYMP_FRAC > 0.325 ) G_CLO_SYMP_FRAC = 0.325;
            //if( G_CLO_SYMP_FRAC < 0.0 ) G_CLO_SYMP_FRAC = 0.0;
        }
        else if( str == "-symp-frac-davies" ) // use the Davies et al (Nature Medicine) numbers
        {
            G_B_SYMP_FRAC_DAVIES_20200616 = true;
            G_B_SYMP_FRAC_SIMPLEAVERAGE = false;
        }
        else if( str == "-symp-frac-equal" ) // here set the symp fraction to be equal for all age groups
        {
            G_B_SYMP_FRAC_DAVIES_20200616 = false;
            G_B_SYMP_FRAC_SIMPLEAVERAGE = false;
            G_CLO_SYMP_FRAC_EQUAL = atof( argv[++i] ); 
        }
        else if( str == "-tf" )                     {            G_CLO_TF = atof( argv[++i] );                      }
        else if( str == "-time-symp-to-hosp" )      {            G_CLO_TIME_SYMPTOHOSP = atof( argv[++i] );         }


        else if( str == "-normalcy-beginday" )   {     G_CLO_NORMALCY_BEGINDAY = atoi( argv[++i] );       }

        // natural immunity
        else if( str == "-len-nat-immunity" )       {        G_CLO_LEN_NAT_IMMU = atof( argv[++i] );      }

        // reporting rate, to be used in distributing vaccines to R-class (among previously asymptomatic cases)
        else if( str == "-rr" )                     {        G_CLO_REPORTING_RATE = atof( argv[++i] );      }
        
        /////////////////////////
        // vaccination
        /////////////////////////
        else if( str == "-no-test-before-vaccinate" ){        G_B_TEST_BEFORE_VACCINATE = false;                   }
        
        else if( str == "-vac1-foi")                 {        G_CLO_VAC_1_FOI = atof( argv[++i] );      }
        else if( str == "-vac1-protect-duration" )   {        G_CLO_VAC_1_PROTECT_DURATION = atof( argv[++i] );      }
        else if( str == "-vac1-efficacy-halflife" )  {        G_CLO_VAC_1_EFF_HALFLIFE = atof( argv[++i] );          }
        else if( str == "-vac1-efficacy-slope" )     {        G_CLO_VAC_1_EFF_SLOPE = atof( argv[++i] );             }

        else if( str == "-vac1-phase1-beginday" ){        
            // G_CLO_VAC_1_PHASE1_BEGINDAY = atoi( argv[++i] );       
            ppc->v_begin_days_phase1_vac1.clear();
            i++;

            // printf("\nread vac1-phase1-beginday\n");

            //BEGIN LOOPING THROUGH THE LIST OF vac1-phase1-beginday
            while(i<argc)
            {
                int d = -1;
                // if the current string is convertible to integer
                if( isdigit( argv[i][0] ) ){
                    d = atoi( argv[i] );
                    if( d < 61 ) { fprintf(stderr,"\n\tvac1-phase1-beginday must be > 60\n"); } // vac1-phase1-beginday must be > 60
                    ppc->v_begin_days_phase1_vac1.push_back( d ); // push to the end of v_begin_days_phase1_vac1 vector
                    i++; // increment and move on in this sub-loop

                    // printf("%d\t", d);
                }
                else {

                    // printf("\n");

                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--; 
                    break;
                }   
            } 
            //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->v_begin_days_phase1_vac1.size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-phase1-beginday' option must be followed by at least one integer number.\n\n");
            }
            
            ppc->index_current_day_phase1_vac1 = 0;
        }

        else if( str == "-vac1-phase1-endday" ){        
            // G_CLO_VAC_1_PHASE1_ENDDAY = atoi( argv[++i] );         
            ppc->v_end_days_phase1_vac1.clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac1-phase1-endday
            while(i<argc)
            {
                int d = -1;
                // if the current string is convertible to integer
                if( isdigit( argv[i][0] ) ){
                    d = atoi( argv[i] );
                    if( d < 62 ) { fprintf(stderr,"\n\tvac1-phase1-endday must be > 61\n"); } // vac1-phase1-endday must be > 61
                    ppc->v_end_days_phase1_vac1.push_back( d ); // push to the end of v_end_days_phase1_vac1 vector
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } 
            //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac1-phase1-endday VALUE WAS READ IN 
            if( ppc->v_end_days_phase1_vac1.size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-phase1-endday' option must be followed by at least one integer number.\n\n");
            }
        }

        else if( str == "-vac1-phase1-dpd" ){        
            // G_CLO_VAC_1_PHASE1_DPD = atoi( argv[++i] );            
            ppc->v_dpd_phase1_vac1.clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac1-phase1-dpd
            while(i<argc)
            {
                int d = -1;
                // if the current string is convertible to integer
                if( isdigit( argv[i][0] ) ){
                    d = atoi( argv[i] );
                    if( d < 0 ) { fprintf(stderr,"\n\tvac1-phase1-dpd must be non-negative\n"); } // vac1-phase1-dpd must be >= 0 
                    ppc->v_dpd_phase1_vac1.push_back( d ); // push to the end of v_dpd_phase1_vac1 vector
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } 
            //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac1-phase1-dpd VALUE WAS READ IN 
            if( ppc->v_dpd_phase1_vac1.size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-phase1-dpd' option must be followed by at least one integer number.\n\n");
            }
        }
        // else if( str == "-vac1-ratio-00" )   {           G_CLO_VAC_1_FRAC_00 = atof( argv[++i] );              }
        else if( str == "-vac1-ratio-10" )   {           
            // G_CLO_VAC_1_FRAC_10 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[1].clear();
            i++;

            // printf("read -vac1-ratio-10\n");
            
            //BEGIN LOOPING THROUGH THE LIST OF vac1-ratio-10
            while(i<argc)
            {
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-10 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-10 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[1].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[1] vector
                    i++; // increment and move on in this sub-loop

                    // printf("%.8f\t", vf);
                }
                else {
                    // printf("\n");

                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } 
            //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac1-ratio-10 VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[1].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-10' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-ratio-20" )   {           
            // G_CLO_VAC_1_FRAC_20 = atof( argv[++i] ); 
            ppc->v_vac_ratios_phase1_vac1[2].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-20
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-20 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-20 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[2].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[2] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[2].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-20' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-ratio-30" )   {           
            // G_CLO_VAC_1_FRAC_30 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[3].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-30
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-30 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-30 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[3].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[3] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[3].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-30' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-ratio-40" )   {           
            // G_CLO_VAC_1_FRAC_40 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[4].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-40
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-40 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-40 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[4].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[4] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[4].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-40' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-ratio-50" )   {           
            // G_CLO_VAC_1_FRAC_50 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[5].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-50
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-50 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-50 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[5].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[5] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[5].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-50' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-ratio-60" )   {           
            // G_CLO_VAC_1_FRAC_60 = atof( argv[++i] );             
            ppc->v_vac_ratios_phase1_vac1[6].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-60 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-60 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[6].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[6] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[6].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-60' option must be followed by at least one integer number.\n\n");
            } 
        }
        else if( str == "-vac1-ratio-70" )   {           
            // G_CLO_VAC_1_FRAC_70 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[7].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-70 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-70 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[7].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[7] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[7].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-70' option must be followed by at least one integer number.\n\n");
            } 
        }
        else if( str == "-vac1-ratio-80" )   {           
            // G_CLO_VAC_1_FRAC_80 = atof( argv[++i] );              
            ppc->v_vac_ratios_phase1_vac1[8].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-ratio-80 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-ratio-80 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_ratios_phase1_vac1[8].push_back( vf ); // push to the end of v_vac_ratios_phase1_vac1[8] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_ratios_phase1_vac1[8].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-ratio-80' option must be followed by at least one integer number.\n\n");
            } 
        }






        /////////////////////
        ///////////////////// real vac frac (to be multiplied by dose-per-day)
        /////////////////////
        else if( str == "-vac1-frac-10" )   {           
            // G_CLO_VAC_1_FRAC_10 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[1].clear();
            i++;

            // printf("read -vac1-frac-10\n");
            
            //BEGIN LOOPING THROUGH THE LIST OF vac1-frac-10
            while(i<argc)
            {
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-10 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-10 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[1].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[1] vector
                    i++; // increment and move on in this sub-loop

                    // printf("%.8f\t", vf);
                }
                else {
                    // printf("\n");

                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } 
            //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac1-frac-10 VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[1].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-10' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-frac-20" )   {           
            // G_CLO_VAC_1_FRAC_20 = atof( argv[++i] ); 
            ppc->v_vac_fracs_phase1_vac1[2].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-20
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-20 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-20 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[2].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[2] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[2].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-20' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-frac-30" )   {           
            // G_CLO_VAC_1_FRAC_30 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[3].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-30
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-30 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-30 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[3].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[3] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[3].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-30' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-frac-40" )   {           
            // G_CLO_VAC_1_FRAC_40 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[4].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-40
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-40 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-40 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[4].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[4] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[4].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-40' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-frac-50" )   {           
            // G_CLO_VAC_1_FRAC_50 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[5].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-50
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-50 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-50 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[5].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[5] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[5].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-50' option must be followed by at least one integer number.\n\n");
            }
        }
        else if( str == "-vac1-frac-60" )   {           
            // G_CLO_VAC_1_FRAC_60 = atof( argv[++i] );             
            ppc->v_vac_fracs_phase1_vac1[6].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-60 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-60 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[6].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[6] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[6].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-60' option must be followed by at least one integer number.\n\n");
            } 
        }
        else if( str == "-vac1-frac-70" )   {           
            // G_CLO_VAC_1_FRAC_70 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[7].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-70 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-70 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[7].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[7] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[7].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-70' option must be followed by at least one integer number.\n\n");
            } 
        }
        else if( str == "-vac1-frac-80" )   {           
            // G_CLO_VAC_1_FRAC_80 = atof( argv[++i] );              
            ppc->v_vac_fracs_phase1_vac1[8].clear();
            i++;
            
            //BEGIN LOOPING THROUGH THE LIST OF vac-frac-60
            while(i<argc){
                double vf = -1.0;
                // if the current string is convertible to float
                if( isdigit( argv[i][0] ) ){
                    vf = atof( argv[i] );
                    if( vf < 0.0 ) { fprintf(stderr,"\n\tvac1-frac-80 must be non-negative\n"); } // vac frac must be >= 0 
                    if( vf > 1.0 ) { fprintf(stderr,"\n\tvac1-frac-80 must not exceed 1.\n"); } // vac frac must be <= 1 
                    ppc->v_vac_fracs_phase1_vac1[8].push_back( vf ); // push to the end of v_vac_fracs_phase1_vac1[8] vector
                    i++; // increment and move on in this sub-loop
                } else {
                    // if the current string is NOT a convertible to float, reduce the counter to parse the next option and break out of the loop
                    i--;    break;
                }   
            }  //END OF WHILE LOOP
             
            // MAKE SURE AT LEAST ONE vac frac VALUE WAS READ IN 
            if( ppc->v_vac_fracs_phase1_vac1[8].size() == 0 ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-vac1-frac-80' option must be followed by at least one integer number.\n\n");
            } 
        }





        else if( str == "-vac1-rel-eff-00" )   {          G_CLO_VAC_1_REL_EFF_00 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-10" )   {          G_CLO_VAC_1_REL_EFF_10 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-20" )   {          G_CLO_VAC_1_REL_EFF_20 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-30" )   {          G_CLO_VAC_1_REL_EFF_30 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-40" )   {          G_CLO_VAC_1_REL_EFF_40 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-50" )   {          G_CLO_VAC_1_REL_EFF_50 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-60" )   {          G_CLO_VAC_1_REL_EFF_60 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-70" )   {          G_CLO_VAC_1_REL_EFF_70 = atof( argv[++i] );              }
        else if( str == "-vac1-rel-eff-80" )   {          G_CLO_VAC_1_REL_EFF_80 = atof( argv[++i] );              }

        else if( str == "-vac1-phase2-beginday" )    {        G_CLO_VAC_1_PHASE2_BEGINDAY = atoi( argv[++i] );       }
        else if( str == "-vac1-phase2-endday" )      {        G_CLO_VAC_1_PHASE2_ENDDAY = atoi( argv[++i] );         }


        /* else if( str == "-vac2-foi")                 {        G_CLO_VAC_2_FOI = atof( argv[++i] );      }
        else if( str == "-vac2-protect-duration" )   {        G_CLO_VAC_2_PROTECT_DURATION = atof( argv[++i] );      }
        else if( str == "-vac2-efficacy-halflife" )  {        G_CLO_VAC_2_EFF_HALFLIFE = atof( argv[++i] );          }
        else if( str == "-vac2-efficacy-slope" )     {        G_CLO_VAC_2_EFF_SLOPE = atof( argv[++i] );             }
        else if( str == "-vac2-phase1-beginday" )    {        G_CLO_VAC_2_PHASE1_BEGINDAY = atoi( argv[++i] );       }
        else if( str == "-vac2-phase1-endday" )      {        G_CLO_VAC_2_PHASE1_ENDDAY = atoi( argv[++i] );         }
        else if( str == "-vac2-phase1-dpd" )         {        G_CLO_VAC_2_PHASE1_DPD = atoi( argv[++i] );            }
        else if( str == "-vac2-phase2-beginday" )    {        G_CLO_VAC_2_PHASE2_BEGINDAY = atoi( argv[++i] );       }
        else if( str == "-vac2-phase2-endday" )      {        G_CLO_VAC_2_PHASE2_ENDDAY = atoi( argv[++i] );         }
        // else if( str == "-vac2-frac-00" )   {           G_CLO_VAC_2_FRAC_00 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-10" )   {           G_CLO_VAC_2_FRAC_10 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-20" )   {           G_CLO_VAC_2_FRAC_20 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-30" )   {           G_CLO_VAC_2_FRAC_30 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-40" )   {           G_CLO_VAC_2_FRAC_40 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-50" )   {           G_CLO_VAC_2_FRAC_50 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-60" )   {           G_CLO_VAC_2_FRAC_60 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-70" )   {           G_CLO_VAC_2_FRAC_70 = atof( argv[++i] );              }
        else if( str == "-vac2-frac-80" )   {           G_CLO_VAC_2_FRAC_80 = atof( argv[++i] );              }
        // relative efficacies
        else if( str == "-vac2-rel-eff-00" )   {          G_CLO_VAC_2_REL_EFF_00 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-10" )   {          G_CLO_VAC_2_REL_EFF_10 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-20" )   {          G_CLO_VAC_2_REL_EFF_20 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-30" )   {          G_CLO_VAC_2_REL_EFF_30 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-40" )   {          G_CLO_VAC_2_REL_EFF_40 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-50" )   {          G_CLO_VAC_2_REL_EFF_50 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-60" )   {          G_CLO_VAC_2_REL_EFF_60 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-70" )   {          G_CLO_VAC_2_REL_EFF_70 = atof( argv[++i] );              }
        else if( str == "-vac2-rel-eff-80" )   {          G_CLO_VAC_2_REL_EFF_80 = atof( argv[++i] );              } */

        // coefficient for lockdown contact matrix
        else if( str == "-contact-coeff-00" )   {           G_CLO_CONTACT_COEFF_00 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-10" )   {           G_CLO_CONTACT_COEFF_10 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-20" )   {           G_CLO_CONTACT_COEFF_20 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-30" )   {           G_CLO_CONTACT_COEFF_30 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-40" )   {           G_CLO_CONTACT_COEFF_40 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-50" )   {           G_CLO_CONTACT_COEFF_50 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-60" )   {           G_CLO_CONTACT_COEFF_60 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-70" )   {           G_CLO_CONTACT_COEFF_70 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-80" )   {           G_CLO_CONTACT_COEFF_80 = atof( argv[++i] );              }

        // coefficient for post-lockdown contact matrix
        else if( str == "-contact-coeff-postld-00" )   {    G_CLO_CONTACT_COEFF_POSTLD_00 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-10" )   {    G_CLO_CONTACT_COEFF_POSTLD_10 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-20" )   {    G_CLO_CONTACT_COEFF_POSTLD_20 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-30" )   {    G_CLO_CONTACT_COEFF_POSTLD_30 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-40" )   {    G_CLO_CONTACT_COEFF_POSTLD_40 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-50" )   {    G_CLO_CONTACT_COEFF_POSTLD_50 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-60" )   {    G_CLO_CONTACT_COEFF_POSTLD_60 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-70" )   {    G_CLO_CONTACT_COEFF_POSTLD_70 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }
        else if( str == "-contact-coeff-postld-80" )   {    G_CLO_CONTACT_COEFF_POSTLD_80 = atof( argv[++i] );  G_CLO_POSTLD_MIXING_SET=true;            }

        // coefficient for POLYMOD contact matrix
        else if( str == "-contact-coeff-norm-00" )   {           G_CLO_CONTACT_COEFF_NORM_00 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-10" )   {           G_CLO_CONTACT_COEFF_NORM_10 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-20" )   {           G_CLO_CONTACT_COEFF_NORM_20 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-30" )   {           G_CLO_CONTACT_COEFF_NORM_30 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-40" )   {           G_CLO_CONTACT_COEFF_NORM_40 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-50" )   {           G_CLO_CONTACT_COEFF_NORM_50 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-60" )   {           G_CLO_CONTACT_COEFF_NORM_60 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-70" )   {           G_CLO_CONTACT_COEFF_NORM_70 = atof( argv[++i] );              }
        else if( str == "-contact-coeff-norm-80" )   {           G_CLO_CONTACT_COEFF_NORM_80 = atof( argv[++i] );              }

        // ### FINAL ### IF BLOCK FOR AN UNKNOWN COMMAND-LINE OPTIONS            
        else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
        
        
        
        //END OF MAIN WHILE-LOOP BLOCK; INCREMENT AND MOVE ON
        i++;
    }
    
    // after parsing all of the arguments, do the following
    
    // ### 1 ### ERROR CHECK
    if( num_postld_contact_rates_set==0 || num_postld_contact_rates_set==8 )
    {
        // do nothing
    }
    else
    {
        fprintf(stderr, "\n\tIf you are setting the -contact-rate-postld-nn parameters, you need to set all eight of them (%d were set).\n\n", num_postld_contact_rates_set);
        exit(-1);
    }
    
    // ### 2 ### set the postld-contact data if it needs to be set
    if( !G_CLO_POSTLD_MIXING_SET )
    {
        // if the post-ld mixing paramters were not set, just set them to the regular (pre-ld) mixing params   
        // for(int ac=1; ac<NUMAC; ac++){
        //     G_CLO_MIXINGLEVEL_POSTLD[ac] = G_CLO_MIXINGLEVEL[ac];
        // }
        G_CLO_CONTACT_COEFF_POSTLD_00 = G_CLO_CONTACT_COEFF_00;
        G_CLO_CONTACT_COEFF_POSTLD_10 = G_CLO_CONTACT_COEFF_10;
        G_CLO_CONTACT_COEFF_POSTLD_20 = G_CLO_CONTACT_COEFF_20;
        G_CLO_CONTACT_COEFF_POSTLD_30 = G_CLO_CONTACT_COEFF_30;
        G_CLO_CONTACT_COEFF_POSTLD_40 = G_CLO_CONTACT_COEFF_40;
        G_CLO_CONTACT_COEFF_POSTLD_50 = G_CLO_CONTACT_COEFF_50;
        G_CLO_CONTACT_COEFF_POSTLD_60 = G_CLO_CONTACT_COEFF_60;
        G_CLO_CONTACT_COEFF_POSTLD_70 = G_CLO_CONTACT_COEFF_70;
        G_CLO_CONTACT_COEFF_POSTLD_80 = G_CLO_CONTACT_COEFF_80;
    }
    
    // ### 3 ### set the specific age structure, population, and introduction time for this location
    SetLocationData( G_CLO_LOCATION );


    // ### 4 ### check vaccination begin and end day 
    /* if ( G_CLO_VAC_1_PHASE1_BEGINDAY <= G_CLO_VAC_1_PHASE1_ENDDAY && 
         G_CLO_VAC_1_PHASE1_ENDDAY <= G_CLO_VAC_1_PHASE2_BEGINDAY && 
         G_CLO_VAC_1_PHASE2_BEGINDAY <= G_CLO_VAC_1_PHASE2_ENDDAY){
            // passed, do nothing
    } else
    {
        fprintf(stderr, "\n\tCheck begin and end days for the two phases of deployment of low-efficacy vaccine.\n\n");
        exit(-1);
    } */

    /* if ( G_CLO_VAC_2_PHASE1_BEGINDAY <= G_CLO_VAC_2_PHASE1_ENDDAY && 
         G_CLO_VAC_2_PHASE1_ENDDAY <= G_CLO_VAC_2_PHASE2_BEGINDAY && 
         G_CLO_VAC_2_PHASE2_BEGINDAY <= G_CLO_VAC_2_PHASE2_ENDDAY){
            // passed, do nothing
    } else
    {
        fprintf(stderr, "\n\tCheck begin and end days for the two phases of deployment of high-efficacy vaccine.\n\n");
        exit(-1);
    } */
    if ( G_CLO_REPORTING_RATE > 0.0 && G_CLO_REPORTING_RATE <= 1.0  ){
             // do nothing
    } else
    {
        fprintf(stderr, "\n\tCheck reporting rate. It must be within (0.0, 1.0].\n\n");
        exit(-1);
    }

    // ### 5a ### check equal length v_begin_days_phase1_vac1, v_end_days_phase1_vac1, v_dpd_phase1_vac1
    if ( ppc->v_begin_days_phase1_vac1.size() == ppc->v_end_days_phase1_vac1.size() && 
         ppc->v_begin_days_phase1_vac1.size() == ppc->v_dpd_phase1_vac1.size() ){
             // do nothing
    } else
    {
        fprintf(stderr, "\n\tCheck begin days, end days, and dose-per-day inputs for the vaccine rollout campaign. The length of each input type must be identical.\n\n");
        exit(-1);
    }
    // ### 5b ### check at least 1 set of vac-ratios or vac-fracs if a v_begin_days_phase1_vac1 was set
    bool b_vacratio_equal_length = true;
    bool b_vacfrac_equal_length = true;
    for (size_t ac = 2; ac < NUMAC; ac++){
        b_vacratio_equal_length = b_vacratio_equal_length && (ppc->v_vac_ratios_phase1_vac1[1].size() == ppc->v_vac_ratios_phase1_vac1[ac].size() );
    }
    for (size_t ac = 2; ac < NUMAC; ac++){
        b_vacfrac_equal_length = b_vacfrac_equal_length && (ppc->v_vac_fracs_phase1_vac1[1].size() == ppc->v_vac_fracs_phase1_vac1[ac].size() );
    }
    if ( ppc->v_begin_days_phase1_vac1.size() > 0 ){
        if( ppc->v_vac_ratios_phase1_vac1[1].size() == 0 && ppc->v_vac_fracs_phase1_vac1[1].size() == 0){
            fprintf(stderr, "\n\tCheck fractions and/or ratios of the population to be vaccinated in the vaccine rollout campaign.\nThe campaign needs at least one set of vaccine fractions or ratios.\n\n");
            exit(-1);
        } 
        if ( ppc->v_vac_ratios_phase1_vac1[1].size() > 0 && b_vacratio_equal_length     ) {
            G_B_USE_VAC_RATIO = true;
        } 
        if( ppc->v_vac_fracs_phase1_vac1[1].size() > 0 && b_vacfrac_equal_length     ) {
            G_B_USE_VAC_FRAC = true;
        } 
        if ((!G_B_USE_VAC_RATIO) && (!G_B_USE_VAC_FRAC)) {
            fprintf(stderr, "\n\tCheck fractions and/or ratios of the population to be vaccinated in the vaccine rollout campaign.\nAll age groups must have identical length of vaccine fractions and/or ratios.\n");
            exit(-1);
        }
        // printf("\nG_B_USE_VAC_RATIO %d G_B_USE_VAC_FRAC %d\n", G_B_USE_VAC_RATIO, G_B_USE_VAC_FRAC);
    }
       

    // ### 5 ### check the relationship between vaccine duration of protection and efficacy half-life ???
    // TODO: discuss with Maciek
    // if( G_CLO_VAC_1_PROTECT_DURATION > G_CLO_VAC_1_EFF_HALFLIFE * 2.0){
    //     G_CLO_VAC_1_EFF_HALFLIFE = G_CLO_VAC_1_PROTECT_DURATION / 2.0;
    // } else {
    //     G_CLO_VAC_1_PROTECT_DURATION = G_CLO_VAC_1_EFF_HALFLIFE * 2.0;
    // }
    // if( G_CLO_VAC_2_PROTECT_DURATION > G_CLO_VAC_2_EFF_HALFLIFE * 2.0){
    //     G_CLO_VAC_2_EFF_HALFLIFE = G_CLO_VAC_2_PROTECT_DURATION / 2.0;
    // } else {
    //     G_CLO_VAC_2_PROTECT_DURATION = G_CLO_VAC_2_EFF_HALFLIFE * 2.0;
    // }
    


    return;
}

void SetLocationData( string loc )
{
    for(int i=0;i<DIMENSION;i++) yic[i]=0.0; // zero everything out
    
    if( loc=="RI" )
    {
        // 2019 est pop is 1,059,361
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 1059361.0;
        yic[0]  = ppc->v[i_N]   *   .105;   //  0-9
        yic[1]  = ppc->v[i_N]   *   .123;   //  10-19
        yic[2]  = ppc->v[i_N]   *   .140;   //  20-29
        yic[3]  = ppc->v[i_N]   *   .127;   //  30-39
        yic[4]  = ppc->v[i_N]   *   .124;   //  40-49
        yic[5]  = ppc->v[i_N]   *   .135;   //  50-59
        yic[6]  = ppc->v[i_N]   *   .120;   //  60-69
        yic[7]  = ppc->v[i_N]   *   .074;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        ppc->v_pop_frac[0] = 0.105; ppc->v_pop_frac[1] = 0.123; ppc->v_pop_frac[2] = 0.140;
        ppc->v_pop_frac[3] = 0.127; ppc->v_pop_frac[4] = 0.124; ppc->v_pop_frac[5] = 0.135;
        ppc->v_pop_frac[6] = 0.120; ppc->v_pop_frac[7] = 0.074; ppc->v_pop_frac[8] = 0.052;

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;
        //G_CLO_EARLYMARCH_ENDDAY = 84.5; 
        
    }
    else if( loc=="MA" )
    {
        // 2018 est pop is 6,897,212
        

        // data from https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-total.html
        ppc->v[i_N] = 6897212.0;
        yic[0]  = ppc->v[i_N]   *   0.10586466;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12243686;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.14498786;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.13384234;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.12230812;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.14064248;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.11801015;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.06958116;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        ppc->v_pop_frac[0] = 0.10586466; ppc->v_pop_frac[1] = 0.12243686; ppc->v_pop_frac[2] = 0.14498786;
        ppc->v_pop_frac[3] = 0.13384234; ppc->v_pop_frac[4] = 0.12230812; ppc->v_pop_frac[5] = 0.14064248;
        ppc->v_pop_frac[6] = 0.11801015; ppc->v_pop_frac[7] = 0.06958116; ppc->v_pop_frac[8] = 0.04232637;

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0; // 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;        
        
    }
    else if( loc=="PA" )
    {
        // 2019 est pop is 12,800,721
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 12800721.0;
        yic[0]  = ppc->v[i_N]   *   0.11160395;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12229803;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.13156525;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.12581869;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.11809624;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.13878546;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.1270166;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.07657303;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        ppc->v_pop_frac[0] = 0.11160395; ppc->v_pop_frac[1] = 0.12229803; ppc->v_pop_frac[2] = 0.13156525;
        ppc->v_pop_frac[3] = 0.12581869; ppc->v_pop_frac[4] = 0.11809624; ppc->v_pop_frac[5] = 0.13878546;
        ppc->v_pop_frac[6] = 0.12701660; ppc->v_pop_frac[7] = 0.07657303; ppc->v_pop_frac[8] = 0.04824275;

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 62.0;
        
        G_CLO_INTRODUCTION_COUNT = 2;        
        
    }
    else if( loc=="FL" )
    {
        // 2019 est pop is 21.5M
        

        // numbers from https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-detail.html
        ppc->v[i_N] = 21477737.0;
        yic[0]  = ppc->v[i_N]   *   0.107371275;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.112138583;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.124684365;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.126472682;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.121251461;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.133090744;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.125931750;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.094975183;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 2;        
        
    }
    else
    {
        fprintf(stderr, "\n\tUnknown location [%s] entered after -loc on command line.\n\n", loc.c_str() );
        exit(-1);
    }
    
}





