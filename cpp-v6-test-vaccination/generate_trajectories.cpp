#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <float.h>
//#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "essentials.h"
#include "prms.h"
#include "assert.h"


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern FILE* OutFile;
extern double* yic;
extern prms* ppc;

extern double G_C_COMM[NUMAC][NUMAC];

extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;

extern int G_CLO_FIRSTLOCKDOWN_ENDDAY;
extern bool G_CLO_POSTLD_MIXING_SET;

extern bool G_B_DIAGNOSTIC_MODE;
extern bool G_B_CHECKPOP_MODE;
extern bool G_B_REMOVE_99COLUMNS;
extern bool G_B_BINARY_OUTPUT;

extern double G_CLO_ICUFRAC_DEV;
extern double G_CLO_ICUFRAC_DEV_SECONDPHASE;
extern int G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY;

extern int G_CLO_STEPS_PER_DAY;

extern int G_CLO_NORMALCY_BEGINDAY;

extern bool G_B_TEST_BEFORE_VACCINATE;
extern bool G_B_USE_VAC_FRAC;
extern bool G_B_USE_VAC_RATIO;

// vaccination
extern double G_CLO_VAC_1_PHASE1_BEGINDAY;
extern double G_CLO_VAC_1_PHASE1_ENDDAY;
extern double G_CLO_VAC_1_PHASE2_BEGINDAY;
extern double G_CLO_VAC_1_PHASE2_ENDDAY;

// extern double G_CLO_VAC_2_PROTECT_DURATION;
extern double G_CLO_VAC_2_PHASE1_BEGINDAY;
extern double G_CLO_VAC_2_PHASE1_ENDDAY;
extern double G_CLO_VAC_2_PHASE2_BEGINDAY;
extern double G_CLO_VAC_2_PHASE2_ENDDAY;

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

void write(double tt, double * yic, size_t dim)
{
    if( G_B_BINARY_OUTPUT){
        write_bin(tt, yic, dim);
        return;
    }
    int startcol = 0;
    if( G_B_REMOVE_99COLUMNS ){ startcol=99; }

    if( OutFile == NULL )
    {
        printf("%1.3f", tt);
        for(int i=startcol; i<dim; i++) printf("\t%1.4f", yic[i]);
        //for(i=0;i<NUMAC+3;i++) printf("\t%1.4f", yic[i]);
        printf("\n");
    }
    else
    {
        fprintf(OutFile, "%1.3f", tt);
        for(int i=startcol; i<dim; i++) fprintf(OutFile, "\t%1.4f", yic[i]);
        //for(i=0;i<NUMAC+3;i++) printf("\t%1.4f", yic[i]);
        fprintf(OutFile, "\n");
    }
}

void write_bin(double tt, double * yic, size_t dim)
{
    int startcol = 0;
    if( G_B_REMOVE_99COLUMNS ){ startcol=99; }

    if( OutFile == NULL )
    {
        fwrite(&tt, sizeof(tt), 1, stdout);
        fwrite(yic+startcol, sizeof(double), dim-startcol, stdout);
    }
    else
    {
        fwrite(&tt, sizeof(tt), 1, OutFile);
        fwrite(yic+startcol, sizeof(double), dim-startcol, OutFile);
    }
}

// [[Rcpp::export]]
void generate_trajectories( double inital_number_of_cases, double param_beta, double t0, double tf, double h )
{
    assert(ppc);
    double tt,rkstep,ttstop,eps,h1,hmin,ttbeforestop;
    int nvar,nok,nbad; //rkqs();
    int i, j;
    int dim = DIMENSION; // STARTK+NUMAC; //STARTZ_2 + NUMAC*NUMZ_2; // this is the dimension of the dynamical system, i.e. the number of equations (August 27 2020: this should be 351)
    
    int startcol=0;         // this is the starting column of the output; since the first 99 columns (indices 0 to 98) are the S, E, A classes,
                            // these may not need to be output since they don't have likelihoods associated with them
                            
    if( G_B_REMOVE_99COLUMNS ) startcol=99;
    
    double totalpop=0.0;
    if( G_B_CHECKPOP_MODE )
    {
        for(int k=0; k<STARTJ; k++) totalpop += yic[k];
        for(int k=STARTZ_1; k<DIMENSION; k++) totalpop += yic[k];
        printf("\n\t totalpop=%1.16f", totalpop);
    }
    
    
    // NOTE 2020-04-04 : this is still very fast and it prevents an off-by-one-day pseudo-error when the step size is close to 0.5 or 1.0
    //                 : just to be clear, a step-size of 1.0 is perfectly fine, but you have to mentally correct for the fact that the 
    //                 : diff-eqs will be about one day late
    int steps_per_day = G_CLO_STEPS_PER_DAY;
    
    // set some values in the RKF integrator
    tt = (double)((int)t0); // chop off any decimals, and just start at an integer time
    ttstop = (double)((int)tf); // chop off any decimals, and just start at an integer time
    rkstep = 1.0 / ((double)steps_per_day);
    eps=0.00000001;
    h1 = 0.1;
    hmin = 1e-13;
    //hmin =    0.00000000001;
    nvar = dim;
    
    int counter = 1;
    bool bIntroductionOccurred = false;
    bool bContactRatesUpdatedAfterLockdown = false;
    bool bPhase2BetterClinicalManagementBegun = false;
    double NextTimeForBetaUpdate=1000000.0; assert( tf < 999999.0 );
    
    if( ppc->v_betatimes.size() > 1 )
    {
        NextTimeForBetaUpdate = ppc->v_betatimes[1];
    }
    ppc->assign_new_beta(); // called here for the first time, so it just puts the first beta into use
    
    // before integration begins, assign the initial beta value

    
    bool b_normalcy_contact_rates_updated = false;

    // vaccination housekeeping
    bool b_vac1_phase1_began = false;
    double next_beginday_vac1_phase1 = (ppc->v_begin_days_phase1_vac1.size() > 0) ? ((double)ppc->v_begin_days_phase1_vac1[0]) : tf+1;
    double next_endday_vac1_phase1 = (ppc->v_end_days_phase1_vac1.size() > 0) ? ((double)ppc->v_end_days_phase1_vac1[0]) : tf+1;
    bool b_vac1_phase2_began = false;
    double min_vac_frac_1 = 1.0;
    double sum_weighted_vac_frac_1 = 0.0;
    double each_step_total_S_Z_1_phase1 = 0.0; // sum of v_number_to_Z_1_phase1 for each integration step
    // bool b_done_vac_frac_1 = true;
    
    
    // bool b_vac2_phase1_began = false;
    // double next_update_vac2_phase1 = G_CLO_VAC_2_PHASE1_BEGINDAY;
    // bool b_vac_high_phase2_began = false;
    // double min_vac_frac_2 = 1.0;
    // double sum_weighted_vac_frac_2 = 0.0;


    
    //
    //BEGIN OF LOOP INTEGRATING ODEs
    //
    while( tt < (ttstop-0.000000001) )      // MFB NOTE 2020-06-09: subtract this tiny amount since sometimes we shouldn't enter the while loop
    {                                       // but we do anyway because tt is just barely below ttstop by some tiny fraction
        // introduce the first infections
        if( !bIntroductionOccurred && tt > G_CLO_INTRODUCTION_TIME )
        {
            yic[STARTI+4]=((double)G_CLO_INTRODUCTION_COUNT); // some number of infected individuals in their forties are introduced 
            bIntroductionOccurred = true;
        }
        
        // post-lockdown contact rates
        if( !bContactRatesUpdatedAfterLockdown && tt > G_CLO_FIRSTLOCKDOWN_ENDDAY )
        {
            /* for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
            {
                for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
                {
                    G_C_COMM[acs][aci] = ppc->v_mixing_level_postld[acs] * ppc->v_mixing_level_postld[aci];
                }
            }     */

            // CoMix - Belgium wave 5 (late June 2020) matrix from `socialmixr` R package;
            G_C_COMM[0][0] = 0.2;
            G_C_COMM[0][1] = 0.136363636363636;     G_C_COMM[1][0] = 0.136363636363636; 
            G_C_COMM[0][2] = 0.165178571428571;     G_C_COMM[2][0] = 0.165178571428571;
            G_C_COMM[0][3] = 1.11034482758621;      G_C_COMM[3][0] = 1.11034482758621;
            G_C_COMM[0][4] = 0.313218390804598;     G_C_COMM[4][0] = 0.313218390804598;
            G_C_COMM[0][5] = 0.072769953051643;     G_C_COMM[5][0] = 0.072769953051643;
            G_C_COMM[0][6] = 0.135514018691589;     G_C_COMM[6][0] = 0.135514018691589;
            G_C_COMM[0][7] = 0.048076923076923;     G_C_COMM[7][0] = 0.048076923076923;
            G_C_COMM[0][8] = 0.001;                 G_C_COMM[8][0] = 0.001;

            G_C_COMM[1][1] = 1.13636363636364;
            G_C_COMM[1][2] = 0.001;                 G_C_COMM[2][1] = 0.59375;
            G_C_COMM[1][3] = 1.13636363636364;      G_C_COMM[3][1] = 0.293103448275862;
            G_C_COMM[1][4] = 0.863636363636364;     G_C_COMM[4][1] = 0.844827586206897;
            G_C_COMM[1][5] = 0.001;                 G_C_COMM[5][1] = 0.375586854460094;
            G_C_COMM[1][6] = 0.001;                 G_C_COMM[6][1] = 0.226635514018692;
            G_C_COMM[1][7] = 0.001;                 G_C_COMM[7][1] = 0.173076923076923;
            G_C_COMM[1][8] = 0.001;                 G_C_COMM[8][1] = 0.001;

            G_C_COMM[2][2] = 2.61607142857143;
            G_C_COMM[2][3] = 1.08035714285714;      G_C_COMM[3][2] = 1.16896551724138;
            G_C_COMM[2][4] = 1.29464285714286;      G_C_COMM[4][2] = 0.985632183908046;
            G_C_COMM[2][5] = 0.763392857142857;     G_C_COMM[5][2] = 1.22300469483568;
            G_C_COMM[2][6] = 0.044642857142857;     G_C_COMM[6][2] = 0.628504672897196;
            G_C_COMM[2][7] = 0.0625;                G_C_COMM[7][2] = 0.336538461538462;
            G_C_COMM[2][8] = 0.15625;               G_C_COMM[8][2] = 0.001;

            G_C_COMM[3][3] = 1.3448275862069;
            G_C_COMM[3][4] = 0.389655172413793;     G_C_COMM[4][3] = 1.12931034482759;
            G_C_COMM[3][5] = 0.513793103448276;     G_C_COMM[5][3] = 0.434272300469484;
            G_C_COMM[3][6] = 0.16551724137931;      G_C_COMM[6][3] = 0.836448598130841;
            G_C_COMM[3][7] = 0.037931034482759;     G_C_COMM[7][3] = 0.403846153846154;
            G_C_COMM[3][8] = 0.010344827586207;     G_C_COMM[8][3] = 0.001;

            G_C_COMM[4][4] = 1.02298850574713;
            G_C_COMM[4][5] = 0.195402298850575;     G_C_COMM[5][4] = 0.910798122065728;
            G_C_COMM[4][6] = 0.086206896551724;     G_C_COMM[6][4] = 0.380841121495327;
            G_C_COMM[4][7] = 0.21551724137931;      G_C_COMM[7][4] = 0.322115384615385;
            G_C_COMM[4][8] = 0.0001;                G_C_COMM[8][4] = 0.6;

            G_C_COMM[5][5] = 0.711267605633803;
            G_C_COMM[5][6] = 0.119718309859155;     G_C_COMM[6][5] = 0.647196261682243;
            G_C_COMM[5][7] = 0.150234741784038;     G_C_COMM[7][5] = 0.173076923076923;
            G_C_COMM[5][8] = 0.096244131455399;     G_C_COMM[8][5] = 0.001;

            G_C_COMM[6][6] = 0.457943925233645;
            G_C_COMM[6][7] = 0.294392523364486;     G_C_COMM[7][6] = 0.605769230769231;
            G_C_COMM[6][8] = 0.135514018691589;     G_C_COMM[8][6] = 0.6;

            G_C_COMM[7][7] = 1.27884615384615;
            G_C_COMM[7][8] = 0.197115384615385;     G_C_COMM[8][7] = 0.2;

            G_C_COMM[8][8] = 0.8;  

            for(int aci=0; aci<NUMAC; aci++){
                G_C_COMM[0][aci] *= G_CLO_CONTACT_COEFF_POSTLD_00;
                G_C_COMM[1][aci] *= G_CLO_CONTACT_COEFF_POSTLD_10;
                G_C_COMM[2][aci] *= G_CLO_CONTACT_COEFF_POSTLD_20;
                G_C_COMM[3][aci] *= G_CLO_CONTACT_COEFF_POSTLD_30;
                G_C_COMM[4][aci] *= G_CLO_CONTACT_COEFF_POSTLD_40;
                G_C_COMM[5][aci] *= G_CLO_CONTACT_COEFF_POSTLD_50;
                G_C_COMM[6][aci] *= G_CLO_CONTACT_COEFF_POSTLD_60;
                G_C_COMM[7][aci] *= G_CLO_CONTACT_COEFF_POSTLD_70;
                G_C_COMM[8][aci] *= G_CLO_CONTACT_COEFF_POSTLD_80;
            }  
            
            bContactRatesUpdatedAfterLockdown = true;
        }

        // normal contact rates from pre-COVID contact surveys
        if( !b_normalcy_contact_rates_updated && tt > G_CLO_NORMALCY_BEGINDAY )
        {
            // POLYMOD - Belgium matrix from `socialmixr` R package;
            G_C_COMM[0][0] = 3.03333333333333;
            G_C_COMM[0][1] = 0.786666666666667;     G_C_COMM[1][0] = 1.07746478873239; 
            G_C_COMM[0][2] = 1.09333333333333;      G_C_COMM[2][0] = 0.615384615384615;
            G_C_COMM[0][3] = 2.36;                  G_C_COMM[3][0] = 1.4578313253012;
            G_C_COMM[0][4] = 0.933333333333333;     G_C_COMM[4][0] = 0.329545454545455;
            G_C_COMM[0][5] = 0.7;                   G_C_COMM[5][0] = 0.735849056603774;
            G_C_COMM[0][6] = 0.426666666666667;     G_C_COMM[6][0] = 0.587301587301587;
            G_C_COMM[0][7] = 0.18;                  G_C_COMM[7][0] = 0.304347826086957;
            G_C_COMM[0][8] = 0.0333333333333333;    G_C_COMM[8][0] = 0.001;

            G_C_COMM[1][1] = 7.07746478873239;
            G_C_COMM[1][2] = 1.72535211267606;      G_C_COMM[2][1] = 1.50549450549451;
            G_C_COMM[1][3] = 1.3943661971831;       G_C_COMM[3][1] = 1.20481927710843;
            G_C_COMM[1][4] = 2.23943661971831;      G_C_COMM[4][1] = 1.36363636363636;
            G_C_COMM[1][5] = 0.633802816901408;     G_C_COMM[5][1] = 0.528301886792453;
            G_C_COMM[1][6] = 0.359154929577465;     G_C_COMM[6][1] = 0.206349206349206;
            G_C_COMM[1][7] = 0.183098591549296;     G_C_COMM[7][1] = 0.695652173913043;
            G_C_COMM[1][8] = 0.0352112676056338;    G_C_COMM[8][1] = 0.25;

            G_C_COMM[2][2] = 4.75824175824176;
            G_C_COMM[2][3] = 2.14285714285714;      G_C_COMM[3][2] = 1.81927710843373;
            G_C_COMM[2][4] = 1.89010989010989;      G_C_COMM[4][2] = 1.69318181818182;
            G_C_COMM[2][5] = 1.14285714285714;      G_C_COMM[5][2] = 2.0377358490566;
            G_C_COMM[2][6] = 0.252747252747253;     G_C_COMM[6][2] = 0.698412698412698;
            G_C_COMM[2][7] = 0.21978021978022;      G_C_COMM[7][2] = 0.434782608695652;
            G_C_COMM[2][8] = 0.0879120879120879;    G_C_COMM[8][2] = 0.001;

            G_C_COMM[3][3] = 3.72289156626506;
            G_C_COMM[3][4] = 2.39759036144578;      G_C_COMM[4][3] = 2.23863636363636;
            G_C_COMM[3][5] = 1.2289156626506;       G_C_COMM[5][3] = 2.29245283018868;
            G_C_COMM[3][6] = 0.734939759036145;     G_C_COMM[6][3] = 1.73015873015873;
            G_C_COMM[3][7] = 0.204819277108434;     G_C_COMM[7][3] = 0.826086956521739;
            G_C_COMM[3][8] = 0.180722891566265;     G_C_COMM[8][3] = 0.75;

            G_C_COMM[4][4] = 2.92045454545455;
            G_C_COMM[4][5] = 1.35227272727273;      G_C_COMM[5][4] = 2.32075471698113;
            G_C_COMM[4][6] = 0.545454545454545;     G_C_COMM[6][4] = 1.58730158730159;
            G_C_COMM[4][7] = 0.397727272727273;     G_C_COMM[7][4] = 1.69565217391304;
            G_C_COMM[4][8] = 0.465909090909091;     G_C_COMM[8][4] = 0.25;

            G_C_COMM[5][5] = 2.77358490566038;
            G_C_COMM[5][6] = 1.16981132075472;      G_C_COMM[6][5] = 1.15873015873016;
            G_C_COMM[5][7] = 0.622641509433962;     G_C_COMM[7][5] = 1.04347826086957;
            G_C_COMM[5][8] = 0.30188679245283;      G_C_COMM[8][5] = 0.75;

            G_C_COMM[6][6] = 2.01587301587302;
            G_C_COMM[6][7] = 0.952380952380952;     G_C_COMM[7][6] = 1.65217391304348;
            G_C_COMM[6][8] = 0.253968253968254;     G_C_COMM[8][6] = 0.75;

            G_C_COMM[7][7] = 1.43478260869565;
            G_C_COMM[7][8] = 0.173913043478261;     G_C_COMM[8][7] = 0.001;

            G_C_COMM[8][8] = 0.5;

            for(int aci=0; aci<NUMAC; aci++){
                G_C_COMM[0][aci] *= G_CLO_CONTACT_COEFF_NORM_00;
                G_C_COMM[1][aci] *= G_CLO_CONTACT_COEFF_NORM_10;
                G_C_COMM[2][aci] *= G_CLO_CONTACT_COEFF_NORM_20;
                G_C_COMM[3][aci] *= G_CLO_CONTACT_COEFF_NORM_30;
                G_C_COMM[4][aci] *= G_CLO_CONTACT_COEFF_NORM_40;
                G_C_COMM[5][aci] *= G_CLO_CONTACT_COEFF_NORM_50;
                G_C_COMM[6][aci] *= G_CLO_CONTACT_COEFF_NORM_60;
                G_C_COMM[7][aci] *= G_CLO_CONTACT_COEFF_NORM_70;
                G_C_COMM[8][aci] *= G_CLO_CONTACT_COEFF_NORM_80;
            } 

            /* for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
            {
                for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
                {
                    if (acs != aci){
                        G_C_COMM[acs][aci] = G_C_COMM[acs][acs] * G_C_COMM[aci][aci];
                    }
                }
            }
            for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
            {
                G_C_COMM[acs][acs] *= G_C_COMM[acs][acs];
            } */


            b_normalcy_contact_rates_updated = true;
        }


        if( !bPhase2BetterClinicalManagementBegun && tt >= G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY )
        {
            for(int ac=0; ac<NUMAC; ac++)
            {
                ppc->v_prob_HA_CA[ac] *= (G_CLO_ICUFRAC_DEV_SECONDPHASE / G_CLO_ICUFRAC_DEV);
            }            

            bPhase2BetterClinicalManagementBegun = true;
        }
        

        //
        // begin vaccination phase 1 (rollout)
        //
        // end the current stage of vaccination campaign
        // reset b_vac1_phase1_began, esp needed when changing to a new stage (so that we won't miss the moving S->Z step)
        // update next_beginday_vac1_phase1, next_endday_vac1_phase1
        if (b_vac1_phase1_began && (tt + DBL_EPSILON) >= (next_endday_vac1_phase1 + 1.0)){
            if (ppc->Increase_index_current_day_phase1_vac1() > -1){
                next_beginday_vac1_phase1 = ((double)ppc->v_begin_days_phase1_vac1[ppc->index_current_day_phase1_vac1]) ;
                next_endday_vac1_phase1 = ((double)ppc->v_end_days_phase1_vac1[ppc->index_current_day_phase1_vac1]) ;
                b_vac1_phase1_began = false; 
                // if (G_B_USE_VAC_FRAC) { b_done_vac_frac_1 = true; }
            }
        }
        // calculate the number of vaccinees per day
        // if (!b_vac1_phase1_began && tt > G_CLO_VAC_1_PHASE1_BEGINDAY){
        if ((tt + DBL_EPSILON) >= next_beginday_vac1_phase1 && tt < (next_endday_vac1_phase1 + 1.0)){
        // if (!b_vac1_phase1_began && (tt + DBL_EPSILON) >= next_beginday_vac1_phase1 && tt < (next_endday_vac1_phase1 + 1.0)){
            b_vac1_phase1_began = true;

            // keep track of whether there are still susceptibles in each age group to be vaccinated
            each_step_total_S_Z_1_phase1 = 0.0; 
            
            int ac = 0;

            if (G_B_USE_VAC_FRAC){
                // index for current day in v_begin_days_phase1_vac1 vector
                // if v_vac_ratios_phase1_vac1 vector is shorter than v_begin_days_phase1_vac1 vector, use the last elements in v_vac_ratios_phase1_vac1
                int icd = (ppc->index_current_day_phase1_vac1 < ppc->v_vac_fracs_phase1_vac1[0].size()) ? ppc->index_current_day_phase1_vac1 : (ppc->v_vac_fracs_phase1_vac1[0].size() - 1);

                for (ac = 0; ac < NUMAC; ac++){
                    // reset v_already_moved_to_Z_1_phase1
                    ppc->v_already_moved_to_Z_1_phase1[ac] = 0.0;
                    // ppc->v_number_to_Z_1_phase1[ac] = 0.0;
                    // absolute number of doses
                    ppc->v_number_to_Z_1_phase1[ac] = ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] * ppc->v_vac_fracs_phase1_vac1[ac][icd] /  steps_per_day; 
                    ppc->v_number_to_Z_1_phase1[ac] = floor(ppc->v_number_to_Z_1_phase1[ac]);
                    // keep track of total S_Z_1 in the whole population, to re-distribute unused doses later
                    each_step_total_S_Z_1_phase1 += ppc->v_number_to_Z_1_phase1[ac];
                    
                    // printf("\n%d calc number_to_Z FRAC %f\n", ac, ppc->v_number_to_Z_1_phase1[ac]);
                }

                /* if (G_B_USE_VAC_RATIO){
                    min_vac_frac_1 = 1.0;
                    sum_weighted_vac_frac_1 = 0.0;

                    // index for current day in v_begin_days_phase1_vac1 vector
                    // if v_vac_ratios_phase1_vac1 vector is shorter than v_begin_days_phase1_vac1 vector, use the last elements in v_vac_ratios_phase1_vac1
                    int icd = (ppc->index_current_day_phase1_vac1 < ppc->v_vac_ratios_phase1_vac1[0].size()) ? ppc->index_current_day_phase1_vac1 : (ppc->v_vac_ratios_phase1_vac1[0].size() - 1);
                    
                    // find min vac_frac
                    for (ac = 0; ac < NUMAC; ac++) {
                        if (ppc->v_vac_ratios_phase1_vac1[ac][icd] > 0.0) { min_vac_frac_1 = Min(min_vac_frac_1, ppc->v_vac_ratios_phase1_vac1[ac][icd]); }
                    }
                    // sum of weighted vac frac * pop frac
                    for (ac =0; ac < NUMAC; ac++){
                        sum_weighted_vac_frac_1 += ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac] / min_vac_frac_1;
                    }
                    for (ac = 0; ac < NUMAC; ac++){
                        if (ppc->v_number_to_Z_1_phase1[ac] < 1.0){
                            ppc->v_number_to_Z_1_phase1[ac] = (ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] * ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac]) / (min_vac_frac_1 * sum_weighted_vac_frac_1 *  steps_per_day); // proportional to pop frac
                        }
                    }
                } */
            }

            if (!G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO){
                min_vac_frac_1 = 1.0;
                sum_weighted_vac_frac_1 = 0.0;
                
                // index for current day in v_begin_days_phase1_vac1 vector
                // if v_vac_ratios_phase1_vac1 vector is shorter than v_begin_days_phase1_vac1 vector, use the last elements in v_vac_ratios_phase1_vac1
                int icd = (ppc->index_current_day_phase1_vac1 < ppc->v_vac_ratios_phase1_vac1[0].size()) ? ppc->index_current_day_phase1_vac1 : (ppc->v_vac_ratios_phase1_vac1[0].size() - 1);
                
                // find min vac_frac
                for (ac = 0; ac < NUMAC; ac++) {
                    if (ppc->v_vac_ratios_phase1_vac1[ac][icd] > 0.0) { min_vac_frac_1 = Min(min_vac_frac_1, ppc->v_vac_ratios_phase1_vac1[ac][icd]); }
                }
                // sum of weighted vac frac * pop frac
                for (ac =0; ac < NUMAC; ac++){
                    sum_weighted_vac_frac_1 += ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac] / min_vac_frac_1;
                }
                for (ac = 0; ac < NUMAC; ac++){
                    // reset
                    ppc->v_already_moved_to_Z_1_phase1[ac] = 0.0;
                    // ppc->v_number_to_Z_1_phase1[ac] = 0.0;
                    
                    ppc->v_number_to_Z_1_phase1[ac] = (ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] * ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac]) / (min_vac_frac_1 * sum_weighted_vac_frac_1 *  steps_per_day); // proportional to pop frac
                    ppc->v_number_to_Z_1_phase1[ac] = ppc->v_number_to_Z_1_phase1[ac];
                    // keep track of total S_Z_1 in the whole population, to re-distribute unused doses later
                    each_step_total_S_Z_1_phase1 += ppc->v_number_to_Z_1_phase1[ac];

                    // printf("\n%d calc number_to_Z RATIO %f\n", ac, ppc->v_number_to_Z_1_phase1[ac]);
                }
            }
            
            // ppc->Increase_index_current_day_phase1_vac1();

            if (each_step_total_S_Z_1_phase1 < ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] / steps_per_day){
                for (ac = NUMAC - 1; ac >= 0; --ac){
                    if (ppc->v_number_to_Z_1_phase1[ac] > 0){
                        ppc->v_number_to_Z_1_phase1[ac] += (ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] / steps_per_day) - each_step_total_S_Z_1_phase1;
                        each_step_total_S_Z_1_phase1 = ppc->v_dpd_phase1_vac1[ppc->index_current_day_phase1_vac1] / steps_per_day;
                        break;
                    }
                }
            }
            // for( ac=0; ac < NUMAC; ac++){
            //     printf("\n%d end calc number_to_Z %f\n", ac, ppc->v_number_to_Z_1_phase1[ac]);
            // }
            
        }


        
        // vaccinate people, i.e. moving S to Z
        // if (b_vac1_phase1_began && tt > G_CLO_VAC_1_PHASE1_BEGINDAY && tt <= G_CLO_VAC_1_PHASE1_ENDDAY){
        // if (b_vac1_phase1_began && tt >= next_beginday_vac1_phase1 && tt < (next_endday_vac1_phase1 + 1.0)){
        if (b_vac1_phase1_began && (tt + DBL_EPSILON) >= next_beginday_vac1_phase1 && tt < (next_endday_vac1_phase1 + 1.0)){
            // printf("%f - %f\n", tt, each_step_total_S_Z_1_phase1);
            double total_S_Z = 0.0, total_to_Z_next = each_step_total_S_Z_1_phase1;
            double tmp = 0.0, frac_R = 1.0; 
            double sub_pop = 0.0, sub_pop_next = 0.0; // to determine vac_frac for the unused doses
            for (int ac = 0; ac < NUMAC; ac++){
                if ( G_B_TEST_BEFORE_VACCINATE ){ frac_R = 1.0 - (1.0 - ppc->v_prob_E_A[ac])*ppc->v[ i_reporting_rate ]; } 
                else { frac_R = 1.0; }

                if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                    double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                    for (size_t i = 0; i < NUME; i++) {         sum_SEAR += floor(yic[STARTE  + i*NUMAC + ac]);  }
                    for (size_t i = 0; i < NUMA; i++) {         sum_SEAR += floor(yic[STARTA  + i*NUMAC + ac]);  }
                    for (size_t i = 0; i < NUMRR; i++){         sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                    // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                    // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                    // sum_SEAR = floor(sum_SEAR);

                    if (sum_SEAR < ppc->v_number_to_Z_1_phase1[ac]){
                        tmp = 1.0;
                    } else {
                        tmp = ppc->v_number_to_Z_1_phase1[ac] / sum_SEAR;
                        // sub_pop_next += ppc->v_pop_frac[ac];
                        // b_done_vac_frac_1 = false;
                    }

                    // if (tmp * sum_SEAR >= 1.0){
                        yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                     // count all vaccine doses being given out

                        yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R)); // but only S and R move to Z 
                        yic[ac]             -= tmp *  floor(yic[ac]         ) ;
                        yic[STARTR + ac]    -= tmp *  floor(yic[STARTR + ac] *frac_R) ;
                        // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                        // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                        for (size_t i = 0; i < NUMRR; i++){   
                            yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                            yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;  
                        }
                    // }
                    total_to_Z_next      -= tmp * sum_SEAR;
                    ppc->v_already_moved_to_Z_1_phase1[ac] = tmp * sum_SEAR;
                }
                
                /* if ( G_B_TEST_BEFORE_VACCINATE ) // test before vaccinating --> equal chance for S, E, A to be vaccinated. S will move to Z; E and A will stay in their classes
                {
                    if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                        double sum_SEA = floor(yic[ac]);
                        for (size_t i = 0; i < NUME; i++){          sum_SEA += floor(yic[STARTE + i*NUMAC + ac]);  }
                        for (size_t i = 0; i < NUMA; i++){          sum_SEA += floor(yic[STARTA + i*NUMAC + ac]);  }
                        // sum_SEA = floor(sum_SEA);

                        if (sum_SEA < ppc->v_number_to_Z_1_phase1[ac]){
                            tmp = 1.0;
                        } else {
                            tmp = ppc->v_number_to_Z_1_phase1[ac] / sum_SEA;
                            // sub_pop_next += ppc->v_pop_frac[ac];
                            // b_done_vac_frac_1 = false;
                        }

                        // if (tmp * sum_SEA >= 1.0){
                            yic[STARTZ_1 + ac]  += tmp * floor(yic[ac]); // only S moves to Z
                            yic[STARTJZ_1 + ac] += tmp * sum_SEA; // but count all vaccine doses being given out
                            // printf("%d - %f\n", ac, yic[STARTJZ_1 + ac]);

                            yic[ac]              -= floor(yic[ac])              * tmp ;
                            // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                            // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }

                        // }
                        
                        total_to_Z_next      -= tmp * sum_SEA;
                        ppc->v_already_moved_to_Z_1_phase1[ac] = tmp * sum_SEA;
                    }
                }
                else // no test before vaccinating --> equal chance for S, E, A, R to be vaccinated. S and R will move to Z; E and A will stay in their classes
                {
                    if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                        double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                        for (size_t i = 0; i < NUME; i++){          sum_SEAR += floor(yic[STARTE + i*NUMAC + ac]);  }
                        for (size_t i = 0; i < NUMA; i++){          sum_SEAR += floor(yic[STARTA + i*NUMAC + ac]);  }
                        for (size_t i = 0; i < NUMRR; i++){         sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                        // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                        // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                        // sum_SEAR = floor(sum_SEAR);

                        if (sum_SEAR < ppc->v_number_to_Z_1_phase1[ac]){
                            tmp = 1.0;
                        } else {
                            tmp = ppc->v_number_to_Z_1_phase1[ac] / sum_SEAR;
                            // sub_pop_next += ppc->v_pop_frac[ac];
                            // b_done_vac_frac_1 = false;
                        }

                        // if (tmp * sum_SEAR >= 1.0){
                            yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                     // count all vaccine doses being given out

                            yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R)); // but only S and R move to Z 
                            yic[ac]             -= tmp *  floor(yic[ac]         ) ;
                            yic[STARTR + ac]    -= tmp *  floor(yic[STARTR + ac] *frac_R) ;
                            // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                            // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                            for (size_t i = 0; i < NUMRR; i++){   
                                yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                                yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;  
                            }
                        // }
                        total_to_Z_next      -= tmp * sum_SEAR;
                        ppc->v_already_moved_to_Z_1_phase1[ac] = tmp * sum_SEAR;
                    }
                } */
            } 
            
            //
            // distribute leftover doses within "vac frac" groups
            //
            sub_pop_next = 0.0; 
            if (G_B_USE_VAC_FRAC && total_to_Z_next > (1.0 - DBL_EPSILON)){
                int ac = 0;
                for (ac = 0; ac < NUMAC; ac++){
                    if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                        if (ppc->v_already_moved_to_Z_1_phase1[ac] < (ppc->v_number_to_Z_1_phase1[ac] - 1) ) {
                            ppc->v_number_to_Z_1_phase1[ac] = -1;
                        } else {
                            // sub_pop_next += ppc->v_pop_frac[ac];
                            sub_pop_next += 1;
                        }
                    } else {
                        ppc->v_number_to_Z_1_phase1[ac] = -1;
                    }
                    
                    
                }
                // while (! b_done_vac_frac_1 && total_to_Z_next >= 1.0 ){
                while (total_to_Z_next > (1.0 - DBL_EPSILON) && sub_pop_next > DBL_EPSILON){
                    sub_pop = sub_pop_next;
                    total_S_Z = total_to_Z_next;
                    for (ac = 0; ac < NUMAC; ac++){
                        if ( G_B_TEST_BEFORE_VACCINATE ){ frac_R = 1.0 - (1.0 - ppc->v_prob_E_A[ac])*ppc->v[ i_reporting_rate ]; } 
                        else { frac_R = 1.0; }

                        if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                        // if (ppc->v_number_to_Z_1_phase1[ac] > 0.9 && ppc->v_already_moved_to_Z_1_phase1[ac] == ppc->v_number_to_Z_1_phase1[ac] ){
                            // tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group
                            tmp = total_S_Z / sub_pop; // max vaccines to be randomly distributed to this age group

                            double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                            for (size_t i = 0; i < NUME; i++) {          sum_SEAR += floor(yic[STARTE + i*NUMAC + ac]);  }
                            for (size_t i = 0; i < NUMA; i++) {          sum_SEAR += floor(yic[STARTA + i*NUMAC + ac]);  }
                            for (size_t i = 0; i < NUMRR; i++){          sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                            // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                            // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                            // sum_SEAR = floor(sum_SEAR);

                            if (sum_SEAR < tmp){
                                tmp = 1.0;
                                // sub_pop_next -= ppc->v_pop_frac[ac];
                            } else {
                                tmp = tmp / sum_SEAR;
                            }

                            if (tmp * sum_SEAR < (1.0 - DBL_EPSILON)){
                                // sub_pop_next -= ppc->v_pop_frac[ac];
                                sub_pop_next -= 1;
                            } else {
                                ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEAR;

                                yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                      // count all vaccine doses being given out

                                yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R));  // but only S and R move to Z  
                                yic[ac]             -= tmp * floor(yic[ac]         ) ;
                                yic[STARTR + ac]    -= tmp * floor(yic[STARTR + ac] *frac_R) ;
                                // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                                // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                                for (size_t i = 0; i < NUMRR; i++){   
                                    yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                                    yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                                }
                                // yic[STARTRHOSP + ac] -= yic[STARTRHOSP + ac] * tmp ;
                                // for (size_t i = 0; i < NUMRRHOSP; i++){     yic[STARTRRHOSP + i*NUMAC + ac] -= yic[STARTRRHOSP + i*NUMAC + ac]  * tmp ;  }
                                
                                total_to_Z_next      -= tmp * sum_SEAR;
                            //     ppc->v_number_to_Z_1_phase1[ac] = ppc->v_already_moved_to_Z_1_phase1[ac];
                            }
                        }

                        // printf("\n%d distribute leftover doses within vac frac ppc->v_number_to_Z_1_phase1[ac] %f ppc->v_already_moved_to_Z_1_phase1[ac] %f\n", ac, ppc->v_number_to_Z_1_phase1[ac], ppc->v_already_moved_to_Z_1_phase1[ac]);
                        /* if ( G_B_TEST_BEFORE_VACCINATE ) // test before vaccinating --> equal chance for S, E, A to be vaccinated. S will move to Z; E and A will stay in their classes
                        {
                            if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                            // if (ppc->v_number_to_Z_1_phase1[ac] > 0.9 && ppc->v_already_moved_to_Z_1_phase1[ac] == ppc->v_number_to_Z_1_phase1[ac] ){
                                // tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group
                                tmp = total_S_Z / sub_pop; // max vaccines to be randomly distributed to this age group

                                double sum_SEA = floor(yic[ac]) ;
                                for (size_t i = 0; i < NUME; i++){          sum_SEA += floor(yic[STARTE + i*NUMAC + ac]);  }
                                for (size_t i = 0; i < NUMA; i++){          sum_SEA += floor(yic[STARTA + i*NUMAC + ac]);  }
                                // sum_SEA = floor(sum_SEA);

                                if (sum_SEA < tmp){
                                    tmp = 1.0;
                                    // sub_pop_next -= ppc->v_pop_frac[ac];
                                } else {
                                    tmp = tmp / sum_SEA;
                                }

                                if (tmp * sum_SEA < (1.0 - DBL_EPSILON)){
                                    // sub_pop_next -= ppc->v_pop_frac[ac];
                                    sub_pop_next -= 1;
                                } else {
                                    ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEA;

                                    yic[STARTZ_1 + ac]  += tmp * floor(yic[ac]); // only S moves to Z
                                    yic[STARTJZ_1 + ac] += tmp * sum_SEA; // but count all vaccine doses being given out

                                    yic[ac]              -= tmp * floor(yic[ac]) ;
                                    // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                                    // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                                    
                                    total_to_Z_next      -= tmp * sum_SEA;
                                    // ppc->v_number_to_Z_1_phase1[ac] = ppc->v_already_moved_to_Z_1_phase1[ac];
                                }
                                
                            }
                        }
                        else // no test before vaccinating --> equal chance for S, E, A, R to be vaccinated. S and R will move to Z; E and A will stay in their classes
                        {
                            if (ppc->v_number_to_Z_1_phase1[ac] > (1.0 - DBL_EPSILON)){
                            // if (ppc->v_number_to_Z_1_phase1[ac] > 0.9 && ppc->v_already_moved_to_Z_1_phase1[ac] == ppc->v_number_to_Z_1_phase1[ac] ){
                                // tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group
                                tmp = total_S_Z / sub_pop; // max vaccines to be randomly distributed to this age group

                                double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                                for (size_t i = 0; i < NUME; i++) {          sum_SEAR += floor(yic[STARTE + i*NUMAC + ac]);  }
                                for (size_t i = 0; i < NUMA; i++) {          sum_SEAR += floor(yic[STARTA + i*NUMAC + ac]);  }
                                for (size_t i = 0; i < NUMRR; i++){          sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                                // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                                // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                                // sum_SEAR = floor(sum_SEAR);

                                if (sum_SEAR < tmp){
                                    tmp = 1.0;
                                    // sub_pop_next -= ppc->v_pop_frac[ac];
                                } else {
                                    tmp = tmp / sum_SEAR;
                                }

                                if (tmp * sum_SEAR < (1.0 - DBL_EPSILON)){
                                    // sub_pop_next -= ppc->v_pop_frac[ac];
                                    sub_pop_next -= 1;
                                } else {
                                    ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEAR;

                                    yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                      // count all vaccine doses being given out

                                    yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R));  // but only S and R move to Z  
                                    yic[ac]             -= tmp * floor(yic[ac]         ) ;
                                    yic[STARTR + ac]    -= tmp * floor(yic[STARTR + ac] *frac_R) ;
                                    // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                                    // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                                    for (size_t i = 0; i < NUMRR; i++){   
                                        yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                                        yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R)  ;
                                    }
                                    // yic[STARTRHOSP + ac] -= yic[STARTRHOSP + ac] * tmp ;
                                    // for (size_t i = 0; i < NUMRRHOSP; i++){     yic[STARTRRHOSP + i*NUMAC + ac] -= yic[STARTRRHOSP + i*NUMAC + ac]  * tmp ;  }
                                    
                                    total_to_Z_next      -= tmp * sum_SEAR;
                                //     ppc->v_number_to_Z_1_phase1[ac] = ppc->v_already_moved_to_Z_1_phase1[ac];
                                }
                            }
                        } */
                        // printf("%d distribute leftover doses within vac frac ppc->v_number_to_Z_1_phase1[ac] %f ppc->v_already_moved_to_Z_1_phase1[ac] %f\n", ac, ppc->v_number_to_Z_1_phase1[ac], ppc->v_already_moved_to_Z_1_phase1[ac]);
                    }
                }

                if (G_B_USE_VAC_RATIO){
                    min_vac_frac_1 = 1.0;
                    sum_weighted_vac_frac_1 = 0.0;
                    double subsum = 0.0;

                    // index for current day in v_begin_days_phase1_vac1 vector
                    // if v_vac_ratios_phase1_vac1 vector is shorter than v_begin_days_phase1_vac1 vector, use the last elements in v_vac_ratios_phase1_vac1
                    int icd = (ppc->index_current_day_phase1_vac1 < ppc->v_vac_ratios_phase1_vac1[0].size()) ? ppc->index_current_day_phase1_vac1 : (ppc->v_vac_ratios_phase1_vac1[0].size() - 1);
                    
                    // find min vac_frac
                    for (ac = 0; ac < NUMAC; ac++) {
                        if (ppc->v_vac_ratios_phase1_vac1[ac][icd] > 0.0 && ppc->v_already_moved_to_Z_1_phase1[ac] < DBL_EPSILON ) { 
                            min_vac_frac_1 = Min(min_vac_frac_1, ppc->v_vac_ratios_phase1_vac1[ac][icd]); 
                        }
                    }
                    // sum of weighted vac frac * pop frac
                    for (ac =0; ac < NUMAC; ac++){
                        if (ppc->v_vac_ratios_phase1_vac1[ac][icd] > 0.0 && ppc->v_already_moved_to_Z_1_phase1[ac] < DBL_EPSILON ) { 
                            sum_weighted_vac_frac_1 += ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac] / min_vac_frac_1;
                        }
                    }
                    for (ac = 0; ac < NUMAC; ac++){
                        if (ppc->v_vac_ratios_phase1_vac1[ac][icd] > 0.0){
                            if ( ppc->v_already_moved_to_Z_1_phase1[ac] < (1.0 - DBL_EPSILON) ) { 
                                ppc->v_number_to_Z_1_phase1[ac] = (total_to_Z_next * ppc->v_vac_ratios_phase1_vac1[ac][icd] * ppc->v_pop_frac[ac]) / (min_vac_frac_1 * sum_weighted_vac_frac_1 *  steps_per_day); // proportional to pop frac
                                ppc->v_number_to_Z_1_phase1[ac] = floor(ppc->v_number_to_Z_1_phase1[ac]);
                                subsum += ppc->v_number_to_Z_1_phase1[ac];
                                // printf("\n%d distribute leftover doses within vac frac USE RATIO %f\n", ac, ppc->v_number_to_Z_1_phase1[ac]);
                            } else {
                                ppc->v_number_to_Z_1_phase1[ac] = -1;// ppc->v_already_moved_to_Z_1_phase1[ac] ;
                            }
                            
                        }
                    }

                    if (subsum < total_to_Z_next / steps_per_day){
                        for (ac = NUMAC - 1; ac >= 0; --ac){
                            if (ppc->v_number_to_Z_1_phase1[ac] > 0){
                                ppc->v_number_to_Z_1_phase1[ac] += (total_to_Z_next / steps_per_day) - subsum;
                                break;
                            }
                        }
                    }
                }
            }
            sub_pop_next = 0.0; 
            // printf("%f - %f\n", tt, total_to_Z_next);
            // distribute the unused doses to available (susceptible) age groups
            if (total_to_Z_next > (1.0 - DBL_EPSILON)){
                for (int ac = 0; ac < NUMAC; ac++){
                    if ( (!G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 && ppc->v_number_to_Z_1_phase1[ac] < 1.0) ||
                         (G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 ) ){
                        sub_pop_next += ppc->v_pop_frac[ac];
                    }
                }
            }
            while (total_to_Z_next > (1.0 - DBL_EPSILON) && sub_pop_next > DBL_EPSILON){
                sub_pop = sub_pop_next;
                total_S_Z = total_to_Z_next;
                for (int ac = 0; ac < NUMAC; ac++){
                    if ( G_B_TEST_BEFORE_VACCINATE ){ frac_R = 1.0 - (1.0 - ppc->v_prob_E_A[ac])*ppc->v[ i_reporting_rate ]; } 
                    else { frac_R = 1.0; }

                    if ( (!G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 && ppc->v_number_to_Z_1_phase1[ac] < 1.0) ||
                            (G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 ) ){
                        tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group

                        double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                        for (size_t i = 0; i < NUME; i++) {         sum_SEAR += floor(yic[STARTE + i*NUMAC + ac]);  }
                        for (size_t i = 0; i < NUMA; i++) {         sum_SEAR += floor(yic[STARTA + i*NUMAC + ac]);  }
                        for (size_t i = 0; i < NUMRR; i++){         sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                        // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                        // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                        // sum_SEAR = floor(sum_SEAR);

                        // sum_SEAR = floor(sum_SEAR);

                        if (sum_SEAR < tmp){
                            tmp = 1.0;
                            sub_pop_next -= ppc->v_pop_frac[ac];
                        } else {
                            tmp = tmp / sum_SEAR;
                        }

                        // if (tmp * sum_SEAR >= 1.0){
                            yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                      // count all vaccine doses being given out

                            yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R));  // but only S and R move to Z  
                            yic[ac]             -= tmp *  floor(yic[ac]         ) ;
                            yic[STARTR + ac]    -= tmp *  floor(yic[STARTR + ac] *frac_R) ;
                            // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                            // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                            for (size_t i = 0; i < NUMRR; i++){   
                                yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R) ;
                                yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R) ;
                            }
                            // yic[STARTRHOSP + ac] -= yic[STARTRHOSP + ac] * tmp ;
                            // for (size_t i = 0; i < NUMRRHOSP; i++){     yic[STARTRRHOSP + i*NUMAC + ac] -= yic[STARTRRHOSP + i*NUMAC + ac]  * tmp ;  }
                        // }

                        // total_to_Z_next      -= ceil(tmp * sum_SEAR);
                        total_to_Z_next      -= tmp * sum_SEAR;
                        ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEAR;
                    }

                    // printf("%d distribute leftover doses within vac ratio ppc->v_number_to_Z_1_phase1[ac] %f ppc->v_already_moved_to_Z_1_phase1[ac] %f\n", ac, ppc->v_number_to_Z_1_phase1[ac], ppc->v_already_moved_to_Z_1_phase1[ac]);
                    /* if ( G_B_TEST_BEFORE_VACCINATE ) // test before vaccinating --> equal chance for S, E, A to be vaccinated. S will move to Z; E and A will stay in their classes
                    {
                        if ( (!G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 && ppc->v_number_to_Z_1_phase1[ac] < 1.0) ||
                             (G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 ) ){
                            tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group

                            double sum_SEA = floor(yic[ac]) ;
                            for (size_t i = 0; i < NUME; i++){          sum_SEA += floor(yic[STARTE + i*NUMAC + ac]);  }
                            for (size_t i = 0; i < NUMA; i++){          sum_SEA += floor(yic[STARTA + i*NUMAC + ac]);  }
                            // sum_SEA = floor(sum_SEA);
                            
                            // sum_SEA = floor(sum_SEA);

                            if (sum_SEA < tmp){
                                tmp = 1.0;
                                sub_pop_next -= ppc->v_pop_frac[ac];
                            } else {
                                tmp = tmp / sum_SEA;
                            }

                            // if (tmp * sum_SEA >= 1.0){
                                yic[STARTZ_1 + ac]  += tmp * floor(yic[ac]); // only S moves to Z
                                yic[STARTJZ_1 + ac] += tmp * sum_SEA; // but count all vaccine doses being given out

                                yic[ac]              -= tmp * floor(yic[ac])  ;
                                // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                                // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                            // }                            
                            // total_to_Z_next      -= ceil(tmp * sum_SEA);
                            total_to_Z_next      -= tmp * sum_SEA;
                            ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEA;
                        }
                    }
                    else // no test before vaccinating --> equal chance for S, E, A, R to be vaccinated. S and R will move to Z; E and A will stay in their classes
                    {
                        if ( (!G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 && ppc->v_number_to_Z_1_phase1[ac] < 1.0) ||
                             (G_B_USE_VAC_FRAC && G_B_USE_VAC_RATIO && ppc->v_number_to_Z_1_phase1[ac] > 0.0 ) ){
                            tmp = total_S_Z * ppc->v_pop_frac[ac] / sub_pop; // max vaccines to be randomly distributed to this age group

                            double sum_SEAR = floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R) ;
                            for (size_t i = 0; i < NUME; i++) {         sum_SEAR += floor(yic[STARTE + i*NUMAC + ac]);  }
                            for (size_t i = 0; i < NUMA; i++) {         sum_SEAR += floor(yic[STARTA + i*NUMAC + ac]);  }
                            for (size_t i = 0; i < NUMRR; i++){         sum_SEAR += floor(yic[STARTRR + i*NUMAC + ac] *frac_R);  }
                            // double sum_SEARRh = yic[ac] + yic[STARTR + ac] + yic[STARTRHOSP + ac];
                            // for (size_t i = 0; i < NUMRRHOSP; i++){     sum_SEARRh += yic[STARTRRHOSP + i*NUMAC + ac];  }
                            // sum_SEAR = floor(sum_SEAR);

                            // sum_SEAR = floor(sum_SEAR);

                            if (sum_SEAR < tmp){
                                tmp = 1.0;
                                sub_pop_next -= ppc->v_pop_frac[ac];
                            } else {
                                tmp = tmp / sum_SEAR;
                            }

                            // if (tmp * sum_SEAR >= 1.0){
                                yic[STARTJZ_1 + ac] += tmp * sum_SEAR;                      // count all vaccine doses being given out

                                yic[STARTZ_1 + ac]  += tmp * (floor(yic[ac]) + floor(yic[STARTR + ac] *frac_R));  // but only S and R move to Z  
                                yic[ac]             -= tmp *  floor(yic[ac]         ) ;
                                yic[STARTR + ac]    -= tmp *  floor(yic[STARTR + ac] *frac_R) ;
                                // for (size_t i = 0; i < NUME; i++){          yic[STARTE + i*NUMAC + ac]      -= yic[STARTE + i*NUMAC + ac]       * tmp ;  } // E and A stay in their classes
                                // for (size_t i = 0; i < NUMA; i++){          yic[STARTA + i*NUMAC + ac]      -= yic[STARTA + i*NUMAC + ac]       * tmp ;  }
                                for (size_t i = 0; i < NUMRR; i++){   
                                    yic[STARTZ_1 + ac]              += tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R) ;
                                    yic[STARTRR + i*NUMAC + ac]     -= tmp * floor(yic[STARTRR + i*NUMAC + ac] *frac_R) ;
                                }
                                // yic[STARTRHOSP + ac] -= yic[STARTRHOSP + ac] * tmp ;
                                // for (size_t i = 0; i < NUMRRHOSP; i++){     yic[STARTRRHOSP + i*NUMAC + ac] -= yic[STARTRRHOSP + i*NUMAC + ac]  * tmp ;  }
                            // }

                            // total_to_Z_next      -= ceil(tmp * sum_SEAR);
                            total_to_Z_next      -= tmp * sum_SEAR;
                            ppc->v_already_moved_to_Z_1_phase1[ac] += tmp * sum_SEAR;
                        }
                    } */
                    // printf("%d distribute leftover doses within vac ratio ppc->v_number_to_Z_1_phase1[ac] %f ppc->v_already_moved_to_Z_1_phase1[ac] %f\n", ac, ppc->v_number_to_Z_1_phase1[ac], ppc->v_already_moved_to_Z_1_phase1[ac]);   
                }
            }

        }
        



        // 
        // high-efficacy vaccine
        //
        // TODO: number_S_Z_2_phase1 must be proportional to pop_frac
        //
        /* if (!b_vac2_phase1_began && tt > G_CLO_VAC_2_PHASE1_BEGINDAY){
            // double phase1_duration = G_CLO_VAC_2_PHASE1_ENDDAY - G_CLO_VAC_2_PHASE1_BEGINDAY + 1;  
            int ac = 0;
            // find min vac_frac
            for (ac = 0; ac < NUMAC; ac++) {
                if (ppc->v_prob_S_Z_2[ac] > 0.0) { min_vac_frac_2 = Min(min_vac_frac_2, ppc->v_prob_S_Z_2[ac]); }
            }  
            // sum of weighted vac frac * pop frac
            for (ac =0; ac < NUMAC; ac++){
                sum_weighted_vac_frac_2 += ppc->v_prob_S_Z_2[ac] * ppc->v_pop_frac[ac] / min_vac_frac_2;
            }
            for (ac = 0; ac < NUMAC; ac++)
            {
                if (ppc->v_prob_S_Z_2[ac] > 0.0){
                    // ppc->v_number_S_Z_2_phase1[ac] = yic[ac] * ppc->v_prob_S_Z_2[ac] / (phase1_duration * steps_per_day) ;
                    // ppc->v_number_S_Z_2_phase1[ac] = ppc->v[i_vac2_phase1_dpd] * ppc->v_prob_S_Z_2[ac] /  steps_per_day; // absolute number of doses
                    ppc->v_number_S_Z_2_phase1[ac] = (ppc->v[i_vac2_phase1_dpd] * ppc->v_prob_S_Z_2[ac] * ppc->v_pop_frac[ac]) / (min_vac_frac_2 * sum_weighted_vac_frac_2 *  steps_per_day); // proportional to pop frac
                }
            }       
            b_vac2_phase1_began = true;
        }
        // vaccinate people, i.e. moving S to Z
        if (b_vac2_phase1_began && tt > G_CLO_VAC_2_PHASE1_BEGINDAY && tt <= G_CLO_VAC_2_PHASE1_ENDDAY){
            double tmp = 0.0;
            for (int ac = 0; ac < NUMAC; ac++){
                tmp = Min(yic[ac], ppc->v_number_S_Z_2_phase1[ac]);
                yic[ac] -= tmp;
                yic[STARTZ_2 + ac] += tmp;
                yic[STARTJZ_2 + ac] += tmp;
            }
            
        } */
        
        
        
        
        // check if the hospitalization rates need to be updated (because we are out of the early period)
        /*if( ppc->earlymarch_highhosp_period )
        {
            if( tt > ppc->earlymarch_highhosp_endday )
            {
                ppc->end_earlymarch_hosprates();
            }
        }*/
        
        // check if the beta value needs to be updated
        if( tt > NextTimeForBetaUpdate )
        {
            ppc->assign_new_beta();
            NextTimeForBetaUpdate = ppc->get_new_update_time();
        }
        



        // #### main ODE integration ####
        if ( odeint(yic, nvar, tt, tt+rkstep, eps, &h1, hmin, &nok, &nbad, derivs, rkqs) != 1)
      	{
	        fprintf(stderr, "\n\nEXITING BC OF ERROR IN ODEINT\n\n");
            exit(-1);
        }




        // #### print state of the system ####
        if( counter%steps_per_day==0 && !G_B_DIAGNOSTIC_MODE && !G_B_CHECKPOP_MODE && tt>49.5 )
        {
            write(tt, yic, dim);
        }
        
        
        // ### BEGIN only enter block below if you're checking the total population size
        if( counter%(50*steps_per_day)==0 && G_B_CHECKPOP_MODE )
        {
            totalpop=0.0;
            for(int k=0; k<STARTJ; k++) totalpop += yic[k];
            for(int k=STARTZ_1; k<STARTJZ_1; k++) totalpop += yic[k];
            printf("\n\t totalpop=%1.16f", totalpop);
            
        }
        // ### END only enter block if you're checking the total population size
        
        tt += rkstep;
        counter++;
    }
    //
    //END OF LOOP INTEGRATING ODEs
    //

    // print out results of the last step
    if( !G_B_DIAGNOSTIC_MODE && !G_B_CHECKPOP_MODE )
    {
	    write(tt, yic, dim);
    }

    if( G_B_CHECKPOP_MODE ) printf("\n\n");
    //printf("\n\n    %d \n\n", Q );
    
    //delete[] yic;    
    
    //return vv;
}

