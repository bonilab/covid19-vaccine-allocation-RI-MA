#ifndef ESSENTIALS
#define ESSENTIALS


#include <math.h>
#include <stdio.h>

// This is the number of age classes in the model, 9 age classes currently: 0-9,10-19, etc up to 80+
#define NUMAC 9


// number of stages for the exposed class (E-class) - these are individuals who have recently become exposed (and thus infected)
// but they are still in their incubation period, they are not showing symptoms, and they are likely not infectious,
// except maybe in the late stages of the incubation period
//
// it is currently set to 6 to give the incubation period a coefficiant of variation of 1/sqrt(6) which is about 0.4
#define NUME 6
#define STARTE 9    // this is the start-index of the exposed (E) individuals;
                    // this has to be NUMAC, because you have the 9 S-classes in the index before you start indexing the E-classes

                    
// number of stages for the asymptomatic class (A-class); these are individuals who are either (1) fully asymptomatic, or (2) have symptoms that
// are mild enough that they would not report to a doctor/clinic under normal circumstances
#define NUMA 4
#define STARTA 63   // this is the start-index of for the asymptomatic individuals


// number of infected classes (I-classes) - these individuals are symptomatic and infectious (when you transition from E6 to I1, that
// transition coincides with the arrival of symptoms
// - these individuals are not hospitalized
// - this class includes individuals that have no dyspnea, no hypoxia, and no other underlying conditions that would cause a hospitalization
// - this class includes individuals that have dyspnea, no hypoxia, and no other underlying conditions that would cause a hospitalization; these individuals would be 
//   treated as out-patient at a hospital or clinic
#define NUMI 4
#define STARTI 99   // this is the start-index of for the infected individuals
                    // this has to be NUMAC + NUMAC*NUME + NUMAC*NUMA because the order of the indexing goes S, then E, then A, then I

                    
// number of Hospitalized individuals in Acute phase (i.e. not recovering and getting better, but just arrived and getting worse)
#define NUMHA 4
#define STARTHA 135


// number of individuals in Critical Care (ICU) in acute phase (i.e. not recovering)
#define STARTCA 171                    
                    

// V = individuals on mechanical ventilation; NUMV is the number of stages in class V
#define NUMV 6
#define STARTV 180


// number of individuals in Critical Care (ICU) in recovering phase (post ventilation)
#define STARTCR 234                    


// number of hospitalized individuals in recovering phase (post ICU)
#define STARTHR 243                    


// D = dead individuals who died at HOME or at their primary place of residence (e.g. an elderly care facility); there will be NUMAC dead classes.  
#define STARTD 252

// DHOSP = dead individuals who died at a HOSPITAL; there will be NUMAC dead classes.  
#define STARTDHOSP 261

                    
// there will be NUMAC recovered classes.  These are individuals who recovered at home or at their primary residence.
// NOTE that the start index of the R class was changed from 261 to 270 on May 20 2020, as the STARTDHOSP class was inserted above
#define STARTR 270

// RHOSP = these are individuals who were discharged from a hospital
#define STARTRHOSP 279


// ### CLASSES THAT KEEP TRACK OF CUMULATIVE SYMPTOMATIC INCIDENCE
// number of J-classes - these are the cumulative symptomatic incidence classes
// NOTE that the start index of the J classes was changed from 270 to 288 on May 20 2020
#define STARTJ 288

// ### CLASSES THAT KEEP TRACK OF CUMULATIVE HOSPITALIZATION INCIDENCE
// these are the K-classes 
// NOTE that the start index of the K classes was changed from 279 to 297 on May 20 2020
#define STARTK 297


// August 24, 2020
// CLASSES TO COUNT TESTED INDIVIDUALS
#define STARTJTEST 306


// CLASSES FOR VACCINATED INDIVIDUALS
// low-efficacy vaccine
// 2020-11-27: change NUMZ_1 from 4 to 24
#define STARTZ_1 315  // probably more efficient to move Z closer to S and E (spatial locality)
#define NUMZ_1 24 // 4 
// high-efficacy vaccines
#define STARTZ_2 531 // 351 
#define NUMZ_2 4


// patch R classes with 3 more stages
#define STARTRR 567 // 387 
#define NUMRR 3

#define STARTRRHOSP 594 // 414
#define NUMRRHOSP 3

// cumulative vaccinees Z_1
#define STARTJZ_1 621
// cumulative vaccinees Z_2
#define STARTJZ_2 630

// #### NOTE:
// whenever DIMENSION changes, please change the Q constant in rkf.cpp accordingly
// check the calculation for population sum in generate_trajectories.cpp (i.e. G_B_CHECKPOP_MODE)
#define DIMENSION 639 // 441 
// number of equations 9 S classes, 54 E-classes, 36 A-classes, 36 I-classes, 9 R-classes, 9 J-classes, 
//                     9 JTEST-classes, 36 Z_LOW-classes, 36 Z_HIGH-classes,
//                     27 RR-classes, 27 RRHOSP-classes (patch for R-classes)
//
// NOTE NOTE NOTE TOTAL NUMBER OF EQUATIONS IS 441


//#define HAVE_INLINE 	// this will make GSL vector operations faster

using namespace std;

//
//
// this function contains the ode system
//
int func(double t, const double y[], double f[], void *params);



//void* jac;	// do this for C-compilation
//
// for C++ compilation we are replacing the void* declaration above with
// the inline dummy declaration below
inline int jac(double a1, const double* a2, double* a3, double* a4, void* a5)
{
    return 0;	
};



inline double agesum( double* p, int numac=NUMAC )
{
    double sum=0;
    for( int i=0; i < numac; i++ ) sum += p[i];
    return sum;
}


inline double Min(const double &a, const double &b){
    return( a < b ? a : b); 
}

inline double Max(const double &a, const double &b){
    return( a > b ? a : b); 
}

inline double Get_Efficacy(const double &time, const double &half_life, const double &slope){
    return( pow(half_life, slope) / (pow(half_life, slope) + pow(time, slope)) );
}

inline double Get_Rel_Sus_Z(const double &force_infection, const double &time, const double &half_life, const double &slope){
    // double tmp = force_infection * time;
    // double e_ntmp = exp(-tmp);
    // double half_life_slope = pow(half_life, slope);
    // return log( ( half_life_slope / (half_life_slope + pow(time, slope)) ) * (1.0 - e_ntmp) + e_ntmp ) / tmp;
    // fprintf(stderr, "efficacy %f\t exp(-force_infection * time) %f\t force_infection %f\n", Get_Efficacy(time, half_life, slope), exp(-force_infection * time), force_infection);
    return -log( ( pow(half_life, slope) / (pow(half_life, slope) + pow(time, slope)) ) * (1.0 - exp(-force_infection * time )) + exp(-force_infection * time ) ) / (force_infection * (time + 1.0e-30));
    // return 0.0;
}

inline double Get_Rel_Sus_Z_pow(const double &force_infection, const double &time, const double &half_life_pow_slope,  const double &slope){
    return -log( ( half_life_pow_slope / (half_life_pow_slope + pow(time, slope)) ) * (1.0 - exp(-force_infection * time )) + exp(-force_infection * time ) ) / (force_infection * (time + 1.0e-30));
}



#endif // ESSENTIALS
