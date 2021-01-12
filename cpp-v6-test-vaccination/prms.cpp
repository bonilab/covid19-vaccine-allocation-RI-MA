//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "essentials.h"
#include "assert.h"
#include "prms.h"


// constructor
prms::prms()
{
    v.insert( v.begin(), num_params, 0.0 );
    assert( v.size()==num_params );
    
    v_rel_susc.insert( v_rel_susc.begin(), NUMAC, 0.0 );
    // v_mixing_level.insert( v_mixing_level.begin(), NUMAC, 0.0 );
    // v_mixing_level_postld.insert( v_mixing_level_postld.begin(), NUMAC, 0.0 );
    v_prob_E_A.insert( v_prob_E_A.begin(), NUMAC, 0.0 );
    v_prob_I2_H.insert(  v_prob_I2_H.begin(),  NUMAC, 0.0 );
    v_prob_I4_D.insert(  v_prob_I4_D.begin(),  NUMAC, 0.0 );
    v_prob_HA4_D.insert( v_prob_HA4_D.begin(), NUMAC, 0.0 );
    v_prob_HA_CA.insert( v_prob_HA_CA.begin(), NUMAC, 0.0 );
    v_prob_V_D.insert( v_prob_V_D.begin(), NUMAC, 0.0 );
    v_prob_CA_D.insert( v_prob_CA_D.begin(), NUMAC, 0.0 );
    v_prob_CR_D.insert( v_prob_CR_D.begin(), NUMAC, 0.0 );
    v_prob_CA_V.insert( v_prob_CA_V.begin(), NUMAC, 0.0 );
    
    // v_prob_S_Z_1.insert( v_prob_S_Z_1.begin(), NUMAC, 0.0 );
    v_number_to_Z_1_phase1.insert( v_number_to_Z_1_phase1.begin(), NUMAC, 0.0 );
    v_already_moved_to_Z_1_phase1.insert( v_already_moved_to_Z_1_phase1.begin(), NUMAC, 0.0 );
    v_efficacy_Z_1.insert( v_efficacy_Z_1.begin(), NUMZ_1, 0.0 );
    v_rel_sus_Z_1.insert( v_rel_sus_Z_1.begin(), NUMAC, std::vector<double>{} );
    v_rel_eff_Z_1.insert( v_rel_eff_Z_1.begin(), NUMAC, 1.0 );

    // time-varying vaccine rollout campaign 
    v_begin_days_phase1_vac1.clear();
    v_end_days_phase1_vac1.clear();
    v_dpd_phase1_vac1.clear();
    v_vac_ratios_phase1_vac1.insert( v_vac_ratios_phase1_vac1.begin(), NUMAC,  std::vector<double>{} );
    v_vac_fracs_phase1_vac1.insert( v_vac_fracs_phase1_vac1.begin(), NUMAC,  std::vector<double>{} );
    index_current_day_phase1_vac1 = -1;
    
    // v_prob_S_Z_2.insert( v_prob_S_Z_2.begin(), NUMAC, 0.0 );
    // v_number_S_Z_2_phase1.insert( v_number_S_Z_2_phase1.begin(), NUMAC, 0.0 );
    // v_efficacy_Z_2.insert( v_efficacy_Z_2.begin(), NUMZ_2, 0.0 );
    // v_rel_sus_Z_2.insert( v_rel_sus_Z_2.begin(), NUMAC, std::vector<double>{}  );
    // v_rel_eff_Z_2.insert( v_rel_eff_Z_2.begin(), NUMAC, 1.0 );

    v_pop_frac.insert( v_pop_frac.begin(), NUMAC, 0.0 );
    
    assert( v_rel_susc.size()==NUMAC );  
    // assert( v_mixing_level.size()==NUMAC );  
    // assert( v_mixing_level_postld.size()==NUMAC );  
    assert( v_prob_E_A.size()==NUMAC );
    assert( v_prob_I2_H.size() ==NUMAC );
    assert( v_prob_I4_D.size() ==NUMAC );
    assert( v_prob_HA4_D.size()==NUMAC );
    assert( v_prob_HA_CA.size()==NUMAC );
    assert( v_prob_V_D.size()==NUMAC );
    assert( v_prob_CA_D.size()==NUMAC );
    assert( v_prob_CR_D.size()==NUMAC );
    assert( v_prob_CA_V.size()==NUMAC );
    
    // assert( v_prob_S_Z_1.size()==NUMAC );
    assert( v_number_to_Z_1_phase1.size()==NUMAC );
    assert( v_already_moved_to_Z_1_phase1.size()==NUMAC );
    assert( v_efficacy_Z_1.size()==NUMZ_1 );
    assert( v_rel_sus_Z_1.size()==NUMAC );
    for (size_t i = 0; i < NUMAC; i++){  assert(v_rel_sus_Z_1[i].size()==0);   }
    assert( v_rel_eff_Z_1.size()==NUMAC );

    assert( v_begin_days_phase1_vac1.size()==0 );
    assert( v_end_days_phase1_vac1.size()==0 );
    assert( v_dpd_phase1_vac1.size()==0 );
    assert( v_vac_ratios_phase1_vac1.size()==NUMAC );
    for (size_t i = 0; i < NUMAC; i++){  assert(v_vac_ratios_phase1_vac1[i].size()==0);   }
    assert( v_vac_fracs_phase1_vac1.size()==NUMAC );
    for (size_t i = 0; i < NUMAC; i++){  assert(v_vac_fracs_phase1_vac1[i].size()==0);   }
    
    
    // assert( v_prob_S_Z_2.size()==NUMAC );
    // assert( v_number_S_Z_2_phase1.size()==NUMAC );
    // assert( v_efficacy_Z_2.size()==NUMZ_2 );
    // assert( v_rel_sus_Z_2.size()==NUMAC );
    // for (size_t i = 0; i < NUMAC; i++){  assert(v_rel_sus_Z_2[i].size()==0);   }
    // assert( v_rel_eff_Z_2.size()==NUMAC );

    assert( v_pop_frac.size()==NUMAC );
    
    v_betas.clear();
    v_betatimes.clear();
    assert( v_betas.size() == 0 );
    assert( v_betatimes.size() == 0 );
    
    earlymarch_highhosp_period = false;
    earlymarch_highhosp_factor = 1.0;
    earlymarch_highhosp_endday = -1.0;
    
    index_current_beta=-1;
    
}

// destructor
prms::~prms()
{
}


void prms::assign_new_beta( void )
{
    assert( v_betas.size() == v_betatimes.size() );
    assert( v_betas.size() > 0 );
    
    // if it's invalid, just set it to zero
    // this is what happens at initialization as well
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        index_current_beta = 0;
    }
    else // if it's valid, increment it if you can
    {
        if( index_current_beta != v_betas.size() - 1 ) // if it's NOT the last element, then increment it
        {
            index_current_beta++;
        }
    }
    
    // assign the possibly new beta value to the main parameter vector "v"
    v[i_beta] = v_betas[ index_current_beta ];
}


double prms::get_new_update_time( void )
{
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        return 1000000.0;
    }
    else // if the index is within range, return the next update time
    {
        if( index_current_beta != v_betas.size() - 1 ) // if the index does NOT point to the last element, then return the next element
        {
            return v_betatimes[ index_current_beta+1 ];
        }
        else // if it is the last element in the array, just return 1000000.0, as in the next update time is a million days from now
        {
            return 1000000.0;
        }
    }
    
    
}

void prms::apply_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) / earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = true;
}

void prms::end_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) * earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = false;
}


int prms::Increase_index_current_day_phase1_vac1(){
    return ( v_begin_days_phase1_vac1.size() > (index_current_day_phase1_vac1 + 1)) ? ++index_current_day_phase1_vac1 : -1;
}
