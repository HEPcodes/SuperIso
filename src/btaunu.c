#include "include.h"


float btaunu(struct parameters* param)
{

	float Vub=3.55e-3; 
	float m_B=5.28;
	float f_B=0.200;
	float life_B=1.638e-12;
	float hbar=6.58211889e-25; 

#ifdef SMONLY
	return life_B/hbar*m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/m_B/m_B),2.);
#endif
	
	return life_B/hbar*m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/m_B/m_B)	
	*(1.-m_B*m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_b(param)*param->tan_beta)),2.);
}

/*--------------------------------------------------------------------*/

float btaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return btaunu(&param);
}

