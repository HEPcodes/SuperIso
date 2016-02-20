#include "include.h"

float Btaunu(struct parameters* param)
/* computes the branching ratio of B-> tau nu */
{
	float Vub=3.95e-3;

#ifdef SMONLY
	return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*param->f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_B/param->m_B),2.);
#endif
	
	return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*param->f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_B/param->m_B)	
	*(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta)),2.);
}

/*--------------------------------------------------------------------*/

float Btaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Btaunu(&param);
}

/*--------------------------------------------------------------------*/

float RBtaunu(struct parameters* param)
/* computes the ratio of BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{

#ifdef SMONLY	
	return 1.;
#endif
	
	return 
	pow(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta),2.);
}

/*--------------------------------------------------------------------*/

float RBtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return RBtaunu(&param);
}
