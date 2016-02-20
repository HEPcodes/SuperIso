#include "include.h"

double Dstaunu(struct parameters* param)
/* computes the branching ratio of Ds -> tau nu */
{
	double Vcs=0.97334;

#ifdef SMONLY
	return 
param->m_Ds/8./pi*pow(param->Gfermi*Vcs*param->mass_tau_pole*param->f_Ds*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_Ds/param->m_Ds),2.)*param->life_Ds/hbar;
#endif
	
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*
H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->McR_Q*param->McR_Q/param->mass_gluino/param->mass_gluino);

	
	return param->m_Ds/8./pi*pow(param->Gfermi*Vcs*param->mass_tau_pole*param->f_Ds*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_Ds/param->m_Ds)		*(1.+param->m_Ds*param->m_Ds/param->mass_H/param->mass_H*(param->mass_c-param->mass_s*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta))/(param->mass_s+param->mass_c)),2.)*param->life_Ds/hbar;
}

/*--------------------------------------------------------------------*/

double Dstaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Ds -> tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Dstaunu(&param);
}

/*--------------------------------------------------------------------*/

double Dsmunu(struct parameters* param)
/* computes the branching ratio of Ds -> mu nu */
{
	double Vcs=0.97334;

#ifdef SMONLY
	return 
param->m_Ds/8./pi*pow(param->Gfermi*Vcs*param->mass_mu*param->f_Ds*(1.-param->mass_mu*param->mass_mu/param->m_Ds/param->m_Ds),2.)*param->life_Ds/hbar;
#endif
	
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*
H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->McR_Q*param->McR_Q/param->mass_gluino/param->mass_gluino);

	
	return param->m_Ds/8./pi*pow(param->Gfermi*Vcs*param->mass_mu*param->f_Ds*(1.-param->mass_mu*param->mass_mu/param->m_Ds/param->m_Ds)		*(1.+param->m_Ds*param->m_Ds/param->mass_H/param->mass_H*(param->mass_c-param->mass_s*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta))/(param->mass_s+param->mass_c)),2.)*param->life_Ds/hbar;
}

/*--------------------------------------------------------------------*/

double Dsmunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Ds -> mu nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Dsmunu(&param);
}
