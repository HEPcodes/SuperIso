#include "include.h"


float Bkmunu_Bpimunu(struct parameters* param)
{

	float Vus_Vud=0.2321; 
	float m_K=0.4937;
	float m_pi=0.1396;
	float fK_fpi=1.189;
	float delta_em=-0.0070;
	float life_pi=2.6033e-8; 
	float life_K=1.2385e-8;  

#ifdef SMONLY
	return m_K/m_pi*(1.+delta_em)*pow(Vus_Vud*fK_fpi*(1.-param->mass_mu*param->mass_mu/m_K/m_K)/(1.-param->mass_mu*param->mass_mu/m_pi/m_pi),2.)*life_K/life_pi;
#endif
	
	return m_K/m_pi*(1.+delta_em)*pow(Vus_Vud*fK_fpi*(1.-param->mass_mu*param->mass_mu/m_K/m_K)/(1.-param->mass_mu*param->mass_mu/m_pi/m_pi)	*(1.-m_K*m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->tan_beta*param->tan_beta/(1.+epsilon_b(param)*param->tan_beta)),2.)*life_K/life_pi;
}

/*--------------------------------------------------------------------*/

float Rl23(struct parameters* param)
{
/* Rl23 from arXiv:0801.1817 */

	float m_K=0.4937;

#ifdef SMONLY
	return 1.;
#endif
	
	return fabs(1.-m_K*m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->tan_beta*param->tan_beta/(1.+epsilon_b(param)*param->tan_beta));
}

/*--------------------------------------------------------------------*/

float Bkmunu_Bpimunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(K-> mu nu)/BR(pi-> mu nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Bkmunu_Bpimunu(&param);
}

/*--------------------------------------------------------------------*/

float Rl23_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the ratio Rl23 (Flavianet, arXiv:0801.1817) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Rl23(&param);
}
