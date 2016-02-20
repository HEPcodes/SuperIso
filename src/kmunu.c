#include "include.h"

double Kmunu_pimunu(struct parameters* param)
/* computes the ratio BR(K-> mu nu)/BR(pi-> mu nu) */
{
	double Vus_Vud=0.2321;
	double m_pi=0.1396;
	double fK_fpi=1.189;
	double delta_em=-0.0070;
	double life_pi=2.6033e-8;
	double life_K=1.2380e-8;  

	if(param->SM==1) return param->m_K/m_pi*(1.+delta_em)*pow(Vus_Vud*fK_fpi*(1.-param->mass_mu*param->mass_mu/param->m_K/param->m_K)/(1.-param->mass_mu*param->mass_mu/m_pi/m_pi),2.)*life_K/life_pi;
	
	if(param->THDM_model>0) return param->m_K/m_pi*(1.+delta_em)*pow(Vus_Vud*fK_fpi*(1.-param->mass_mu*param->mass_mu/param->m_K/param->m_K)/(1.-param->mass_mu*param->mass_mu/m_pi/m_pi)*(1.-param->m_K*param->m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->lambda_d[2][2]*param->lambda_l[2][2]),2.)*life_K/life_pi;

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*
H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->MsR_Q*param->MsR_Q/param->mass_gluino/param->mass_gluino);
	
	return param->m_K/m_pi*(1.+delta_em)*pow(Vus_Vud*fK_fpi*(1.-param->mass_mu*param->mass_mu/param->m_K/param->m_K)/(1.-param->mass_mu*param->mass_mu/m_pi/m_pi)*(1.-param->m_K*param->m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta)),2.)*life_K/life_pi;
}

/*--------------------------------------------------------------------*/

double Rl23(struct parameters* param)
/* computes the ratio Rl23 */
{
	if(param->SM==1) return 1.;

	if(param->THDM_model>0) return fabs(1.-param->m_K*param->m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->lambda_d[2][2]*param->lambda_l[2][2]);

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*
H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->MsR_Q*param->MsR_Q/param->mass_gluino/param->mass_gluino);

	return fabs(1.-param->m_K*param->m_K/param->mass_H/param->mass_H*(1.-param->mass_d/param->mass_s)*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta));
}

/*--------------------------------------------------------------------*/

double Kmunu_pimunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(K-> mu nu)/BR(pi-> mu nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Kmunu_pimunu(&param);
}

/*--------------------------------------------------------------------*/

double Rl23_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the ratio Rl23 */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Rl23(&param);
}
