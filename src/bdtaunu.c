#include "include.h"


float GBDlnu(float w)
{
	float G1=1.03;
	float rho2=1.17;
	float zw=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
	
	return G1*(1.-8.*rho2*zw+(51.*rho2-10.)*zw*zw-(252.*rho2-84.)*zw*zw*zw);
}

/*--------------------------------------------------------------------*/

float tBDlnu(float w, float m_B, float m_D)
{
	return m_B*m_B+m_D*m_D-2.*w*m_D*m_B;
}

/*--------------------------------------------------------------------*/

float rhoV(float w, float ml, float m_B, float m_D)
{
	return 4.*pow(1.+m_D/m_B,2.)*pow(m_D/m_B,3.)*pow(w*w-1.,1.5)*pow(1.-ml*ml/tBDlnu(w,m_B,m_D),2.)*(1.+ml*ml/2./tBDlnu(w,m_B,m_D))*pow(GBDlnu(w),2.);
}

/*--------------------------------------------------------------------*/

float rhoS(float w, float ml, float m_B, float m_D)
{
	float Deltaw=0.46;
	
	return 1.5*m_B*m_B/tBDlnu(w,m_B,m_D)/(1.+ml*ml/2./tBDlnu(w,m_B,m_D))*(1.+w)/(1.-w)*Deltaw*Deltaw;
}

/*--------------------------------------------------------------------*/

float dGammaBDlnu_dw(float w, float ml, struct parameters* param)
{
	float Vcb=4.17e-2;

#ifdef SMONLY
	return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(param->m_B,5.)/192./pow(pi,3.)*rhoV(w,ml,param->m_B,param->m_D)*(1.-ml*ml/param->m_B/param->m_B*rhoS(w,ml,param->m_B,param->m_D));
#endif
	
	return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(param->m_B,5.)/192./pow(pi,3.)*rhoV(w,ml,param->m_B,param->m_D)*	(1.-ml*ml/param->m_B/param->m_B*pow(1.-tBDlnu(w,param->m_B,param->m_D)/(param->mass_b-param->mass_c)*param->mass_b/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta),2.)*rhoS(w,ml,param->m_B,param->m_D));
}

/*--------------------------------------------------------------------*/

float GammaBDlnu(float ml, struct parameters* param)
{
	int ie;
	int nmax=100.;
	float Gamma=0.;
	float w;
	float wmin=1.;
	float wmax=(1.+param->m_D*param->m_D/param->m_B/param->m_B-ml*ml/param->m_B/param->m_B)/2./(param->m_D/param->m_B);
	
	for(ie=1;ie<=nmax;ie++)
	{
		w=wmin+(wmax-wmin)*ie/nmax;
		Gamma+=dGammaBDlnu_dw(w,ml,param);
	}
	Gamma*=(wmax-wmin)/nmax;

	return Gamma;
}

/*--------------------------------------------------------------------*/

float BDtaunu(struct parameters* param)
/* computes the branching ratio of B-> D0 tau nu */
{
	return param->life_B/hbar*GammaBDlnu(param->mass_tau_pole,param);
}

/*--------------------------------------------------------------------*/

float BDtaunu_BDenu(struct parameters* param)
/* computes the ratio BR(B-> D0 tau nu)/BR(B-> D0 e nu) */
{

	return GammaBDlnu(param->mass_tau_pole,param)/GammaBDlnu(param->mass_e,param);
	
}

/*--------------------------------------------------------------------*/

float BDtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDtaunu(&param);
}

/*--------------------------------------------------------------------*/

float BDtaunu_BDenu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu)/BR(B-> D0 e nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDtaunu_BDenu(&param);
}

