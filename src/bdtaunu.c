#include "include.h"


float GBDlnu(float w)  
{
	float G1=1.03;
	float rho2=1.12; 
	float zw=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
	
	return G1*(1.-8.*rho2*zw+(51.*rho2-10.)*zw*zw-(252.*rho2-84.)*zw*zw*zw);
}

/*--------------------------------------------------------------------*/
float tBDlnu(float w)
{
	float m_D=1.87;
	float m_B=5.28;
	return m_B*m_B+m_D*m_D-2.*w*m_D*m_B;
}

/*--------------------------------------------------------------------*/
float rhoV(float w, float ml) 
{
	float m_D=1.87;
	float m_B=5.28;
	
	return 4.*pow(1.+m_D/m_B,2.)*pow(m_D/m_B,3.)*pow(w*w-1.,1.5)*pow(1.-ml*ml/tBDlnu(w),2.)*(1.+ml*ml/2./tBDlnu(w))*pow(GBDlnu(w),2.);
}

/*--------------------------------------------------------------------*/
float rhoS(float w, float ml) 
{
	float m_B=5.28;
	float Deltaw=0.46;
	
	return 1.5*m_B*m_B/tBDlnu(w)/(1.+ml*ml/2./tBDlnu(w))*(1.+w)/(1.-w)*Deltaw*Deltaw;
}

/*--------------------------------------------------------------------*/
float dGammaBDlnu_dw(float w, float ml, struct parameters* param) 
{
	float Vcb=4.17e-2;
	float m_B=5.28;

#ifdef SMONLY
	return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(m_B,5.)/192./pow(pi,3.)*rhoV(w,ml)*(1.-ml*ml/m_B/m_B*rhoS(w,ml));
#endif
	
	return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(m_B,5.)/192./pow(pi,3.)*rhoV(w,ml)*	(1.-ml*ml/m_B/m_B*pow(1.-tBDlnu(w)/(param->mass_b-param->mass_c)*param->mass_b/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_b(param)*param->tan_beta),2.)*rhoS(w,ml));
}

/*--------------------------------------------------------------------*/
float GammaBDlnu(float ml, struct parameters* param)
{
	int ie;
	int nmax=100.;
	float Gamma=0.;
	float m_D=1.87;
	float m_B=5.28;
	float w;
	float wmin=1.;
	float wmax=(1.+m_D*m_D/m_B/m_B-ml*ml/m_B/m_B)/2./(m_D/m_B); 
	
	for(ie=1;ie<=nmax;ie++)
	{
		w=wmin+(wmax-wmin)*ie/nmax;
		Gamma+=dGammaBDlnu_dw(w,ml,param);
	}
	Gamma*=(wmax-wmin)/nmax;

	return Gamma;
}

/*--------------------------------------------------------------------*/
float Bbdtaunu(struct parameters* param)
{
	
	float life_B=1.638e-12; 
	float hbar=6.58211889e-25; 

	return life_B/hbar*GammaBDlnu(param->mass_tau_pole,param);
	
}

/*--------------------------------------------------------------------*/
float Bbdtaunu_Bbdenu(struct parameters* param)
{

	return GammaBDlnu(param->mass_tau_pole,param)/GammaBDlnu(param->mass_e,param);
	
}

/*--------------------------------------------------------------------*/

float Bbdtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Bbdtaunu(&param);
}

/*--------------------------------------------------------------------*/

float Bbdtaunu_Bbdenu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu)/BR(B-> D0 e nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Bbdtaunu_Bbdenu(&param);
}

