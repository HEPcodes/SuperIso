#include "include.h"


float delta0m(float C0[],float C0_spec[],float C1[],float C1_spec[],struct parameters* param,float mub,float muspec, float lambda_h)
/* computes the isospin asymmetry of B -> K* gamma */
/* C0 and C1: respectively LO and NLO contributions of the Wilson coefficients at scale mu=O(mb) */
/* C0_spec and C1_spec: respectively LO and NLO contributions of the Wilson coefficients at scale mu_spec=O(sqrt(lambda_h*mb)) */
{
	float C[9],C_spec[9];
	int ie;
	
	float pi=4.*atan(1.);
	float T1=0.3;  
	float lambda_B=0.35;
	float f_K=0.226;
	float f_K_orth=0.175;
	float m_B=5.28;
	float m_K=0.892;
	float F_orth=1.21;
	complex float G_orth=2.82+I*0.81;
	complex float H_orth= 2.32+I*0.50;

	float alpha_s_muspec;
	float alpha_s_mub;

	float alpha_s_mb=alpha_s_running(param->mass_b,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);	

 	float mass_b_pole=param->mass_b*(1+alpha_s_mb/pi*(4./3.+alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*((param->mass_u+param->mass_d+param->mass_s+param->mass_c)/param->mass_b))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));

	float mass_b_mub=running_mass(param->mass_b,param->mass_b,mub,param->mass_top_pole,param->mass_b,param->alpha_s_MZ,param->mass_Z);
	alpha_s_mub=alpha_s_running(mub,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z);
	alpha_s_muspec=alpha_s_running(muspec,param->mass_top_pole,mass_b_pole,param->alpha_s_MZ,param->mass_Z);	
	 
	for (ie=1;ie<=8;ie++) 
	{
		C[ie]=C0[ie] + alpha_s_mub/4./pi*C1[ie];
		C_spec[ie]=C0_spec[ie] + alpha_s_mub/4./pi*C1_spec[ie];
	}
	
	float f_B=0.2;
	int nf=5;
	int N=3;
	float mu0=mub;
	float r1=(8./3.*C[3]+4./3.*nf*(C[4]+C[6])-8.*(N*C[6]+C[5]))*F_orth*log(mub/mu0);
	float r2=(-44./3.*C[3]-4./3.*nf*(C[4]+C[6]))*log(mub/mu0);
	complex float lambda_u_lambda_c=0.011;

	float zeta3=1.2020569031595942855;
	float sc=pow(param->mass_c/mass_b_mub,2.);
	
	float a_orth1=0.19;
	float ubarm1_orth=3.7;
	complex float h_spec=4.8+I*1.5;

	float um1_orth = ubarm1_orth - 6.*a_orth1;
	float H8a7= 4./9.*pi*pi*f_B*f_K_orth/T1/m_B/lambda_B*um1_orth;
	
	complex float H1a7= 4./9.*pi*pi*f_B*f_K_orth/T1/m_B/lambda_B*(ubarm1_orth-h_spec);
	
	complex float G8a7= 11./3.-2.*pi*pi/9.+2.*I*pi/3.;
	
	complex float G1a7= (-833./162. - 20.*I*pi/27. + 8.*pi*pi/9.*pow(sc,3./2.) + 2./9.*(48.+30.*I*pi-5.*pi*pi-2.*I*pi*pi*pi-36.*zeta3+(36.+6.*I*pi-9.*pi*pi)*log(sc) + (3.+6.*I*pi)*pow(log(sc),2.) + pow(log(sc),3.))*sc + 2./9.*(18.+2.*pi*pi-2.*I*pi*pi*pi+(12.-6.*pi*pi)*log(sc)+6.*I*pi*pow(log(sc),2.)+pow(log(sc),3.))*sc*sc + 1./27.*(-9.+112.*I*pi-14.*pi*pi+(182.-48.*I*pi)*log(sc)-126.*pow(log(sc),2.))*sc*sc*sc);

	complex float a7c=C[7] + alpha_s_mub/4./pi*4./3.*(C[1]*G1a7+C[8]*G8a7) + alpha_s_muspec/4./pi*4./3.*(C_spec[8]*H8a7+C_spec[1]*H1a7);

	float rho=0.;
	float phi=0.;
	
	complex float  X_orth=-3.91+3.44*log(m_B/lambda_h)*(1.+rho*(cos(phi)+I*sin(phi)));

	complex float K1=-(C[6]+C[5]/3.)*F_orth + 4./9.*alpha_s_mub/4./pi*(pow(mass_b_mub/m_B,2.)*C[8]*X_orth-C[1]*(2./3.*F_orth-G_orth)+r1); 	

	complex float K2u=lambda_u_lambda_c*(C[1]+C[2]/3.)+(C[4]+C[3]/3.)+4./9.*alpha_s_mub/4./pi*(C[1]*(2./3.-H_orth)+r2);
	
	complex float K2d=(C[4]+C[3]/3.)+4./9.*alpha_s_mub/4./pi*(C[1]*(2./3.-H_orth)+r2); 
		
	complex float b_d=12.*pi*pi*f_B*(-1./3.)/mass_b_mub/T1/a7c*(f_K_orth/mass_b_mub*K1+f_K*m_K/6./lambda_B/m_B*K2d);
	
	complex float b_u=12.*pi*pi*f_B*(2./3.)/mass_b_mub/T1/a7c*(f_K_orth/mass_b_mub*K1+f_K*m_K/6./lambda_B/m_B*K2u);
	
#ifdef DEBUG
	for(ie=1;ie<=8;ie++) printf("C%d_0=%f\t C%d_1=%f\t C%d_0_spec=%f\n",ie,C0[ie],ie,C1[ie],ie,C0_spec[ie]);
	printf("a7=%f+I*%f\n",creal(a7c),cimag(a7c));
#endif

	return creal(b_d-b_u);
}

/*--------------------------------------------------------------------*/

float delta0_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the isospin asymmetry */
{
	float C0[9],C0spec[9],C1[9],C1spec[9];
	struct parameters param;
	
	float pi=4.*atan(1.);
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	float alpha_s_mb=alpha_s_running(param.mass_b,param.mass_top_pole,param.mass_b,param.alpha_s_MZ,param.mass_Z);	

	float mass_b_pole=param.mass_b*(1+alpha_s_mb/pi*(4./3.+alpha_s_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*(1.61/4.62))+alpha_s_mb/pi*(190.595-4.*(26.655-4.*0.6527)))));
	
	float mu_b=mass_b_pole;
	
	float lambda_h=0.5;
	float muspec=sqrt(lambda_h*mu_b);

	C_calculator(C0,C1,mu_b,&param);
	C_calculator(C0spec,C1spec,muspec,&param);

	return delta0m(C0,C0spec,C1,C1spec,&param,mu_b,muspec,lambda_h);
}
