#include "include.h"


float F_orth(float a1_orth, float a2_orth)
{
	int n1=500;
	int ie;

	float t=0.;
	float int1=6.*t*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1 +=6.*t*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	}
	int1 *= 1./n1;
	
	return int1/3.;
}

/*--------------------------------------------------------------------*/

float X_orth1(float a1_orth, float a2_orth)
{
	int n1=500;
	int ie;

	float t=0.;
	float int1 = 6.*t*(-6.*a1_orth-a2_orth*30.*t)+6.*(t-1.)*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1 += 6.*t*(-6.*a1_orth-a2_orth*30.*t)
		+6.*(t-1.)*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	}
	int1 *= 1./n1;
	
	return int1/3.;
}

/*--------------------------------------------------------------------*/

float X_orth2(float a1_orth, float a2_orth)
{
	return 2.*(1.+3.*a1_orth+6.*a2_orth);
}

/*--------------------------------------------------------------------*/

float complex G(float s, float x)
{	
	int n1=50;
	int ie;

	float t=0.;
	float complex int1=t*(1.-t)*clog(s-t*(1.-t)*x-I*1.e-10);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1 +=t*(1.-t)*clog(s-t*(1.-t)*x-I*1.e-10);
	}
	int1 *= 1./n1;
	
	return -4.*int1;
}

/*--------------------------------------------------------------------*/

float complex G_orth(float s, float a1_orth, float a2_orth)
{
	int n1=100.;
	int ie;

	float t=0.;
	float complex int1=6.*t*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.))*G(s,1.-t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1 +=6.*t*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.))*G(s,1.-t);
	}
	int1 *= 1./n1;
	
	return int1/3.;
}

/*--------------------------------------------------------------------*/

float complex H_orth(float s, float a1_par, float a2_par)
{
	float zeta3A=0.032;
	float zeta3V=0.013;
	float zeta3T=0.024;
	float wA10=-2.1;

	float deltatp=0.16;
	float deltatm=-0.16;
		
	int n1=100;
	int ie;

 	float t=0.;
	float complex int1=((3./4.*(1.+pow(2.*t-1.,2.))+a1_par*3./2.*pow(2.*t-1.,3.)+(3./7.*a2_par+5.*zeta3A)*(3.*pow(2.*t-1.,2.)-1.)
	+(9./122.*a2_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*t-1.,2.)+35.*pow(2.*t-1.,4.))
	+3.*deltatp+3.*deltatm*(2.*t-1.))
	-1./4.*(6.*(1.-2.*t)*(1.+a1_par*(2.*t-1.)+(a2_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*t-1.,2.)-1.))
	+6.*t*(1.-t)*(2.*a1_par*t+(a2_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*t*(2.*t-1.)))
	+18.*deltatp*(1.-2.*t)-12.*deltatm))*G(s,1.-t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1+=
		((3./4.*(1.+pow(2.*t-1.,2.))+a1_par*3./2.*pow(2.*t-1.,3.)+(3./7.*a2_par+5.*zeta3A)*(3.*pow(2.*t-1.,2.)-1.)
		+(9./122.*a2_par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*t-1.,2.)+35.*pow(2.*t-1.,4.))
		+3.*deltatp+3.*deltatm*(2.*t-1.))
		-1./4.*(6.*(1.-2.*t)*(1.+a1_par*(2.*t-1.)+(a2_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*t-1.,2.)-1.))
		+6.*t*(1.-t)*(2.*a1_par*t+(a2_par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*t*(2.*t-1.)))
		+18.*deltatp*(1.-2.*t)-12.*deltatm)
		)*G(s,1.-t);
	}
	int1 *= 1./n1;
	
	return int1;
}

/*--------------------------------------------------------------------*/

float H8_orth(float a1_orth, float a2_orth)
{
	int n1=500;
	int ie;

	float t=0.;
	float int1=6.*(1.-t)*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=1./n1;
		int1 +=6.*(1.-t)*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	}
	int1 *= 1./n1;
	
	return int1;
}

/*--------------------------------------------------------------------*/

float complex h(float u,float s)
{
	return 4.*s/u/u*(CLi2(2./(1.-csqrt((u-4.*s+I*1.e-10)/u)))+CLi2(2./(1.+csqrt((u-4.*s+I*1.e-10)/u))))-2./u;
}

/*--------------------------------------------------------------------*/

float complex H2_orth(float s, float a1_orth, float a2_orth)
{
	int n1=100;
	int ie;

	float t=0.;
	float complex int1=0.;
	
	for(ie=1;ie<=n1-1.;ie++)
	{
		t+=1./n1;
		int1 +=h(1.-t,s)*6.*t*(1.-t)*(1.+3.*a1_orth*(2.*t-1.)+a2_orth*3./2.*(5.*pow(2.*t-1.,2.)-1.));
	}
	int1 *= 1./n1;
	
	return int1;
}

/*--------------------------------------------------------------------*/

float delta0(float C0[],float C0_spec[],float C1[],float C1_spec[],struct parameters* param,float mub,float muspec, float lambda_h)
/* computes the isospin asymmetry of B -> K* gamma */
/* C0 and C1: respectively LO and NLO contributions of the Wilson coefficients at scale mu=O(mb) */
/* C0_spec and C1_spec: respectively LO and NLO contributions of the Wilson coefficients at scale mu_spec=O(sqrt(lambda_h*mb)) */
{
	float C[9],C_spec[9];
	int ie;	
	
	float mass_b_mub=running_mass(param->mass_b,param->mass_b,mub,param->mass_top_pole,param->mass_b,param);
	float alphas_mub=alphas_running(mub,param->mass_top_pole,param->mass_b_pole,param);
	float alphas_muspec=alphas_running(muspec,param->mass_top_pole,param->mass_b_pole,param);	
	
	float T1=0.31;
	float lambda_B=0.46;
	
	float f_K=0.220;
	float f_K_orth=0.185;

	float a1_orth=0.04;
	float a2_orth=0.15;
	
	float a1_par=0.03;
	float a2_par=0.11;

	float eta=alphas_mub /
	alphas_running(sqrt(1.),param->mass_top_pole,param->mass_b_1S,param);

	int nf=5;
	f_K_orth *= pow(eta,4./3./(11.-2./3.*nf));

	lambda_B /= 1.+alphas_mub/3./pi*log(pow(mub,2.))*(1.-2.*1.4);

	a1_orth*=pow(eta,4./(11.-2./3.*nf));
	a2_orth*=pow(eta,4./3.*(1.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	a1_par*=pow(eta,4./3.*(1.-1./3.+2.)*(11.-2./3.*nf));
	a2_par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));
		 
	for (ie=1;ie<=8;ie++) 
	{
		C[ie]=C0[ie] + alphas_mub/4./pi*C1[ie];
		C_spec[ie]=C0_spec[ie] + alphas_muspec/4./pi*C1_spec[ie];
	}
	
	int N=3;
	
	float mu0=mub;
	float r1=(8./3.*C[3]+4./3.*nf*(C[4]+C[6])-8.*(N*C[6]+C[5]))*F_orth(a1_orth,a2_orth)*log(mub/mu0);
	float r2=(-44./3.*C[3]-4./3.*nf*(C[4]+C[6]))*log(mub/mu0);
	
	float lambda_u_lambda_c=0.011;

	float mass_c_mub=running_mass(param->mass_c,param->mass_c,mub,param->mass_top_pole,param->mass_b_pole,param);
 	float sc=pow(mass_c_mub/mass_b_mub,2.);
	
	complex float Horth= H_orth(sc,a1_par,a2_par);
	float Forth=F_orth(a1_orth,a2_orth);
	complex float Gorth=G_orth(sc,a1_orth,a2_orth);

	float H8a7= 4./3./N*pi*pi*param->f_B*f_K_orth/T1/param->m_B/lambda_B*H8_orth(a1_orth,a2_orth);
	
	complex float H2a7= -2./3./N*pi*pi*param->f_B*f_K_orth/T1/param->m_B/lambda_B*H2_orth(sc,a1_orth,a2_orth);
	
	complex float G8a7= -104./27.*log(mub/param->mass_b_1S)+11./3.-2.*pi*pi/9.+2.*I*pi/3.;
	
	complex float G2a7= 8./3.*log(mub/param->mass_b_1S)+(-833./162. - 20.*I*pi/27. + 8.*pi*pi/9.*pow(sc,3./2.) + 2./9.*(48.+30.*I*pi-5.*pi*pi-2.*I*pi*pi*pi-36.*zeta3+(36.+6.*I*pi-9.*pi*pi)*log(sc) + (3.+6.*I*pi)*pow(log(sc),2.) + pow(log(sc),3.))*sc + 2./9.*(18.+2.*pi*pi-2.*I*pi*pi*pi+(12.-6.*pi*pi)*log(sc)+6.*I*pi*pow(log(sc),2.)+pow(log(sc),3.))*sc*sc + 1./27.*(-9.+112.*I*pi-14.*pi*pi+(182.-48.*I*pi)*log(sc)-126.*pow(log(sc),2.))*sc*sc*sc);

	complex float a7c=C[7] + alphas_mub/4./pi*4./3.*(C[2]*G2a7+C[8]*G8a7) + alphas_muspec/4./pi*4./3.*(C_spec[8]*H8a7+C_spec[2]*H2a7);

	float rho=0.;
	float phi=0.;
	
	complex float  Xorth=X_orth1(a1_orth,a2_orth)+X_orth2(a1_orth,a2_orth)*log(param->m_B/lambda_h)*(1.+rho*(cos(phi)+I*sin(phi)));

	complex float K1=-(C[6]+C[5]/3.)*Forth + 4./9.*alphas_mub/4./pi*(pow(mass_b_mub/param->m_B,2.)*C[8]*Xorth-C[2]*((4./3.*log(param->mass_b_1S/mub)+2./3.)*Forth-Gorth)+r1); 	

	complex float K2u=lambda_u_lambda_c*(C[2]+C[1]/3.)+(C[4]+C[3]/3.)+4./9.*alphas_mub/4./pi*(C[2]*(4./3.*log(param->mass_b_1S/mub)+2./3.-Horth)+r2);
	
	complex float K2d=(C[4]+C[3]/3.)+4./9.*alphas_mub/4./pi*(C[2]*(4./3.*log(param->mass_b_1S/mub)+2./3.-Horth)+r2);
		
	complex float b_d=12.*pi*pi*param->f_B*(-1./3.)/mass_b_mub/T1/a7c*(f_K_orth/mass_b_mub*K1+f_K*param->m_Kstar/6./lambda_B/param->m_B*K2d);
	
	complex float b_u=12.*pi*pi*param->f_B*(2./3.)/mass_b_mub/T1/a7c*(f_K_orth/mass_b_mub*K1+f_K*param->m_Kstar/6./lambda_B/param->m_B*K2u);
	
#ifdef DEBUG
	printf("-----------------\n");
	printf("Isospin Asymmetry\n");
	printf("-----------------\n");
	for(ie=1;ie<=8;ie++) printf("C0[%d]=%f\t C1[%d]=%f\t C0spec[%d]=%f\n",ie,C0[ie],ie,C1[ie],ie,C0_spec[ie]);
	printf("a7=%f+I*%f\n\n",creal(a7c),cimag(a7c));
#endif

	return creal(b_d-b_u);
}

/*--------------------------------------------------------------------*/

float delta0_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the isospin asymmetry */
{
	float C0b[9],C0spec[9],C1b[9],C1spec[9],C0w[9],C1w[9],C2w[9];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	float mu_W=2.*param.mass_W;
	
	float mu_b=param.mass_b_1S;
	
	float lambda_h=0.5;
	float mu_spec=sqrt(lambda_h*param.mass_b_1S);
			
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);

	return delta0(C0b,C0spec,C1b,C1spec,&param,mu_b,mu_spec,lambda_h);
}
