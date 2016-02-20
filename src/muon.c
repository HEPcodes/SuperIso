#include "include.h"


/*--------------------------------------------------------------------*/

float F1N(float x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(1.-6.*x+3.*x*x+2.*pow(x,3.)-6.*x*x*log(x));
}

/*--------------------------------------------------------------------*/

float F2N(float x)
{
	if(x==1.) return 1.;
	return 3./pow(1.-x,3.)*(1.-x*x+2.*x*log(x));
}

/*--------------------------------------------------------------------*/

float F1C(float x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(2.+3.*x-6.*x*x+pow(x,3.)+6.*x*log(x));
}

/*--------------------------------------------------------------------*/

float F2C(float x)
{
	if(x==1.) return 1.;
	return -3./2./pow(1.-x,3.)*(3.-4.*x+x*x+2.*log(x));
}

/*---------------------------*/
/* Calculation of muon (g-2) */
/*---------------------------*/

float muon_gm2(struct parameters* param)
{

#ifdef SMONLY
	return 0.;
#endif

#ifdef SM_ChargedHiggs
	return 0.;
#endif	

	float mass_smu[3],mass_charg[3];
	float n_L[5][3],n_R[5][3],xim[5][3],xk[3],c_L[3],c_R[3],X[3][3];
	int ie,me,ke;

	mass_charg[1]=param->mass_cha1;
	mass_charg[2]=param->mass_cha2;

	float M2smu11=pow(param->mass_mul,2.);
	float M2smu22=pow(param->mass_mur,2.);
	float M2smu12=(param->A_mu-param->mu_Q*param->tan_beta)*param->mass_mu;

mass_smu[1]=sqrt((M2smu11+M2smu22-sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);
	mass_smu[2]=sqrt((M2smu11+M2smu22+sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);

	X[1][1]=cos(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[1][2]=sin(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[2][1]=-X[1][2];
	X[2][2]=X[1][1];

	
	float ymu=param->g2*param->mass_mu/sqrt(2.)/param->mass_W/cos(atan(param->tan_beta));
	
	for(ie=1;ie<=4;ie++) for(me=1;me<=2;me++)
	{
		n_L[ie][me]=1./sqrt(2.)*(param->gp*param->neut_mix[ie][1]+param->g2*param->neut_mix[ie][2])*X[me][1]
		-ymu*param->neut_mix[ie][3]*X[me][2];
		n_R[ie][me]=sqrt(2.)*param->gp*param->neut_mix[ie][1]*X[me][2]+ymu*param->neut_mix[ie][3]*X[me][1];
		xim[ie][me]=pow(param->mass_neut[ie]/mass_smu[me],2.);
	}
	
	for(ke=1;ke<=2;ke++)
	{
		c_L[ke]=-param->g2*param->charg_Vmix[ke][1];
		c_R[ke]=ymu*param->charg_Umix[ke][2];
		xk[ke]=pow(mass_charg[ke]/param->mass_numl,2.);
	}
	

	float pi=4.*atan(1.);

	float a_neut=0.;
	for(ie=1;ie<=4;ie++) for(me=1;me<=2;me++)
	{
		a_neut+=
		-param->mass_mu/12./pow(mass_smu[me],2.)*(pow(n_L[ie][me],2.)+pow(n_R[ie][me],2.))*F1N(xim[ie][me])
		+param->mass_neut[ie]/3./pow(mass_smu[me],2.)*n_L[ie][me]*n_R[ie][me]*F2N(xim[ie][me]);
	}
	a_neut*=param->mass_mu/16./pi/pi;
	
	
	
	float a_charg=0.;
	for(ke=1;ke<=2;ke++)
	{
		a_charg+=
		param->mass_mu/12./pow(param->mass_numl,2.)*(pow(c_L[ke],2.)+pow(c_R[ke],2.))*F1C(xk[ke])
		+2.*mass_charg[ke]/3./pow(param->mass_numl,2.)*c_L[ke]*c_R[ke]*F2C(xk[ke]);
	}
	a_charg*=param->mass_mu/16./pi/pi;

	return a_neut+a_charg;
}

/*--------------------------------------------------------------------*/

float muon_gm2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating muon (g-2) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return muon_gm2(&param);
}

