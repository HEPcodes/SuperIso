#include "include.h"

/*--------------------------------------------------------------------*/

float alpha_s_running(float Q, float mtop, float mbot, float alphas_MZ, float MZ)
/* computes the QCD coupling constant alpha_s at the energy Q */
/* valid for at least 4 active flavors */

{
	float beta0,beta1,beta2,alpha_s_running,Lambda4,Lambda5,Lambda6,pi,Lambda_min,Lambda_max,Lambda_moy,alphas_min,alphas_max,alphas_moy;
	int nf;

	pi=4.*atan(1.);

	nf=5;
	beta0 = 11.-2./3.*nf;
	beta1=51.-19./3.*nf;
	beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

	Lambda_min=1.e-3;
	Lambda_max=1.;
	alphas_min=0.;

	while(fabs(1.-alphas_min/alphas_MZ)>=1.e-4)
	{
		alphas_min=4.*pi/beta0/log(pow(MZ/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_min,2.)))/log(pow(MZ/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_min,2.)),2.)*(pow(log(log(pow(MZ/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		alphas_max=4.*pi/beta0/log(pow(MZ/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_max,2.)))/log(pow(MZ/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_max,2.)),2.)*(pow(log(log(pow(MZ/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		Lambda_moy=(Lambda_min+Lambda_max)/2.;

		alphas_moy=4.*pi/beta0/log(pow(MZ/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_moy,2.)))/log(pow(MZ/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_moy,2.)),2.)*(pow(log(log(pow(MZ/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		if((alphas_MZ>=alphas_min)&&(alphas_MZ<=alphas_moy))
			Lambda_max=Lambda_moy;
		else Lambda_min=Lambda_moy;
	}

	Lambda5=Lambda_min;

	if((Q<=mtop)&&(Q>=mbot))
	/* 5 active flavors */
	{
		alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda5,2.)))/log(pow(Q/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda5,2.)),2.)*(pow(log(log(pow(Q/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
		return alpha_s_running;
	}
	else if((Q>mtop))
	/* 6 active flavors */
	{
		alpha_s_running=4.*pi/beta0/log(pow(mtop/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda5,2.)))/log(pow(mtop/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda5,2.)),2.)*(pow(log(log(pow(mtop/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		nf=6;
		beta0 = 11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while(fabs(1.-alphas_min/alpha_s_running)>=1.e-4)
		{
			alphas_min=4.*pi/beta0/log(pow(mtop/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_min,2.)))/log(pow(mtop/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_min,2.)),2.)*(pow(log(log(pow(mtop/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			alphas_max=4.*pi/beta0/log(pow(mtop/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_max,2.)))/log(pow(mtop/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_max,2.)),2.)*(pow(log(log(pow(mtop/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;

			alphas_moy=4.*pi/beta0/log(pow(mtop/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_moy,2.)))/log(pow(mtop/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_moy,2.)),2.)*(pow(log(log(pow(mtop/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			if((alpha_s_running>=alphas_min)&&(alpha_s_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}

		Lambda6=Lambda_min;

		alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda6,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda6,2.)))/log(pow(Q/Lambda6,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda6,2.)),2.)*(pow(log(log(pow(Q/Lambda6,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		return alpha_s_running;
	}
	else
	/* 4 active flavors */
	{
		alpha_s_running=4.*pi/beta0/log(pow(mbot/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda5,2.)))/log(pow(mbot/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda5,2.)),2.)*(pow(log(log(pow(mbot/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		nf=4;
		beta0 = 11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while(fabs(1.-alphas_min/alpha_s_running)>=1.e-4)
		{
			alphas_min=4.*pi/beta0/log(pow(mbot/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_min,2.)))/log(pow(mbot/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_min,2.)),2.)*(pow(log(log(pow(mbot/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			alphas_max=4.*pi/beta0/log(pow(mbot/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_max,2.)))/log(pow(mbot/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_max,2.)),2.)*(pow(log(log(pow(mbot/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;

			alphas_moy=4.*pi/beta0/log(pow(mbot/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_moy,2.)))/log(pow(mbot/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_moy,2.)),2.)*(pow(log(log(pow(mbot/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			if((alpha_s_running>=alphas_min)&&(alpha_s_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}
		Lambda4=Lambda_min;

		alpha_s_running=4.*pi/beta0/log(pow(Q/Lambda4,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda4,2.)))/log(pow(Q/Lambda4,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda4,2.)),2.)*(pow(log(log(pow(Q/Lambda4,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		return alpha_s_running;
	}
}
