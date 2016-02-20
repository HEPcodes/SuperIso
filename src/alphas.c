#include "include.h"

float alpha_s_running(float Q, float mtop, float mbot, float alphas_MZ, float MZ)
/* computes the QCD coupling constant alpha_s at the energy Q */
/* based on the Particle Data Group reviews */
/* valid for at least 4 active flavors */

{
	float beta0,beta1,beta2,alpha_s_running,Lambda4,Lambda5,Lambda6,pi,Lambda_min,Lambda_max,Lambda_moy,alphas_min,alphas_max,alphas_moy;
	int nf;

	pi=4.*atan(1.);

	nf=5;
	beta0 = 11.-2./3.*nf;
	beta1=51.-19./3.*nf;
	beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

	Lambda_min=0.1;
	Lambda_max=0.3;
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

		Lambda_min=0.01;
		Lambda_max=0.2;
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

		Lambda_min=0.2;
		Lambda_max=0.4;
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


float running_mass(float quark_mass, float Qinit, float Qfin,  float mtop, float mbot, float alphas_MZ, float MZ)

/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
/* valid for at least 4 active flavors */

{

	float alphas_Qinit,alphas_Qfin,running_mass;
	float beta0,beta1,beta2,gamma0,gamma1,gamma2;
	int nf;

	float pi=4.*atan(1.);
	float zeta3=1.2020569031595942855;

	alphas_Qinit=alpha_s_running(Qinit,mtop,mbot,alphas_MZ,MZ);
	
	float R_Qinit,R_Qfin;
		
	if(Qinit <= mbot) 
	/* 4 active flavors at Qinit */
	{
		nf=4.;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin <= mbot) 
		{
		/* if Qinit and Qfin are in the same range, just calculate R_Qfin */
		/* 4 active flavors at Qinit and Qfin */
			nf=4.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else if(Qfin <= mtop)
		{
		/* if Qinit and Qfin are NOT in the same range, evolve first the running mass towards the range limit(s), then from the limit(s) towards Qfin */
		/* 4 active flavors at Qinit, 5 active flavors at Qfin */
			nf=4.;

			alphas_Qfin=alpha_s_running(mbot,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;


			alphas_Qinit=alphas_Qfin;

			nf=5.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else
		/* 4 active flavors at Qinit, 6 active flavors at Qfin */
		{
			nf=4.;

			alphas_Qfin=alpha_s_running(mbot,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;


			alphas_Qinit=alphas_Qfin;

			nf=5.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5.;

			alphas_Qfin=alpha_s_running(mtop,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
	
	
			alphas_Qinit=alphas_Qfin;

			nf=6.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=6.;
	
			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
	}
	else if(Qinit <= mtop) 
	/* 5 active flavors at Qinit */
	{
		nf=5.;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin <= mbot) 
		/* 5 active flavors at Qinit, 4 active flavors at Qfin */
		{
			nf=5.;

			alphas_Qfin=alpha_s_running(mbot,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;

			alphas_Qinit=alphas_Qfin;

			nf=4.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=4.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else if(Qfin <= mtop)
		/* 5 active flavors at Qinit and Qfin */
		{
			nf=5.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else
		/* 5 active flavors at Qinit, 6 active flavors at Qfin */
		{
			nf=5.;

			alphas_Qfin=alpha_s_running(mtop,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=6.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=6.;
	
			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}	
	}
	else
	/* 6 active flavors at Qinit */
	{
		nf=6.;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin >= mtop) 
		/* 6 active flavors at Qinit and Qfin */
		{
			nf=6.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else if(Qfin >= mbot)
		/* 6 active flavors at Qinit, 5 active flavors at Qfin */
		{
			nf=6.;

			alphas_Qfin=alpha_s_running(mtop,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=5.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else
		/* 6 active flavors at Qinit, 4 active flavors at Qfin */
		{
			nf=6.;

			alphas_Qfin=alpha_s_running(mtop,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=5.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			alphas_Qfin=alpha_s_running(mbot,mtop,mbot,alphas_MZ,MZ);

			nf=5.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=4.;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=4.;

			alphas_Qfin=alpha_s_running(Qfin,mtop,mbot,alphas_MZ,MZ);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
	}
	
	return running_mass;
}
