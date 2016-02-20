#include "include.h"

/*---------------------------------------------------------------------*/

double Bsmumu(double C0b[], double C1b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* computes the inclusive branching ratio of Bs -> mu+ mu- */
{
	double alpha_em=1./137.036;
	double VtbVts=cabs(param->Vtb*param->Vts);
	
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	double C10=C0b[10]+alphas_mub/4./pi*C1b[10];
	double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
	double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
	double C10p=Cpb[10];
	double complex CQp1=CQpb[1];
	double complex CQp2=CQpb[2];
	
		
#ifdef DEBUG
	printf("-----------------\n");
	printf("BR(Bs -> mu+ mu-)\n");
	printf("-----------------\n");
	printf("C10=%.5e\t CQ1=%.5e\t CQ2=%.5e\n",C10,creal(CQ1),creal(CQ2));
#endif	
	
	double BRmumu=param->Gfermi*param->Gfermi*alpha_em*alpha_em*pow(param->m_Bs,3.)*param->f_Bs*param->f_Bs*param->life_Bs/hbar/64./pi/pi/pi*VtbVts*VtbVts*sqrt(1.-4.*param->mass_mu*param->mass_mu/param->m_Bs/param->m_Bs)
	*((1.-4.*param->mass_mu*param->mass_mu/param->m_Bs/param->m_Bs)*pow(cabs(CQ1-CQp1),2.) + pow(cabs((CQ2-CQp2)+2.*(C10-C10p)*param->mass_mu/param->m_Bs),2.));
	
	return BRmumu;
}

/*--------------------------------------------------------------------*/

double Bsmumu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Bs-> mu+ mu-) */
{	
	double C0b[11],C1b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;
	
	double mu_W=2.*param.mass_W;
	
	double mu_b=param.mass_b;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return Bsmumu(C0b,C1b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}
