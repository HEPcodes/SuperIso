#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int main()
{ 
	struct parameters param;
		
	Init_param(&param);
	slha_adjust(&param);
	param.SM=1;

	double C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11];
	double complex CQpb[3],CQ0b[3],CQ1b[3];
	double obs[2];
	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b;
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b);		
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	printf("delta0=%.3e\n",delta0(C0b,C0spec,C1b,C1spec,&param,mu_b,mu_spec,lambda_h));

	mu_W=2.*param.mass_W;
	mu_b=param.mass_b_1S/2.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	printf("BR_bsgamma=%.3e\n",bsgamma(C0b,C1b,C2b,mu_b,mu_W,&param));

	mu_b=param.mass_b;
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=0.;
	Cpb[10]=CQpb[1]=CQpb[2]=0.;
	printf("BR_Bsmumu=%.3e\n",Bsmumu(C0b,C1b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR_Bdmumu=%.3e\n",Bdmumu(C0b,C1b,CQ0b,CQ1b,&param,mu_b));

      	printf("BR_Btaunu=%.3e\n",Btaunu(&param));
      	printf("Rtaunu=%.3e\n",RBtaunu(&param));
      	printf("BR_BDtaunu=%.3e\n",BDtaunu(&param));
      	printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu(&param));
     	printf("BR_Dstaunu=%.3e\n",Dstaunu(&param));
     	printf("BR_Dsmunu=%.3e\n",Dsmunu(&param));
     	printf("BR_Dmunu=%.3e\n",Dmunu(&param));
      	printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu(&param));
     	printf("Rmu23=%.3e\n\n",Rmu23(&param));

	return 1;
}
