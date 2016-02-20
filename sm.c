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

	double C0b[9],C0spec[9],C1b[9],C1spec[9],C0w[9],C1w[9],C2w[9],C2b[9];
	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_1S;
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b_1S);		
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	printf("delta0=%.3e\n",delta0(C0b,C0spec,C1b,C1spec,&param,mu_b,mu_spec,lambda_h));

	mu_W=2.*param.mass_W;
	mu_b=param.mass_b_1S/2.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	printf("BR_bsgamma=%.3e\n",bsgamma(C0b,C1b,C2b,mu_b,mu_W,&param));

      	printf("BR_Btaunu=%.3e\n",Btaunu(&param));
      	printf("Rtaunu=%.3e\n",RBtaunu(&param));
      	printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu(&param));
     	printf("Rl23=%.3e\n",Rl23(&param));
      	printf("BR_BDtaunu=%.3e\n",BDtaunu(&param));
      	printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu(&param));
	printf("BR_Bsmumu=%.3e\n",Bsmumu(&param));
     	printf("BR_Dstaunu=%.3e\n",Dstaunu(&param));
     	printf("BR_Dsmunu=%.3e\n",Dsmunu(&param));
	
	return 1;
}
