#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables using a given SLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	int test;

  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the SLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}
  
	test=test_slha(name);
	
	if(test==1)
	{
		printf("delta0=%.3e\n",delta0_calculator(name));
       		printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      		printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      		printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      		printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     		printf("Rl23=%.3e\n",Rl23_calculator(name));
      		printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      		printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
		printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
		printf("a_muon=%.3e\n",muon_gm2_calculator(name));
 		printf("excluded_mass=%d\n\n",excluded_mass_calculator(name));
	}
	else if(test==-1) printf("Invalid point\n\n");
	else if(test==-2) printf("Model not yet implemented\n\n");
	else if(test==-3) printf("Invalid SLHA file\n\n");
	
	return 1;
}
