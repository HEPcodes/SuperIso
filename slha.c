#include "src/include.h"

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma, and some other observables using a given SLHA file */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float delta0;

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
  
	delta0=delta0_calculator(name);
	if(delta0 !=0.)
	{
		printf("delta0=%.3e\n",delta0);
       		printf("BR_bsgamma=%.3e\n",BRbsgamma_calculator(name));
      		printf("BR_Btaunu=%.3e\n",btaunu_calculator(name));
      		printf("BR_Kmunu/BR_pimunu=%.3e\n",Bkmunu_Bpimunu_calculator(name));
     		printf("Rl23=%.3e\n",Rl23_calculator(name));
      		printf("BR_BDtaunu=%.3e\n",Bbdtaunu_calculator(name));
      		printf("BR_BDtaunu/BR_BDenu=%.3e\n",Bbdtaunu_Bbdenu_calculator(name));
		printf("a_muon=%.3e\n",muon_gm2_calculator(name));
		printf("excluded_mass=%d\n\n",excluded_mass_calculator(name));
	}
	else printf("Invalid point\n\n");

	return 1;
}
