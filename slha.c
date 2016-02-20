#include "src/include.h"

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to a SLHA file */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float delta0,BR;

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
	BR=BRbsgamma_calculator(name);
	if((delta0 !=0.)&&(BR >0.))
	{
		printf("delta0=%.3e\n",delta0);
       		printf("BR=%.3e\n",BR);
		printf("excluded_mass=%d\n",excluded_mass_calculator(name));
		printf("(g-2)=%.3e\n\n",muon_gm2_calculator(name));
	}
	else printf("Invalid point\n\n");

	return 1;
}
