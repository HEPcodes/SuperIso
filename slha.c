#include "src/include.h"

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to a SLHA file */
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
		printf("delta0=%f\n",delta0);
       		printf("BR=%f\n",BRbsgamma_calculator(name));
		printf("excluded_mass=%d\n\n",excluded_mass_calculator(name));
	}
	else printf("Invalid point or invalid SLHA file\n\n");

	return 1;
}
