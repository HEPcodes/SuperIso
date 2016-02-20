#include "src/include.h"

#define USE_SOFTSUSY

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to a GMSB point generated with SOFTSUSY */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float Lambda,Mmess,tanb,cGrav,sgnmu,mtop,mbot,alphas_mz,delta0;
	int N5;

  	if(argc<5) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   Lambda  scale of the SUSY breaking (10000-100000 GeV)\n"
           	"   Mmess   messenger mass scale > Lambda\n"
	   	"   N5      equivalent number of 5+5bar messenger fields\n"
           	"   tanb    tan(beta) \n");
    		printf(" Auxiliary parameters are:\n"
           	"   sgnmu   +/-1,  sign of Higgsino mass term (default 1)\n" 
	   	"   cGrav   >=1, ratio of the gravitino mass to its value for a breaking scale of F_m \n"
           	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass\n"
           	"   alphas_mz  strong coupling at MZ\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%f",&Lambda);
     		sscanf(argv[2],"%f",&Mmess);
     		sscanf(argv[3],"%d",&N5);
     		sscanf(argv[4],"%f",&tanb);
     		if(argc>5) sscanf(argv[5],"%f",&sgnmu); else sgnmu=1.;
     		if(argc>6) sscanf(argv[6],"%f",&cGrav); else cGrav=1.;   
     		if(argc>7) sscanf(argv[7],"%f",&mtop); else mtop=172.5;   
     		if(argc>8) sscanf(argv[8],"%f",&mbot); else mbot=4.2;
     		if(argc>9) sscanf(argv[9],"%f",&alphas_mz); else alphas_mz=0.1172;
  	}	

	if(Lambda>Mmess) printf("Lambda=%.0f must be smaller than Mmess=%.0f\n\n",Lambda,Mmess);
	else
	{
		sprintf(name,"gmsb.lha");

#ifdef USE_SOFTSUSY
		softsusy_gmsb(Lambda, Mmess, tanb, N5, cGrav, sgnmu, mtop, mbot, alphas_mz, name);

		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{
			printf("delta0=%f\n",delta0);
       			printf("BR=%f\n",BRbsgamma_calculator(name));
			printf("excluded_mass=%d\n\n",excluded_mass_calculator(name));
		}
		else printf("Invalid point\n\n");
#endif

		sprintf(name,"rm gmsb.lha");
		system(name);		
		return 1;
	}
}
