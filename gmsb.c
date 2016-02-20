#include "src/include.h"

#define USE_SOFTSUSY

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to a GMSB point generated by SOFTSUSY */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	double Lambda,Mmess,tanb,cGrav,sgnmu,mtop,mbot,alphas_mz,delta0;
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
		sscanf(argv[1],"%lf",&Lambda);
     		sscanf(argv[2],"%lf",&Mmess);
     		sscanf(argv[3],"%d",&N5);
     		sscanf(argv[4],"%lf",&tanb);
     		if(argc>5) sscanf(argv[5],"%lf",&sgnmu); else sgnmu=1.;
     		if(argc>6) sscanf(argv[6],"%lf",&cGrav); else cGrav=1.;   
     		if(argc>7) sscanf(argv[7],"%lf",&mtop); else mtop=172.4;   
     		if(argc>8) sscanf(argv[8],"%lf",&mbot); else mbot=4.2;
     		if(argc>9) sscanf(argv[9],"%lf",&alphas_mz); else alphas_mz=0.1176;
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
			printf("delta0=%.3e\n",delta0);
       			printf("BR_bsgamma=%.3e\n",bsgamma_calculator(name));
      			printf("BR_Btaunu=%.3e\n",Btaunu_calculator(name));
      			printf("Rtaunu=%.3e\n",RBtaunu_calculator(name));
      			printf("BR_Kmunu/BR_pimunu=%.3e\n",Kmunu_pimunu_calculator(name));
     			printf("Rl23=%.3e\n",Rl23_calculator(name));
      			printf("BR_BDtaunu=%.3e\n",BDtaunu_calculator(name));
      			printf("BR_BDtaunu/BR_BDenu=%.3e\n",BDtaunu_BDenu_calculator(name));
			printf("BR_Bsmumu=%.3e\n",Bsmumu_calculator(name));
     			printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     			printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
			printf("a_muon=%.3e\n\n",muon_gm2_calculator(name));
			printf("excluded_mass=%d\n",excluded_mass_calculator(name));
		}
		else printf("Invalid point\n\n");
#endif

		sprintf(name,"rm gmsb.lha");
		system(name);		
		return 1;
	}
}
