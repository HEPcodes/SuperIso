#include "src/include.h"

/*#define USE_ISAJET*/
#define USE_SOFTSUSY

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to an AMSB point generated with SOFTSUSY */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	float m0,m32,tanb,A0,sgnmu,mtop,mbot,alphas_mz,delta0,BR;

	if(argc<4) 
	{ 
		printf(" This program needs 3 parameters:\n"
		"   m0      universal scalar mass at GUT scale\n"
		"   m3/2    gravitino mass\n"
		"   tanb      tan(beta) \n");
		printf(" Auxiliary parameters are:\n"
		"   sgnmu    +/-1,  sign of Higgsino mass term (default 1)\n"    
		"   mtop     top quark pole mass\n"
		"   mbot     Mb(Mb) scale independent b-quark mass\n"
		"   alphas_mz  strong coupling at MZ\n");
		exit(1); 
  	} 
	else  
  	{
		sscanf(argv[1],"%f",&m0);
     		sscanf(argv[2],"%f",&m32);
     		sscanf(argv[3],"%f",&tanb);
     		if(argc>4) sscanf(argv[4],"%f",&sgnmu); else sgnmu=1;
     		if(argc>5) sscanf(argv[5],"%f",&mtop); else mtop=172.5;   
     		if(argc>6) sscanf(argv[6],"%f",&mbot); else mbot=4.2;
     		if(argc>7) sscanf(argv[7],"%f",&alphas_mz); else alphas_mz=0.1172;
  	}	
	sprintf(name,"amsb.lha");

#ifdef USE_SOFTSUSY
	softsusy_amsb(m0, m32, tanb, sgnmu, mtop, mbot, alphas_mz, name);
	
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
#endif

	sprintf(name,"rm amsb.lha");
	system(name);
	return 1;
}
