#include "src/include.h"
#include "src/2hdmc.h"

#define USE_2HDMC

/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to a 2HDM point generated by 2HDMC */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	double lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,m12_2,tanb,mA,delta0;
	int type;

	if(argc<4) 
	{ 
		printf(" This program needs 3 parameters:\n"
		"   type    Yukawa type (1-4)\n"
		"   tanb    tan(beta)\n"
		"   mA      CP-odd Higgs mass\n");
		printf(" Auxiliary parameters are:\n"
		"   lambda1...7  Higgs potential parameters (by default MSSM like)\n"
		"   m12^2        if set, mA is not used (otherwise MSSM like)\n");                        
		exit(1); 
	} 
	else  
	{
		sscanf(argv[1],"%d",&type);
		sscanf(argv[2],"%lf",&tanb);
		sscanf(argv[3],"%lf",&mA);

		double g=6.48408288e-1;
		double gp=3.58051564e-1;

		if(argc>4) sscanf(argv[4],"%lf",&lambda1); else lambda1=(g*g+gp*gp)/4.;   
		if(argc>5) sscanf(argv[5],"%lf",&lambda2); else lambda2=lambda1;
		if(argc>6) sscanf(argv[6],"%lf",&lambda3); else lambda3=(g*g-gp*gp)/4.;
		if(argc>7) sscanf(argv[7],"%lf",&lambda4); else lambda4=-g*g/2.;
		if(argc>8) sscanf(argv[8],"%lf",&lambda5); else lambda5=0.;
		if(argc>9) sscanf(argv[9],"%lf",&lambda6); else lambda6=0.;
		if(argc>10) sscanf(argv[10],"%lf",&lambda7); else lambda7=0.;
		if(argc>11) sscanf(argv[11],"%lf",&m12_2); else m12_2=mA*mA*sin(2.*atan(tanb))/2.;
	}	
	sprintf(name,"thdm%d.lha",getpid());


	int filesOK=1;
#ifdef USE_2HDMC
	if(!test_file(THDMC)) 
	{
		printf("\"%s\" absent. Please check the 2HDMC path or comment \"#define USE_2HDMC\" in thdm.c\n",THDMC);
		filesOK=0;
	};
#endif
	if(!filesOK) return 1;


	if(!test_file("tmp")) system("mkdir tmp");
	chdir("tmp");

#ifdef USE_2HDMC
	thdmc_types(lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,m12_2,tanb,type,name);
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
     		printf("BR_Dmunu=%.3e\n",Dmunu_calculator(name));
     		printf("BR_Dstaunu=%.3e\n",Dstaunu_calculator(name));
     		printf("BR_Dsmunu=%.3e\n",Dsmunu_calculator(name));
		printf("a_muon=%.3e\n",muon_gm2_calculator(name));
		flha_generator(name,"output.flha");
		printf("output.flha generated\n\n");	
	}
	else printf("Invalid point\n\n");
#endif

	sprintf(name,"rm -f thdm%d.lha",getpid());
	system(name);
	return 1;
}
