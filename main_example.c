#include "src/include.h"

/*#define USE_ISAJET*/ /* to be commented if ISAJET is unavailable */
/*#define USE_SOFTSUSY*/ /* to be commented if SOFTSUSY is unavailable */

/*--------------------------------------------------------------------*/
/* Calculation of the isospin asymmetry and inclusive branching ratio of b->s gamma corresponding to several MSSM test-points */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	
	sprintf(name,"test.lha");

#ifdef USE_ISAJET
	isajet_sugra(500., 500., 50., 0., 1., 172., name);
	printf("delta0_isajet_sugra=%f\n",delta0_calculator(name));
       	printf("BR_isajet_sugra=%f\n\n",BRbsgamma_calculator(name));
	system("rm test.lha");
#endif

#ifdef USE_SOFTSUSY
	softsusy_sugra(500., 500., 50., 0., 1., 172.5, 4.2, 0.1172, name);
	printf("delta0_softsusy_sugra=%f\n",delta0_calculator(name));
       	printf("BR_softsusy_sugra=%f\n\n",BRbsgamma_calculator(name));
 	
        softsusy_amsb(500., 40000., 60., -1., 172., 4.2, 0.1172, name);
	printf("delta0_softsusy_amsb=%f\n",delta0_calculator(name));
       	printf("BR_softsusy_amsb=%f\n\n",BRbsgamma_calculator(name));

        softsusy_gmsb(50000., 70000., 50., 5, 1., 1., 172., 4.2, 0.1172, name);
	printf("delta0_softsusy_gmsb=%f\n",delta0_calculator(name));
       	printf("BR_softsusy_gmsb=%f\n\n",BRbsgamma_calculator(name));
	system("rm test.lha");
#endif

	sprintf(name,"example.lha");
	printf("delta0_example_slha=%f\n",delta0_calculator(name));
       	printf("BR_example_slha=%f\n\n",BRbsgamma_calculator(name));
	
	
	/* example of direct call of the delta0m function */
	
	/* necessary parameters */
	struct parameters param;
	param.mass_s = 0.1;
	param.mass_c = 1.3;
	param.mass_b = 4.2;
	param.mass_top_pole = 172.5;
	param.mass_Z=91.1876;
	param.alpha_s_MZ=0.1172;
	float mass_b_pole=4.9;
	
	float C0[9],C0spec[9],C1[9],C1spec[9];
	int ie;
	
	for(ie=1;ie<=8;ie++) C0[ie]=C1[ie]=C0spec[ie]=C1spec[ie]=0.;
	
	/* necessary Wilson coefficients at mb scale and spectator scale (for the SM) */
	C0[1]=1.10;
	C1[1]=-2.01;
	C0spec[1]=1.24;

	C0[2]=-0.24;
	C1[2]=4.41;
		
	C0[3]=0.011;
	C1[3]=0.086;
	
	C0[4]=-0.025;
	C1[4]=-0.45;
	
	C0[5]=0.007;
	C1[5]=0.051;
	
	C0[6]=-0.030;
	C1[6]=-0.41;
	
	C0[7]=-0.31;
	C1[7]=0.48;
	
	C0[8]=-0.15;
	C0spec[8]=-0.18;
	
	printf("delta0_direct_SM=%f\n", delta0m(C0,C0spec,C1,C1spec,&param,mass_b_pole,sqrt(0.5*mass_b_pole),0.5));

	return 1;
}
