#include "include.h"

int excluded_masses(struct parameters* param)
/* tests whether the SUSY point is excluded by the collider contraints */
/* contraints taken from PDG 2006 */
/* if excluded, return 1, otherwise 0 */
{
	int excluded=0;
	
	if(param->mass_h0!=0.) excluded=(excluded||(fabs(param->mass_h0)<111.)); /* Higgs */
	if(param->mass_neut[1]!=0.) excluded=(excluded||(fabs(param->mass_neut[1])<46.)); /* neutralino */
	if(param->mass_cha1!=0.) excluded=(excluded||(fabs(param->mass_cha1)<67.7)); /* chargino */
	if(param->mass_er!=0.) excluded=(excluded||(fabs(param->mass_er)<88.)); /* slepton R */
	if(param->mass_mur!=0.) excluded=(excluded||(fabs(param->mass_mur)<88.)); /* slepton R */
	if(param->mass_nuel!=0.) excluded=(excluded||(fabs(param->mass_nuel)<43.7)); /* sneutrino */
	if(param->mass_numl!=0.) excluded=(excluded||(fabs(param->mass_numl)<43.7)); /* sneutrino */
	if(param->mass_dnr!=0.) excluded=(excluded||(fabs(param->mass_dnr)<250.)); /* squark R */
	if(param->mass_upr!=0.) excluded=(excluded||(fabs(param->mass_upr)<250.)); /* squark R */
	if(param->mass_str!=0.) excluded=(excluded||(fabs(param->mass_str)<250.)); /* squark R */
	if(param->mass_chr!=0.) excluded=(excluded||(fabs(param->mass_chr)<250.)); /* squark R */
	if(param->mass_t1!=0.) excluded=(excluded||(fabs(param->mass_t1)<92.6)); /* stop */	
	if(param->mass_gluino!=0.) excluded=(excluded||(fabs(param->mass_gluino)<195.)); /* gluino */
	if(param->mass_b1!=0.) excluded=(excluded||(fabs(param->mass_b1)<89.)); /* sbottom */
	if(param->mass_tau1!=0.) excluded=(excluded||(fabs(param->mass_tau1)<81.9)); /* stau */

	return excluded;
}

/*--------------------------------------------------------------------*/

int excluded_mass_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point is excluded by the mass contraints */
{
	float C0[9],C0spec[9],C1[9],C1spec[9];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return excluded_masses(&param);
}

/*--------------------------------------------------------------------*/

int charged_LSP(struct parameters* param)
/* tests whether the SUSY point corresponds to a charged LSP (NLSP if the LSP is a gravitino) */
/* if the LSP is charged, return 1, otherwise 0 */
{
	float mass_non_chargee;
	int charged_LSP=0;

	mass_non_chargee=fabs(param->mass_neut[1]);
	
	if(param->mass_nuel!=0.) mass_non_chargee=min(fabs(param->mass_nuel),mass_non_chargee);
	if(param->mass_numl!=0.) mass_non_chargee=min(fabs(param->mass_numl),mass_non_chargee);
	if(param->mass_nutl!=0.) mass_non_chargee=min(fabs(param->mass_nutl),mass_non_chargee);
	if(param->mass_gluino!=0.) mass_non_chargee=min(fabs(param->mass_gluino),mass_non_chargee);
	
	if(param->mass_tau1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_tau1)<mass_non_chargee));
	if(param->mass_t1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_t1)<mass_non_chargee));
	if(param->mass_cha1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_cha1)<mass_non_chargee));
	
	return charged_LSP;
}

/*--------------------------------------------------------------------*/

int charged_LSP_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point corresponds to a charged LSP */
{
	float C0[9],C0spec[9],C1[9],C1spec[9];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return charged_LSP(&param);
}

