#include "include.h"
#include "higgsbounds.h"

int higgsbounds_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calling HiggsBounds */
{
	struct parameters param;
	int ie;
	char tmp_char[500],namemod[200],dummy[200];
	FILE *input,*output;
	char *curdir;
	curdir=getcwd(NULL, 500);

	sprintf(namemod,"%s.hbtmp",name);

	Init_param(&param);

	if(!Les_Houches_Reader(name,&param)) return -1;
	
	if((param.THDM_model==0)&&(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.)) Hdecay(name,&param);
	
 	param.mass_hSM[1]=param.mass_h0;
	param.mass_hSM[2]=param.mass_H0;
	param.mass_hSM[3]=param.mass_A0;
	for(ie=1;ie<=3;ie++) HdecaySM(name,ie,&param);
	
	sprintf(tmp_char,"rm -rf %s",namemod);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namemod);
 	system(tmp_char);
	chdir(namemod);
	
	sprintf(tmp_char,"%s_MH_GammaTot.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e  %.5e  %.5e  ",param.mass_h0,param.mass_H0,param.mass_A0);
		
	fprintf(output,"%.5e  %.5e  %.5e  ",param.width_h0,param.width_H0,param.width_A0);
	fprintf(output,"\n");
	fclose(output);
	
	sprintf(tmp_char,"%s_MHplus_GammaTot.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e  ",param.mass_H);
	fprintf(output,"%.5e\n",param.width_H);
	fclose(output);
	
	sprintf(tmp_char,"%s_effC.dat",namemod);
	output=fopen(tmp_char,"w");
	
	fprintf(output,"1  %.5e  %.5e  %.5e  ",param.g2h0ss/param.g2hss_SM[1],param.g2H0ss/param.g2hss_SM[2],0.);
	fprintf(output,"%.5e  %.5e  %.5e  ",0.,0.,param.g2A0ss/param.g2hss_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0cc/param.g2hcc_SM[1],param.g2H0cc/param.g2hcc_SM[2],0.);
	fprintf(output,"%.5e  %.5e  %.5e  ",0.,0.,param.g2A0cc/param.g2hcc_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0bb/param.g2hbb_SM[1],param.g2H0bb/param.g2hbb_SM[2],0.);
	fprintf(output,"%.5e  %.5e  %.5e  ",0.,0.,param.g2A0bb/param.g2hbb_SM[3]);

	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0toptop/param.g2htoptop_SM[1],param.g2H0toptop/param.g2htoptop_SM[2],0.);
	fprintf(output,"%.5e  %.5e  %.5e  ",0.,0.,param.g2A0toptop/param.g2htoptop_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0tautau/param.g2htautau_SM[1],param.g2H0tautau/param.g2htautau_SM[2],0.);
	fprintf(output,"%.5e  %.5e  %.5e  ",0.,0.,param.g2A0tautau/param.g2htautau_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0WW/param.g2hWW_SM[1],param.g2H0WW/param.g2hWW_SM[2],param.g2A0WW/param.g2hWW_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0ZZ/param.g2hZZ_SM[1],param.g2H0ZZ/param.g2hZZ_SM[2],param.g2A0ZZ/param.g2hZZ_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0Zga/param.g2hZga_SM[1],param.g2H0Zga/param.g2hZga_SM[2],param.g2A0Zga/param.g2hZga_SM[3]);

	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0gaga/param.g2hgaga_SM[1],param.g2H0gaga/param.g2hgaga_SM[2],param.g2A0gaga/param.g2hgaga_SM[3]);
	
	fprintf(output,"%.5e  %.5e  %.5e  ",param.g2h0gg/param.g2hgg_SM[1],param.g2H0gg/param.g2hgg_SM[2],param.g2A0gg/param.g2hgg_SM[3]);
	
	double cw2=param.mass_W/param.mass_Z*param.mass_W/param.mass_Z;
	double norm=1./param.inv_alpha_em*pi/cw2/(1.-cw2);
 	
	fprintf(output,"%.5e  %.5e  %.5e  %.5e  %.5e  %.5e",param.g2h0Zh0/norm,param.g2H0Zh0/norm,param.g2H0ZH0/norm,param.g2A0Zh0/norm,param.g2A0ZH0/norm,param.g2A0ZA0/norm);

	fprintf(output,"\n");
	fclose(output);
	
	sprintf(tmp_char,"%s_BR_H_NP.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n",param.BRh0invisible,param.BRH0invisible,param.BRA0invisible,param.BRh0H0H0,param.BRh0A0A0,param.BRH0h0h0,param.BRH0A0A0,param.BRA0h0h0,param.BRA0H0H0);

	fclose(output);

	sprintf(tmp_char,"%s_BR_Hplus.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e  %.5e  %.5e\n",param.BRHcs,param.BRHcb,param.BRHtaunu);
	fclose(output);

	sprintf(tmp_char,"%s_BR_t.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e  %.5e\n",param.BRtWb,param.BRtHb);
	fclose(output);

	sprintf(tmp_char,"%s_LEP_HpHm_CS_ratios.dat",namemod);
	output=fopen(tmp_char,"w");
	fprintf(output,"1  %.5e\n",1.);
	fclose(output);

	sprintf(tmp_char,"%s LandT effC 3 1 %s_ > %s_tmp",HIGGSBOUNDS,namemod,namemod);
	system(tmp_char);
	

	sprintf(tmp_char,"%s_HiggsBounds_results.dat",namemod);
	if(!test_file(tmp_char)) return -1;
	output=fopen(tmp_char,"r");
	while(EOF != fscanf(output,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(output,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
	;

		if((atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))==1)
		{
			for(ie=1;ie<=5;ie++) fscanf(output,"%s",dummy);
			break;
		}
	}
	fclose(output);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namemod);
 	system(tmp_char);
	
	return (atoi(dummy)==0);
}
