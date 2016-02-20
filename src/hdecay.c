#include "include.h"
#include "hdecay.h"

int Hdecay(char name[], struct parameters* param)
{
	FILE *tmp;
	char tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.hdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);
	
	tmp=fopen("hdecay.in","w");
	
	fprintf(tmp,"SLHAIN   = 0\n");
	fprintf(tmp,"SLHAOUT  = 1\n");
	fprintf(tmp,"HIGGS    = 5\n");
	fprintf(tmp,"MODEL    = 1\n");
	fprintf(tmp,"TGBET    = %.5e\n",param->tan_beta);
	fprintf(tmp,"MABEG    = %.5e\n",param->mass_A0);
	fprintf(tmp,"MAEND    = %.5e\n",param->mass_A0);
	fprintf(tmp,"NMA      = 1\n");
	fprintf(tmp,"ALS(MZ)  = %.5e\n",param->alphas_MZ);
	fprintf(tmp,"MSBAR(1) = %.5e\n",param->mass_s);
	fprintf(tmp,"MC       = %.5e\n",param->mass_c);
	fprintf(tmp,"MB       = %.5e\n",param->mass_b);
	fprintf(tmp,"MT       = %.5e\n",param->mass_top_pole);
	fprintf(tmp,"MTAU     = %.5e\n",param->mass_tau_pole);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GF       = %.5e\n",param->Gfermi);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"MZ       = %.5e\n",param->mass_Z);
	fprintf(tmp,"MW       = %.5e\n",param->mass_W);
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VUB/VCB  = %.5e\n",cabs(param->Vub/param->Vcb));
	fprintf(tmp,"MU       = %.5e\n",param->mu_Q);
	fprintf(tmp,"M2       = %.5e\n",param->M2_Q);
	fprintf(tmp,"MGLUINO  = %.5e\n",param->mass_gluino);
	fprintf(tmp,"MSL1     = %.5e\n",param->MeL_Q);
	fprintf(tmp,"MER1     = %.5e\n",param->MeR_Q);
	fprintf(tmp,"MQL1     = %.5e\n",param->MqL1_Q);
	fprintf(tmp,"MUR1     = %.5e\n",param->MuR_Q);
	fprintf(tmp,"MDR1     = %.5e\n",param->MdR_Q);
	fprintf(tmp,"MSL      = %.5e\n",param->MtauL_Q);
	fprintf(tmp,"MER      = %.5e\n",param->MtauR_Q);
	fprintf(tmp,"MSQ      = %.5e\n",param->MqL3_Q);
	fprintf(tmp,"MUR      = %.5e\n",param->MtR_Q);
	fprintf(tmp,"MDR      = %.5e\n",param->MbR_Q);
	fprintf(tmp,"AL       = %.5e\n",param->A_tau);
	fprintf(tmp,"AU       = %.5e\n",param->A_t);
	fprintf(tmp,"AD       = %.5e\n",param->A_b);
	fprintf(tmp,"NNLO (M) = 1\n");
	fprintf(tmp,"ON-SHELL = 0\n");
	fprintf(tmp,"ON-SH-WZ = 0\n");
	fprintf(tmp,"IPOLE    = 0\n");
	fprintf(tmp,"OFF-SUSY = 0\n");
	fprintf(tmp,"INDIDEC  = 0\n");
	fprintf(tmp,"NF-GG    = 5\n");
	fprintf(tmp,"IGOLD    = 0\n");
	fprintf(tmp,"MPLANCK  = %.5e\n",2.4e18);
	fprintf(tmp,"MGOLD    = %.5e\n",1.e-13);
	
	fclose(tmp);
	
	sprintf(tmp_char,"%s",HDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		printf("Problem with Hdecay, using Tree level...\n");
		return HdecayTree(name,param);
	}
	
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500];
			
	param->BRh0invisible=param->BRH0invisible=param->BRA0invisible=0.;
	
	lecture = fopen("slha.out","r");
	
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 25: fscanf(lecture,"%lf",&param->mass_h0); break;
					case 35: fscanf(lecture,"%lf",&param->mass_H0); break;
					case 36: fscanf(lecture,"%lf",&param->mass_A0); break;
					case 37: fscanf(lecture,"%lf",&param->mass_H); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
	
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
			{
				case 25:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{	
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
														if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0ss=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0cc=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2h0bb=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2h0toptop=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2h0tautau=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2h0WW=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2h0gg=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2h0gaga=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2h0Zga=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2h0ZZ=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2h0ZH0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2h0ZA0=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRh0H0H0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRh0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRh0invisible+=atof(dummy);

						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 26:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0ss=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0cc=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2h0bb=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2h0toptop=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2h0tautau=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2h0WW=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2h0gg=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2h0gaga=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2h0Zga=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2h0ZZ=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2h0ZH0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2h0ZA0=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRh0H0H0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRh0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRh0invisible+=atof(dummy);

						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 35:
				{
					fscanf(lecture,"%lf",&param->width_H0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);					if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2H0ss=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2H0cc=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2H0bb=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2H0toptop=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2H0tautau=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2H0WW=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2H0gg=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2H0gaga=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2H0Zga=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2H0ZZ=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2H0Zh0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2H0Zh0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2H0ZH0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2H0ZA0=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRH0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRH0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRH0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRH0invisible+=atof(dummy);
						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 36:
				{
					fscanf(lecture,"%lf",&param->width_A0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2A0ss=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2A0cc=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2A0bb=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2A0toptop=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2A0tautau=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2A0WW=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2A0gg=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2A0gaga=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2A0Zga=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2A0ZZ=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2A0Zh0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2A0Zh0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2A0ZH0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2A0ZA0=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRA0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRA0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRA0H0H0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRA0invisible+=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 37: 
				{
					fscanf(lecture,"%lf",&param->width_H);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"4")))) param->BRHcs=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"4")))) param->BRHcb=atof(dummy);
							else if(((!strcasecmp(id1,"16"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"16")))) param->BRHtaunu=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.) 
	{
		printf("Problem with Hdecay, using Tree level...\n");
		return HdecayTree(name,param);
	} 
	else  return 1;
}

/* -------------------------------- */

int HdecayTree(char name[], struct parameters* param)
{
	FILE *tmp;
	char tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.hdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);
	
	tmp=fopen("hdecay.in","w");
	
	fprintf(tmp,"SLHAIN   = 0\n");
	fprintf(tmp,"SLHAOUT  = 1\n");
	fprintf(tmp,"HIGGS    = 5\n");
	fprintf(tmp,"MODEL    = 1\n");
	fprintf(tmp,"TGBET    = %.5e\n",param->tan_beta);
	fprintf(tmp,"MABEG    = %.5e\n",param->mass_A0);
	fprintf(tmp,"MAEND    = %.5e\n",param->mass_A0);
	fprintf(tmp,"NMA      = 1\n");
	fprintf(tmp,"ALS(MZ)  = %.5e\n",param->alphas_MZ);
	fprintf(tmp,"MSBAR(1) = %.5e\n",param->mass_s);
	fprintf(tmp,"MC       = %.5e\n",param->mass_c);
	fprintf(tmp,"MB       = %.5e\n",param->mass_b);
	fprintf(tmp,"MT       = %.5e\n",param->mass_top_pole);
	fprintf(tmp,"MTAU     = %.5e\n",param->mass_tau_pole);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GF       = %.5e\n",param->Gfermi);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"MZ       = %.5e\n",param->mass_Z);
	fprintf(tmp,"MW       = %.5e\n",param->mass_W);
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VUB/VCB  = %.5e\n",cabs(param->Vub/param->Vcb));
	fprintf(tmp,"MU       = %.5e\n",param->mu_Q);
	fprintf(tmp,"M2       = %.5e\n",param->M2_Q);
	fprintf(tmp,"MGLUINO  = %.5e\n",param->mass_gluino);
	fprintf(tmp,"MSL1     = %.5e\n",param->MeL_Q);
	fprintf(tmp,"MER1     = %.5e\n",param->MeR_Q);
	fprintf(tmp,"MQL1     = %.5e\n",param->MqL1_Q);
	fprintf(tmp,"MUR1     = %.5e\n",param->MuR_Q);
	fprintf(tmp,"MDR1     = %.5e\n",param->MdR_Q);
	fprintf(tmp,"MSL      = %.5e\n",param->MtauL_Q);
	fprintf(tmp,"MER      = %.5e\n",param->MtauR_Q);
	fprintf(tmp,"MSQ      = %.5e\n",param->MqL3_Q);
	fprintf(tmp,"MUR      = %.5e\n",param->MtR_Q);
	fprintf(tmp,"MDR      = %.5e\n",param->MbR_Q);
	fprintf(tmp,"AL       = %.5e\n",param->A_tau);
	fprintf(tmp,"AU       = %.5e\n",param->A_t);
	fprintf(tmp,"AD       = %.5e\n",param->A_b);
	fprintf(tmp,"NNLO (M) = 0\n");
	fprintf(tmp,"ON-SHELL = 0\n");
	fprintf(tmp,"ON-SH-WZ = 0\n");
	fprintf(tmp,"IPOLE    = 0\n");
	fprintf(tmp,"OFF-SUSY = 0\n");
	fprintf(tmp,"INDIDEC  = 0\n");
	fprintf(tmp,"NF-GG    = 5\n");
	fprintf(tmp,"IGOLD    = 0\n");
	fprintf(tmp,"MPLANCK  = %.5e\n",2.4e18);
	fprintf(tmp,"MGOLD    = %.5e\n",1.e-13);
	
	fclose(tmp);
	
	sprintf(tmp_char,"%s",HDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		return 0;
	}
	
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500];
			
	param->BRh0invisible=param->BRH0invisible=param->BRA0invisible=0.;
	
	lecture = fopen("slha.out","r");
	
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 25: fscanf(lecture,"%lf",&param->mass_h0); break;
					case 35: fscanf(lecture,"%lf",&param->mass_H0); break;
					case 36: fscanf(lecture,"%lf",&param->mass_A0); break;
					case 37: fscanf(lecture,"%lf",&param->mass_H); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
	
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
			{
				case 25:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{	
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
														if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0ss=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2h0cc=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2h0bb=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2h0toptop=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2h0tautau=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2h0WW=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2h0gg=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2h0gaga=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2h0Zga=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2h0ZZ=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2h0ZH0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2h0ZA0=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRh0H0H0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRh0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRh0invisible+=atof(dummy);

						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 26:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0ss=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2h0cc=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2h0bb=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2h0toptop=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2h0tautau=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2h0WW=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2h0gg=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2h0gaga=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2h0Zga=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2h0ZZ=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2h0Zh0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2h0ZH0=atof(dummy)*param->width_h0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2h0ZA0=atof(dummy)*param->width_h0;
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRh0H0H0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRh0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRh0invisible+=atof(dummy);

						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 35:
				{
					fscanf(lecture,"%lf",&param->width_H0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);					if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2H0ss=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2H0cc=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2H0bb=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2H0toptop=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2H0tautau=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2H0WW=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2H0gg=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2H0gaga=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2H0Zga=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2H0ZZ=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2H0Zh0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2H0Zh0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2H0ZH0=atof(dummy)*param->width_H0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2H0ZA0=atof(dummy)*param->width_H0;
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRH0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRH0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"36"))) param->BRH0A0A0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRH0invisible+=atof(dummy);
						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 36:
				{
					fscanf(lecture,"%lf",&param->width_A0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2A0ss=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2A0cc=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2A0bb=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2A0toptop=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2A0tautau=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2A0WW=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2A0gg=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2A0gaga=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2A0Zga=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2A0ZZ=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))) param->g2A0Zh0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))) param->g2A0Zh0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"35")))||((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"23")))) param->g2A0ZH0=atof(dummy)*param->width_A0;
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"23")))) param->g2A0ZA0=atof(dummy)*param->width_A0;
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRA0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRA0h0h0=atof(dummy);
							else if((!strcasecmp(id1,"35"))&&(!strcasecmp(id2,"35"))) param->BRA0H0H0=atof(dummy);
							else if((abs(atoi(id1))>1000000)&&(abs(atoi(id2))>1000000)) param->BRA0invisible+=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 37: 
				{
					fscanf(lecture,"%lf",&param->width_H);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"4")))) param->BRHcs=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"4")))) param->BRHcb=atof(dummy);
							else if(((!strcasecmp(id1,"16"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"16")))) param->BRHtaunu=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.) 
	return 0; else return 2;
}

/* -------------------------------- */

int HdecaySM(char name[], int ie, struct parameters* param)
{
	FILE *tmp;
	char tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.hdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);
	
	tmp=fopen("hdecay.in","w");
	
	fprintf(tmp,"SLHAIN   = 0\n");
	fprintf(tmp,"SLHAOUT  = 1\n");
	fprintf(tmp,"HIGGS    = 0\n");
	fprintf(tmp,"MODEL    = 1\n");
	fprintf(tmp,"TGBET    = %.5e\n",0.);
	fprintf(tmp,"MABEG    = %.5e\n",param->mass_hSM[ie]);
	fprintf(tmp,"MAEND    = %.5e\n",param->mass_hSM[ie]);
	fprintf(tmp,"NMA      = 1\n");
	fprintf(tmp,"ALS(MZ)  = %.5e\n",param->alphas_MZ);
	fprintf(tmp,"MSBAR(1) = %.5e\n",param->mass_s);
	fprintf(tmp,"MC       = %.5e\n",param->mass_c);
	fprintf(tmp,"MB       = %.5e\n",param->mass_b);
	fprintf(tmp,"MT       = %.5e\n",param->mass_top_pole);
	fprintf(tmp,"MTAU     = %.5e\n",param->mass_tau_pole);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GF       = %.5e\n",param->Gfermi);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"MZ       = %.5e\n",param->mass_Z);
	fprintf(tmp,"MW       = %.5e\n",param->mass_W);
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VUB/VCB  = %.5e\n",cabs(param->Vub/param->Vcb));
	fprintf(tmp,"MU       = %.5e\n",0.);
	fprintf(tmp,"M2       = %.5e\n",0.);
	fprintf(tmp,"MGLUINO  = %.5e\n",0.);
	fprintf(tmp,"MSL1     = %.5e\n",0.);
	fprintf(tmp,"MER1     = %.5e\n",0.);
	fprintf(tmp,"MQL1     = %.5e\n",0.);
	fprintf(tmp,"MUR1     = %.5e\n",0.);
	fprintf(tmp,"MDR1     = %.5e\n",0.);
	fprintf(tmp,"MSL      = %.5e\n",0.);
	fprintf(tmp,"MER      = %.5e\n",0.);
	fprintf(tmp,"MSQ      = %.5e\n",0.);
	fprintf(tmp,"MUR      = %.5e\n",0.);
	fprintf(tmp,"MDR      = %.5e\n",0.);
	fprintf(tmp,"AL       = %.5e\n",0.);
	fprintf(tmp,"AU       = %.5e\n",0.);
	fprintf(tmp,"AD       = %.5e\n",0.);
	fprintf(tmp,"NNLO (M) = 1\n");
	fprintf(tmp,"ON-SHELL = 0\n");
	fprintf(tmp,"ON-SH-WZ = 0\n");
	fprintf(tmp,"IPOLE    = 0\n");
	fprintf(tmp,"OFF-SUSY = 0\n");
	fprintf(tmp,"INDIDEC  = 0\n");
	fprintf(tmp,"NF-GG    = 5\n");
	fprintf(tmp,"IGOLD    = 0\n");
	fprintf(tmp,"MPLANCK  = %.5e\n",2.4e18);
	fprintf(tmp,"MGOLD    = %.5e\n",1.e-13);
	
	fclose(tmp);
	
	sprintf(tmp_char,"%s",HDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		return 0;
	}
	
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500];
	
	lecture = fopen("slha.out","r");
	
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 36: fscanf(lecture,"%lf",&param->mass_hSM[ie]); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);	
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(abs(atoi(dummy)))
			{
				case 25: 
				{
					fscanf(lecture,"%lf",&param->width_hSM[ie]);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(atof(dummy)>0.) 
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->g2hss_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->g2hcc_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->g2hbb_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->g2htoptop_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->g2htautau_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->g2hWW_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->g2hgaga_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->g2hgg_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if(((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))||((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))) param->g2hZga_SM[ie]=atof(dummy)*param->width_hSM[ie];
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->g2hZZ_SM[ie]=atof(dummy)*param->width_hSM[ie];
						}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					}	
					break;
				}
				
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
