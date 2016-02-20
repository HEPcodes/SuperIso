#include "include.h" 
#include "isajet.h"


int isajet_sugra(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a mSUGRA parameter space point using ISAJET */
{
	FILE *tmp;
	char tmp_char[300];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"is1_%s",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"is2_%s\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"1\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n",m0,m12,A0,tanb,sgnmu,mtop);
	fprintf(tmp,"0\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < is1_%s > is3_%s",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f is1_%s is2_%s is3_%s ISALHD.out",name,name,name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a GMSB parameter space point using ISAJET */
{
	FILE *tmp;
	char tmp_char[300];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"is1_%s",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"is2_%s\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"2\n");
	fprintf(tmp,"%.5e,%.5e,%d,%.5e,%.5e,%.5e,%.5e\n", Lambda, Mmess, N5, tanb, sgnmu, mtop, cGrav);
	fprintf(tmp,"0\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < is1_%s > is3_%s",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f is1_%s is2_%s is3_%s ISALHD.out",name,name,name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a AMSB parameter space point using ISAJET */
{
 	FILE *tmp;
	char tmp_char[300];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"is1_%s",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"is2_%s\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"7\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e\n", m0, m32, tanb, sgnmu, mtop);
	fprintf(tmp,"0\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < is1_%s > is3_%s",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f is1_%s is2_%s is3_%s ISALHD.out",name,name,name);
	system(tmp_char);
	return 1;
}
/*--------------------------------------------------------------------*/

int isajet_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, char name[])
/* generates a SLHA file for a NUHM parameter space point using ISAJET */
{
	FILE *tmp;
	char tmp_char[300];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"is1_%s",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"is2_%s\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"3\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n",m0,m12,A0,tanb,mu/fabs(mu),mtop);
	fprintf(tmp,"8\n");
	fprintf(tmp,"%.5e,%.5e\n",mu,mA);
	fprintf(tmp,"0\n");
	fprintf(tmp,"0\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < is1_%s > is3_%s",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f is1_%s is2_%s is3_%s ISALHD.out",name,name,name);
	system(tmp_char);
	return 1;
}
