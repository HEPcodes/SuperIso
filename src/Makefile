
.KEEP_STATE:

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libisospin.a

include FlagsForMake

.c.o:
	$(CC) -c $(CFLAGS) $<
.c.a:
	$(CC) -c $(CFLAGS) $<
	ar r $@ $*.o; rm $*.o


all: libisospin.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac
	   
clean:
	rm *.a isajet.h softsusy.h FlagsForMake
	
distclean: 
	rm *.a *.o *.x isajet.h softsusy.h FlagsForMake
	
libisospin.a: libisospin.a(softsusy.o) libisospin.a(leshouches.o)\
	  libisospin.a(isospin.o) libisospin.a(isajet.o) libisospin.a(wilson.o)\
	  libisospin.a(bsgamma.o) libisospin.a(alphas.o) libisospin.a(general.o)\
	  libisospin.a(excluded_masses.o) libisospin.a(quark_masses.o) libisospin.a(muon.o)
	$(RANL)

libisospin.a(softsusy.o): softsusy.c include.h softsusy.h
libisospin.a(leshouches.o): leshouches.c include.h
libisospin.a(isospin.o): isospin.c include.h
libisospin.a(isajet.o): isajet.c include.h isajet.h
libisospin.a(wilson.o): wilson.c include.h
libisospin.a(bsgamma.o): bsgamma.c include.h
libisospin.a(alphas.o): alphas.c include.h
libisospin.a(quark_masses.o): quark_masses.c include.h
libisospin.a(general.o): general.c include.h
libisospin.a(excluded_masses.o): excluded_masses.c include.h
libisospin.a(muon.o): muon.c include.h

softsusy.h:
	@echo \#define SOFTSUSY \"$(SOFTSUSY)\" > softsusy.h;

isajet.h:
	@echo \#define ISAJET \"$(ISAJET)\" > isajet.h;


