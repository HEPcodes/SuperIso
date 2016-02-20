.KEEP_STATE:

#
VERSION = v2.8

# Choose your compilers here (in general gcc/gfortran on Linux systems):
CC = gcc
CFLAGS= -O3 -pipe -fomit-frame-pointer

#CC = icc
#CFLAGS = -O3

MAKE = make
AR = ar

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libisospin.a librelic.a

# Add the link to Softsusy, Isajet, SuSpect and SPheno main programs, if available.
# Otherwise, comment them in the main programs */
# If you change the path, you must "make distclean" first.
SOFTSUSY = ~/softsusy/softpoint.x
ISAJET = ~/isajet/isasugra.x
SUSPECT = ~/suspect/suspect2
SPHENO = ~/spheno/bin/SPheno
# Add the link to 2HDMC directory, if available.
THDMC = ~/2HDMC
# Add the links to Hdecay and HiggsBounds, if available.
HDECAY = ~/hdecay
HIGGSBOUNDS = ~/higgsbounds/HiggsBounds-f77/HiggsBounds

CINCLUDE= -I./src -L./src

all: libisospin.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.c" into "name.x"';\
	echo ' '

%.c:: %.c libisospin.a
	$(CC) -c $(CFLAGS) $@;
	$(CC) -o $*.x $(CFLAGS) $(CINCLUDE) $*.o -lisospin -lm;
	@rm -f $*.o;
	@touch $*.x

%:: %.c libisospin.a
	$(CC) -c $(CFLAGS) $*.c;
	$(CC) -o $*.x $(CFLAGS) $(CINCLUDE) $*.o -lisospin -lm;
	@rm -f $*.o;
	@touch $*.x

clean:
	rm -rf tmp *.x *.bin .*.bin *.tmplha *.hbtmp *.fhtmp *.hdtmp *.fh *.is? *.ss *.sptmp *.sutmp *.fhmod;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ clean
	
distclean: 
	rm -rf tmp *.x *.bin .*.bin *.tmplha *.hbtmp *.fhtmp *.hdtmp *.fh *.is? *.ss *.sptmp *.sutmp *.fhmod;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ distclean
	
libisospin.a: 
	@echo;
	@echo SuperIso $(VERSION) - F.N. Mahmoudi 2010;
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo MAKE = $(MAKE) >> src/FlagsForMake;\
	echo AR = $(AR) >> src/FlagsForMake;\
	echo SOFTSUSY = $(SOFTSUSY) >> src/FlagsForMake;\
	echo ISAJET = $(ISAJET) >> src/FlagsForMake;\
	echo SPHENO = $(SPHENO) >> src/FlagsForMake;\
	echo SUSPECT = $(SUSPECT) >> src/FlagsForMake;\
	echo THDMC = $(THDMC) >> src/FlagsForMake;\
        echo FEYNHIGGS = $(FEYNHIGGS)/build/FeynHiggs >> src/FlagsForMake;\
 	echo HDECAY = $(HDECAY)/run >> src/FlagsForMake;\
 	echo HIGGSBOUNDS = $(HIGGSBOUNDS) >> src/FlagsForMake;
	$(MAKE) -C src/ libisospin.a

save: 
	rm -f superiso_$(VERSION).tgz;\
	mkdir superiso_$(VERSION);\
	cp -p README superiso_$(VERSION)/;\
	cp -p example.lha superiso_$(VERSION)/;\
	cp -p amsb.c superiso_$(VERSION)/;\
	cp -p hcamsb.c superiso_$(VERSION)/;\
	cp -p mmamsb.c superiso_$(VERSION)/;\
	cp -p gmsb.c superiso_$(VERSION)/;\
	cp -p msugra.c superiso_$(VERSION)/;\
	cp -p nuhm.c superiso_$(VERSION)/;\
	cp -p slha.c superiso_$(VERSION)/;\
	cp -p sm.c superiso_$(VERSION)/;\
	cp -p thdm.c superiso_$(VERSION)/;\
	cp -p Makefile superiso_$(VERSION)/;\
	mkdir superiso_$(VERSION)/src;\
	cp -p src/*.h superiso_$(VERSION)/src/;\
	cp -p src/*.c superiso_$(VERSION)/src/;\
	rm -f superiso_$(VERSION)/src/softsusy.h;\
	rm -f superiso_$(VERSION)/src/isajet.h;\
	cp -p src/Makefile superiso_$(VERSION)/src/;\
	tar czvf superiso_$(VERSION).tgz superiso_$(VERSION);\
	rm -rf superiso_$(VERSION)
