
.KEEP_STATE:

#
VERSION = v2.0

# Choose your C compiler here (in general gcc on Linux systems):
CC = gcc
#CC = icc

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libisospin.a

#Optimisation level, eg: -O3
OPT= -O3
#OR debug level: -g(n=1,2,3)
DEBUG= 

CFLAGS= -I./src -L./src $(DEBUG) $(OPT)

# Add the link to Softsusy and Isajet, if available.
# Otherwise, comment them in the main programs */
SOFTSUSY = ~/softsusy/softpoint.x
ISAJET = ~/isajet/isasugra.x

all: libisospin.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.c" into "name.x"';\
	echo ' '

%:: %.c libisospin.a
	$(CC) -c $(CFLAGS) $*.c;
	$(CC) -o $@.x $(CFLAGS) $*.o -lm -lisospin;
	@rm $*.o;
	@touch $*.x

clean:
	rm *.x;
	@echo > src/FlagsForMake;
	make -C src/ clean
	
distclean: 
	rm *.a *.o *.x;
	@echo > src/FlagsForMake;
	make -C src/ distclean
	
libisospin.a: 
	@echo;
	@echo SuperIso $(VERSION) - F.N. Mahmoudi 2008;
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo SOFTSUSY = $(SOFTSUSY) >> src/FlagsForMake;\
	echo ISAJET = $(ISAJET) >> src/FlagsForMake;\
	make -C src/

save: 
	rm superiso_$(VERSION).tgz;\
	mkdir superiso_$(VERSION);\
	cp -p README superiso_$(VERSION)/;\
	cp -p example.lha superiso_$(VERSION)/;\
	cp -p amsb.c superiso_$(VERSION)/;\
	cp -p gmsb.c superiso_$(VERSION)/;\
	cp -p msugra.c superiso_$(VERSION)/;\
	cp -p nuhm.c superiso_$(VERSION)/;\
	cp -p slha.c superiso_$(VERSION)/;\
	cp -p main_example.c superiso_$(VERSION)/;\
	cp -p Makefile superiso_$(VERSION)/;\
	mkdir superiso_$(VERSION)/src;\
	cp -p src/*.h superiso_$(VERSION)/src/;\
	cp -p src/*.c superiso_$(VERSION)/src/;\
	rm superiso_$(VERSION)/src/softsusy.h;\
	rm superiso_$(VERSION)/src/isajet.h;\
	cp -p src/Makefile superiso_$(VERSION)/src/;\
	tar czvf superiso_$(VERSION).tgz superiso_$(VERSION);\
	rm -r superiso_$(VERSION)
