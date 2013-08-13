#Make file for shear_grh.  G. Helffrich/U. Bristol/3 Oct. 2008
#
#Edit SACLIB as appropriate for your SAC library location.
#
#There is no automatic install make target.  Place splstacksac and shear_grh
#   in a directory containing executable code referenced by your search path.
FC = g77
FC = gfortran

FFLAGS  = -O
FFLAGS  = -g

LIBDIR  = ../lib
SACLIBDIR= /usr/local/lib
#Alternative for g77 SAC io libraries; see INSTALL.
SACLIBG77  = $(SACLIBDIR)/libg77-sacio.a
#Default SAC io libraries (gfortran)
SACLIBDEF  = $(SACLIBDIR)/libsacio.a
SACLIB = `test $(FC) = g77 && echo $(SACLIBG77) || echo $(SACLIBDEF)`

FOBJS = \
	wigint.o \
	rsple.o \
	rspln.o \
	ftest_mse_2.o \
	f_test2.o

OBJECTS = \
	ndf_fun2.o \
	shear_min.o \
	shear_bars_ev.o \
	e_t.o \
	e_ev.o \
	crosscorr2.o \
	spline_int.o \
	rotcomp.o \
	aniput.o \
	$(FOBJS) \
	ndf_spect.o \
	rotsub.o \
	addtime.o \
	disaz.o

LIBS =	$(LIBDIR)/tslib.a 

OPT = wigint.o crosscorr2.o rotcomp.o rotsub.o rsple.o rspln.o

all: shear_grh splstacksac
default: shear_grh

${OPT}:
	$(FC) -O ${FFLAGS} -c $*.f

rsach.o: rsach.sachdr.com

$(LIBDIR)/tslib.a:
	(cd ../tslib; make FC=$(FC) install)

shear_grh: shear_grh.o $(OBJECTS) $(LIBS)
	$(FC) ${FFLAGS} shear_grh.o $(OBJECTS) $(LIBS) $(SACLIB) -o shear_grh
	chmod g=rwx,o=rx shear_grh

splstacksac: stacksac.o $(FOBJS) $(LIBS)
	$(FC) ${FFLAGS} stacksac.o $(FOBJS) $(LIBS) $(SACLIB) -o splstacksac
	chmod g=rwx,o=rx splstacksac

clean:
	(cd ../tslib; make clean)
	/bin/rm -f core *.o

distclean:
	(cd ../tslib; make clean)
	/bin/rm -f core *.o ../lib/tslib.a shear_grh splstacksac

dist: distclean
	-rmdir ../lib
	/bin/rm -f ../*.var ../*.imea ../*.split ../.sachist
	(cd ../..; tar cfz /tmp/split.tgz split)
