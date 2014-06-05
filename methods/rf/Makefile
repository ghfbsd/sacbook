# Receiver function tools -- G. Helffrich/U Bristol
#    Last updated 8 Aug. 2013
# 
# Configuration:
#  FC = your Fortran compiler
#  TPTT = location of the Buland-Kennett tau-P subroutine interface
#         See http://www1.gly.bris.ac.uk/~george/sac-bugs.html#ttimes for code
#  SACLIB = location of the SAC I/O library
#         Distributed with SAC

FC= gfortran

OPT= \
	brent.o \
	dcool.o \
	deconsubs.o \
	dilog.o \
	four1.o \
	four1d.o \
	gcdist.o \
	ludcmp.o \
	lubksb.o \
	medsq.o \
	meigen.o \
	minv.o \
	modps.o \
	mtrefps.o \
	mtaper.o \
	rayps.o \
	realftalt.o \
	rkqc.o \
	rk4.o \
	sort2.o \
	sphcap.o \
	spline.o \
	savgol.o \
	tpray.o \
	tsenv.o \
	wigint.o \
	zbrent.o

SLINK= realft.f four1.f

SPECS=
DEBUG=	-O2
DEBUG=	-g -fbounds-check
FFLAGS= ${DEBUG} ${SPECS}
PFLAGS= ${DEBUG}
CFLAGS= -g

TPTT= /usr/local/lib/tpttsub-64.o
TPTT= /usr/local/lib/tpttsub-32.o
TPTT= /usr/local/lib/tpttsub.o
SACLIB= /usr/local/lib/libg77-sacio.a
SACLIB= /usr/local/lib/libsacio.a

all: hk mtdecon sacsetomarker setrfslow

${OPT}:
	$(FC) -O2 -c $*.f

realft.o: realft.f
four1.o: four1.f

hk: hk.o wigint.o gnum.o tokens.o hilbert.o realftalt.o
	$(FC) ${FFLAGS} -o hk \
	   hk.o wigint.o gnum.o tokens.o hilbert.o realftalt.o ${SACLIB}

setrfslow: setrfslow.o gcdist.o tokens.o ${SACLIB} ${TPTT}
	$(FC) ${FFLAGS} -o setrfslow setrfslow.o tokens.o gcdist.o \
	   ${SACLIB} ${TPTT}

mtdecon: mtdecon.o realft.o four1.o mtaper.o gnum.o tokens.o
	$(FC) ${FFLAGS} -o mtdecon mtdecon.o realft.o four1.o mtaper.o gnum.o \
	   tokens.o ${SACLIB}

sacsetomarker: sacsetomarker.o gnum.o tokens.o gcdist.o
	$(FC) ${FFLAGS} -o sacsetomarker sacsetomarker.o gnum.o gcdist.o \
	   tokens.o ${SACLIB} ${TPTT}

clean:
	/bin/rm -f *.o

distclean:
	/bin/rm -f *.o hk setrfslow mtdecon sacsetomarker four1.f realft.f
	/bin/rm -f *.[ns]rf[rtz] *.dpar *.bh[lqt]

$(SLINK): links
links:
	for f in $(SLINK) ; do \
	    [ -r $$f ] || ln -sf ~/seismo/src/$$f $$f ; \
	done
