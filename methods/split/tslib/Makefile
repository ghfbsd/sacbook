# Make file for time series library.  Cut into separate modules from tslib.f

DEST = ../lib
OWNER = -o george
GROUP = -g geophys
OWNER =
GROUP =

FFLAGS = -O
FFLAGS = 

MODULES = autocorr.o \
	convolv.o convolvnp.o cosfilt.o costap.o crosscorr.o crosspec.o \
	dcm.o demean.o dtr.o dtrend.o \
	fft.o fftl.o filt.o ftaper.o \
	hamwin.o \
	linpol.o \
	ntwo.o \
	realtr.o rmean.o \
	scalar.o scdcm.o sincfilt.o \
	trpfil.o tsadd.o tscopy.o tscostaper.o tscut.o tsdemean.o tserr.o \
	tsgauss.o tsline.o tsmult.o tsnorm.o tsplotset.o tsraise.o \
	tsread.o tsreverse.o tsscale.o tswrite.o tszero.o two.o \
	zero.o zeroarray.o zerochk.o

tslib: tslib.a

install: tslib.a
	[ -d ${DEST} ] || mkdir ${DEST}
	install ${OWNER} ${GROUP} tslib.a ${DEST}/tslib.a
	ranlib ${DEST}/tslib.a

tslib.a: ${MODULES}
	ar cq tslib.a ${MODULES}

clean:
	/bin/rm -f core tslib.a *.o
