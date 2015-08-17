### Makefile for surface waves
### Gustavo Catao Alves
### Stanford Exploration Project
### July 2015
include ${SEPINC}/SEP.top

SRCDIR = ./Src_working
OBJDIR = ./Obj
MODDIR = ./Mod
BINDIR = ./Bin
PARDIR = ./Par
RESDIR = ./Fig
DATDIR = ./Dat
DATDIR2 = ~/models/marmousi2
BKPDIR = ./Fig-backup

S = ${SRCDIR}
O = ${OBJDIR}
M = ${MODDIR}
B = ${BINDIR}
P = ${PARDIR}
R = ${RESDIR}
D = ${DATDIR}
E = ${DATDIR2}
K = ${BKPDIR}



CUDIR = /usr/local/cuda-7.0
CULIB = ${CUDIR}/lib64
CUINC = ${CUDIR}/include
#UPSTEXFLAGS = invras=n color=y

UF90INCLUDES = -I${CUINC}
UF90LIBS =${GEELIB} -L/usr/lib64 -lgomp -L/usr/local/cuda-7.0/lib64 -lcuda -lcudart -L/usr/lib -lstdc++ 	#-L/opt/FFTW/lib -L./Lib -lfftw3 -lfftw3f
#UF90FLAGS = -fpenmp -nocheck -C -module ${M}
UF90FLAGS = -nocheck ${M}
CC = gcc
CFLAG = -fopenmp=0 -lstdc++

#UF90FLAGS = -fast -override-limits -module ${M}

################################################################################
./Obj/rtmlib-mod.o: $O/rtmlib-modC.o $O/rtmlib-modC-gpu.o
 

./Bin/rtm.x: $O/rtmlib-modC.o $O/rtmlib-modC-gpu.o

#$O/rtmlib-modC.o:rtmlib-modC.c
#	gcc -c $O/rtmlib-modC.o $S/rtmlib-modC.c

$O/rtmlib-modC-gpu.o: $S/rtmlib-modC-gpu.cu
	nvcc -c $O/rtmlib-modC-gpu.o -gencode arch=compute_35,code=sm_35 -Xptxas=-v -fmad=false -use_fast_math -maxrregcount=64 $S/rtmlib-modC-gpu.cu -I${CUINC}
	mv rtmlib-modC-gpu.o $O/rtmlib-modC-gpu.o
	
#nvcc -c $O/rtmlib-modC-gpu.o -gencode arch=compute_35,code=sm_35 -Xptxas=-v -fmad=false -use_fast_math $S/rtmlib-modC-gpu.cu -I${CUINC}
	
# Rules for dot product tests

dptest-der: $B/DPtest-der.x
	$< n1=20 n2=20 >/dev/null

$D/dptest-pml.H: $B/DPtest-pml.x
	$< n1=500 n2=500 vel=2000.0 abc=50 dt=0.0004 dx=5.0 >$@

# New dot product test. Rembember to change the Src path to Src_dotproduct2
dptest-image: $B/DPtest-imaging.x
	$< nz=100 nx=100 nt=100 lame=0 >/dev/null

test.H: $B/src_sgn.x
	$< >$@
#
dptest-rtm: $B/rtm-dp.x
	$< lame=1 >/dev/null
################################################################################

# Rules for building models
# model build
$D/lambda-gauss.H $D/lambda-gauss2.H $D/lambda-cte.H $D/lambda-vz.H $D/lambda-vz2.H: $D/lambda-%.H: $B/create-model-lame.x $P/%-lambda.p
	$B/create-model-lame.x par=$P/$*-lambda.p mu=${D}/mu-$*.H rho=${D}/rho-$*.H > $@
	Smooth < $@          rect1=8 rect2=8 > $D/lambdaS-$*.H
	Smooth < $D/mu-$*.H  rect1=8 rect2=8 > $D/muS-$*.H
	Smooth < $D/rho-$*.H rect1=8 rect2=8 > $D/rhoS-$*.H

$D/lambda-layer.H: $D/lambda-%.H: $B/create-model-lame.x $P/%-lambda.p
	$B/create-model-lame.x par=$P/$*-lambda.p mu=${D}/mu-$*.H rho=${D}/rho-$*.H > $@
	Smooth < $@          rect1=20 rect2=20 > $D/lambdaS-$*.H
	Smooth < $D/mu-$*.H  rect1=20 rect2=20 > $D/muS-$*.H
	Smooth < $D/rho-$*.H rect1=20 rect2=20 > $D/rhoS-$*.H
$D/vp-layer.H: $D/vp-%.H: $B/create-model-vp.x $P/%-vp.p
	$B/create-model-vp.x par=$P/$*-vp.p vs=${D}/vs-$*.H rho=${D}/rho-$*.H > $@
	Smooth < $@          rect1=20 rect2=20 > $D/vpS-$*.H
	Smooth < $D/vs-$*.H  rect1=20 rect2=20 > $D/vsS-$*.H
	Smooth < $D/rho-$*.H rect1=20 rect2=20 > $D/rhoS-$*.H

$D/dlambda-gauss.H $D/dlambda-gauss2.H $D/dlambda-cte.H $D/dlambda-layer.H: $D/dlambda-%.H: $D/lambda-%.H
	Add scale=1,-1 $< $D/lambdaS-$*.H > $@
	Add scale=1,-1 $D/mu-$*.H $D/muS-$*.H > $D/dmu-$*.H
	Add scale=1,-1 $D/rho-$*.H $D/rhoS-$*.H > $D/drho-$*.H
$D/dvp-gauss.H $D/dvp-gauss2.H $D/dvp-cte.H $D/dvp-layer.H: $D/dvp-%.H: $D/vp-%.H
	Add scale=1,-1 $< $D/vpS-$*.H > $@
	Add scale=1,-1 $D/vs-$*.H $D/vsS-$*.H > $D/dvs-$*.H
	Add scale=1,-1 $D/rho-$*.H $D/rhoS-$*.H > $D/drho-$*.H

$D/dlambda-marmousi.H: $D/dlambda-%.H: $D/lambda-%.H
	Add scale=1,-1 $D/lambda-$*-small.H $D/lambdaS-$*.H > $@
	Add scale=1,-1 $D/mu-$*-small.H $D/muS-$*.H > $D/dmu-$*.H
	Add scale=1,-1 $D/rho-$*-small.H $D/rhoS-$*.H > $D/drho-$*.H

$D/lambda-marmousi-small.H: $D/lambda-marmousi.H
	Interp < $< lsinc=12 type=2 o1out=0 n1out=701 d1out=5.0 o2out=0 n2out=3401 d2out=5.0 maxsize=162 >$@

$D/mu-marmousi-small.H: $D/mu-marmousi.H
	Interp < $< lsinc=12 type=2 o1out=0 n1out=701 d1out=5.0 o2out=0 n2out=3401 d2out=5.0 maxsize=162 >$@

$D/rho-marmousi-small.H: $D/rho-marmousi.H
	Interp < $< lsinc=12 type=2 o1out=0 n1out=701 d1out=5.0 o2out=0 n2out=3401 d2out=5.0 maxsize=162 >$@

$D/lambdaS-marmousi.H: $D/lambdaS-%.H: $D/lambda-marmousi-small.H
	Smooth < $<          maxmem=2000000000 repeat=10 rect1=20 rect2=20 > $D/lambdaS-$*.H
	Smooth < $D/mu-$*-small.H  maxmem=2000000000 repeat=10 rect1=20 rect2=20 > $D/muS-$*.H
	Smooth < $D/rho-$*-small.H maxmem=2000000000 repeat=10 rect1=20 rect2=20 > $D/rhoS-$*.H
	
################################################################################
# Rules for running elastic modeling and migration
#
# Nonlinear elastic modeling
elastic-vz elastic-vz2: elastic-%: $B/elastic.x ${D}/lambda-%.H ${P}/elastic-%.p
	${B}/elastic.x <${D}/lambda-$*.H mu=${D}/mu-$*.H rho=${D}/rho-$*.H par=${P}/elastic-$*.p Vx_rec=${D}/seis-Vx-$*.H Vz_rec=${D}/seis-Vz-$*.H sigmaxx_rec=${D}/seis-sigmaxx-$*.H sigmazz_rec=${D}/seis-sigmazz-$*.H sigmaxz_rec=${D}/seis-sigmaxz-$*.H >/dev/null

# Born modeling
born-gauss born-gauss2 born-cte born-layer: born-%: $B/rtm.x $D/lambda-%.H $D/dlambda-%.H $P/rtm-%.p
	$< <${D}/lambdaS-$*.H mu=${D}/muS-$*.H rho=${D}/rhoS-$*.H dlambda=${D}/dlambda-$*.H dmu=${D}/dmu-$*.H drho=${D}/drho-$*.H Vx_rec=${D}/seis-Vx-$*.H Vz_rec=${D}/seis-Vz-$*.H sigmaxx_rec=${D}/seis-sigmaxx-$*.H sigmazz_rec=${D}/seis-sigmazz-$*.H sigmaxz_rec=${D}/seis-sigmaxz-$*.H par=$P/rtm-$*.p >/dev/null adj=0 lame=1

born-vp-layer: born-vp-%: $B/rtm.x $D/vp-%.H $D/dvp-%.H $P/rtm-%.p
	$< <${D}/vpS-$*.H vs=${D}/vsS-$*.H rho=${D}/rhoS-$*.H dvp=${D}/dvp-$*.H dvs=${D}/dvs-$*.H drho=${D}/drho-$*.H Vx_rec=${D}/seis-Vx-$*.H Vz_rec=${D}/seis-Vz-$*.H sigmaxx_rec=${D}/seis-sigmaxx-$*.H sigmazz_rec=${D}/seis-sigmazz-$*.H sigmaxz_rec=${D}/seis-sigmaxz-$*.H par=$P/rtm-$*.p >/dev/null adj=0 lame=0

# RTM
rtm-gauss rtm-gauss2 rtm-cte rtm-layer: rtm-%: $B/rtm.x $D/lambda-%.H $D/dlambda-%.H $P/rtm-%.p
	$< <${D}/lambdaS-$*.H mu=${D}/muS-$*.H rho=${D}/rhoS-$*.H dlambda=${D}/dlambda-$*.H dmu=${D}/dmu-$*.H drho=${D}/drho-$*.H Vx_rec=${D}/seis-Vx-$*.H Vz_rec=${D}/seis-Vz-$*.H sigmaxx_rec=${D}/seis-sigmaxx-$*.H sigmazz_rec=${D}/seis-sigmazz-$*.H sigmaxz_rec=${D}/seis-sigmaxz-$*.H par=$P/rtm-$*.p >/dev/null adj=1 lame=1

rtm-vp-layer: rtm-vp-%: $B/rtm.x $D/vp-%.H $D/dvp-%.H $P/rtm-%.p
	$< <${D}/vpS-$*.H vs=${D}/vsS-$*.H rho=${D}/rhoS-$*.H dvp=${D}/dvp-$*.H dvs=${D}/dvs-$*.H drho=${D}/drho-$*.H Vx_rec=${D}/seis-Vx-$*.H Vz_rec=${D}/seis-Vz-$*.H sigmaxx_rec=${D}/seis-sigmaxx-$*.H sigmazz_rec=${D}/seis-sigmazz-$*.H sigmaxz_rec=${D}/seis-sigmaxz-$*.H par=$P/rtm-$*.p >/dev/null adj=1 lame=0

#rtm-marmousi-data: $B/rtm.x ${D}/lambda-marmousi.H ${D}/dlambda-marmousi.H $P/marmousi.p
#	$< <${D}/lambdaS-marmousi.H par=$P/marmousi.p mu=${D}/muS-marmousi.H rho=${D}/rhoS-marmousi.H dlambda=${D}/dlambda-marmousi.H dmu=${D}/dmu-marmousi.H drho=${D}/drho-marmousi.H >/dev/null

################################################################################
# Rules for visualizing models and data
view.gauss view.gauss2 view.cte view.layer view.vz view.vz2: view.%: $D/lambda-%.H
	Grey <Dat/lambda-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="lambda-$*" | Tube &
	Grey <Dat/mu-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="mu-$*" | Tube &
	Grey <Dat/rho-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="rho-$*" | Tube &
	Grey <Dat/lambdaS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="lambdaS-$*" | Tube &
	Grey <Dat/muS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="muS-$*" | Tube &
	Grey <Dat/rhoS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="rhoS-$*" | Tube &

view.vp.layer: view.vp.%: $D/vp-%.H
	Grey <Dat/vp-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="vp-$*" | Tube &
	Grey <Dat/vs-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="vs-$*" | Tube &
	Grey <Dat/rho-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="rho-$*" | Tube &
	Grey <Dat/vpS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="vpS-$*" | Tube &
	Grey <Dat/vsS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="vsS-$*" | Tube &
	Grey <Dat/rhoS-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="rhoS-$*" | Tube &

view.dgauss view.dgauss2 view.dlayer: view.d%: $D/dlambda-%.H
	Grey <Dat/dlambda-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="dlambda-$*" | Tube &
	Grey <Dat/dmu-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="dmu-$*" | Tube &
	Grey <Dat/drho-$*.H color=j newclip=1 epclip=100 wantscalebar=yes title="drho-$*" | Tube &

seis.layer seis.cte: seis.%:
	Grey <Dat/seis-Vx-$*.H | Tube &
	Grey <Dat/seis-Vz-$*.H | Tube &
	Grey <Dat/seis-sigmaxx-$*.H | Tube &
	Grey <Dat/seis-sigmazz-$*.H | Tube &
	Grey <Dat/seis-sigmaxz-$*.H | Tube &

image.layer: image.%: $D/dlambda-%.H
	Grey <Dat/dlambda-$*.H | Tube &
	Grey <Dat/dmu-$*.H | Tube &
	Grey <Dat/drho-$*.H | Tube &

image.vp.layer: image.vp.%: $D/dvp-%.H
	Grey <Dat/dvp-$*.H | Tube &
	Grey <Dat/dvs-$*.H | Tube &
	Grey <Dat/drho-vp-$*.H | Tube &

################################################################################
# Rules for building ER figures
INVRAS = y
G0 = Grey title=' ' label1="z(m)" label2= > /dev/null
G1 = Grey title=' ' label1="z(m)" label2="x(m)" plotfat=1 titlefat=2 labelfat=2 > /dev/null
G2 = Grey title=' ' label1="t(s)" label2="x(m)" >/dev/null
G3 = Grey title=' ' label1="z(m)" label2="x(m)" newclip=1 color=j plotfat=1 titlefat=2 labelfat=2 > /dev/null

rho.gif: $D/r_rho_time.H
	<$< ${G0} gainpanel=a screenwd=50.0 screenht=10.0 title="Cross-correlation of wave fields for the density image" out=$R/temp.v
	./vplot2gif $R/temp.v Fig/$@
	rm -f $R/temp.v

lambda.gif: $D/r_l_time.H
	<$< ${G0} gainpanel=a screenwd=50.0 screenht=10.0 title="Cross-correlation of wave fields for the lambda image" out=$R/temp.v
	./vplot2gif $R/temp.v Fig/$@
	rm -f $R/temp.v

mu.gif: $D/r_mu_time.H
	<$< ${G0} gainpanel=a screenwd=49.0 screenht=7.0 title="Cross-correlation of wave fields for the mu image" out=$R/temp.v
	./vplot2gif $R/temp.v Fig/$@
	rm -f $R/temp.v

$R/lambda-layer.v: $D/lambda-layer.H
	Math file1=$< exp='file1/1000000000.0' | ${G3} epclip=100 title="1 layer model for lambda" barlabel="GPa" wantscalebar=yes out=$@ 

$R/dlambda-layer.v: $D/dlambda-layer.H
	Math file1=$< exp='file1/1000000000.0' | ${G3} epclip=100 title="1 layer model difference for lambda" barlabel="GPa" wantscalebar=yes out=$@ 

$R/lambdaS-layer.v: $D/lambdaS-layer.H
	Math file1=$< exp='file1/1000000000.0' | ${G3} epclip=100 title="1 layer smoothed model for lambda" barlabel="GPa" wantscalebar=yes out=$@ 

#$R/r-rho.v: $D/r-rho.H
#	<$< ${G1} title="RTM image for the density model" out=$@

#$R/r-lambda.v: $D/r-l.H
#	<$< ${G1} title="RTM image for the lambda model" out=$@

#$R/r-mu.v: $D/r-mu.H
#	<$< ${G1} title="RTM image for the mu model" out=$@

$R/r-rho-vp.v:
	<$D/drho-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the density (vel) model" out=$@

$R/r-vp.v:
	<$D/dvp-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the Vp model" out=$@

$R/r-vs.v:
	<$D/dvs-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the Vs model" out=$@

$R/r-rho-lambda.v:
	< $D/drho-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the density model" out=$@

$R/r-l.v:
	<$D/dlambda-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the lambda model" out=$@

$R/r-mu.v:
	<$D/dmu-layer.H Window f1=40 n1=220 f2=40 n2=220 | ${G1} title="RTM image for the mu model" out=$@

$R/seis-P-vz.v $R/seis-P-vz2.v: $R/seis-P-%.v:
	Math file1="$D/seis-sigmaxx-$*.H" file2="$D/seis-sigmazz-$*.H" exp="file1+file2" > temp.H
	<temp.H ${G2} >/dev/null title="P-seismogram" pclip=98 out=$@
	rm -f temp.H

$R/seis-Vx-vz.v $R/seis-Vx-vz2.v: $R/seis-Vx-%.v: $D/seis-Vx-%.H
	<$< ${G2} >/dev/null title="Vx-seismogram" pclip=98 out=$@

$R/seis-Vz-vz.v $R/seis-Vz-vz2.v: $R/seis-Vz-%.v: $D/seis-Vz-%.H
	<$< ${G2} >/dev/null title="Vz-seismogram" pclip=98 out=$@

$R/vp-vz.v $R/vp-vz2.v: $R/vp-%.v: $D/lambda-%.H
	Math file1="$<" file2="$D/mu-$*.H" file3="$D/rho-$*.H" exp="@SQRT((file1+file2)/file3)" >temp.H
	<temp.H ${G3} title="Vp velocity model" wantscalebar=yes out=$@
	rm -f temp.H

$R/vs-vz.v $R/vs-vz2.v: $R/vs-%.v: $D/lambda-%.H
	Math file1="$D/mu-$*.H" file2="$D/rho-$*.H" exp="@SQRT(file1/file2)" >temp.H
	<temp.H ${G3} title="Vp velocity model" wantscalebar=yes out=$@
	rm -f temp.H

################################################################################

# Rules for building CR figures

################################################################################

# Rules for building NR figures

################################################################################

# Rules for building .gif animations (run from a Tesla machine)
snap-vz-vz.gif snap-vz-vz2.gif: snap-vz-%.gif:
	<$D/Vz0-$*.H ${G0} gainpanel=a screenwd=30.0 screenht=7.5 title="Elastic modeling in a V(z) layer model - z-velocity" out=$R/temp.v
	./vplot2gif $R/temp.v Fig/$@
	rm -f $R/temp.v

snap-vx-vz.gif snap-vx-vz2.gif: snap-vx-%.gif:
	<$D/Vx0-$*.H ${G0} gainpanel=a screenwd=30.0 screenht=7.5 title="Elastic modeling in a V(z) layer model - x-velocity" out=$R/temp.v
	./vplot2gif $R/temp.v Fig/$@
	rm -f $R/temp.v

################################################################################

# Rules for cleaning

clean: jclean

burn: clean
	rm -f ${B}/*
	rm -f ${O}/*
	rm -f ${M}/*
	rm -f ${R}/*.pdf
	rm -f ${R}/*.ps
	rm -f ${R}/*.v*
	rm -f ${D}/*

clean.v:
	rm -f ${R}/*.v*

################################################################################

	#cp $< $@	 
	#-fopenmp -I/opt/SEP/SEP7.0/include -I/opt/SEP/SEP7.0/include -g -O2 -c -DLINUX86_64 $< -DSOURCE='"'/net/server/homes/sep/lingan27/elastic/2d_rtm_kernel/Src_cuda/rtmlib-modC.c'"' -o $@

# Rules for compiling C objects (currently not being used)

#./Bin/elastic.x: ./Obj/clib-modCispc.o

#./Obj/clib-modCispc.o: ./Obj_bkp/clib-modCispc.o
#	cp $< $@



#./Obj/rtmlib-mod.o: ./Obj/rtmlib-modC.o

#./Bin/rtm.x: ./Obj/rtmlib-modC.o
	

################################################################################

# default rules
default: $R/lambda-layer.pdf $R/dlambda-layer.pdf $R/lambdaS-layer.pdf $R/r-rho.pdf $R/r-lambda.pdf $R/r-mu.pdf


SOURCES= $(notdir $(wildcard ${SRCDIR}/*.r90)  $(wildcard ${SRCDIR}/*.F90) $(wildcard ${SRCDIR}/*.rs90) $(wildcard ${SRCDIR}/*.f90) $(wildcard ${SRCDIR}/*.lop) $(wildcard ${SRCDIR}/*.dclop) $(wildcard ${SRCDIR}/*.clop) $(wildcard ${SRCDIR}/*.F89))

check:
	Add scale=1,-1 Dat/dlambda-layer.H Dat_check/dlambda-layer-check.H > result-dlambda
	Add scale=1,-1 Dat/dmu-layer.H Dat_check/dmu-layer-check.H > result-dmu
	Add scale=1,-1 Dat/drho-layer.H Dat_check/drho-layer-check.H > result-drho
	echo dlambda
	Attr < result-dlambda  
	echo dmu
	Attr < result-dmu
	echo drho
	Attr < result-drho

check-view:
	Add scale=1,-1 Dat/dlambda-layer.H Dat_fortran_no_BC/dlambda-layer-no-BC.H > result-dlambda
	Add scale=1,-1 Dat/dmu-layer.H Dat_fortran_no_BC/dmu-layer-no-BC.H > result-dmu
	Add scale=1,-1 Dat/drho-layer.H Dat_fortran_no_BC/drho-layer-no-BC.H > result-drho
	echo dlambda
	Attr < result-dlambda  
	echo dmu
	Attr < result-dmu
	echo drho
	Attr < result-drho
	Grey < result-dlambda pclip=100 | Tube
	Grey < result-dmu pclip=100 | Tube
	Grey < result-drho pclip=100 | Tube

cmp-l:
	Grey < Dat_check/dlambda-layer-check.H pclip=100 > /dev/null out=x.v
	Grey < Dat/dlambda-layer.H pclip=100 > /dev/null out=y.v
	tube x.v y.v


cmp-m:
	Grey < Dat_check/dmu-layer-check.H pclip=100 > /dev/null out=x.v
	Grey < Dat/dmu-layer.H pclip=100 > /dev/null out=y.v
	tube x.v y.v
cmp-r:
	Grey < Dat_check/drho-layer-check.H pclip=100 > /dev/null out=x.v
	Grey < Dat/drho-layer.H pclip=100 > /dev/null out=y.v
	tube x.v y.v

diff:
	Add scale=1,-1 Dat/dlambda-layer.H Dat_check/dlambda-layer-check.H > result-dlambda
	Add scale=1,-1 Dat/dmu-layer.H Dat_check/dmu-layer-check.H > result-dmu
	Add scale=1,-1 Dat/drho-layer.H Dat_check/drho-layer-check.H > result-drho
	Grey<result-dlambda pclip=100> /dev/null out=diff-dlambda.v
	Grey<result-dmu pclip=100> /dev/null out=diff-dmu.v
	Grey<result-drho pclip=100> /dev/null out=diff-drho.v
	tube diff-dlambda.v
	tube diff-dmu.v
	tube diff-drho.v


cmp:
	make cmp-l
	make cmp-m
	make cmp-r



	

include ${SEPINC}/SEP.bottom
