CUDIR = /usr/local/cuda-7.0
CULIB = ${CUDIR}/lib64 -lcudart
CUINC = ${CUDIR}/include


all: a.out

a.out:	rtmlib-modC.o rtmlib-modC-gpu.o
	g++ -fopenmp rtmlib-modC.o rtmlib-modC-gpu.o -L${CULIB} -o a.out  



rtmlib-modC.o:	rtmlib-modC.c
	g++ -c -fopenmp rtmlib-modC.c



rtmlib-modC-gpu.o: rtmlib-modC-gpu.cu
	nvcc -c rtmlib-modC-gpu.o -gencode arch=compute_35,code=sm_35 -Xptxas=-v -fmad=false -use_fast_math -maxrregcount=56 rtmlib-modC-gpu.cu -I${CUINC}



run:	
	./a.out

clean:
	rm -f *.o
	rm *.out
