compiler=gcc
flags=-std=gnu11 -march=native
densenauty32=-DWORDSIZE=32 -DMAXN=WORDSIZE nautyW1.a
densenauty64=-DWORDSIZE=64 -DMAXN=WORDSIZE nautyL1.a

all: genRGLambda genRGLambdaImbalance


genRGLambda: generateRGLambdaGraphs.c bitset.h read_graph/readGraph6.c
	$(compiler) -o generateRGLambdaGraphsExecutable generateRGLambdaGraphs.c read_graph/readGraph6.c $(densenauty64) $(flags) -O4

genRGLambdaImbalance: generateRGLambdaGraphsAlsoConsiderImbalance.c bitset.h read_graph/readGraph6.c
	$(compiler) -o generateRGLambdaGraphsAlsoConsiderImbalanceExecutable generateRGLambdaGraphsAlsoConsiderImbalance.c read_graph/readGraph6.c $(densenauty64) $(flags) -O4

.PHONY: clean
clean:
	rm -f generateRGLambdaGraphsExecutable generateRGLambdaGraphsAlsoConsiderImbalanceExecutable
