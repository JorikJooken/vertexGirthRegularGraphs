# On vertex-girth-regular graphs: (Non-)existence, bounds and enumeration

This repository contains code and data related to the paper "On vertex-girth-regular graphs: (Non-)existence, bounds and enumeration". All code can be found in the directory "Code", whereas all data can be found in the directory "Data".

For integers $v$, $k$, $g$ and $\lambda$ a $vgr(v,k,g,\lambda)$-graph is a $k$-regular graph with girth $g$ on $v$ vertices such that every vertex is contained in exactly $\lambda$ cycles of length $g$.

Below, we briefly describe the different programs and data.

## DATA
The files "vertexGirthRegularGraphs.g6" contain the vertex-girth-regular graphs in graph6 format that were found in the paper. The file "commentsPerGraph.txt" indicates for each graph the parameters $v$, $k$, $g$ and $\lambda$. Finally, the file "someSmallVertexGirthRegularGraphs.g6" contains a number of small vertex-girth-regular graphs that we will use later to illustrate the code.

## CODE

This folder contains multiple files containing code for the exhaustive generation of vertex-girth-regular graphs. There are also several files containing helper functions (e.g. to support reading graphs in graph6 format, efficiently representing sets as bitsets and so on). We will only briefly discuss the code that is directly related to vertex-girth-regular graphs.

All code for generating vertex-girth-regular graphs supports graphs that have at most 64 vertices. The code in the folder "veryLargeGraphs" can handle arbitrary large graphs (as large as your RAM is able to handle).

### generateRGLambdaGraphs.c, generateRGLambdaGraphsAlsoConsiderImbalance.c

These programs are different implementations of an algorithm that can be used to generate all connected $vgr(v,k,g,\lambda)$-graphs for fixed parameters $v$, $k$, $g$ and $\lambda$.

The programs can be compiled by executing the makefile:
```bash
make
```

These programs expect the four parameters $v$, $k$, $g$ and $\lambda$ as input and will output a list of all $vgr(v,k,g,\lambda)$-graphs (one graph in graph6 format per line). For example, executing the following command:

```bash
./generateRGLambdaGraphsAlsoConsiderImbalanceExecutable 20 3 5 1
```
 will output a list of all $vgr(20,3,5,1)$ graphs:

```bash
Computation started for v, k, g and lambda: 20 3 5 1
SsP@P?OC?P?O?E?W??GA?@@G?OoGC??Gc
SsP@P?OC?P?O?E?W??GA?@@G?OoG?_?KC
SsP@P?OC?P?O?E?W??GA?@@G?AoG_??Gc
SsP@P?OC?P?O?E?W??GA?@@G?AoG?_?gC
SsP@P?OC?P?O?E?W??GA??CK?C_GG_GCG
SsP@P?OC?P?O?E?G?GG?A?OcC?_?gGGG_
Computation finished for v, k, g and lambda: 20 3 5 1
```

The other programs can be invoked by changing the name of the executable. For example, invoking the command:

```bash
./generateRGLambdaGraphsExecutable 10 4 4 12
```

will produce an exhaustive list of all $vgr(10,4,4,12)$ graphs (there is in fact one such graph):
```bash
Computation started for v, k, g and lambda: 10 4 4 12
Is`rQow@w
Computation finished for v, k, g and lambda: 10 4 4 12
```

The different programs can have different speeds, depending on the precise values of $v$, $k$, $g$ and $\lambda$. For large orders, we advice you to first empirically check which program is the fastest one for smaller orders and then use this program for the large order as well.

### veryLargeGraphs/calcLambdaRegular.cpp

This program can be compiled by executing the following command:
```bash
g++ -g -std=c++11 -O3 calcLambdaRegular.cpp -o calcLambdaRegularExecutable
```

This program expects as input a list of connected regular graphs in graph6 format (one graph per line). It will output for each graph the parameters $v$, $k$, $g$ and $\lambda$ (where $\lambda$ will be indicated as "-1" if the graph is not vertex-girth-regular). For example, executing the following command:

```bash
./calcLambdaRegularExecutable < ../../Data/someSmallVertexGirthRegularGraphs.g6
```

will produce the following output:
 ```bash
H{_yqgj
9 4 3 1
Os?GOO??AHGQGSDCAK@D?
16 3 6 9
Js`@IStU`w?
11 4 4 8
```

This indicates that for example the graph 

```bash
H{_yqgj
```

is a 9-vertex graph with girth 3 in which every vertex has degree 4 and each vertex is contained in 1 cycle of length 3 (i.e. a $vgr(9,4,3,1)$ graph).
