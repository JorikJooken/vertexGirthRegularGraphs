/**
 * generateKVertexCriticalGraphs.c
 * 
 * Author: Jorik Jooken (jorik.jooken@kuleuven.be)
 * 
 */

// This program contains parts of the generator for K2-hypohamiltonian graphs https://github.com/JarneRenders/GenK2Hypohamiltonian

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "bitset.h"
#include "nauty.h"
#include "read_graph/readGraph6.h"

#define USAGE \
"\nUsage:./generateVertexCriticalGraphsC5Coloring [-l] [-k] [-h] n\n\n"

#define HELPTEXT \
"Generate all vertex-critical graphs for C5-coloring of order at most n.\n\
Graphs are sent to stdout in graph6 format. For more information on the format,\n\
see http://users.cecs.anu.edu.au/~bdm/data/formats.txt.\n\
\n\
\n\
Arguments.\n\
    -l, --length                number of vertices in forbidden path\n\
    -k, --coloring              (k+1)-vertex-critical graphs \n\
    -h, --help                  print help message\n"

//  Macro's for nauty representation
#define FOREACH(element,nautySet)\
 for(int element = nextelement((nautySet),MAXM,-1); (element) >= 0;\
 (element) = nextelement((nautySet),MAXM,(element)))
#define REMOVEONEEDGE(g,i,j,MAXM)\
 DELELEMENT(GRAPHROW(g,i,MAXM),j); DELELEMENT(GRAPHROW(g,j,MAXM),i)

// Graph structure containing nauty and bitset representations.
struct graph {
    graph nautyGraph[MAXN];
    int numberOfVertices;
    int numberOfEdges;
    bitset* adjacencyList;
    bitset* verticesOfDeg;
}; 

//  Struct for passing options.
struct options {
    int pathLength;
    int modulo;
    int remainder;
};

// Data structures for algorithm

int n, r, girthRequired, lambda;
int level[MAXVERTICES];
int vertexQ[MAXVERTICES];
int qStart;
int qFinish;
int degree[MAXVERTICES];
int nextIsolatedVertex;
int startLeaf;
bitset verticesWithCorrectDegree;
int dist[MAXVERTICES][MAXVERTICES];
bitset validNeighbors[MAXVERTICES];

int distBFS1[MAXVERTICES];
int distBFS2[MAXVERTICES];

int numberShortestPaths[MAXVERTICES];

#define MAXR 15
bitset oldValidNeighbors1[MAXR*MAXVERTICES][MAXVERTICES]; //[iteration][vertex]
bitset oldValidNeighbors2[MAXR*MAXVERTICES][MAXVERTICES];

int oldDist[MAXR*MAXVERTICES][MAXVERTICES][MAXVERTICES];
int endpoint1[MAXR*MAXVERTICES][MAXR*MAXVERTICES];
int endpoint2[MAXR*MAXVERTICES][MAXR*MAXVERTICES];

int tempArrayStoringCounts[MAXVERTICES];

/************************************************************************************
 * 
 * 
 *                      Macro's for dealing with graphs
 * 
 * 
 ************************************************************************************/

//  Initializer for empty graph.
#define emptyGraph(g) EMPTYGRAPH((g)->nautyGraph, (g)->numberOfVertices, MAXM);\
 (g)->verticesOfDeg[0] = compl(EMPTY,(g)->numberOfVertices);\
 (g)->numberOfEdges = 0;\
 for(int i = 0; i < (g)->numberOfVertices; i++) {(g)->adjacencyList[i] = EMPTY;} 

//  Add one edge. Assume that every edge has already deg2 or more.
#define addEdge(g,i,j) {ADDONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);\
 (g)->numberOfEdges++;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) - 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) - 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
 degree[(i)]++;\
 degree[(j)]++;\
}

//  Remove one edge.
#define removeEdge(g,i,j) {REMOVEONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 removeElement((g)->adjacencyList[i], j); removeElement((g)->adjacencyList[j],i);\
 (g)->numberOfEdges--;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) + 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) + 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
 degree[(i)]--;\
 degree[(j)]--;\
}

#define areNeighbours(g,i,j) contains((g)->adjacencyList[(i)], j)


/************************************************************************************
 * 
 * 
 *                      Definitions for Nauty's splay tree
 * 
 * 
 ************************************************************************************/
typedef struct SPLAYNODE {
    graph* canonForm;
    struct SPLAYNODE *left, *right, *parent;
} SPLAYNODE;

SPLAYNODE *splayTreeArray[MAXN*(MAXN-1)/2] = {NULL};

#define SCAN_ARGS

#define ACTION(p) 

#define INSERT_ARGS , graph gCan[], int numberOfVertices, bool *isPresent

int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices);

#define COMPARE(p) compareSplayNodeToGraph(p, gCan, numberOfVertices);

#define PRESENT(p) {(*isPresent) = true;}

#define NOT_PRESENT(p) {p->canonForm = gCan; (*isPresent) = false;}

#define LOOKUP_ARGS , graph gCan[], int numberOfVertices

#include "splay.c"


/************************************************************************************
 * 
 * 
 *                          Functions for printing data 
 * 
 * 
 ************************************************************************************/

void printGraph(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g->nautyGraph, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printNautyGraph(graph g[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printAdjacencyList(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        forEach(neighbour, g->adjacencyList[i]) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

//  Print gCan to stdout in graph6 format. 
void writeToG6(graph gCan[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's 
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(ISELEMENT(GRAPHROW(gCan, i, MAXM), j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}

/************************************************************************************/

//  Uses Nauty to get canonical form.
void createCanonicalForm(graph g[], graph gCan[], int numberOfVertices) {

    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    statsblk stats;

    densenauty(g, lab, ptn, orbits, &options, &stats, MAXM,
     numberOfVertices, gCan);
}

//  Splaynode contains the canonical form of a graph checked earlier.
int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices) {
    return memcmp(p->canonForm, gCan, numberOfVertices * sizeof(graph));
}

/*
void someFunction(struct graph* g, int numberOfVertices)
{
    bool isPresent=false;
    graph* gCan = malloc(sizeof(graph)*(numberOfVertices));
    if(gCan == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
    createCanonicalForm(g->nautyGraph, gCan, numberOfVertices);
    SPLAYNODE *sp=splay_lookup(&splayTreeArray[g->numberOfEdges], gCan, numberOfVertices);

    sp=splay_lookup(&splayTreeArray[g->numberOfEdges], gCan, numberOfVertices);
    if(sp==NULL)
    {
        fprintf(stderr,"It was not inside!\n");
    }
    else
    {
        fprintf(stderr,"It was inside!\n");
    }

    splay_insert(&splayTreeArray[g->numberOfEdges], gCan, numberOfVertices, &isPresent);
	//if(isPresent)
	//{
	//    free(gCan);
	//}
    fprintf(stderr,"is Present: %d\n",isPresent);
    sp=splay_lookup(&splayTreeArray[g->numberOfEdges], gCan, numberOfVertices);
    if(sp==NULL)
    {
        fprintf(stderr,"It was not inside!\n");
    }
    else
    {
        fprintf(stderr,"It was inside!\n");
    }    
}*/

void buildMooreTree(struct graph* g)
{
    qStart=0;
    qFinish=-1;
    vertexQ[qFinish+1]=0;
    qFinish++;
    level[0]=0;
    nextIsolatedVertex=1;
    if(girthRequired%2==0)
    {
        addEdge(g,0,1);
        level[1]=0;
        vertexQ[qFinish+1]=nextIsolatedVertex;
        qFinish++;
        nextIsolatedVertex++;
    }
    while(true)
    {
        int now=vertexQ[qStart];
        qStart++;
        if(level[now]>=(girthRequired-1)/2) 
        {
            startLeaf=now;
            break;
        }
        while(degree[now]<r)
        {
            addEdge(g,now,nextIsolatedVertex);
            level[nextIsolatedVertex]=level[now]+1;
            vertexQ[qFinish+1]=nextIsolatedVertex;
            qFinish++;
            nextIsolatedVertex++;
        }
    }    
}

void calcDistances(struct graph* g, int start)
{
    qStart=0;
    qFinish=0;
    vertexQ[qStart]=start;
    dist[start][start]=0;
    while(qStart<=qFinish)
    {
        int now=vertexQ[qStart];
        qStart++;
        forEach(neigh,g->adjacencyList[now])
        {
            if(dist[start][neigh]>n)
            {
                dist[start][neigh]=dist[start][now]+1;
                vertexQ[qFinish+1]=neigh;
                qFinish++;
            }
        }
    }
}

void calcNumberShortestPaths(struct graph *g, int u, int v, int* length, int* numPaths)
{
    //fprintf(stderr,"function called with u and v: %d %d\n",u,v);
    removeEdge(g,u,v);
    qStart=0;
    qFinish=0;
    vertexQ[qStart]=u;
    for(int i=0; i<n; i++) distBFS1[i]=n+100;
    numberShortestPaths[u]=1;
    distBFS1[u]=0;
    while(qStart<=qFinish)
    {
        int now=vertexQ[qStart];
        //fprintf(stderr,"now: %d\n",now);
        if(now==v) break;
        qStart++;
        forEach(neigh,g->adjacencyList[now])
        {
            //fprintf(stderr,"now and neigh: %d %d\n",now,neigh);
            if(distBFS1[neigh]>n) // new
            {
                distBFS1[neigh]=distBFS1[now]+1;
                numberShortestPaths[neigh]=numberShortestPaths[now];
                qFinish++;
                vertexQ[qFinish]=neigh;
            }
            else if(distBFS1[neigh]==distBFS1[now]+1) // tie
            {
                numberShortestPaths[neigh]+=numberShortestPaths[now];
            }
        }
    }
    addEdge(g,u,v);
    //fprintf(stderr,"Distance between u and v: %d\n",distBFS1[v]);
    (*length)=distBFS1[v];
    (*numPaths)=numberShortestPaths[v];
}

int nbInSplayTree=0;
// this parameter should be set differently depending on the system and depending on MAXN
// I experimentally set it to 5 Million, which seems OK on my machine (around 15% of available memory used)
#define MAXNBINSPLAYTREE 5000000
bool warningAlreadySent=false;

SPLAYNODE *sp=NULL;

void recursivelyAddEdges(struct graph* g, int iteration)
{
    /*fprintf(stderr,"Entered at iteration %d\n",iteration);
    for(int i=0; i<n; i++)
    {
        fprintf(stderr,"Neighs of %d:\n",i);
        forEach(neigh,g->adjacencyList[i])
        {
            fprintf(stderr," %d",neigh);
        }
        fprintf(stderr,"\n");
    }
    for(int i=0; i<n; i++)
    {
        fprintf(stderr,"Valid neighs of %d:\n",i);
        forEach(neigh,validNeighbors[i])
        {
            fprintf(stderr," %d",neigh);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"Next isolated vertex: %d\n",nextIsolatedVertex);*/

    graph* gCan = malloc(sizeof(graph)*(nextIsolatedVertex));
    createCanonicalForm(g->nautyGraph, gCan, nextIsolatedVertex);
    // careful, do not free this pointer because it will either point to NULL or to an existing node in the splay tree
    sp=splay_lookup(&splayTreeArray[g->numberOfEdges], gCan, nextIsolatedVertex);
    // this graph was already seen    
    if(sp!=NULL)
    {
        //fprintf(stderr,"Was already seen!\n");
        free(gCan);
        return;
    }
    if(nbInSplayTree<MAXNBINSPLAYTREE || isEmpty(compl(verticesWithCorrectDegree,n)))
    {
        bool isPresent;
        nbInSplayTree++;
        splay_insert(&splayTreeArray[g->numberOfEdges], gCan, nextIsolatedVertex, &isPresent);
    }
    else
    {
        free(gCan);
        if(!warningAlreadySent)
        {
            warningAlreadySent=true;
            fprintf(stderr,"Not adding anything more to the splay tree, memory is almost full!\n");
        }
    }

    // all vertices have the correct degree, write graph to output and return
    // for biregular case: careful about return, maybe not allowed
    if(isEmpty(compl(verticesWithCorrectDegree,n)))
    {
        bool shouldOutput=true;
        for(int i=0; i<n && shouldOutput; i++)
        {
            long long sum=0;
            forEach(neigh,g->adjacencyList[i])
            {
                int lengthPath;
                int numPaths;
                calcNumberShortestPaths(g,i,neigh,&lengthPath,&numPaths);
                //fprintf(stderr,"i, neigh, lengthPath and numPaths: %d %d %d %d\n",i,neigh,lengthPath,numPaths);
                if(1+lengthPath==girthRequired)
                {
                    sum+=numPaths;
                    if(sum>2*lambda)
                    {
                        shouldOutput=false;
                        break;
                    }
                }
            }
            if(sum!=2*lambda)
            {
                shouldOutput=false;
                break;
            }
        }
        if(shouldOutput)
        {
            //fprintf(stderr,"Writing graph6 string\n");
            writeToG6(g->nautyGraph,n);
        }
        return;
    }
    
    int leastNumberOptions=n+5;
    int argMinn=-1;
    // check which vertex has least number of options
    forEach(i,compl(verticesWithCorrectDegree,n))
    {
        if(i>nextIsolatedVertex) break;
        int nbOptions=size(validNeighbors[i]);
        if(nbOptions<leastNumberOptions)
        {
            leastNumberOptions=nbOptions;
            argMinn=i;
        }
    }
  
    oldValidNeighbors1[iteration][argMinn]=validNeighbors[argMinn];

    bitset valChanged1=singleton(argMinn); // for which vertices x did we change validNeighbors[x]
    // add edge between argMinn and neigh
    bitset setToIterateOver=validNeighbors[argMinn];
    forEach(neigh,setToIterateOver)
    {
        if(neigh>nextIsolatedVertex) break;
        bitset valChanged2=singleton(argMinn);
        oldValidNeighbors2[iteration][argMinn]=validNeighbors[argMinn];
        oldValidNeighbors2[iteration][neigh]=validNeighbors[neigh];
        if(!contains(valChanged1,neigh))
        {
            oldValidNeighbors1[iteration][neigh]=validNeighbors[neigh];
        }
        add(valChanged1,neigh);
        add(valChanged2,neigh);
        addEdge(g,argMinn,neigh);
        //fprintf(stderr,"Added edge between %d and %d\n",argMinn,neigh);
        removeElement(validNeighbors[argMinn],neigh);
        removeElement(validNeighbors[neigh],argMinn);
        
        // update nextIsolatedVertex
        bool nextIncremented=false;
        if(neigh==nextIsolatedVertex || argMinn==nextIsolatedVertex)
        {
            nextIsolatedVertex++;
            nextIncremented=true;
        }

        bool shouldRecurse=true;
        // update verticesWithCorrectDegree and validNeighbors
        if(degree[argMinn]==r)
        {
            add(verticesWithCorrectDegree,argMinn);
            forEach(neigh2,validNeighbors[argMinn])
            {
                if(!contains(valChanged1,neigh2))
                {
                    oldValidNeighbors1[iteration][neigh2]=validNeighbors[neigh2];
                }
                if(!contains(valChanged2,neigh2))
                {
                    oldValidNeighbors2[iteration][neigh2]=validNeighbors[neigh2];
                }
                add(valChanged1,neigh2);
                add(valChanged2,neigh2);
                removeElement(validNeighbors[neigh2],argMinn);
                if((degree[neigh2]+size(validNeighbors[neigh2])<r)) shouldRecurse=false;
            }
            if(!contains(valChanged1,argMinn))
            {
                oldValidNeighbors1[iteration][argMinn]=validNeighbors[argMinn];
            }
            if(!contains(valChanged2,argMinn))
            {
                oldValidNeighbors2[iteration][argMinn]=validNeighbors[argMinn];
            }
            add(valChanged1,argMinn);
            add(valChanged2,argMinn);
            validNeighbors[argMinn]=EMPTY;
        }   
        if(degree[neigh]==r)
        {
            add(verticesWithCorrectDegree,neigh);
            forEach(neigh2,validNeighbors[neigh])
            {
                if(!contains(valChanged1,neigh2))
                {
                    oldValidNeighbors1[iteration][neigh2]=validNeighbors[neigh2];
                }
                if(!contains(valChanged2,neigh2))
                {
                    oldValidNeighbors2[iteration][neigh2]=validNeighbors[neigh2];
                }
                add(valChanged1,neigh2);
                add(valChanged2,neigh2);
                removeElement(validNeighbors[neigh2],neigh);
                if((degree[neigh2]+size(validNeighbors[neigh2])<r)) shouldRecurse=false;
            }
            if(!contains(valChanged1,neigh))
            {
                oldValidNeighbors1[iteration][neigh]=validNeighbors[neigh];
            }
            if(!contains(valChanged2,neigh))
            {
                oldValidNeighbors2[iteration][neigh]=validNeighbors[neigh];
            }
            add(valChanged1,neigh);
            add(valChanged2,neigh);
            validNeighbors[neigh]=EMPTY;
        }
 
        int endpointPtr=-1;
        // update other stuff, e.g. dist and validNeighbors
        if(shouldRecurse)
        {
            for(int i=0; i<n; i++)
            {
                distBFS1[i]=distBFS2[i]=n+5;
            }
            //fprintf(stderr,"ok1\n");
            // BFS from argMinn
            qStart=0;
            qFinish=0;
            vertexQ[0]=argMinn;
            distBFS1[argMinn]=0;
            while(qStart<=qFinish)
            {
                int now=vertexQ[qStart];
                //fprintf(stderr,"now: %d\n",now);
                qStart++;
                forEach(neigh2,g->adjacencyList[now])
                {
                    if(distBFS1[neigh2]>=n)
                    {
                        distBFS1[neigh2]=distBFS1[now]+1;
                        vertexQ[qFinish+1]=neigh2;
                        qFinish++;
                    }
                }
            }   
            //fprintf(stderr,"ok2\n");
            // BFS from neigh
            qStart=0;
            qFinish=0;
            vertexQ[0]=neigh;
            distBFS2[neigh]=0;
            while(qStart<=qFinish)
            {
                int now=vertexQ[qStart];
                qStart++;
                forEach(neigh2,g->adjacencyList[now])
                {
                    if(distBFS2[neigh2]>=n)
                    {
                        distBFS2[neigh2]=distBFS2[now]+1;
                        vertexQ[qFinish+1]=neigh2;
                        qFinish++;
                    }
                }
            }

            //fprintf(stderr,"ok3\n");
            forEach(u,compl(verticesWithCorrectDegree,n))
            {
                forEach(v,validNeighbors[u])
                {
                    if(u==v) continue;
                    int newDist=dist[u][v];
                    if(1+distBFS1[u]+distBFS2[v]<newDist) newDist=1+distBFS1[u]+distBFS2[v];
                    if(1+distBFS2[u]+distBFS1[v]<newDist) newDist=1+distBFS2[u]+distBFS1[v];
                    if(newDist<dist[u][v])
                    {
                        endpointPtr++;
                        endpoint1[iteration][endpointPtr]=u;
                        endpoint2[iteration][endpointPtr]=v;
                        oldDist[iteration][u][v]=dist[u][v];
                        dist[u][v]=newDist;
                        if(1+dist[u][v]<girthRequired)
                        {
                            if(!contains(valChanged1,u))
                            {
                                oldValidNeighbors1[iteration][u]=validNeighbors[u];
                            }
                            if(!contains(valChanged2,u))
                            {
                                oldValidNeighbors2[iteration][u]=validNeighbors[u];
                            }
                            add(valChanged1,u);
                            add(valChanged2,u);
                            removeElement(validNeighbors[u],v);
                            if((degree[u]+size(validNeighbors[u])<r)) shouldRecurse=false;
                        }
                    }
                }
            }
        }

        // After the edge is added, for many vertices the number of girth cycles on which they lie can change
        // In this if, I do a cheap check and only check if it is violated for vertex "argMinn" or vertex "neigh"
        // Alternatively, we could check it for all vertices, or some vertices (trade-off between pruning and cost to check pruning condition)
        //
        // In any case, before the graph is output, this condition is checked for all vertices to ensure
        // that the graph is vertex-girth-regular
        if(shouldRecurse)
        {
            int sumForArgMinn=0;
            forEach(neigh2,g->adjacencyList[argMinn])
            {
                int lengthPath;
                int numPaths;
                calcNumberShortestPaths(g,argMinn,neigh2,&lengthPath,&numPaths);
                if(1+lengthPath==girthRequired)
                {
                    sumForArgMinn+=numPaths; 
                }
                if(sumForArgMinn>2*lambda)
                {
                    shouldRecurse=false;
                    break;
                }
            }
            if(shouldRecurse)
            {
                int sumForNeigh=0;
                forEach(neigh2,g->adjacencyList[neigh])
                {
                    int lengthPath;
                    int numPaths;
                    calcNumberShortestPaths(g,neigh,neigh2,&lengthPath,&numPaths);
                    if(1+lengthPath==girthRequired)
                    {
                        sumForNeigh+=numPaths; 
                    }
                    if(sumForNeigh>2*lambda)
                    {
                        shouldRecurse=false;
                        break;
                    }
                }
            }
        }

        // check whether there exists some vertex which is "imbalanced": if it still needs to lie on x additional girth cycles, some conditions
        // must hold on the number of girth cycles that its neighbors lie on, because each of these x additional girth cycles contains at least
        // two of these neighbors

        // Remark that we can only check this condition for vertices that already have the right degree!
        if(shouldRecurse)
        {
            for(int i=0; i<n; i++) tempArrayStoringCounts[i]=0;
            for(int i=0; i<n; i++)
            {
                forEach(neigh2,g->adjacencyList[i])
                {
                    if(i>neigh2) continue;

                    int lengthPath;
                    int numPaths;
                    calcNumberShortestPaths(g,i,neigh2,&lengthPath,&numPaths);
                    if(1+lengthPath==girthRequired)
                    {
                        tempArrayStoringCounts[i]+=numPaths;
                        tempArrayStoringCounts[neigh2]+=numPaths;  
                    }
                }
            }
            for(int i=0; i<n; i++) tempArrayStoringCounts[i]/=2;


            for(int i=0; i<n && shouldRecurse; i++)
            {
                tempArrayStoringCounts[i]=lambda-tempArrayStoringCounts[i];
                if(tempArrayStoringCounts[i]<0) shouldRecurse=false;
            }
            
            // check some sum requirement on neighbors
            forEach(i,verticesWithCorrectDegree) 
            {
                //fprintf(stderr,"i and should be on how many more cycles: %d %d\n",i,tempArrayStoringCounts[i]);
                int sumOfNeighbors=0;
                forEach(neigh2,g->adjacencyList[i])
                {
                    sumOfNeighbors+=tempArrayStoringCounts[neigh2];
                }
                if(sumOfNeighbors<2*tempArrayStoringCounts[i])
                {
                    shouldRecurse=false;
                    break;
                }
            }
            // check some balance requirement on neighbors
            // I think this condition is easy in fact: if the neighbors still need x_1, x_2, ..., x_k additional cycles
            // there are only 2 cases: 
            // Case 1: max(x_1,x_2,...,x_k)>=sum(x_1,x_2,...,x_k)-max(x_1,x_2,...,x_k). In this case the maximum number of moves we can do is
            // clearly equal to sum(x_1,x_2,...,x_k)-max(x_1,x_2,...,x_k).
            // Case 2: max(x_1,x_2,...,x_k)<sum(x_1,x_2,...,x_k)-max(x_1,x_2,...,x_k). In this case one can use induction to show that the maximum
            // number of moves one can do is equal to floor(sum(x_1,x_2,...,x_k)/2) (this is an ovious upper bound and one can show that it is also a lower bound by induction; by subtracting 1 from the 2 highest numbers we again end up with Case 2 or Case 1 where all numbers are 0, except for perhaps 1 which can be 1).

            // (Every time i lies on an additional girth cycle, we need to subtract two times 1 from x_j and x_k for j not equal to k)
            forEach(i,verticesWithCorrectDegree) 
            {
                int maxAmongNeighbors=-1;
                int sumOfNeighbors=0;
                forEach(neigh2,g->adjacencyList[i])
                {
                    sumOfNeighbors+=tempArrayStoringCounts[neigh2];
                    if(tempArrayStoringCounts[neigh2]>maxAmongNeighbors) maxAmongNeighbors=tempArrayStoringCounts[neigh2];
                }
                int nb_moves_possible=0;
                if(maxAmongNeighbors>=sumOfNeighbors-maxAmongNeighbors) nb_moves_possible=sumOfNeighbors-maxAmongNeighbors;
                else
                {
                    nb_moves_possible=(sumOfNeighbors/2);
                }
                if(nb_moves_possible<tempArrayStoringCounts[i])
                {
                    shouldRecurse=false;
                    break;
                }
            }
        }
        // recurse
        if(shouldRecurse)
        {
            //fprintf(stderr,"Added edge between %d and %d\n",argMinn,neigh);
            recursivelyAddEdges(g,iteration+1);
        }

        // undo state change
        // verticesWithCorrectDegree
        if(degree[argMinn]==r)
        {
            removeElement(verticesWithCorrectDegree,argMinn); 
        }    
        if(degree[neigh]==r)
        {
            removeElement(verticesWithCorrectDegree,neigh);
        }

        // validNeighbors
        forEach(u,valChanged2)
        {
            validNeighbors[u]=oldValidNeighbors2[iteration][u];
        }
        removeElement(validNeighbors[argMinn],neigh);
        removeElement(validNeighbors[neigh],argMinn);
        removeEdge(g,argMinn,neigh);
        if(nextIncremented) nextIsolatedVertex--;

        // dist
        for(int i=0; i<=endpointPtr; i++)
        {
            int u=endpoint1[iteration][i];
            int v=endpoint2[iteration][i];
            dist[u][v]=oldDist[iteration][u][v];
        }
        // note that validNeighbors[argMinn] and validNeighbors[neigh] are not changed back again to its original form, because
        // we demand that edges are added in sequential order: if we decide not to add an edge between x and y, then
        // we never want to add that edge in deeper recursion levels either

        // too few edges left
        if((degree[argMinn]+size(validNeighbors[argMinn])<r) ||(degree[neigh]+size(validNeighbors[neigh])<r))
        {
            forEach(u,valChanged1)
            {
                validNeighbors[u]=oldValidNeighbors1[iteration][u];
            }
            return;
        }
    }
    forEach(u,valChanged1)
    {
        validNeighbors[u]=oldValidNeighbors1[iteration][u];
    }
}

// command line arguments: n, r, girthRequired and lambda
int main(int argc, char ** argv) {

    if(argc!=5)
    {
        fprintf(stderr,"Wrong number of command line arguments!\nExpected n, r and girthRequired\n");
    }
    n=atoi(argv[1]);
    r=atoi(argv[2]);
    girthRequired=atoi(argv[3]);
    lambda=atoi(argv[4]);
    fprintf(stderr,"Computation started for v, k, g and lambda: %d %d %d %d\n",n,r,girthRequired,lambda);
    if(r>MAXR)
    {
        fprintf(stderr,"r is too big!\n");
        exit(0);
    }
    //fprintf(stderr,"n, r, girthRequired and lambda: %d %d %d %d\n",n,r,girthRequired,lambda);
    struct graph g = {.numberOfVertices = n};
    g.adjacencyList = malloc(sizeof(bitset)*n);
    g.verticesOfDeg = malloc(sizeof(bitset)*n);
    for(int i = 0; i < g.numberOfVertices; i++) {
        g.verticesOfDeg[i] = EMPTY;
    }
    emptyGraph(&g);
    for(int i=0; i<n; i++) degree[i]=0;
    buildMooreTree(&g);
    verticesWithCorrectDegree=EMPTY;
    for(int i=0; i<startLeaf; i++)
    {
        add(verticesWithCorrectDegree,i);
    }
    /*for(int i=0; i<n; i++)
    {
        fprintf(stderr,"Neighbors of %d:\n",i);
        forEach(neigh,g.adjacencyList[i])
        {
            fprintf(stderr," %d",neigh);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"StartLeaf: %d\n",startLeaf);*/
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            dist[i][j]=n+100;
        }
    }
    for(int i=0; i<n; i++)
    {
        calcDistances(&g, i);
    }
    /*for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
                fprintf(stderr,"dist between %d and %d: %d\n",i,j,dist[i][j]);
        }
    }*/
    for(int i=0; i<n; i++)
    {
        validNeighbors[i]=EMPTY;
        if(i<startLeaf) continue;
        for(int j=startLeaf; j<n; j++)
        {
            if(i==j) continue;
            if(contains(g.adjacencyList[i],j)) continue;
            if(1+dist[i][j]>=girthRequired) add(validNeighbors[i],j);
        }
    }
    recursivelyAddEdges(&g,0);
    /*for(int i=0; i<n;i++)
    {
        if(contains(verticesWithCorrectDegree,i))
        {
            fprintf(stderr,"Correct degree: %d\n",i);
        }
    }*/
    /*for(int i=0; i<n; i++)
    {
        fprintf(stderr,"Valid neighs of %d:\n",i);
        forEach(neigh,validNeighbors[i])
        {
            fprintf(stderr," %d",neigh);
        }
        fprintf(stderr,"\n");
    }*/
    fprintf(stderr,"Computation finished for v, k, g and lambda: %d %d %d %d\n",n,r,girthRequired,lambda);
    free(g.adjacencyList);
    free(g.verticesOfDeg);
    return 0;
}
