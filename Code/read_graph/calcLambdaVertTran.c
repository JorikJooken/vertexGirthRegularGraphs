#include "bitset.h"
#include <stdio.h>
#include "read_graph/readGraph6.h"

#define MAXGORDER 256
bitset adjacencyListG[MAXGORDER];
int nVerticesG;


char * graphString = NULL;
size_t size;

int vertexQ[MAXVERTICES];
int qStart;
int qFinish;

int distBFS1[MAXVERTICES];
int numberShortestPaths[MAXVERTICES];

int n;

void calcNumberShortestPaths(int u, int v, int* length, int* numPaths)
{
    removeElement(adjacencyListG[u],v);
    removeElement(adjacencyListG[v],u);
    //fprintf(stderr,"function called with u and v: %d %d\n",u,v);
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
        forEach(neigh,adjacencyListG[now])
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
    add(adjacencyListG[u],v);
    add(adjacencyListG[v],u);
    //fprintf(stderr,"Distance between u and v: %d\n",distBFS1[v]);
    (*length)=distBFS1[v];
    (*numPaths)=numberShortestPaths[v];
}
// reads graphs F and calculates the number of girth cycles through each edge
// if this number is the same for every edge, it will output this number and otherwise -1
int main(int argc, char ** argv) {
    while(getline(&graphString, &size, stdin) != -1) {
        nVerticesG = getNumberOfVertices(graphString);
        n=nVerticesG;
        loadGraph(graphString, nVerticesG, adjacencyListG);

        int lambda=-2;
        int girth=-2;
        for(int u=0; u<=0 && lambda!=-1; u++)
        {
            forEach(v,adjacencyListG[u])
            {
                int newLen=-1;
                int newNumPaths=-1;
                calcNumberShortestPaths(u, v, &newLen, &newNumPaths);
                if(lambda==-2)
                {
                    lambda=newNumPaths;
                    girth=newLen+1;
                    continue;
                }
                if(!(newLen+1==girth && newNumPaths==lambda))
                {
                    lambda=-1;
                    break;
                }
            }
        }
        printf("%s",graphString);
        printf("%d %d %d %d\n",n,size(adjacencyListG[0]),girth,lambda);
    }
    return 0;
}
