#include <bits/stdc++.h>

// Unsafe because no defined behaviour if character = 0. ctz and clz work with 32 bit numbers.
#define unsafePrev(character, current) (__builtin_ctz(character) - current >= 0 ? -1 : current -__builtin_clz((character) << (32 - current)) - 1)

#define prev(character,current) (character ? unsafePrev(character,current) : -1)

using namespace std;

const int nb_bits=4096;
#define MAXVERTICES 4096

int vertexQ[MAXVERTICES];
int qStart;
int qFinish;

int distBFS1[MAXVERTICES];
int numberShortestPaths[MAXVERTICES];

int n;
vector< vector<int> > graph;

void calcNumberShortestPaths(int u, int v, int* length, int* numPaths)
{
    /*auto it = find(graph[u].begin(), graph[u].end(), v);
    graph[u].erase(it);
    it = find(graph[v].begin(), graph[v].end(), u);
    graph[v].erase(it);*/
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
        for(int neigh : graph[now])
        {
            if(now==u && neigh==v) continue;
            if(now==v && neigh==u) continue;
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
    /*graph[u].push_back(v);
    graph[v].push_back(u);*/
    //fprintf(stderr,"Distance between u and v: %d\n",distBFS1[v]);
    (*length)=distBFS1[v];
    (*numPaths)=numberShortestPaths[v];
}

int getNumberOfVertices(string graphString) 
{
	if(graphString.size() == 0){
        printf("Error: String is empty.\n");
        abort();
    }
    else if((graphString[0] < 63 || graphString[0] > 126) && graphString[0] != '>') {
    	printf("Error: Invalid start of graphstring.\n");
    	abort();
    }

	int index = 0;
	if (graphString[index] == '>') { // Skip >>graph6<< header.
		index += 10;
	}

	if(graphString[index] < 126) { // 0 <= n <= 62
		return (int) graphString[index] - 63;
	}

	else if(graphString[++index] < 126) { 
		int number = 0;
		for(int i = 2; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	else if (graphString[++index] < 126) {
		int number = 0;
		for (int i = 5; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	else {
		printf("Error: Format only works for graphs up to 68719476735 vertices.\n");
		abort();
	}
}

void loadGraph(string graphString, int numberOfVertices) {
    vector<int> emp;
    graph.assign(n,emp);
	int startIndex = 0;
	if (graphString[startIndex] == '>') { // Skip >>graph6<< header.
		startIndex += 10;
	}
	if (numberOfVertices <= 62) {
		startIndex += 1;
	}
	else if (numberOfVertices <= MAXVERTICES) {
		startIndex += 4;
	}
	else {
		printf("Error: Program can only handle graphs with %d vertices or fewer.\n",MAXVERTICES);
		abort();
	}

	int currentVertex = 1;
	int sum = 0; 
	for (int index = startIndex; index<graphString.size(); index++) {
		int i;
		for (i = prev(graphString[index] - 63, 6); i != -1; i = prev(graphString[index] - 63, i)) {
			while(5-i+(index-startIndex)*6 - sum >= 0) {
				sum += currentVertex;
				currentVertex++;
			}
			sum -= --currentVertex;
			int neighbour = 5-i+(index - startIndex)*6 - sum;
            graph[currentVertex].push_back(neighbour);
            graph[neighbour].push_back(currentVertex);
		}
	}
}

void printGraph()
{
    cerr << "Start printing graph" << endl;
    for(int u=0; u<n; u++)
    {
        cerr << "neighbors of " << u << ":" << endl;
        for(int v: graph[u])
        {
            cerr << v << " ";
        }
        cerr << endl;
    }
    cerr << "End printing graph" << endl;
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    long long nb_graphs_read_from_input=0;
    string line;
    while(getline(cin,line))
    {
        //line+="\n";
        nb_graphs_read_from_input++;
        n = getNumberOfVertices(line);
        loadGraph(line,n);
        int girth=1e9;
        for(int u=0; u<n; u++)
        {
            for(int v : graph[u])
            {
                int newLen=-1;
                int newNumPaths=-1;
                calcNumberShortestPaths(u, v, &newLen, &newNumPaths);
                girth=min(girth,newLen+1);
            }
        }
        int lambda=-2;
        for(int u=0; u<n && lambda!=-1; u++)
        {
            long long sum=0;
            vector<int> neighsOfU=graph[u]; // necessary to make explicit copy, because graph[u] will be altered in "calcNumberShortestPaths"
            for(int v : neighsOfU)
            {
                int newLen=-1;
                int newNumPaths=-1;
                calcNumberShortestPaths(u, v, &newLen, &newNumPaths);
                if(newLen+1==girth)
                {
                    sum+=newNumPaths;
                }
            }
            if(sum%2!=0)
            {
                cout << "LOGICAL ERROR! The sum in the program calcLambdaRegular should have been even!" << endl;
                exit(0);
            }
            sum/=2;
            if(lambda==-2)
            {
                lambda=sum;
            }
            else if(sum != lambda)
            {
                lambda=-1;
            }
        }
        printf("%s\n",line.c_str());
        printf("%d %d %d %d\n",n,(int)graph[0].size(),girth,lambda);
    }
    return 0;
}
