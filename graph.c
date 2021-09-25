/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include "graph.h"
#include "utils.h"
#include "pq.h"
#include "list.h"

#define INITIALEDGES 32
#define INF 99
#define NIL -1
#define WEIGHT 1

struct edge;

/* Definition of a graph. */
struct graph {
  int numVertices;
  int numEdges;
  int allocedEdges;
  struct edge **edgeList;
};

/* Definition of an edge. */
struct edge {
  int start;
  int end;
};

struct graph *newGraph(int numVertices){
  struct graph *g = (struct graph *) malloc(sizeof(struct graph));
  assert(g);
  /* Initialise edges. */
  g->numVertices = numVertices;
  g->numEdges = 0;
  g->allocedEdges = 0;
  g->edgeList = NULL;
  return g;
}

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end){
  assert(g);
  struct edge *newEdge = NULL;
  /* Check we have enough space for the new edge. */
  if((g->numEdges + 1) > g->allocedEdges){
    if(g->allocedEdges == 0){
      g->allocedEdges = INITIALEDGES;
    } else {
      (g->allocedEdges) *= 2;
    }
    g->edgeList = (struct edge **) realloc(g->edgeList,
      sizeof(struct edge *) * g->allocedEdges);
    assert(g->edgeList);
  }

  /* Create the edge */
  newEdge = (struct edge *) malloc(sizeof(struct edge));
  assert(newEdge);
  newEdge->start = start;
  newEdge->end = end;

  /* Add the edge to the list of edges. */
  g->edgeList[g->numEdges] = newEdge;
  (g->numEdges)++;
}

/* Frees all memory used by graph. */
void freeGraph(struct graph *g){
  int i;
  for(i = 0; i < g->numEdges; i++){
    free((g->edgeList)[i]);
  }
  if(g->edgeList){
    free(g->edgeList);
  }
  free(g);
}

/* Finds:
  - Number of connected subnetworks (before outage) (Task 2)
  - Number of servers in largest subnetwork (before outage) (Task 3)
  - SIDs of servers in largest subnetwork (before outage) (Task 3)
  - Diameter of largest subnetworks (after outage) (Task 4)
  - Number of servers in path with largest diameter - should be one more than
    Diameter if a path exists (after outage) (Task 4)
  - SIDs in path with largest diameter (after outage) (Task 4)
  - Number of critical servers (before outage) (Task 7)
  - SIDs of critical servers (before outage) (Task 7)
*/
struct solution *graphSolve(struct graph *g, enum problemPart part,
                        int numServers, int numOutages, int *outages) {
  struct solution *solution = (struct solution *)
                                malloc(sizeof(struct solution));
  assert(solution);
  /* Initialise solution values */
  initaliseSolution(solution);

  int i;

  /********* BEFORE OUTAGE *************/ 
  // 2-d array containing all vertices, row: subnet number, column: vertices
  int** subVertices_b = (int**)malloc(sizeof(int*) * numServers);
  assert(subVertices_b);

  // allocate menmory to each row of the array 
  for (i = 0; i < numServers; i++) {
      subVertices_b[i] = (int*)malloc(sizeof(int) * numServers);
      assert(subVertices_b[i]);
  }

  // number of vertices in each row of subVertices
  int* numSubVert_b = (int*)calloc(numServers, sizeof(int));
  assert(numSubVert_b);

  // the number of subnets
  int numSub_b = getSubnets(g, subVertices_b, numSubVert_b, numServers, 0, NULL);


  /************ AFTER OUTAGE ***************/
  int** subVertices_o = (int**)malloc(sizeof(int*) * numServers);
  assert(subVertices_o);
  for (i = 0; i < numServers; i++) {
      subVertices_o[i] = (int*)malloc(sizeof(int) * numServers);
      assert(subVertices_o[i]);
  }
  int* numSubVert_o = (int*)calloc(numServers, sizeof(int));
  assert(numSubVert_o);
  int numSub_o = getSubnets(g, subVertices_o, numSubVert_o, numServers, numOutages, outages);



  if(part == TASK_2) {
    /* IMPLEMENT TASK 2 SOLUTION HERE */
    solution->connectedSubnets = count_connected_subnets(g, numServers);
  } else if(part == TASK_3) {
    /* IMPLEMENT TASK 3 SOLUTION HERE */    
    int* SubnetSIDs = (int*) malloc(numServers * sizeof(int));
    assert(SubnetSIDs);
    solution->largestSubnet = count_largest_subnet_servers(g, numServers, subVertices_b, 
                                                    numSubVert_b, numSub_b, SubnetSIDs);
    solution->largestSubnetSIDs = SubnetSIDs;    
  } else if(part == TASK_4) {
    /* IMPLEMENT TASK 4 SOLUTION HERE */
    int* diamSIDs = (int*) malloc(numServers * sizeof(int));
    int diam = largest_subnet_diameter(g, numServers,  subVertices_o, 
                    numSubVert_o, numSub_o, numOutages, outages, diamSIDs);
    solution->postOutageDiameter = diam;
    solution->postOutageDiameterCount = diam + 1;
    solution->postOutageDiameterSIDs = diamSIDs;
  } else if(part == TASK_7) {
    /* IMPLEMENT TASK 7 SOLUTION HERE */
    int* criticals = malloc(sizeof(int) * numServers);
    assert(criticals);
    solution->criticalServerCount = findCritical(g, numServers, subVertices_b, 
                                            numSubVert_b, numSub_b, criticals);
    solution->criticalServerSIDs = criticals;
  }


  for (i = 0; i < numServers; i++) {
      free(subVertices_b[i]);
      free(subVertices_o[i]);
  }
  free(subVertices_b);
  free(numSubVert_b);
  free(subVertices_o);
  free(numSubVert_o);


  return solution;
}


/* Task2: find the number of connected subnetworks (before outage) */
int count_connected_subnets(struct graph* g, int numServers) {

    // mark all vertices with false
    bool* explored = (bool*)calloc(numServers, sizeof(bool));
    assert(explored);

    int count = 0;
    int i;
    for (i = 0; i < numServers; i++) {
        if (!explored[i]) {
            dfsExplore(g, explored, i, &count, numServers, NULL, NULL, 0, NULL);
        }
    }
    free(explored);

    return count;
}

/* DFS traversal */
void dfsExplore(struct graph* g, bool* explored, int vertex, int* numSub, int numServers,
                 int** subVertices, int* numSubVert, int numOutages, int* outages) {
    int i;
    // terminate the function if the server is affected by teh outage
    if (outages != NULL) {
        if(inArray(vertex, outages, numOutages)) return;
    }
    
    explored[vertex] = true;
    
    // retrive neighbours
    int* neighbours = (int*) malloc(sizeof(int) * numServers);
    assert(neighbours);
    int numNeighbours = get_neighbours(g, vertex, neighbours, numOutages, outages);
 
    // the first vertex of a subnet can be located if there is no neighbour explored
    if (has_explored_neighb(explored, neighbours, numNeighbours)) {
        *numSub += 1;
    }

    // store all the vertices into subVertices array if needed
    if (subVertices != NULL && numSubVert != NULL) {
        int numSubIndex = *numSub - 1;
        subVertices[numSubIndex][numSubVert[numSubIndex]] = vertex;
        numSubVert[numSubIndex] += 1;
    }
      
    // DFS recurrence
    for (i = 0; i < numNeighbours; i++) {
        int neighbour = neighbours[i];
        if (!explored[neighbour]) {
            dfsExplore(g, explored, neighbour, numSub, numServers, 
                        subVertices, numSubVert, numOutages, outages);
        }
    }
    
}


/* find neighbours of a specific vertice */
int get_neighbours(struct graph* g, int vertice, int* neighbours, 
                                    int numOutages, int* outages) {
    assert(g != NULL);

    int numNeighbours = 0;
    int i;
    for (i = 0; i < g->numEdges; i++) {
        struct edge *currEdge = (g->edgeList)[i];
        int start = currEdge->start;
        int end = currEdge->end;
       
        // case with outage
        if (outages != NULL) {
            // ignore the servers with outage 
            if (inArray(start, outages, numOutages) || 
                inArray(end, outages, numOutages)) continue;                                                
        }
        
        if (start == vertice) {
            neighbours[numNeighbours++] = end;
        }
        else if (end == vertice) {
            neighbours[numNeighbours++] = start;
        }
    }  
    
    return numNeighbours;
}

/* chech if an element is in an array */
bool inArray(int val, int* values, int range) {
    int i;
    for (i = 0; i < range; i++) {
        if (values[i] == val) {
            return true;
        }
    }
    return false;
}


/* check if a vertex has neighbours to explore */
bool has_explored_neighb(bool* explored, int* neighbours, int numNeighbours) {
    bool has = true;
    int i;
    for (i = 0; i < numNeighbours; i++) {
        if (explored[neighbours[i]]) {
            has = false;
            break;
        }
    }

    return has;
}



/* Task3: Number of servers in largest subnetwork (before outage) */
int count_largest_subnet_servers(struct graph* g, int numServers, 
            int** subVertices, int* numSubVert, int numSub, int* SubnetSIDs) {   
    int i;
    // initilaize the index of the largest subnet in subVertices array
    int subnetIndex = 0;
    // the number of vertices in hte largest subnet
    int largest = getMax(numSubVert, numSub, &subnetIndex);

    // an array to hold all indices of largest subnet in subVertices array
    int* largestSubnets = malloc(sizeof(int) * numSub);
    assert(largestSubnets);

    // the number of largest subnet
    int numLargest = countLargest(largest, numSubVert, numSub, largestSubnets);

    // sort the first largest subnet
    insertionSort(subVertices[largestSubnets[0]], numSubVert[largestSubnets[0]]);

    // in case the largest subnet is not unique
    if (numLargest > 1) {
        // Compare the first vertice (smallest ID) of each largest subnet after sorting
        // Find the index of the subnet with smallest server ID in subVertices array
        int minVert = subVertices[largestSubnets[0]][0];
        int minIndex = largestSubnets[0];

        for (i = 1; i < numLargest; i++) {
            int currIndex = largestSubnets[i];
            insertionSort(subVertices[currIndex], numSubVert[currIndex]);

            if (subVertices[currIndex][0] < minVert) {
                minIndex = currIndex;
                minVert = subVertices[currIndex][0];
            }
        }
        subnetIndex = minIndex;
    }
    
    // copy the largest subnet from subVertices to SubnetSIDs
    for (i = 0; i < largest; i++) {
        SubnetSIDs[i] = subVertices[subnetIndex][i];
    }
    
    free(largestSubnets);    

    return largest;
}


/* store all subnetworks into a 2-d array, return the number of subnets */
int getSubnets(struct graph* g, int** subVertices, int* numSubVert, 
                        int numServers, int numOutages, int* outages) {
    
    // mark all vertices with false
    bool* explored = (bool*)calloc(numServers, sizeof(bool));
    assert(explored);

    // number of subnetworks
    int numSub = 0;

    int i;
    // traverse the graph by using DFS
    for (i = 0; i < numServers; i++) {
        if (!explored[i]) {
            dfsExplore(g, explored, i, &numSub, numServers, 
                subVertices, numSubVert, numOutages, outages);
        }
    }
    free(explored);

    return numSub;
}





/* find the largest value in an array */
int getMax(int* values, int range, int* index) {
    assert(values != NULL);
    int max = values[0];
    int i;
    for (i=0; i < range; i++) {
        if (values[i] > max) {
            *index = i;
            max = values[i];
        }
    }

    return max;
}

/* count the number of the largest value in an array */
int countLargest(int max, int* values, int range, int* largestSubnets) {
    assert(values != NULL);
    int count = 0;
    int i;
    for (i = 0; i < range; i++) {
        if (values[i] == max) {
            largestSubnets[count] = i;
            count += 1;
        }
    }
    return count;
}

/* insortion sort algorithm */
void insertionSort(int* values, int range) {
    int i, j;
    for (i = 1; i < range; i++) {
        for (j = i - 1; j >= 0; j--) {
            if (values[j + 1] < values[j]) {
                // swap
                int temp = values[j + 1];
                values[j + 1] = values[j];
                values[j] = temp;
            }
        }
    }
}



/* Task4: find diameter of largest subnetworks (after outage) */
int largest_subnet_diameter(struct graph* g, int numServers, int** subVertices, 
        int* numSubVert,  int numSub, int numOutages, int* outages, int* diamSIDs) {

    int i, v, j; // index
    int longestDiameter = 0;
    int longestDiameter_i;

    // by performing Dijkstra, use indices of serverIDs in each subnet instead
    // store the final indices into path array
    int* path = malloc(numServers * sizeof(int));
  
    // apply Dijkstra algorithm to each server
    for (i = 0; i < numSub; i++) {
        for (v = 0; v < numSubVert[i]; v++) {
            
            int currLongest;
            int currLongest_i;
            
            int* dist = malloc(numSubVert[i] * sizeof(int));
            assert(dist);
            int* prev = malloc((numSubVert[i]+1) * sizeof(int));
            assert(prev);
     
            dijkstra(g, subVertices[i], numSubVert[i], subVertices[i][v], 
                                         dist, prev, numOutages, outages);
            
            // get the longest distance 
            currLongest = getMax(dist, numSubVert[i], &currLongest_i);

            if (currLongest > longestDiameter) {
                longestDiameter = currLongest;
                longestDiameter_i = i;

                // retrive the path by looking backwards in the prev array
                path[longestDiameter] = prev[numSubVert[i]];
                for (j = longestDiameter-1; j >=0; j--) {
                    path[j] = prev[path[j+1]];
                }
            }
            free(dist);
            free(prev);
        }
    }

    // change the indices back to serverIDs
    for (j = 0; j <= longestDiameter; j++) {
        int index = path[j];       
        diamSIDs[j] = subVertices[longestDiameter_i][index];
    }

    free(path);
    
    return longestDiameter;
}


/* given an element of an array, find its index in that array */
int getArrayIndex(int var, int* values, int range) {
    int i;
    for (i = 0; i < range; i++) {
        if (values[i] == var) {
            return i;
        }
    }
    return 0;
}



/* Dijkstra Algorithm, find the shortest path from a specific server to all other servers */
void dijkstra(struct graph* g, int* subnet, int numVert,
            int start, int* dist, int* prev, int numOutages, int* outages) {  

    // note that the dist and prev arrays store the 
    // indices of the serverIDs in subnet array
    // but enqueue the serverIDs

    int i;
    // initialize
    for (i = 0; i < numVert; i++) {
        dist[i] = INF;
        prev[i] = NIL;
    }
    *(prev+1) = NIL;

    int start_i = getArrayIndex(start, subnet, numVert);
    dist[start_i] = 0;

    struct pq* pq = newPQ();

    
    for (i = 0; i < numVert; i++) {       
        enqueue(pq, &subnet[i], dist[i]);
    }
    

    while (!empty(pq)) {
        int* u = deletemin(pq);
        int u_i = getArrayIndex(*u, subnet, numVert);       

        // retrive neighbours
        int* neighbours = (int*)malloc(sizeof(int) * numVert);
        assert(neighbours);
        int numNeighbours = get_neighbours(g, *u, neighbours, numOutages, outages);


        for (i = 0; i < numNeighbours; i++) {
            int neighb = neighbours[i];
            int neighb_i = getArrayIndex(neighb, subnet, numVert);

            if (inPQ(pq, &subnet[neighb_i]) && dist[u_i] + WEIGHT < dist[neighb_i]) {
                dist[neighb_i] = dist[u_i] + WEIGHT;
                prev[neighb_i] = u_i;
                updatePQ(pq, &subnet[neighb_i], dist[neighb_i]);
            }
        }

        free(neighbours);

        if (empty(pq)) {
            prev[numVert] = u_i;
        }
    }

   
    freePQ(pq);

}



/* Task7: Find SIDs of critical servers (before outage) */
int findCritical(struct graph* g, int numServers, int** subVertices, 
                        int* numSubVert, int numSub, int* criticals) {
    int i;

    int** children = malloc(sizeof(int*) * numServers);
    assert(children);
    
    // allocate menmory to each row of the array 
    for (i = 0; i < numServers; i++) {        
        children[i] = malloc(sizeof(int) * numServers);
        assert(children[i]);
    }

    int* numChildren = (int*)calloc(numServers, sizeof(int));
    assert(numChildren);

    // mark all vertices with false
    bool* explored = (bool*)calloc(numServers, sizeof(bool));
    assert(explored);
    
    struct list* plist = NULL;  
    int order = 0;
    int num = 0;
    int* porder = malloc(sizeof(int) * numServers);
    assert(porder);
    int* hra = malloc(sizeof(int) * numServers);
    assert(hra);

    for (i = 0; i < numServers; i++) {
        hra[i] = NIL;
    }

    int s;
    // construct DFS tree in each connected component
    for (s = 0; s < numSub; s++) {
        insertionSort(subVertices[s], numSubVert[s]);
        for (i = 0; i < numSubVert[s]; i++) {
            int v = subVertices[s][i];
            if (!explored[v]) {
                dfsFindCriticals(g, numServers, subVertices[s], numSubVert[s], &subVertices[s][i], 
                    &order, porder, hra, explored, children, numChildren, plist);
            }
        }
    }
    
    

    int j;
    // go through each vertex again, find their children and 
    // get critical servers
    for (s = 0; s < numSub; s++) {
        // root of one DFS tree
        int root = subVertices[s][0];
        
        if (numChildren[root] > 1) {
            criticals[num++] = root;
        }

        // go through the rest except for the leaf
        for (i = 1; i < numSubVert[s] - 1; i++) {
            int v = subVertices[s][i];
           
            for (j = 0; j < numChildren[v]; j++) {
                int child = children[v][j];
                if (porder[hra[child]] >= porder[v]) {
                    criticals[num++] = v;
                }
            }
        }
    }

    
    free(explored);
    free(porder);
    free(hra);
    free(children);
    for (i = 0; i < numServers; i++) {
        free(children[i]);
    }
    free(numChildren);
    freeList(plist);


    return num;
}


/* DFS to find critical points */
void dfsFindCriticals(struct graph* g, int numServers, int* subVertices, int numSubVert, int* v, int* order,
        int* porder, int* hra, bool* explored, int** children, int* numChildren, struct list* plist) {

    explored[*v] = true;
    porder[*v] = *order + 1;
    *order += 1;
    
    if (hra[*v] == NIL) {
        hra[*v] = *v;
    }
   

    // retrive neighbours
    int* neighbours = (int*)malloc(sizeof(int) * numServers);
    assert(neighbours);
    int numNeighbours = get_neighbours(g, *v, neighbours, 0, NULL);
    // make sure to explore the neighbour with lower ID first
    insertionSort(neighbours, numNeighbours);
    
    int i;
    // record v’s HRA as the smallest ID of any of its reachable ancestors 
    for (i = 0; i < numNeighbours; i++) {

        int neighb = neighbours[i];
        // list stores the pointer of items in subVertices
        int neighb_i = getArrayIndex(neighb, subVertices, numSubVert);   
        if (explored[neighb]) {            
            if (inList(plist, &subVertices[neighb_i]) && neighb != *peekHead(plist) &&
                                                         neighb < hra[*v]) {                
                hra[*v] = neighb;
                                 
            }
        }
    }
    

    if (plist == NULL) {
        plist = newlist(v);
    }
    else {
        plist = prependList(plist, v);
    }
      

    for (i = 0; i < numNeighbours; i++) {
        int neighb = neighbours[i];
        int neighb_i = getArrayIndex(neighb, subVertices, numSubVert);
        if (!explored[neighb]) {
            children[*v][numChildren[*v]++] = neighb;
            dfsFindCriticals(g, numServers, subVertices, numSubVert, &subVertices[neighb_i], 
                            order, porder, hra, explored, children, numChildren, plist);
        }
      
        // pass HRA of the child of v, if it is smaller, to HRA of v before poping
        int* head = peekHead(plist);
        if (*head != *v) {           
            if (hra[*head] < hra[*v]) {              
                hra[*v] = hra[*head];               
            }
            plist = deleteHead(plist);
            
        }              
    }
  
    free(neighbours);

}




