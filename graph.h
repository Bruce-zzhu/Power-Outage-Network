/*
graph.h

Visible structs and functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021 and
  modified for Assignment 2 2021
*/

#include <stdbool.h>

struct list;

/* Definition of a graph. */
struct graph;

enum problemPart;

struct solution;

/* A particular solution to a graph problem. */
#ifndef SOLUTION_STRUCT
#define SOLUTION_STRUCT
struct solution {
  int connectedSubnets;
  int largestSubnet;
  int *largestSubnetSIDs;
  int postOutageDiameter;
  int postOutageDiameterCount;
  int *postOutageDiameterSIDs;
  int criticalServerCount;
  int *criticalServerSIDs;
};
#endif

/* Which part the program should find a solution for. */
#ifndef PART_ENUM
#define PART_ENUM
enum problemPart {
  TASK_2=0,
  TASK_3=1,
  TASK_4=2,
  TASK_7=3
};
#endif

/* Creates an undirected graph with the given numVertices and no edges and
returns a pointer to it. NumEdges is the number of expected edges. */
struct graph *newGraph(int numVertices);

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end);

/* Finds:
  - Number of connected subnetworks (before outage) (Task 2)
  - Number of servers in largest subnetwork (before outage) (Task 3)
  - SIDs of servers in largest subnetwork (before outage) (Task 3)
  - Diameter of largest subnetworks (after outage) (Task 4)
  - Number of servers in path with largest diameter - should be one more than
    Diameter if a path exists (after outage) (Task 4)
  - SIDs in largest subnetwork (after outage) (Task 4)
  - Number of critical servers (before outage) (Task 7)
  - SIDs of critical servers (before outage) (Task 7)
 */
struct solution *graphSolve(struct graph *g, enum problemPart part,
  int numServers, int numOutages, int *outages);

/* Frees all memory used by graph. */
void freeGraph(struct graph *g);

/* Sets all values to initial values so free can work for all tasks without
  change. */
void initaliseSolution(struct solution *solution);

/* Frees all data used by solution. */
void freeSolution(struct solution *solution);


/* Task2: find the number of connected subnetworks (before outage) */
int count_connected_subnets(struct graph* g, int numServers);

/* find neighbours of a specific vertice */
int get_neighbours(struct graph* g, int vertice, int* neighbours,
                                    int numOutages, int* outages);

/* chech if an element is in an array */
bool inArray(int val, int* values, int range);

/* DFS traversal */
void dfsExplore(struct graph* g, bool *explored, int vertice, int* numSub, 
     int numServers, int** subVertices, int* numSubVert, int numOutages, int* outages);

/* check if a vertex has neighbours to explore */
bool has_explored_neighb(bool* explored, int* neighbours, int numNeighbours);

/* Task3: Number of servers in largest subnetwork (before outage) */
int count_largest_subnet_servers(struct graph* g, int numServers, int** subVertices, 
                                int* numSubVert, int numSub, int* SubnetSIDs);
/* store all subnetworks into a 2-d array, return the number of subnets */
int getSubnets(struct graph* g, int** subVertices, int* numSubVert, 
                        int numServers, int numOutages, int* outages);

/* find the largest value in an array */
int getMax(int* values, int range, int *index);

/* count the number of the largest value in an array */
int countLargest(int max, int* values, int range, int* largestSubnets);

/* insortion sort algorithm */
void insertionSort(int* values, int range);

/* Task4: find diameter of largest subnetworks (after outage) */
int largest_subnet_diameter(struct graph* g, int numServers, int** subVertices, int* numSubVert, int numSub,
              int numOutages, int* outages, int* diamSIDs);

/* Dijkstra Algorithm, find the shortest path from a specific server to all other servers */
void dijkstra(struct graph* g, int* subnet, int numVert, 
                int vertice, int* dist, int* prev, int numOutages, int *outages);

/* given an element of an array, find its index in that array */
int getArrayIndex(int var, int* values, int range);

/* Task7: Find SIDs of critical servers (before outage) */
int findCritical(struct graph* g, int numServers, int** subVertices, int* numSubVert, int numSub, int* criticals);

/* DFS to find critical points */
void dfsFindCriticals(struct graph* g, int numServers, int* subVertices, int numSubVert, int* v, int* order, 
                int* porder, int* hra, bool* explored, int** children, int* numChildren, struct list* plist);



