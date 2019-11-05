/*
John Kolden, SCPD
Standford University CS106B
Filename: trailblazer.cpp
Assignment 7
June 1, 2015
Section Leader: Sarah Spikes

Purpose of Program:
This program tests several different graph traversal algorithms: Depth First Search (DFS), Breadth First Search (BFS),
Djikstra's Algorithm and aStar Algorithm. The program also implements Kruskal's Algorithm for creating a Minimum Spanning Tree.
*/

#include "trailblazer.h"
#include "queue.h"
#include "pqueue.h"
#include "basicgraph.h"
#include "hashset.h"

using namespace std;

bool dfsHelper(BasicGraph&, Vertex*, Vertex*,Vector<Vertex*>&);
void bfsHelper (Vector<Vertex*>&, Vertex*, Vertex*);

Vector<Vertex*> depthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {

    graph.resetData();
    Vector<Vertex*> path;
    path.add(start);
    dfsHelper(graph, start, end, path);
    start->setColor(GREEN);

    return path;
}

//This is a recursive function that implements a depth first search between
//two vertexes on a graph. It uses backtracking to prune unsuccessful searches
//before continuing to find a path.
bool dfsHelper(BasicGraph& graph, Vertex* start, Vertex* end,Vector<Vertex*>& path ) {

    if (graph.compare(start, end) == 0) {//base case
        return true;
    } else {
        if (start->visited) {
            return false;
        }

        start->visited = true;
        start->setColor(GREEN);

        for (Edge* n : start->edges) {
            path.add(n->finish);
            if (n->finish->visited == false) {//if we haven't visited this vertex yet
                if (dfsHelper(graph, n->finish, end, path))
                {   n->finish->setColor(GREEN);
                    return true;
                }
            }
            path.remove(path.size() - 1);//backtrack
            start->setColor(GRAY);
        }
    }
    return false;
}


//function implements a breadth first search between two vertexes on a graph. It enqueues
//the starting vertex's neighbor vertexes and then proceeds to test each successive neighbor
//until a path is found.
Vector<Vertex*> breadthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {

    graph.resetData();
    Queue<Vertex*> queue;
    start->visited = true;
    start->previous = NULL;

    queue.add(start);
    Vector<Vertex*> path;

    while (!queue.isEmpty()) {
        Vertex* v = queue.dequeue();
        v->setColor(GREEN);//mark as visited

        if (v == end) {
            bfsHelper(path, end, start);//success; trace a path from the end back to the beginning
            return path;
        }

        for (Edge* edge : v->edges) {
            if (edge->finish->visited == false) {
                Vertex* n = edge->finish;
                n->visited = true;
                n->previous = v;
                n->setColor(YELLOW);
                queue.enqueue(n);//enqueue neighbors of v
            }
        }
    }
    return path;
}

//this function traces a path from the end to the start based on the value
//of each vertex's previous value
void bfsHelper (Vector<Vertex*>& path, Vertex* v, Vertex* start) {

    if (v == start) {
        path.insert(0, v);
        return;
    } else {
        path.insert(0, v);
        bfsHelper(path, v->previous, start);
    }
}


//this function implements Dijkastra's Algorithm for finding the lowest cost
//path between two vertexes. It is similar to the breadthFirstSearch but every neighbor
//is enqueued with a cost such that only the lowest cost alternatives are selected for exploration.
Vector<Vertex*> dijkstrasAlgorithm(BasicGraph& graph, Vertex* start, Vertex* end) {

    graph.resetData();
    PriorityQueue<Vertex*> pqueue;
    Vector<Vertex*> path;

    Set<Vertex*> vertices = graph.getVertexSet();
    for (Vertex* vertex : vertices) {//mark initial cost of all vertexes as infinity
        vertex->cost = POSITIVE_INFINITY;
    }


    start->visited = true;
    start->previous = NULL;
    start->cost = 0;

    pqueue.enqueue(start, 0); //enqueue the starting point with a cost of zero

    while (!pqueue.isEmpty()) {
        Vertex* v = pqueue.dequeue();
        v->visited = true;
        v->setColor(GREEN);

        if (v == end) {//success, trace the path from end to start
            bfsHelper(path, end, start);
            return path;
        }

        for (Edge* edge : v->edges) {
            if (!edge->finish->visited) {
                Vertex* n = edge->finish;
                double cost = v->cost + edge->cost;

                if (cost < n->cost) {//if the cost to get here is cheaper, update n's data with this new info
                    n->previous = v;
                    n->cost = cost;

                    if (n->getColor() == YELLOW) {//if this is in the queue but not yet officially 'visited'
                        pqueue.changePriority(n, n->cost);
                    } else {
                        n->setColor(YELLOW);
                        pqueue.enqueue(n, n->cost);
                    }
                }
            }
        }
    }

    return path;
}

/*this function implements the aStar Algorithm which is very similar to Dijkstra's, with
  the excepption that the cost to get to a vertext is appended by a heuristic component that
  estimates the remaining cost to the end. This allows the algorithm to exclude expensive searches
  that are not likely to result in successfull lowest cost path
*/
Vector<Vertex*> aStar(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();
    PriorityQueue<Vertex*> pqueue;
    Vector<Vertex*> path;

    Set<Vertex*> vertices = graph.getVertexSet();
        for (Vertex* vertex : vertices) {
            vertex->cost = POSITIVE_INFINITY;
        }


    start->visited = true;
    start->previous = NULL;
    start->cost = 0;

    pqueue.enqueue(start, heuristicFunction(start, end));//enqueue initial cost of 0 + heuristic component

    while (!pqueue.isEmpty()) {
        Vertex* v = pqueue.dequeue();
        v->visited = true;
        v->setColor(GREEN);

        if (v == end) {
            bfsHelper(path, end, start);
            return path;
        }

        for (Edge* edge : v->edges) {
            if (!edge->finish->visited) {
                Vertex* n = edge->finish;
                double cost = v->cost + edge->cost;

                if (cost < n->cost) {
                    n->previous = v;
                    n->cost = cost;

                    if (n->getColor() == YELLOW) {
                     pqueue.changePriority(n, n->cost + heuristicFunction(n, end));
                    } else {
                   n->setColor(YELLOW);
                   pqueue.enqueue(n, n->cost + heuristicFunction(n, end));
                    }
                }
            }
        }
    }

    return path;
}

/* this function implements Kruskal's Algorithm for finding a Minimum Spanning Tree.
 * Initially, each vertex is placed into it's own 'cluster' of vertexes and these clusters
 * are merged as minimum weight paths are found between the vertexes. The data structure chosen
 * for the cluster is a vertex of hashsets.
*/
Set<Edge*> kruskal(BasicGraph& graph) {

    PriorityQueue<Edge*> pq;
    Set<Edge*> mst;
    Vector<HashSet<string> > cluster;

    for (Vertex* v : graph.getVertexSet()) {//add each vertex to it's own hashset
        HashSet<string> newSet;
        newSet.add(v->name);
        cluster.add(newSet);//enqueue each hashset into a vertex
    }

    for (Edge* e : graph.getEdgeSet()) {
        pq.enqueue(e, e->weight); //enqueue all edges into a pqueue by weight

    }

    graph.clearEdges();//remove all edges from graph. We will add the lowest cost edges back

    while (cluster.size() > 1 && !pq.isEmpty()) {

        Edge* edge = pq.dequeue();
        int start = 0;
        int end = 0;

        for (int i = 0; i < cluster.size(); i++) {
            if (cluster[i].contains(edge->start->name)) {
                start = i;
            }

            if (cluster[i].contains(edge->finish->name)) {
                end = i;
            }
        }

        if (start != end) {
            cluster[start] += cluster[end];//merge hashset that contains end into hashset that contains start
            cluster.remove(end);
            mst.add(edge);
        }
    }
    return mst;
}
