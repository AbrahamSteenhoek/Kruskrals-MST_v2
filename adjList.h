#ifndef ADJ_LIST
#define ADJ_LIST

#include <vector>
#include <iostream>
#include <ostream>
#include <queue>
#include <list>
#include <exception>
#include "disjoint.h"

// cd /mnt/c/Users/AJ\ Steenhoek/Documents/Fall_2017/Data_Structures/Prog4_v2

using namespace std;

template <typename Weight, bool Directed=false>
class AdjList {
public:
    typedef int Vertex;
    struct Edge3 { Vertex u; Vertex v; Weight w; };
private:
    std::vector<std::vector<Edge3> > adj;
    unsigned int nEdges;
    DisjointSet dSet;
    list<Edge3> mst;
public:
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // min priority queue that stores Edge3's
    // not a direct heap tree structure, just a linked list
    class Edge3PQ {
        vector<Edge3> edgeList;
    public:
        unsigned int parent(int i) { return (i - 1)/2; }
        unsigned int leftChild(int i) { return 2 * i + 1; }
        unsigned int rightChild(int i) { return 2 * i + 2; }
        unsigned int minChild(int i) { return (edgeList[leftChild(i)].w < edgeList[rightChild(i)].w) ? leftChild(i) : rightChild(i); }
        bool empty() { return edgeList.size() == 0; }

        Edge3 top() { return edgeList[0]; }

        void push(Edge3 edge) {
            edgeList.push_back(edge);
            int i = edgeList.size() - 1;
            while( i > 0 && edgeList[i].w < edgeList[parent(i)].w) {
                swap(edgeList[i], edgeList[parent(i)]);
                i = (i - 1)/2;
            }
        }

        void pop() {
            //cout << "pop" << endl;
            if (edgeList.size() > 0) {
                edgeList[0] = edgeList.back();
                int i = 0;
                edgeList.pop_back();
                while (i < edgeList.size()) {
                    int min = minChild(i);
                    if(rightChild(i) >= edgeList.size()) {
                        if(leftChild(i) >= edgeList.size()) {
                            break;
                        }
                        else {
                            min = leftChild(i);
                        }
                    }
                    
                    //cout << "i = " << i << endl;
                    if (edgeList.at(i).w > edgeList.at(min).w ) {
                        //cout << "swap: " << i << ", " << minChild(i) << endl;
                        swap(edgeList.at(i), edgeList.at(min));
                        
                    }
                    i++;
                }   
                
            }
        }

        void print() {
            for (Edge3 edge : edgeList) {
                cout << "{" << edge.u << ", " << edge.v << ", " << edge.w << "}";
            }
            cout << endl;
        }

    }; // end of class Edge3PQ
    //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    unsigned int EdgeCount() const { return nEdges; }
    unsigned int VertexCount() const { return adj.size(); }
    
    void addVertex(Vertex u) {
        if (VertexCount() <= u) {
            adj.resize(u + 1);
        }     
    }

    void addEdge(Vertex u, Vertex v, Weight w) {
        addVertex(std::max(u,v));
        adj[u].push_back({u,v,w});
        if (!Directed) {
            adj[v].push_back({u,v,w});
        }
        nEdges++;

    }

    AdjList() : nEdges(0) { ; }
    AdjList(std::istream &is) {
        Vertex u,v; Weight w;
        while (is >> u >> v >> w) {
            addEdge(u,v,w);
        }
        for(int i = 0; i < adj.size(); i++) {
            dSet.addSingleton(i);
        }
        //dSet.print();
    }

    void print(std::ostream &os) {
        std::cout << "Graph:" << std::endl;
        for(Vertex u = 0; u < VertexCount(); u++) {
            os << u << ":";
            for(auto e : adj[u]) {
                os << " {" << e.u << ", " << e.v << ", " << e.w << "}";
            }
            os << std::endl;
        }
        //std::cout << "\nParents:" << std::endl;
        
    }

    void printMST() {
        Weight mstWeight = calculateMST();
        for (Edge3 edge : mst) {
            cout << edge.u << " " << edge.v << endl;
        }
        cout << "The total weight is " << mstWeight << endl;
    }

    bool isCycle(Vertex x, Vertex y) {
        return dSet.find(x) == dSet.find(y);
    }

    struct EdgeComp {
        bool operator() (Edge3 edge1, Edge3 edge2) {
            return edge1.w > edge2.w;
        }
    };
    //Kruskrals
    Weight calculateMST() {
        
        Edge3PQ minpq;

        for(auto row : adj) {
            for(Edge3 edge: row) {
                minpq.push(edge);
            }
        }
        //minpq.print();

        // std::cout << std::endl;
        // minpq.pop();
        // minpq.print();
        // std::cout << std::endl;
        // minpq.pop();
        // minpq.print();
        // std::cout << std::endl;
        // minpq.pop();
        // minpq.print();
        // std::cout << std::endl;
        // minpq.pop();
        // minpq.print();
        // std::cout << std::endl;
        // minpq.pop();
        // minpq.print();


        mst.clear();
        Weight totalWeight = 0;
        
        while (!minpq.empty()) {
            Edge3 edge = minpq.top();
            minpq.pop();
            // if next edge creates a cycle
            if (isCycle(edge.u, edge.v))
                continue;
            // else
            dSet.merge(edge.u, edge.v);
            mst.push_back(edge);
            totalWeight += edge.w;
        }
        return totalWeight;
    }

}; // end of class AdjList

#endif