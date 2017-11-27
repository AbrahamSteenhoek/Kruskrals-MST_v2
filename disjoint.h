#include <iostream>
#include <vector>

class DisjointSet {
	unsigned int nSets;
	std::vector<int> parent;
public:
	DisjointSet() : nSets(0) { ; }
	unsigned int NumberOfSets() const { return nSets; }
	
	void addSingleton(int v) {
		if (parent.size() > v) { return; }
		unsigned int oldSize = parent.size();
		parent.resize(v+1, -1);
		nSets += v - oldSize + 1;
	}
	
	int find(int value) {
		// am I at a root or the child of the root
		if (parent[value] < 0) { return value; }
		if (parent[parent[value]] < 0) { return parent[value]; }

		int prnt = find(parent[value]);
		parent[value] = prnt;
		return prnt;
	}
	
	// Merge the set containing u with the set containing v
	void merge(int u, int v) {
		if ((u = find(u)) == (v = find(v))) { return; }

		--nSets;
		if (parent[v] < parent[u]) { std::swap(u,v); }
		parent[u] += parent[v];
		parent[v] = u;
	}

    int setSize(int v) { return -parent[find(v)]; }
    
    void print() {
        for (int i=0; i<parent.size();i++) {
            std::cout << i << ": " << parent[i] << std::endl;
        }
    }
}; 