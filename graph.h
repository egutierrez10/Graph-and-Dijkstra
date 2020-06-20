/*graph.h*/

//
// << Everardo Gutierrez >>
//
// Basic graph class using adjacency list representation.

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <set>

using namespace std;


template<typename VertexT, typename WeightT>
class graph
{
  
  vector<VertexT>  Vertices; // vector of vertices in the graph
  map<VertexT,map<VertexT, WeightT>> adjList; // adjacency list for the graph
  int numVertices; // store number of vertices in the graph
  int numEdges; // store number of edges in the graph

public:
  
  //
  // Default Constructor
  //
  graph(){
	  numVertices = 0;
	  numEdges = 0;
  }
	
  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const
  {
    return numVertices;
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const
  {
    return numEdges;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph and returns true. If the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v)
  {

    //
    // Return false if vertex exist already
    //
    auto itr = adjList.find(v);
	if(itr != adjList.end())
		return false;

    //
    // Vertex does not exist in the graph
    //
	map<VertexT, WeightT> data;
	adjList.insert(make_pair(v, data));
	this->Vertices.push_back(v);
	numVertices++;
	
    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist false is returned.
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  // 
  bool addEdge(VertexT from, VertexT to, WeightT weight)
  {
	auto itr = adjList.find(from);
	auto itr2 = adjList.find(to);
	if(itr == adjList.end() || itr2 == adjList.end())
		return false;
	auto itr3 = itr->second.find(to);
	if(itr3 != itr->second.end()){ //exists, overwrite
 		itr3->second = weight;
		return true;
	}else{
 		adjList[from].insert(make_pair(to, weight));
		numEdges++;
	}
	
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If 
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not 
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const
  {
	auto itr = adjList.find(from);
	if(itr == adjList.end())
		return false;
	auto itr2=(*itr).second.find(to);
	if(itr2 == (*itr).second.end())
		return false;
	weight = (*itr2).second;
	
    return true;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const
  {
    set<VertexT> S;
	
	auto itr = adjList.find(v);
	if(itr == adjList.end())
		return S;
    for( auto itr2=(*itr).second.begin(); itr2 != (*itr).second.end(); itr2++)
		S.insert((*itr2).first);
    
    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const
  {
    return this->Vertices;  // returns a copy:
  }

  //
  // dump
  // 
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G;
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const
  {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << numVertices << endl;
    output << "**Num edges: " << numEdges << endl;

    output << endl;
	WeightT weight = WeightT{};   
    output << "**Vertices:" << endl;
	for ( auto itr = adjList.begin() ; itr != adjList.end(); itr++ ) {
		set<VertexT> neighbors = this->neighbors((*itr).first);
		output << (*itr).first << ": ";
		for (VertexT n : neighbors){ 
		  getWeight((*itr).first, n, weight);
		  output << "(" << (*itr).first << "," << n << "," << weight << ")"; 
		}
		cout << endl;
	}
    output << endl;
    output << "**************************************************" << endl;
  }

};