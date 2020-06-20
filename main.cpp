/*main.cpp*/

// Name: Everardo Gutierrez
// U. of Illinois, Chicago
// CS 251: Spring 2020
// Project #07: open street maps, graphs, and Dijkstra's alg
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:  
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <stack>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h" 

using namespace std;
using namespace tinyxml2;

// Function object class to order our priority queue properly with the smaller 
// values at the front.
class prioritize{
public: 
bool operator()(const pair<long long, double>& p1, const pair<long long, double>& p2) const
{
   if(p1.second == p2.second)
      return p1.first > p2.first;
	return p1.second > p2.second;
}
};

//
// Dijkstra:
//
// Performs Dijkstra's shortest weighted path algorithm from
// the given start vertex.  Returns a vector of vertices in
// the order they were visited, along with a map of (long long,double)
// pairs where the long long is a vertex V and the double is the 
// distance from the start vertex to V; if no such path exists,
// the distance is INF.
//
vector<long long> Dijkstra(graph<long long,double>& G, 
  long long startV, 
  map<long long, double>& distances,
  map<long long, long long>& predecessors)
{
  const double INF = numeric_limits<double>::max();
  vector<long long>  visited;
  priority_queue<
  pair<long long,double>,
  vector<pair<long long,double>>,
  prioritize> unvisitedQueue;
  
  vector<long long> vertices = G.getVertices();
  for(auto& currentV: vertices){
     distances[currentV] = INF;
     unvisitedQueue.push(make_pair(currentV,INF));
	 predecessors[currentV] = 0.0;
  }
	
  distances[startV] = 0.0;
  unvisitedQueue.push(make_pair(startV,0.0));
  set<long long> neighbors;
  double edgeWeight = 0.0; 
  double alternativePathDistance = 0.0;
  while(!unvisitedQueue.empty()){
     pair<long long, double>curr = unvisitedQueue.top();
     long long currentV = curr.first;
     unvisitedQueue.pop();
     auto itr = find(visited.begin(),visited.end(),currentV);
     if(distances[currentV] == INF){
        break;
     }else if(itr != visited.end()){
       continue;
    }else{
       visited.push_back(currentV);
    }
    neighbors = G.neighbors(currentV);
    for(auto& adjV: neighbors){
       G.getWeight(currentV, adjV,edgeWeight);
       alternativePathDistance = distances[currentV] + edgeWeight;
       if(alternativePathDistance < distances[adjV]){
          distances[adjV] = alternativePathDistance;
		  predecessors[adjV] = currentV;
          unvisitedQueue.push(make_pair(adjV,alternativePathDistance));
       }
    }
  }
  return visited;
}

//
// Found:
//
// Searches in the buildings vector for the searchBuilding string parameter either by abbreviation or fullname.
// If the building exist we save the buildings object to searchCords and return true. If the building is not 
// found as either abbreviation or fullname we return false.
//
bool found(vector<BuildingInfo>& Buildings, BuildingInfo& searchCords, string searchBuilding, string position){
	auto itr = find_if( Buildings.begin(), Buildings.end(),[searchBuilding] ( const BuildingInfo& a ) { return a.Abbrev.find(searchBuilding) != string::npos;}); 
	auto itr2 = find_if( Buildings.begin(), Buildings.end(),[searchBuilding] ( const BuildingInfo& a ) { return a.Fullname.find(searchBuilding) != string::npos;}); 
	
	if(itr != Buildings.end()){ // searched by abbreviation
		searchCords = (*itr);
		return true;
	}else if(itr2 != Buildings.end()){ // searched by fullname
		searchCords = (*itr2);
		return true;
	}else{ // building passed does not exist
		return false;
	}
}

//
// searchFootways:
//
// Searches in the footways map for the smallest starting and ending point. This is done 
// by looping through and calling distBetween2Points to calculate the distance between the 
// starting or destination building and the current footway node. 
//
void searchFootways(map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways, 
					Coordinates& c1, Coordinates& c2, long long& fStart, long long& fEnd){
	double distance = 0; 
	double distance2 = 0;
	double minStart = 1.0, minEnd = 1.0;
	vector<long long> nodes; // vector to hold footway nodes to loop through
	long long f1;
	
	// search for nearest start and end point to the starting and destination building
	for(auto f: Footways){
		nodes = f.Nodes;
		f1 = nodes[0]; // set the first footway node as the minimum
		distance = distBetween2Points(c1.Lat, c1.Lon, Nodes[f1].Lat, Nodes[f1].Lon); // calculate distance between start building and footway node
		distance2 = distBetween2Points(Nodes[f1].Lat, Nodes[f1].Lon, c2.Lat, c2.Lon); // calculate distance between footway node and destination building
		if(distance < minStart){ // check to see if distance calculated is lower than our set minStart
			minStart = distance;  // it is so update minStart
			fStart = f1;
		}	
		if(distance2 < minEnd){ // check to see if distance calculated is lower than our set minEnd
			minEnd = distance2; // it is so update minEnd
			fEnd = f1;
		}
		// loop through each footway node that exist and calculate its distance and update minStart and minEnd
		for(unsigned int i = 1; i < nodes.size(); i++){ 
			f1 = nodes[i];
			distance = distBetween2Points(c1.Lat, c1.Lon, Nodes[f1].Lat, Nodes[f1].Lon);
			distance2 = distBetween2Points(Nodes[f1].Lat, Nodes[f1].Lon, c2.Lat, c2.Lon);
			if(distance < minStart){
				minStart = distance;
				fStart = f1;
			}
			if(distance2 < minEnd){
				minEnd = distance2;
				fEnd = f1;
			}  
		 }
	}
	cout << endl;
	cout << "Nearest start node:" <<endl;
	cout << "  " << fStart<< endl;
	cout << "  (" << Nodes[fStart].Lat << ", " << Nodes[fStart].Lon << ")" << endl;
	cout << "Nearest destination node:" <<endl;
	cout << "  " << fEnd<< endl;
	cout << "  (" << Nodes[fEnd].Lat << ", " << Nodes[fEnd].Lon << ")" << endl;
}

//
// outputDijkstra:
//
// Function is to essentially loop through the predecessors map until we hit 0 
// or no previous vertex to the current. While doing this we push each vertex onto 
// the stack to then be printed down below. 
//
void outputDijkstra(long long& fEnd, map<long long, long long>& predecessors){
	stack<long long> s;
	
	while(fEnd != 0){ // keep looping until no previous vector exist
		s.push(fEnd); // push the vertex onto the stack
		fEnd = predecessors[fEnd]; // update fEnd with the previous vector
	}
	
  cout << "Path: ";
	while(!s.empty()){ // pop the top off the stack and print the vertex
		if(s.size() == 1)
			cout << s.top();
		else
			cout << s.top() << "->";
		s.pop();
	}
	cout << endl;
}
//////////////////////////////////////////////////////////////////
//
// main
//
int main()
{
  graph<long long, double>     G;         // graph in which vertices are nodes, weights are distances
  map<long long, Coordinates>  Nodes;     // maps a Node ID to it's coordinates (lat, lon)
  vector<FootwayInfo>          Footways;  // info about each footway, in no particular order
  vector<BuildingInfo>         Buildings; // info about each building, in no particular order
  XMLDocument                  xmldoc;
  
  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "")
  {
    filename = def_filename;
  }

  //
  // Load XML-based map file 
  //
  if (!LoadOpenStreetMap(filename, xmldoc))
  {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }
  
  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert((unsigned)nodeCount == Nodes.size());
  assert((unsigned)footwayCount == Footways.size());
  assert((unsigned)buildingCount == Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


  //
  // TODO: build the graph, output stats:
  //
  
  // Step 5: Add vertices to graph from nodes map
  for(auto n: Nodes){
	G.addVertex(n.first);
  }
  
  // Step 6: Add edges to graph from footways vector of nodes
  double distance;
  long long f1, f2;
  vector<long long> nodes;
  for(auto f: Footways){
	  nodes = f.Nodes;
	  for(unsigned int i = 0; i < nodes.size()-1; i++){
		  f1 = nodes[i];
		  
		  f2 = nodes[i+1];
		  
		  distance = distBetween2Points(Nodes[f1].Lat, Nodes[f1].Lon, Nodes[f2].Lat, Nodes[f2].Lon);
		  
		  G.addEdge(f1, f2, distance);
		  G.addEdge(f2, f1, distance);	
	  }
  }
	
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Navigation from building to building
  //
  string startBuilding, destBuilding;

  cout << "Enter start (partial name or abbreviation), or #> ";
  getline(cin, startBuilding);

  while (startBuilding != "#")
  {
    cout << "Enter destination (partial name or abbreviation)> ";
    getline(cin, destBuilding);

	
    //
    // TODO: lookup buildings, find nearest start and dest nodes,
    // run Dijkstra's alg, output distance and path to destination:
    //
    
	Coordinates c1, c2;
	BuildingInfo b1, b2;
	// Step 7: Lookup start and destination buildings in Buildings vector, output their information
	if(found(Buildings, b1, startBuilding, "Starting")){
		if(found(Buildings, b2, destBuilding, "Destination")){
			cout << "Starting point:" << endl;
			cout << "  " << b1.Fullname << endl;
			c1 = b1.Coords;
			cout << "  (" << c1.Lat << ", " << c1.Lon << ")" << endl;
			
			cout << "Destination point:" << endl;
			cout << "  " << b2.Fullname << endl;
			c2 = b2.Coords;
			cout << "  (" << c2.Lat << ", " << c2.Lon << ")" << endl;
			
			long long fStart = 0, fEnd = 0;
			
			// Step 8: Search footways vector for the nearest start and dest nodes
			searchFootways(Nodes, Footways, c1, c2, fStart, fEnd);
			
			// Step 9: Run Dijkstra's algorithm
			cout << endl;
			cout << "Navigating with Dijkstra..." << endl;
			vector<long long> visited;
			map<long long, double> distances;
			map<long long, long long> predecessors;
			const double INF = numeric_limits<double>::max();
			visited = Dijkstra(G, fStart, distances, predecessors);
			
			// Step 10: Output distance and path to destination
			if(distances[fEnd] == INF){
				cout << "Sorry, destination unreachable" << endl;
			}else{
				cout << "Distance to dest: " << distances[fEnd] << " miles" << endl;
				outputDijkstra(fEnd, predecessors);
			}
			
		}else{
			cout << "Destination building not found" << endl;
		}
	}else{
		cout << "Start building not found" << endl;
	}
	  
    //
    // another navigation?
    //
	  
    cout << endl;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
	  
  }

  //
  // done:
  //
  cout << "** Done **" << endl;

  return 0;
}
