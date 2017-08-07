// NautyGraph.hpp -- C++ interface for nauty dense graphs with dynamic allocation
// Version 0.0.1 -- Last modified 1/27/15
#ifndef NAUTYGRAPH_HPP_
#define NAUTYGRAPH_HPP_

#include "nauty.h"
#include "naututil.h"
#include <iostream>
#include <utilib/PackObject.h> //for packing shit

/* add macros for removing edges */
#define DELONEARC1(g,v,w,m) (g)[v] &= ~BITT[w]
#define DELONEEDGE1(g,v,w,m) { DELONEARC1(g,v,w,m); DELONEARC1(g,w,v,m); }

#define DELONEARC0(g,v,w,m) DELELEMENT0(GRAPHROW0(g,v,m),w)
#define DELONEEDGE0(g,v,w,m) { DELONEARC0(g,v,w,m); DELONEARC0(g,w,v,m); }

#define ISARC1(g,v,w,m) (((g)[v] & BITT[w]) != 0)

#define ISARC0(g,v,w,m) ISELEMENT0(GRAPHROW0(g,v,m),w)

#if  (MAXM==1) && defined(ONE_WORD_SETS)
#define DELONEARC DELONEARC1
#define DELONEEDGE DELONEEDGE1
#define ISARC ISARC1
#else
#define DELONEARC DELONEARC0
#define DELONEEDGE DELONEEDGE0
#define ISARC ISARC0
#endif
/* end macros */

/*namespace DenseGraphHelpers {
  bool setEmpty(set *set1, int m);
  }*/

class DenseGraph : public utilib::PackObject 
{
protected:
  graph *g;
  int n; //number of vertices
  int m; //=SETWORDSNEEDED(n)
private:
  int *deg_g; //for keeping track of degree of each vertex
  int edges; //keep how many edges in graph
public:
  static float packFactor;
  class Row
  {
  private:
    set * row;
    int m;
  public:
    Row() {row = NULL;}
    Row(set * _row, int _m) {row = _row; m = _m;}
    //copying and everything should be done memberwise...
    void nextRow() { row += m; }
    int nextVertex(int i = -1) const {return nextelement(row,m,i);}
  };
  DenseGraph(int _n); //constructor that just takes the number of vertices
  DenseGraph(FILE *fp); //constructor that takes a file pointer to DIMACS file
  DenseGraph(); //default constructor, not very useful
  virtual ~DenseGraph(); //destructor (frees memory)
  DenseGraph & operator=(const DenseGraph &); //copy assignment operator
  DenseGraph(const DenseGraph &); //copy constructor
  void dimacsReader(FILE *fp);
  void addEdge(int i, int j) {ADDONEEDGE(g,i,j,m); edges++; deg_g[i]++; deg_g[j]++;}
  void delEdge(int i, int j) {DELONEEDGE(g,i,j,m); edges--; deg_g[i]--; deg_g[j]--;}
  bool isEdge(int i, int j) const {return ISARC(g,i,j,m);}// only safe for undirected graphs
  //bool isEdge(int i, int j){return (ISARC(g,i,j,m) && ISARC(g,j,i,m));} //or use this
  void putGraph(FILE *f = stdout) const {putgraph(f,g,0,m,n);}
  int numVertices() const {return n;}
  int numEdges() const {return edges;}
  int deg(int i) const {return deg_g[i];}
  int components() const;
  int components(int* lab, int* ptn) const;
  Row row(int i = 0) const {return Row(GRAPHROW(g,i,m),m);}
  bool isIsolated(int i) const;
  virtual void write(utilib::PackBuffer&) const;
  virtual void read(utilib::UnPackBuffer&);
  virtual int packSize(){return (3*sizeof(int) + n*m*sizeof(graph));}
  friend std::ostream & operator<<(std::ostream & os, const DenseGraph & dg);
};

class NautyGraph :public DenseGraph
{
private:
  optionblk options; //default nauty options...may want to change later
//  statsblk stats; //for nauty's stats
public:
  statsblk stats; //for nauty's stats ... clean up, need public for parFindAlmostSub::unpack
  int *orbits; //for orbits (want deep copy so we can access for each graph once done)
  // is making this public a good idea? Nauty treats it as write-only...so we can't do much harm...
  NautyGraph(int _n); //constructor that just takes the number of vertices
  NautyGraph(FILE *fp);
  NautyGraph(); //default constructor, not very useful
  ~NautyGraph(); //destructor (frees memory)
  NautyGraph & operator=(const NautyGraph &); //copy assignment operator
  NautyGraph(const NautyGraph &); //copy constructor
  void callNauty();
  void callNauty(int* _lab, int* _ptn, bool usePtn = TRUE) {
    options.defaultptn = !usePtn;
    densenauty(g,_lab,_ptn,orbits,&options,&stats,m,n,NULL);}
  int numOrbits() const {return stats.numorbits;}
  double groupSize() const {return stats.grpsize1*pow(10., stats.grpsize2);} 
  virtual void write(utilib::PackBuffer&) const; //warning: These only load stats.numorbits, stats.grpsize1, stats.grpsize2
  virtual void read(utilib::UnPackBuffer&);      //IF other stats need to go they need to be added
  virtual int packSize(){return (DenseGraph::packSize() + (n+2)*sizeof(int) + sizeof(double));}
};

struct Edge { //provides at least 64bits of resolution
  static const short bits = sizeof(size_t)*4; //we the numb of bits in half of size_type
  static const size_t mask = (size_t(1) << bits) - 1;
  int u,v; 
  bool operator==(const Edge & e){ return (((u == e.u) && (v == e.v)) //2nd for safety
					   /*|| ((v == e.u) && (u == e.v))*/) ; }
  operator size_t();
  Edge(size_t _edge);
  Edge() {}
  Edge(int _u, int _v) {u = _u; v = _v;}
  friend std::ostream & operator<<(std::ostream & os, const Edge & e);
};

class EdgeList : public utilib::PackObject 
{
private:
  int max; //for the most edges we can fit
  int n; //how many edges we have now
  Edge * edges;
public:
  EdgeList(int _max); //constructor that just takes the max num of edges
  EdgeList(); //default constructor
  ~EdgeList(); //destructor (frees memory)
  EdgeList & operator=(const EdgeList &); //copy assignment operator
  bool operator==(const EdgeList &) const;
  EdgeList(const EdgeList &); //copy constructor
  void add(int _u, int _v) {edges[n].u = _u; edges[n].v = _v; n++;}
  void add(const Edge & e) {edges[n] = e; n++;}
  int num() const {return n;}
  int capacity() const {return max;}
  int left() const {return (max - n);}
  Edge get(int i) const {return edges[i];}
  int getu(int i) const {return edges[i].u;}
  int getv(int i) const {return edges[i].v;}
  size_t hashValue();
  virtual void write(utilib::PackBuffer&) const;
  virtual void read(utilib::UnPackBuffer&);
  virtual int packSize(){return (max + max + 2)*sizeof(int);}
  friend std::ostream & operator<<(std::ostream & os, const EdgeList & e);
};

#endif
