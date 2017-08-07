// NautyGraph.cpp -- C++ interface for dense nauty graphs with dynamic allocation
// Version 0.0.1 -- Last modified 1/27/15
#include "NautyGraph.hpp"

float DenseGraph::packFactor = 1.0;

//good constructor
NautyGraph::NautyGraph(int _n)
  : DenseGraph(_n)
{
  //dynamically allocate memory
  orbits = new int [n];
  //use default options for nauty
  options = {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL};
  //end of constructor (for now)
}

NautyGraph::NautyGraph(FILE *fp)
  : DenseGraph(fp)
{
  //dynamically allocate memory
  orbits = new int [n];
  //use default options for nauty
  options = {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL};
  //end of constructor (for now)
}

//default constructor ... not very useful
NautyGraph::NautyGraph()
  : DenseGraph()
{
  orbits = NULL;
  options = {0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH,NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL};
}

//destructor
NautyGraph::~NautyGraph() {
  //free memory
  delete [] orbits;
}

//copy assignment operator
NautyGraph & NautyGraph::operator=(const NautyGraph & ng){
  if (this == &ng) return *this; //if object assigned to itself, done
  //free memory
  DenseGraph::operator=(ng);
  delete [] orbits;
  //copy everything as in copy constructor
  options = ng.options;
  stats = ng.stats;
  orbits = new int [n];
  std::copy(&ng.orbits[0], &ng.orbits[n], orbits);
  //return reference to invoking object
  return *this;
}
    

//copy constructor
NautyGraph::NautyGraph(const NautyGraph & ng)
  : DenseGraph(ng)
{
  options = ng.options;
  stats = ng.stats;
  orbits = new int [n];
  std::copy(&ng.orbits[0], &ng.orbits[n], orbits);
}

//for calling nauty, simple version
void NautyGraph::callNauty() {
 options.defaultptn = TRUE; 
 int * lab = new int [n];
 int * ptn = new int [n];
 densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
 delete [] lab;
 delete [] ptn;
}

//DenseGraph methods

//good constructor
DenseGraph::DenseGraph(int _n){
  n = _n;
  m = SETWORDSNEEDED(n); //get m
  //dynamically allocate memory
  g = new graph [m*n];
  deg_g = new int [n];
  for(int i = 0; i < n; i++) deg_g[i] = 0; //initialize degrees to 0
  edges = 0; //initialize edges to 0
  //initialize g to empty graph
  EMPTYGRAPH(g,m,n);
  //end of constructor (for now)
}

//file constructor
DenseGraph::DenseGraph(FILE *fp){
  g = NULL;
  deg_g = NULL;
  dimacsReader(fp);
}

//default constructor ... not very useful
DenseGraph::DenseGraph(){
  n = 0;
  m = 0;
  deg_g = NULL;
  g = NULL;
  edges = 0;
}

//destructor
DenseGraph::~DenseGraph() {
  //free memory
  delete [] deg_g;
  delete [] g;
}

//copy assignment operator
DenseGraph & DenseGraph::operator=(const DenseGraph & dg){
  if (this == &dg) return *this; //if object assigned to itself, done
  //free memory
  delete [] deg_g;
  delete [] g;
  //copy everything as in copy constructor
  n = dg.n;
  m = dg.m;
  edges = dg.edges;
  g = new graph [m*n];
  deg_g = new int [n];
  std::copy(&dg.g[0], &dg.g[m*n], g);
  std::copy(&dg.deg_g[0], &dg.deg_g[n], deg_g);
  //return reference to invoking object
  return *this;
}

void DenseGraph::dimacsReader(FILE *f){
  int c, ch, i, j, fileEdges;
  bool allocated = false;
  char name[10];
  //free memory if already allocated
  delete [] deg_g;
  delete [] g;
  //read in file
  while(( ch = getc(f)) != EOF ){
    switch(ch) {
    case 'p':
      if(!allocated) {
	if (fscanf(f, " %s %d %d", name, &n, &fileEdges) != 3) {
	  fprintf(stderr, "Something went wrong reading in the graph size.\n");
	  exit(EXIT_FAILURE);
	}
	m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	g = new graph [m*n];
	deg_g = new int [n];
	allocated = true;
	EMPTYGRAPH(g,m,n); //initialize g to empty
	for(int i = 0; i < n; i++) deg_g[i] = 0; //initialize degrees to 0
	edges = 0; //initialize edges to 0
      } else {
	fprintf(stderr, "We found more then one line beginning with 'p'!\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'c': //comment line, we'll print it out
      fprintf(stdout, "Comment:");
      //printLine(f); inline this function
      while((c = getc(f)) != EOF && c != '\n') putchar(c);
      if (c == EOF) ungetc(c,f); 
      //back this up so there's something to read if getc is queried again
      fprintf(stdout, "\n");
      break;
    case 'e':
      if(!allocated){
	fprintf(stderr, "We found an edge before the graph's memory was allocated!\n");
	exit(EXIT_FAILURE);
      } else {
	if (fscanf(f, " %d %d", &i, &j) != 2){
	  fprintf(stderr, "Trouble reading edge.\n");
	  exit(EXIT_FAILURE);
	} else {
	  //printf("Adding edge (%d, %d)\n", i-1, j-1);
	  if(!(isEdge(i-1,j-1))) addEdge(i-1,j-1);
	}
      }
      break;
    case '\n':
      break;
    case EOF :
      return;
    default:
      fprintf(stderr, "We found something we weren't expecting reading the file, quitting.\n");
      exit(EXIT_FAILURE);
    }
  }
  if(edges != fileEdges && edges != fileEdges/2){ 
    fprintf(stderr, "Number of edges doesn't match; quitting.\n"); 
    exit(EXIT_FAILURE); 
  }
}
    

//copy constructor
DenseGraph::DenseGraph(const DenseGraph & dg){
  n = dg.n;
  m = dg.m;
  edges = dg.edges;
  g = new graph [m*n];
  deg_g = new int [n];
  std::copy(&dg.g[0], &dg.g[m*n], g);
  std::copy(&dg.deg_g[0], &dg.deg_g[n], deg_g);
}

void DenseGraph::write(utilib::PackBuffer& outBuffer) const {
	//go in order of listing
	outBuffer << n;
	outBuffer << m;
	bool asEdgeList = ((2*edges + 1)*sizeof(int) < packFactor*m*n*sizeof(graph));
	outBuffer << asEdgeList;
	if(asEdgeList) { //pack up the graph as list of adjacentcies, assumes simple graph
	  //std::cout << "Packing as EdgeList!\n";
	  EdgeList gAsEdgeList(edges);
	  Row gi; int i,j;
	  for(i=0, gi=row(i); i<n; i++, gi.nextRow()){
	   for(j=i; (j = gi.nextVertex(j)) >=0; ){
	     gAsEdgeList.add(i,j);
	   }
	  }
	  outBuffer << gAsEdgeList;
	} else {
	  //std::cout << "Packing as packed graph!\n";
	  for(int i=0; i < m*n; i++){
	    outBuffer << g[i];
	  }
	  outBuffer << edges;
	}
}

void DenseGraph::read(utilib::UnPackBuffer& inBuffer) {
	inBuffer >> n;
	inBuffer >> m;
	g = new graph[m*n];
	deg_g = new int[n];
	bool asEdgeList;
	inBuffer >> asEdgeList;
	if(asEdgeList) {
	  //std::cout << "Unpacking as EdgeList!\n";
	  EdgeList gAsEdgeList;
	  inBuffer >> gAsEdgeList;
	  EMPTYGRAPH(g,m,n); //make sure g is empty
	  for(int i = 0; i < n; i++) deg_g[i] = 0; //initialize degrees to 0
	  for(int i=0; i<gAsEdgeList.num(); i++) addEdge(gAsEdgeList.getu(i),gAsEdgeList.getv(i));
	} else {
	  //std::cout << "Unpacking as packed graph!\n";
	  for(int i=0; i < m*n; i++){
	    inBuffer >> g[i];
	  }
	  inBuffer >> edges;
	  //reconsturct deg_g
	  int i; graph* gp;
	  for(i=0, gp = g; i < n; ++i, gp += m) {
		deg_g[i] = setsize(gp,m);
	  }
	}  
}

void NautyGraph::write(utilib::PackBuffer& outBuffer) const {
	DenseGraph::write(outBuffer); //write DenseGraph contents
	for(int i=0; i < n; i++){
		outBuffer << orbits[i];
	}
	outBuffer << stats.numorbits;
	outBuffer << stats.grpsize1;
	outBuffer << stats.grpsize2;
}

void NautyGraph::read(utilib::UnPackBuffer& inBuffer) {
	DenseGraph::read(inBuffer); //read in DenseGraph contents
	orbits = new int[n];
	for(int i=0; i<n; i++){
		inBuffer >> orbits[i];
	}
	inBuffer >> stats.numorbits;
	inBuffer >> stats.grpsize1;
	inBuffer >> stats.grpsize2;
}

int DenseGraph::components() const{
  int i,tail,w,k,compNum;
  set *gw;
  int * queue = new int [n];
  int * visited = new int [n];

  for (i = 0; i < n; ++i) visited[i] = 0;

  compNum = 0;
  tail = -1;
  
  for (i = 0; i < n; ++i) {
    if (!visited[i]){
      compNum++;
      //printf("Component %d: ", compNum);
      visited[i] = 1;
      queue[++tail] = i; //printList((tail+1), queue);
      while(tail >= 0) { //i.e., while stuff is still on the sack 
	w = queue[tail--]; //pop w from stack
	//printf("%d ", w);
	gw = GRAPHROW(g,w,m);
	for (k = -1; (k = nextelement(gw,m,k)) >= 0;) {
	  if (!visited[k]) { //we've found another node in this components
	    visited[k] = 1;
	    queue[++tail] = k;
	  }
	}
      }
      //printf("tail = %d, i = %d\n", tail, i);
    }
  }
  delete [] queue;
  delete [] visited;
  return compNum;
}

int DenseGraph::components(int* lab, int* ptn) const{
  int i,tail,w,k,compNum,ptnIndex;
  set *gw;
  int * queue = new int [n];
  int * visited = new int [n];
  for (i = 0; i < n; ++i) visited[i] = 0;

  compNum = 0;
  tail = -1;
  ptnIndex = -1;
    
  for (i = 0; i < n; ++i) {
    if (!visited[i]){
      compNum++;
      ptnIndex++; //start new partition
      lab[ptnIndex] = i; //add this to the partition
      ptn[ptnIndex] = 1; //assume it's not the end
      //printf("Component %d: ", compNum);
      visited[i] = 1;
      queue[++tail] = i; //printList((tail+1), queue);
      while(tail >= 0) { //i.e., while stuff is still on the sack 
	w = queue[tail--]; //pop w from stack
	//printf("%d ", w);
	gw = GRAPHROW(g,w,m);
	for (k = -1; (k = nextelement(gw,m,k)) >= 0;) {
	  if (!visited[k]) { //we've found another node in this components
	    visited[k] = 1;
	    queue[++tail] = k;
	    ptnIndex++; //add to partition
	    lab[ptnIndex] = k; //add this to the partition
	    ptn[ptnIndex] = 1; //assume it's not the end	      
	  }
	}
      }
      ptn[ptnIndex] = 0; //once we've got here, this is the end of the partition
	//printf("tail = %d, i = %d\n", tail, i);
    }
  }
  delete [] queue;
  delete [] visited;
  //printf("\nlab: "); printList(n, lab);
  //printf("ptn: "); printList(n, ptn);
  return compNum;
}

//says if a vertex is isolated
bool DenseGraph::isIsolated(int i) const{
  //if ((deg(i) == 2) && isEdge(i,i)) return true;
  if (deg(i) == 0) return true;
  return false;
}

//prints a DenseGraph in DIMACS format
std::ostream & operator<<(std::ostream & os, const DenseGraph & dg){
  DenseGraph::Row gi; //for graph row i
  int i,j;
  os << "p edge " << dg.numVertices() << " " << dg.numEdges() << std::endl;
  for(i = 0, gi = dg.row(i); i < dg.numVertices(); i++, gi.nextRow()){
    for(j = i-1; (j = gi.nextVertex(j)) >= 0; ) { //check upper triangle, including self-loops
      os << "e " << i+1 << " " << j+1 << std::endl;
    }
  }
  return os;
}

//Begin EdgeList definitions
Edge::operator size_t(){
  size_t _v = v; //std::cout << std::hex << "_v now: " << _v << '\n';
    _v <<= bits; //sift _v to the left bits bits
    //std::cout << "After bit shift _v = "  << std::hex << _v << '\n';
    //std::cout << "u = " << std::hex << u << '\n';
    //std::cout << "Returning _v + u = " << std::hex << _v + u << '\n';
    return _v + u;
}


Edge::Edge(size_t _edge){
  //this restricts our graph sizes to have of the bits of the computer
  //std::cout << "_edge = "  << std::hex << _edge << '\n';
  v = _edge >> bits;
  //std::cout << "v = " << std::hex << v << '\n';
  u = _edge & mask;
  //std::cout << "u = " << std::hex << u << '\n';
}
std::ostream & operator<<(std::ostream & os, const Edge & e){
  os << '(' << (e.u+1) << ',' << (e.v+1) << ')';
  return os;
}

EdgeList::EdgeList(int _max) {
  max = _max;
  n = 0;
  edges = new Edge [max];
}

EdgeList::EdgeList() {
  max = 0; n = 0; edges = NULL;
}

EdgeList::~EdgeList() {
  delete [] edges;
}

EdgeList::EdgeList(const EdgeList & e){
  max = e.max;
  n = e.n;
  edges = new Edge [max];
  std::copy(&e.edges[0],&e.edges[n],edges);
}

EdgeList & EdgeList::operator=(const EdgeList & e){
  if (this == &e) return *this;
  delete [] edges;
  max = e.max;
  n = e.n;
  edges = new Edge [max];
  std::copy(&e.edges[0],&e.edges[n],edges);
  return *this;
}

void EdgeList::write(utilib::PackBuffer & outBuffer) const {
	outBuffer << max;
	outBuffer << n;
	for(int i=0; i<n; i++) {
		outBuffer << edges[i].u << edges[i].v;
	}
}

void EdgeList::read(utilib::UnPackBuffer & inBuffer){
	inBuffer >> max;
	inBuffer >> n;
	edges = new Edge[max];
	for(int i=0; i<n; i++){
		inBuffer >> edges[i].u >> edges[i].v;
	}
}

bool EdgeList::operator==(const EdgeList & e) const{
  if (n != e.n) return false;
  //else we'll do the check
  int i,j;
  for(i=0; i<n; i++) { //look through this's edges for a match in e
    for(j=0; j<n; j++) {
      if (edges[i] == e.edges[j]) break;
    }
    if (j == n) return false;
  }
  for(i=0; i<n; i++) { //look through e's edges for a match in this
    for(j=0; j<n; j++) {
      if (e.edges[i] == edges[j]) break;
    }
    if (j == n) return false;
  }
  return true;
}

size_t EdgeList::hashValue(){
  size_t hash = 0;
  for(int i = 0; i < n; i++) hash ^= edges[i]; //use bitwise XOR
  return hash;
}

std::ostream & operator<<(std::ostream & os, const EdgeList & e){
  for(int i = 0; i < e.n; i++){
    os << e.edges[i];
  }
  return os;
}
//End EdgeList definitions
