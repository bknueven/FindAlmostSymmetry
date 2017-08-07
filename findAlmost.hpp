//findAlmost.hpp - header for findAlmost.cpp and related files

#ifndef FINDALMOST_HPP_
#define FINDALMOST_HPP_

#include "NautyGraph.hpp" //includes <iostream>
#include <pebbl/branching.h> //pebbl libraries
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace pebbl; //for ease...fix later?
//These macros and this struct are from libhungarian by Cyrill Stachniss #define HUNGARIAN_NOT_ASSIGNED 0 
#define HUNGARIAN_ASSIGNED 1

#define HUNGARIAN_MODE_MINIMIZE_COST   0
#define HUNGARIAN_MODE_MAXIMIZE_UTIL 1

typedef struct {
  int num_rows;
  int num_cols;
  int** cost;
  int** assignment;  
} hungarian_problem_t;

int hungarian_solve(hungarian_problem_t *__restrict__ p);

#define INF (0x7FFFFFFF)

#define hungarian_test_alloc(X) do {if ((void *)(X) == NULL) fprintf(stderr, "Out of memory in %s, (%s, line %d).\n", __FUNCTION__, __FILE__, __LINE__); } while (0)
//end libhungarian code

//forward declarations

//LTIdx allows rempresenting a lower triangular matrix as an array
inline int LTIdx(int bigger, int smaller) { return (bigger*(bigger-1))/2 + smaller; }

void initialRefineByDegrees(const DenseGraph &g, DenseGraph &perg, int B);

void refineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph & perg, int B, int* edge_use, bool always_collect_edge_use = false);

void buildCostMatrix(int i, int j,const DenseGraph & g,const DenseGraph & g_fixed,const DenseGraph & perg, int size, int inf, bool neighbors, int** biGraph, int* Ni, int* Nj);

//prints and nxm array
void printArray(int n, int m, int** array);

//independent set heuristic
int independentSetHeuristic(const DenseGraph & h);

void refineByDegreeDiffFixedEdges(const DenseGraph & g, const DenseGraph & g_fixed, DenseGraph & perg);

void refineByDegreeDiff(const DenseGraph & g, DenseGraph & perg, int B);

//begin PEBBL declarations
//enum BranchDirection { DELETE, FIX };

class findAlmostSub; //forward declaration
class findAlmostSol; //again
class parFindAlmostSub;

class findAlmost : virtual public branching
{
  friend class findAlmostSub;
  friend class parFindAlmostSub;
protected:
  NautyGraph g;
  int budget;
  bool randomBranching;
  bool orbitalBranching;
  bool justDiveLeft;
  bool localBranching;
  bool trackEdges;
  bool testBranchingStrength;
  int disjunctNumber;
  std::vector<int> rankInChild;
  std::vector<int> refinedByFixedBranch;
  int refinedByDefaultBranch;
  std::string file_name;
public:
  findAlmost() : g(), rankInChild(), refinedByFixedBranch(), file_name() { //constructor
    branchingInit(minimization);
    budget = 0;
    disjunctNumber = 1;
    randomBranching = false;
    justDiveLeft = false;
    localBranching = false;
    trackEdges = false;
    testBranchingStrength = false;
    DenseGraph::packFactor = 0.5;
    create_categorized_parameter("budget",
				 budget,
				 "<int>","0",
				 "Number of edge deletions to allow",
				 "Find Almost Symmetry",
				 ParameterNonnegative<int>());				 

    create_categorized_parameter("disjunctNumber",
				 disjunctNumber,
				 "<int>","1",
				 "Number of edges to branch on",
				 "Find Almost Symmetry",
				 ParameterPositive<int>());				 
    
    create_categorized_parameter("packFactor",
		    		 DenseGraph::packFactor,
				 "<double>","0.5",
				 "Cutoff for packing graphs in subproblems.\n0.0 specifies no packing, 1.0 specifies to always pack. Default is 0.5",
				 "Find Almost Symmetry",
				 ParameterBounds<double>(0.0,1.0));

    create_categorized_parameter("justDiveLeft",
		                 justDiveLeft,
				 "<bool>","false",
				 "Heuristic that returns the solution found by diving left",
				 "Find Almost Symmetry");

    create_categorized_parameter("localBranching",
		                 localBranching,
				 "<bool>","false",
				 "Branching decisions are made based on a budget of 1",
				 "Find Almost Symmetry");

    create_categorized_parameter("randomBranching",
		                 randomBranching,
				 "<bool>","false",
				 "Make random branching decisions",
				 "Find Almost Symmetry");

    create_categorized_parameter("trackEdges",
		                 trackEdges,
				 "<bool>","false",
				 "Track 2nd best branch candidate and see where it is in child nodes. Writes summary to file.Only works in serial mode.",
				 "Find Almost Symmetry");

    create_categorized_parameter("testBranchingStrength",
		                 testBranchingStrength,
				 "<bool>","false",
				 "Tests branching strength at root node of tree and terminates. Writes summary to file.",
				 "Find Almost Symmetry");
  
  }
  virtual ~findAlmost() {
    if (trackEdges 
#ifdef MPI_VERSION
		    && (uMPI::rank == uMPI::ioProc)
#endif
       ) {
      std::string tracked_file_name;
      std::ofstream usefile;

      tracked_file_name = file_name + ".trackedEdges.budget"+std::to_string(budget)+ ".txt";
      tracked_file_name = tracked_file_name.substr( tracked_file_name.find_last_of("/\\")+1 );
      usefile.open( tracked_file_name );
      for (auto it = rankInChild.begin(); it < rankInChild.end(); it++) usefile << *it << std::endl;
      usefile.close();
    }
    if (testBranchingStrength 
#ifdef MPI_VERSION
		    && (uMPI::rank == uMPI::ioProc)
#endif
       ) {
      std::string refined_by_fn;
      std::ofstream usefile;

      refined_by_fn = file_name + ".budget"+std::to_string(budget)+ ".refinedByFixedBranch.txt";
      refined_by_fn = refined_by_fn.substr( refined_by_fn.find_last_of("/\\")+1 );
      usefile.open( refined_by_fn );
      for (auto it = refinedByFixedBranch.begin(); it < refinedByFixedBranch.end(); it++) usefile << *it << std::endl;
      usefile.close();

      refined_by_fn = file_name + ".budget"+std::to_string(budget)+ ".refinedByDefaultBranch.txt";
      refined_by_fn = refined_by_fn.substr( refined_by_fn.find_last_of("/\\")+1 );
      usefile.open( refined_by_fn );
      usefile << refinedByDefaultBranch << std::endl; 
      usefile.close();
    }
  }

  branchSub* blankSub();

  bool setupProblem(int& argc,char**& argv) {
    FILE *fp; //file pointer
    if ((fp = fopen(argv[1], "r")) == NULL) {
      printf("Can't open %s\n", argv[1]);
      exit(EXIT_FAILURE);
    }
    //using constructor, reads dimacs files
    g = NautyGraph(fp);
    //close file
    fclose(fp);
    file_name = argv[1];
#ifdef MPI_VERSION
    if (trackEdges && (uMPI::size > 1)) {
      if (uMPI::rank == uMPI::ioProc) ucout << "--trackEges=true can only be used in serial mode. Terminating.\n";
      return false;
    }
#endif
    return true;
  }

  const NautyGraph& gRef() const { return g; }
  int getB() const { return budget; }
  int numVertices() const { return g.numVertices(); }
  int DisjunctNum() const { return disjunctNumber; } 
};

class findAlmostSub : virtual public branchSub
{
  friend class findAlmost;
private:

protected:
  findAlmost* globalPtr; //pointer to global branching object

public:
  NautyGraph g;//bGlobal()->incumbentValue) 
  EdgeList delEdges;
  DenseGraph g_fixed;
  DenseGraph perg;
  bool delNode; //node was just finished deleting
  bool boundCompFirstPass;
  int lowerBound;
  int n;
  int B; //edges deletions left at this node
  EdgeList disjunctEdges;
  Edge prev_second_choice;

  // Return a pointer to global branching object
  findAlmost* global() const { return globalPtr; }

  // Return a pointer to the base class of findAlmost
  branching* bGlobal() const { return global(); }

  // set up problem
  void setGlobalInfo(findAlmost* global_){
    globalPtr = global_;
    n = global_->numVertices();
  }

  findAlmostSub() : g(), delEdges(), g_fixed(), perg(), disjunctEdges(), prev_second_choice() {}  //default constructor
  virtual ~findAlmostSub() {}

  void boundComputation(double* controlParam); //matching problem, nauty call..malloc edge use first time through

  // in serial case this is just a wrapper for refineByMatching
  virtual void runRefineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph & perg, int B, int* edge_use, bool always_collect_edge_use = false) { 
  refineByMatching(g, g_fixed, perg, B, edge_use, always_collect_edge_use);
  }

  // nothing to do for serial case
  virtual void makeBranchingDecisionsAgree() {}

  bool candidateSolution() 
  {return false;} 

  //need to create solution class to implement this....
  solution* extractSolution()
  { 
    return NULL;
  }

  void findBranchEdge(int *); //basically does most of what splitComputation did
  //We don't want to carry edge_use between call to boundComputation() and splitComputation() 
  int splitComputation() { setState(separated); 
	  if (global()->justDiveLeft) return disjunctEdges.num(); 
	  else return (disjunctEdges.num() + 1); } //separates the subproblem, need to return number of children
  //created and setState(separated)

  void findAlmostSubAsChildOf(findAlmostSub* parent, int whichChild);
  
  //makes the child 0 <= whichChild < splitComputation()
  //returns a pointer to it...make sure globalPtr is set an all local data are initialized
  branchSub* makeChild(int whichChild) {
    findAlmostSub * temp = new findAlmostSub;
    DEBUGPR(600, ucout << "new findAlmostSub, address: " << (void*)temp << std::endl);
    temp->findAlmostSubAsChildOf(this, whichChild);
    return temp;
  } 

  void setRootComputation() {  //called when object should be made into a root subproblem
    g = global()->gRef();
    B = global()->getB(); //make delEdges of size B
    g_fixed = DenseGraph(n); //empty graph on n vertices
    perg = DenseGraph(n); //empty graph on n vertices
    delEdges = EdgeList(B);
    disjunctEdges = EdgeList(global()->DisjunctNum());
    for(int i = 0; i < n; i++){
      for(int j = 0; j < i; j++){
	perg.addEdge(i,j);
      }
    } //make perg complete with self-loops
    delNode = true; // we'll do root node refinement in boundComputation
    boundCompFirstPass = true;    
    if (global()->trackEdges){ 
      prev_second_choice.u = -1; prev_second_choice.v = -1;
    }
  }

  //helper function for testing strong branching (also could be used in place of runRefineByMatching to change
  //PEBBL's behavior)
  void refineUntilNone(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph & perg, int B, int* edge_use);

  //helper function for strong branching
  int simulateLeftBranch(const Edge &branch_edge);

  //helper fucntion to write to global variables
  virtual void testBranchingWriteToGlobal(int refinedByDefaultBranch, std::vector<int> & refinedByFixedBranch) {
    global()->refinedByDefaultBranch = refinedByDefaultBranch;
    global()->refinedByFixedBranch = refinedByFixedBranch;
  }

};

class findAlmostSol : public solution
{
private:
  int n; //number of vertices
public:

  //override some basic solution class methods
  const char* typeDescription() const { return "Almost symmetries solution"; }
  void printContents(std::ostream& s);  //implement something to print edge deletions && numOrbits && orbits
  
  EdgeList delEdges;
  int* orbits;
  double group_size;

  //for hashing
  size_type computeHashValue() { return delEdges.hashValue(); }
  bool duplicateOf(findAlmostSol& other) { return (other.delEdges == delEdges); }

  //constuctor
  findAlmostSol(findAlmost* global_, EdgeList& delEdges_, const int* orbits_, int numOrbits_, double group_size_) :
    solution(global_), delEdges(delEdges_)//use EdgeList copy constructor
  {
    value = numOrbits_; //member of parent class, sets interal representation of solution 
    n = global_->numVertices();
    //deep copy orbits
    orbits = new int[n];
    std::copy(&orbits_[0],&orbits_[n],orbits);
    group_size = group_size_;
  }

  //create reference solution for PEBBL
  findAlmostSol(findAlmost* global_) : solution(global_), delEdges(global_->getB()),
       				       n(global_->numVertices()), orbits(NULL), group_size(-1.) {}
  
  //solution to copy
  findAlmostSol(findAlmostSol* toCopy) : solution(toCopy), delEdges(toCopy->delEdges.capacity()), 
					 n(toCopy->n), orbits(NULL), group_size(-1.) {}

  //destructor
  ~findAlmostSol(){delete [] orbits;}

  solution* blankClone() { return new findAlmostSol(this); }

  void packContents(utilib::PackBuffer& outBuf);
  void unpackContents(utilib::UnPackBuffer& inBuf);
  int maxContentsBufSize();

};
#endif
