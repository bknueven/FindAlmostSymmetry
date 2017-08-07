// parFindAlmost.hpp - header file for findAlmostPar.cpp
// and class declarations.

#ifndef PARFINDALMOST_HPP_
#define PARFINDALMOST_HPP_

#include "findAlmost.hpp"
#include <pebbl/parBranching.h>

using namespace pebbl;

void parRefineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph & perg, int B, int* edge_use, bool always_collect_edge_use = false);

class parFindAlmost :
	public parallelBranching,
	public findAlmost
{
public:
	parFindAlmost() {}
	~parFindAlmost() {}
	void reset(bool VBflag = true){
		findAlmost::reset();
		registerFirstSolution(new findAlmostSol(this));
		parallelBranching::reset();
	}
	parallelBranchSub* blankParallelSub();
	void pack(utilib::PackBuffer& outBuffer);
	void unpack(utilib::UnPackBuffer& inBuffer);
	int spPackSize();
	//will automatically call application serial setup routine
	bool setup(int& argc, char**& argv) { return parallelBranching::setup(argc,argv);}
};

class parFindAlmostSub :
	public parallelBranchSub,
	public findAlmostSub
{
protected:
	parFindAlmost* globalPtr;

public:
	//int nodes;
	parFindAlmost* global() const { return globalPtr; }
	parallelBranching* pGlobal() const { return global(); }

	//constructor
	parFindAlmostSub() {}
	~parFindAlmostSub() {
           DEBUGPR(600, ucout << "parFindAlmostSub destructor called for " << this 
			   << " at address " << (void*)this << std::endl); }
	void setGlobalInfo(parFindAlmost* global_){
		globalPtr = global_;
		findAlmostSub::setGlobalInfo(global_);
	}

	void pack(utilib::PackBuffer& outBuffer);
	void unpack(utilib::UnPackBuffer& inBuffer);

	void parFindAlmostSubAsChildOf(parFindAlmostSub* parent, int whichChild){
		globalPtr = parent->globalPtr;
		//initialize other pieces
		findAlmostSubAsChildOf(parent,whichChild);
	};
	
	parallelBranchSub* makeParallelChild(int whichChild) {
		parFindAlmostSub *tmp = new parFindAlmostSub;
		DEBUGPR(600, ucout << "new parFindAlmostSub from parFindAlmostSub::makeParallelChild(int), address: " <<(void*) tmp << std::endl);	
		tmp->parFindAlmostSubAsChildOf(this, whichChild);
		return tmp;
	}

  	// in parallel case this wrapper determines if refineByMatching or parRefineByMatching is called
  	void runRefineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph & perg, int B, int* edge_use, bool always_collect_edge_use = false);

  	// in parallel case this forces branching decisions to agree when
	// the solver is ramping up and there is some randomness in branching
  	void makeBranchingDecisionsAgree(); 

	 
        //helper fucntion to write to global variables
        void testBranchingWriteToGlobal(int refinedByDefaultBranch, std::vector<int> & refinedByFixedBranch) {
		if (uMPI::rank == uMPI::ioProc) {
          		global()->refinedByDefaultBranch = refinedByDefaultBranch;
          		global()->refinedByFixedBranch = refinedByFixedBranch;
		}
        }

};

#endif
