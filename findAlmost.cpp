#include "parFindAlmost.hpp"
#include <unistd.h>

/* being driver */
int main(int argc, char** argv) {
  std::srand( unsigned (std::time(0) + getpid()) );
  return driver<findAlmost,parFindAlmost>(argc, argv);
}
/* end driver */

/* begin findAlmost methods */
branchSub* findAlmost::blankSub() {
  findAlmostSub* newSP = new findAlmostSub;
  newSP->setGlobalInfo(this);
  return newSP;
}
/* end findAlmost methods */

/* begin findAlmostSub methods */
void findAlmostSub::boundComputation(double* controlParam) {  
  DEBUGPR(700, ucout << "findAlmostSub::boundComputation(double*) called for " << this
	                              << " at address " << (void*)this << std::endl);   
  // if the flag testBranchingStrength is set,
  // this only works at the root node, and then terminates
  // by killing the root node
  // TODO: extend to work at other nodes?
  if(global()->testBranchingStrength) {
    // do initial refinement pass
    refineByDegreeDiff(g,perg,B);
    {
      int lab[n];
      int ptn[n];
      perg.components(lab,ptn);
      g.callNauty(lab,ptn);
    }
    
    // do root node refinement
    int * edge_use = new int [LTIdx(n,0)](); //zero-initize new array of ints 
    refineUntilNone(g, g_fixed, perg, B, edge_use);
    findBranchEdge(edge_use);
    delete [] edge_use;


    std::vector<int> refinedByFixedBranch;

    // need an escape in case there's no edge to branch on because the refinement was perfect
    if(disjunctEdges.num() == 0) {
      testBranchingWriteToGlobal(0, refinedByFixedBranch);
      setState(dead);
      return;  
    }

    // collect strong branching statistics for left (fixing) child
    // first for the branch we would have taken based on the branching rule
    // (assuming we just branch on one edge at a time..)
    int refinedDefaultBranch = simulateLeftBranch( disjunctEdges.get(0) );

    // now collect for all other edges in g
    refinedByFixedBranch.reserve(g.numEdges());
    DenseGraph::Row gi;
    // iterate over all edges
    for (int i = 0; i < n; i++){
      gi = g.row(i);
      for (int j = i; (j = gi.nextVertex(j)) >= 0; ){ //check upper triangle only
	Edge sim_branch_edge = Edge(i,j);
	refinedByFixedBranch.push_back(simulateLeftBranch(sim_branch_edge));
      }
    }

    testBranchingWriteToGlobal(refinedDefaultBranch, refinedByFixedBranch);
    setState(dead);	  
    return;
  }

  if(boundCompFirstPass) { //called after node is first deleted
    if(delNode) { 
      refineByDegreeDiff(g,perg,B);
      {
        int lab[n];
        int ptn[n];
        perg.components(lab,ptn);
        g.callNauty(lab,ptn);
      }
      //get new lowerBound 
      lowerBound = independentSetHeuristic(perg);
      DEBUGPR(500, ucout << "Right branch, after refineByDegreeDiff lowerBound = " << lowerBound << std::endl);
      bound = lowerBound;
      DEBUGPR(300, ucout << "Right branch, nauty called. Edge deletions here: " << delEdges << std::endl
	      << "Orbits here: " << g.orbits << std::endl);
      solution * sol = new findAlmostSol(global(),delEdges,g.orbits,g.numOrbits(),g.groupSize());
      global()->foundSolution(sol);
      if( B <= 0 || g_fixed.numEdges() >= g.numEdges() || lowerBound >= g.numOrbits()) {
	DEBUGPR(300, if(B<=0) ucout << "Pruning here because B<=0\n");
	DEBUGPR(300, if(g_fixed.numEdges() >= g.numEdges()) ucout << "Pruning here because we fixed all the edges\n");
	DEBUGPR(300, if(lowerBound >= g.numOrbits()) ucout << "Pruning here because lowerBound >= g.numOrbits()\n");
	DEBUGPR(400, ucout << "g_fixed here: " << g_fixed << std::endl);
	setState(dead); //fathom node
	return;
      }
    } else { //delNode is false, so this was a fixed node.
      //it's possible we've fixed all the edges, in which case we want to prune
      if (g_fixed.numEdges() >= g.numEdges()) {
	DEBUGPR(300, ucout << "Pruning because g_fixed.numEdges() >= g.numEdges()\n");
	DEBUGPR(400, ucout << "g_fixed here: " << g_fixed << std::endl);
	setState(dead);
	return;
      }
      //this is very cheap, so we might as well do it here
      refineByDegreeDiffFixedEdges(g, g_fixed, perg);
      //get new lowerBound 
      lowerBound = independentSetHeuristic(perg);
      DEBUGPR(500, ucout << "Left branch, after refineByDegreeDiffFixedEdges lowerBound = " << lowerBound << std::endl);
      bound = lowerBound;
      if ( lowerBound >= g.numOrbits() /* || fixBranch->lowerBound >= bGlobal()->incumbentValue*/ ){ //PEBBL should do the latter
	DEBUGPR(300, ucout << "Pruning because comp>=g.numOrbits().\n");
	DEBUGPR(400, ucout << "g_fixed here: " << g_fixed << std::endl);
	setState(dead);
	return;
      }     
    }
    boundCompFirstPass = false;
    return; //this will force PEBBL to update these bounds before doing a matching pass
  }
  int * edge_use = new int [LTIdx(n,0)](); //zero-initize new array of ints
  int perg_edges_bf = perg.numEdges();
  // new virtual function
  runRefineByMatching(g, g_fixed, perg, B, edge_use);
  lowerBound = independentSetHeuristic(perg);
  DEBUGPR(500, ucout << "After refineByMatching lowerBound = " << lowerBound << std::endl);
  bound = lowerBound;
  if(lowerBound >= g.numOrbits()){
    delete [] edge_use;	  
    DEBUGPR(300, if(lowerBound >= g.numOrbits()) ucout << "Pruning here because lowerBound >= g.numOrbits()\n");
    DEBUGPR(400, ucout << "g_fixed here: " << g_fixed << std::endl);
    setState(dead);
    return;
  }
  if(perg_edges_bf == perg.numEdges()){
    if (global()->localBranching) {
	DenseGraph perg_copy = perg;
	DenseGraph g_fixed_copy = g_fixed;
	std::fill(edge_use, edge_use+LTIdx(n,0), 0);
	int B_heur = 1;
	//refineByDegreeDiff(g, perg_copy, B_heur);
	// new virtual function
	runRefineByMatching(g, g_fixed_copy, perg_copy, B_heur, edge_use, true);
    }
    findBranchEdge(edge_use);
    delete [] edge_use;

    // new virtual function
    makeBranchingDecisionsAgree();

    if (disjunctEdges.num() == 0) {
	setState(dead);
    } else {
    	setState(bounded);
    }
    return;
  }
  delete [] edge_use;  
}

// returns the number refined on simulated fixing disjuction
int findAlmostSub::simulateLeftBranch(const Edge &branch_edge) {
  // copy perg, g_fixed
  DenseGraph perg_branch = perg;
  DenseGraph g_fixed_branch = g_fixed;
  
  // add edge to g_fixed
  g_fixed_branch.addEdge(branch_edge.u, branch_edge.v);
  int perg_edges_start = perg_branch.numEdges();

  // may as well do this first
  refineByDegreeDiffFixedEdges(g, g_fixed_branch, perg_branch);
  
  int * edge_use = new int [LTIdx(n,0)]();
  refineUntilNone(g, g_fixed_branch, perg_branch, B, edge_use);
  delete [] edge_use;

  return (perg_edges_start - perg_branch.numEdges());
  
}

void findAlmostSub::refineUntilNone(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph &perg, int B, int* edge_use) {
  for(;;) {
    std::fill(edge_use, edge_use+LTIdx(n,0), 0);
    int perg_edges_bf = perg.numEdges();
    // new virtual function
    runRefineByMatching(g, g_fixed, perg, B, edge_use);
    if(perg_edges_bf == perg.numEdges()) return;
  }
}

void findAlmostSub::findBranchEdge(int * edge_use) {
  DEBUGPR(700, ucout << "findAlmostSub::splitComputation() called for " << this
	                              << " at address " << (void*)this << std::endl); 
  int i,j,k,max,best_i,best_j;
  DenseGraph::Row gi;
  if (global()->trackEdges && (prev_second_choice.u != -1)) {
    int place = 0;
    int val = edge_use[LTIdx(prev_second_choice.v, prev_second_choice.u)];
    for (i = 0, gi = g.row(i) ; i < n; ++i, gi.nextRow()){
      for (j = i; (j = gi.nextVertex(j)) >= 0; ){
	if ( edge_use[LTIdx(j,i)] > val ) place++;
      }
    }
    global()->rankInChild.push_back(place);
  }
  if(!(global()->randomBranching)){ 
    for (k = 0; k < global()->DisjunctNum()+int(global()->trackEdges) ; ++k) {
      max=0; best_i = -1; best_j = -1;
      for (i = 0, gi = g.row(i) ; i < n; ++i, gi.nextRow()){
        for (j = i; (j = gi.nextVertex(j)) >= 0; ){
          DEBUGPR(1100, /*if(uMPI::rank == uMPI::ioProc)*/ ucout << "edge_use(" << i << ',' << j << ") = "<< edge_use[LTIdx(j,i)] << std::endl);	
          if (max < edge_use[LTIdx(j,i)]){
            max = edge_use[LTIdx(j,i)];
	    best_i = i;
	    best_j = j;
          }
        }   
      }
      //if we didn't find a better edge on this pass
      if (max == 0) {
	break;
      } else {
	//zero this out so we don't see him again
	if ( k < global()->DisjunctNum() ) {
	  edge_use[LTIdx(best_j, best_i)] = 0;
	  disjunctEdges.add(best_i, best_j);
	} else {
	  prev_second_choice.u = best_i; prev_second_choice.v = best_j;
	}
      }
    }
   
    if (disjunctEdges.num() > 0) { //a max was found, and one end point isn't isolated in perg
      //(if it was, we wouldn't have set up a matching problem for it in the
      // most recent pass. Also, fixed edges won't be deleted because they're
      // not used in any solution <= 2B.
      DEBUGPR(300, /*if (uMPI::rank == uMPI::ioProc)*/ ucout << "EDGE_FOUND, branching on edge(s) " << disjunctEdges << std::endl
	      << "Deleted edges here: " << delEdges << std::endl);
      DEBUGPR(300, /*if (uMPI::rank == uMPI::ioProc)*/ ucout << "Fixed edges here: " << g_fixed << std::endl);
      return;
    } else { //if we got here, that means no edge had to be deleted in order
    // for the matching to work...this probably means the graph has a decent
    // amount of symmetry, so let's try to branch on something that will break
    // it in order to fathom quickly further down (for now keep implementation
    // the same as code on server)
      DEBUGPR(1, ucout << "Warning, no edge had to be deleted in order to make the matching work" << std::endl);
      for (i=0, gi = g.row(i); i<n; ++i, gi.nextRow()) {
        for (j=i; (j = gi.nextVertex(j)) >= 0; ) {
	  //don't delete fixed edges
	  if (g_fixed.isEdge(i,j)) {continue;}
	  //no point in deleting if this is true
	  if (perg.isIsolated(i) && perg.isIsolated(j)) {continue;}
	  // keep adding
	  disjunctEdges.add(i,j);
	  // until this is full
	  if (disjunctEdges.left() == 0) {
	    DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "EDGE_FOUND, branching on edge(s) " << disjunctEdges << std::endl
	      << "Deleted edges here: " << delEdges << std::endl);
	    DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "Fixed edges here: " << g_fixed << std::endl);
	    return;
	  }
        }
      }
      if (disjunctEdges.num() > 0) {
	DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "EDGE_FOUND, branching on edge(s) " << disjunctEdges << std::endl
	      << "Deleted edges here: " << delEdges << std::endl);
	DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "Fixed edges here: " << g_fixed << std::endl);
	return;
      }
    }
  } else { //implement random branching
    std::vector<int> randomIterator (n);   //make a vector for iterator
    std::vector<int> jRandomIterator;      //make vector for neighbors
    for(i=0; i<n; i++) randomIterator[i] = i;
    std::random_shuffle( randomIterator.begin(), randomIterator.end() );
    for (i=0; i<n; i++) {
      jRandomIterator.clear(); //clear out randomIterator[i]'s neighbors
      gi = g.row(randomIterator[i]);
      for (j=0; (j = gi.nextVertex(j)) >= 0;) {
	jRandomIterator.push_back(j);  //load neighbor vector
      }
      std::random_shuffle( jRandomIterator.begin(), jRandomIterator.end() ); //shuffle neighbor vector      
      for (std::vector<int>::iterator it=jRandomIterator.begin(); it!=jRandomIterator.end(); ++it) {
	  //don't delete fixed edges
	  if (g_fixed.isEdge(randomIterator[i],*it)) {continue;}
	  //no point in deleting if this is true
	  if (perg.isIsolated(randomIterator[i]) && perg.isIsolated(*it)) {continue;}
	  disjunctEdges.add(randomIterator[i], *it);
          if (disjunctEdges.left() == 0) {
	    DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "EDGE_FOUND, branching on edges " 
	  	    << disjunctEdges << std::endl
	  	    << "Deleted edges here: " << delEdges << std::endl);
	    DEBUGPR(300,/* if (uMPI::rank == uMPI::ioProc)*/ ucout << "Fixed edges here: " << g_fixed << std::endl);
	    return;
	  }
      }
    }
    if (disjunctEdges.num() > 0) return;
  }
  // if we got here either all the edges are fixed or isolated in perg, so this
  // node can die
  DEBUGPR(400, ucout << "We ran out of edges to fix - deleted edges here:" << delEdges << std::endl);
  DEBUGPR(400, ucout << "g_fixed here: " << g_fixed << std::endl);
  return;
}

void findAlmostSub::findAlmostSubAsChildOf(findAlmostSub* parent, int whichChild){
  //to keep left dive first, i==disjunctEdges.num() ~ fixed child, i<disjunctEdges.num() ~ deleted child
  
  branchSubAsChildOf(parent);
  bound = lowerBound = parent->lowerBound;
  globalPtr = parent->global();
  n = parent->n;
  B = parent->B;
  g = parent->g;
  delEdges = parent->delEdges;
  g_fixed = parent->g_fixed;
  perg = parent->perg;
  disjunctEdges = EdgeList(global()->DisjunctNum());
  boundCompFirstPass = true;  
  if (global()->trackEdges) prev_second_choice = parent->prev_second_choice;
  
  //**************************begin del edge disjunction*******************//
  if (whichChild < parent->disjunctEdges.num()) {
    delNode = true;
    B--; //decrement B
    g.delEdge(parent->disjunctEdges.getu(whichChild), parent->disjunctEdges.getv(whichChild));//delete edge from graph
    DEBUGPR(300, /*if (uMPI::rank == uMPI::ioProc)*/ ucout << "Deleting edge " << parent->disjunctEdges.get(whichChild) << std::endl);
    delEdges.add(parent->disjunctEdges.get(whichChild)); //add to list of deleted edges
  }
  //*************************end del edge disjunction********************//
  
  //********************Edge fixing child**********************//
  if( whichChild == parent->disjunctEdges.num() ){
    delNode = false;
    for (int k = 0; k < parent->disjunctEdges.num(); ++k)  
    	g_fixed.addEdge(parent->disjunctEdges.getu(k),parent->disjunctEdges.getv(k)); //add to graph of fixings
    DEBUGPR(300, /*if (uMPI::rank == uMPI::ioProc)*/ ucout << "Fixing edge(s)" << parent->disjunctEdges << std::endl);
  } 
  //***************************end Fix edge disjunction*********************//

 
}

/* end findAlmostSub methods */

/* begin parFindAlmost methods */
void parFindAlmostSub::runRefineByMatching( const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph &perg, int B, int* edge_use, bool always_collect_edge_use ) {  
  if (rampingUp()) {
    parRefineByMatching(g, g_fixed, perg, B, edge_use, always_collect_edge_use);
    //DEBUGPR(200,
    // if (uMPI::rank == uMPI::ioProc)  ucout << "During syncronous ramp-up, we refined " << perg_edges_bf - perg.numEdges() << std::endl);
  } else {
    refineByMatching(g, g_fixed, perg, B, edge_use, always_collect_edge_use);
  }
}

void parFindAlmostSub::makeBranchingDecisionsAgree() {
    if (global()->randomBranching && rampingUp()) {
      	// make the disjunctions all the same
	// we'll use ioProc's rank
	auto packed_edges = utilib::PackBuffer(disjunctEdges.packSize());
	disjunctEdges.write(packed_edges);
	uMPI::broadcast((char*)packed_edges.buf(),packed_edges.size(),MPI::CHAR,uMPI::ioProc); 
	auto unpacked_edges = utilib::UnPackBuffer(packed_edges.buf(), packed_edges.size());
	disjunctEdges.read(unpacked_edges);
    }
}

void parFindAlmostSub::pack(utilib::PackBuffer& outBuffer) {
	// We can do better...we don't need to pack g, instead let's get it from global()->gRef() and 
	// then delete the appropriate edges when we unpack...but we need to pack g's orbits
	outBuffer << delEdges << g_fixed << perg << delNode << boundCompFirstPass
		  << lowerBound << n << B << disjunctEdges;	
	for(int i=0; i<n; i++){
		outBuffer << g.orbits[i];
	}
	outBuffer << g.stats.numorbits;
}

void parFindAlmostSub::unpack(utilib::UnPackBuffer& inBuffer) {
	inBuffer >> delEdges >> g_fixed >> perg >> delNode >> boundCompFirstPass
		 >> lowerBound >> n >> B >> disjunctEdges;
	g = global()->gRef(); //get g from global branching object
	//fix so it's actually the current node
	for(int i=0; i<delEdges.num(); i++){
		g.delEdge(delEdges.getu(i),delEdges.getv(i));
	}
	for(int i=0; i<n; i++){
		inBuffer >> g.orbits[i];
	}
	inBuffer >> g.stats.numorbits;

}
int parFindAlmost::spPackSize(){
		//  g.orbits + g.stats.numorbits  +  delEdges      +  g_fixed  + perg
	return ( ((g.numVertices()+1) + 2*budget + 2)*sizeof(int) + 2*g.DenseGraph::packSize() 
		+ 3*sizeof(bool) + 5*sizeof(int) );
}

void parFindAlmost::pack(utilib::PackBuffer& outBuffer) {
	outBuffer << g << budget;
}

void parFindAlmost::unpack(utilib::UnPackBuffer& inBuffer) {
	inBuffer >> g >> budget;
}

parallelBranchSub* parFindAlmost::blankParallelSub() {
	parFindAlmostSub* newSP = new parFindAlmostSub;
	DEBUGPR(700, ucout << "new parFindAlmostSub from parFindAlmost::blankParallelSub(), address: " <<(void*) newSP << std::endl);
	newSP->setGlobalInfo(this);
	return newSP;
}

/* end parFindAlmost methodes */

/* being findAlmostSol methods */
void findAlmostSol::packContents(utilib::PackBuffer& outBuf){
	outBuf << delEdges << n;
	for(int i=0; i<n; i++) {
		outBuf << orbits[i];
	}
	outBuf << group_size;
}

void findAlmostSol::unpackContents(utilib::UnPackBuffer& inBuf) {
	inBuf >> delEdges >> n;
	orbits = new int[n];
	for(int i=0; i<n; i++) {
		inBuf >> orbits[i];
	}
	inBuf >> group_size;
}

int findAlmostSol::maxContentsBufSize() {
	return ( delEdges.packSize() + (n+1)*sizeof(int) + sizeof(double) );
}

/* writes the solution to the given ostream object */
void findAlmostSol::printContents(std::ostream& os){
  os << "Edges deleted: " << delEdges << std::endl;
  os << "Nontrival orbits:";

  int * thisOrbit = new int[n];
  int count;

  for( int i=0; i < n ; i++){
	  if (orbits[i] == i) { //this is the first vertex in this orbit
		count = 0;
		thisOrbit[count] = i;
		count++;
		for ( int j=i+1; j < n; j++) { //find the others
			if (orbits[j] == i){
				thisOrbit[count] = j;
				count++;
			}
		}
		if(count > 1) { //if there's more than a singleton we'll put to os
			for( int j=0; j<count; j++) os << " " << (thisOrbit[j]+1);
			//put number of orbits at end
			os << " " << "(" << count << ");";
		}
	  }
  }
  os << '\n';
  os << "Group size: " << group_size << std::endl;
  delete [] thisOrbit;
}
/* end findAlmostSol methods */

//************* begin refinement algorithms ******************//
// This function is from libhungarian by Cyrill Stachniss, with slight modifications for speed when repeatedly calling //
// libhungarian is available at: http://www2.informatik.uni-freiburg.de/~stachnis/misc.html
int hungarian_solve(hungarian_problem_t* p)
{
  int i, j, m, n, k, l, s, t, q, unmatched, cost;

  cost=0;
  m =p->num_rows;
  n =p->num_cols;

  // put these on the stack; they're small 
  int col_mate[m];
  int unchosen_row[m];
  int row_dec[m];
  int slack_row[m];

  int row_mate[n];
  int parent_row[n];
  int col_inc[n];
  int slack[n];

  /* //this should be done already
  for (i=0;i<p->num_rows;i++) {
    col_mate[i]=0;
    unchosen_row[i]=0;
    row_dec[i]=0;
    slack_row[i]=0;
  }
  for (j=0;j<p->num_cols;j++) {
    row_mate[j]=0;
    parent_row[j] = 0;
    col_inc[j]=0;
    slack[j]=0;
  }

  for (i=0;i<p->num_rows;++i)
    for (j=0;j<p->num_cols;++j)
      p->assignment[i][j]=HUNGARIAN_NOT_ASSIGNED; */

  // Begin subtract column minima in order to start with lots of zeroes 12
  for (l=0;l<n;l++)
    {
      s=p->cost[0][l];
      for (k=1;k<m;k++) 
	if (p->cost[k][l]<s)
	  s=p->cost[k][l];
      cost+=s;
      if (s!=0)
	for (k=0;k<m;k++)
	  p->cost[k][l]-=s;
    }
  // End subtract column minima in order to start with lots of zeroes 12

  // Begin initial state 16
  t=0;
   // replace this with calls to fill
  for (l=0;l<n;l++) row_mate[l] = -1;
  for (l=0;l<n;l++) parent_row[l] = -1;
  for (l=0;l<n;l++) col_inc[l] = 0;
  for (l=0;l<n;l++) slack[l] = INF;
  for (k=0;k<m;k++) slack_row[k] = 0;
  for (k=0;k<m;k++) unchosen_row[k] = 0;

  for (k=0;k<m;k++)
    {
      slack_row[k] = 0;
      s=p->cost[k][0];
      for (l=1;l<n;l++)
	if (p->cost[k][l]<s)
	  s=p->cost[k][l];
      row_dec[k]=s;
      for (l=0;l<n;l++)
	if (s==p->cost[k][l] && row_mate[l]<0)
	  {
	    col_mate[k]=l;
	    row_mate[l]=k;
	    goto row_done;
	  }
      col_mate[k]= -1;
      unchosen_row[t++]=k;
    row_done:
      ;
    }
  // End initial state 16
 
  // Begin Hungarian algorithm 18
  if (t==0)
    goto done;
  unmatched=t;
  while (1)
    {
      q=0;
      while (1)
	{
	  while (q<t)
	    {
	      // Begin explore node q of the forest 19
	      {
		k=unchosen_row[q];
		s=row_dec[k];
		for (l=0;l<n;l++)
		  if (slack[l])
		    {
		      int del;
		      del=p->cost[k][l]-s+col_inc[l];
		      if (del<slack[l])
			{
			  if (del==0)
			    {
			      if (row_mate[l]<0)
				goto breakthru;
			      slack[l]=0;
			      parent_row[l]=k;
			      unchosen_row[t++]=row_mate[l];
			    }
			  else
			    {
			      slack[l]=del;
			      slack_row[l]=k;
			    }
			}
		    }
	      }
	      // End explore node q of the forest 19
	      q++;
	    }
 
	  // Begin introduce a new zero into the matrix 21
	  s=INF;
	  for (l=0;l<n;l++)
	    if (slack[l] && slack[l]<s)
	      s=slack[l];
	  for (q=0;q<t;q++)
	    row_dec[unchosen_row[q]]+=s;
	  for (l=0;l<n;l++)
	    if (slack[l])
	      {
		slack[l]-=s;
		if (slack[l]==0)
		  {
		    // Begin look at a new zero 22
		    k=slack_row[l];
		    if (row_mate[l]<0)
		      {
			for (j=l+1;j<n;j++)
			  if (slack[j]==0)
			    col_inc[j]+=s;
			goto breakthru;
		      }
		    else
		      {
			parent_row[l]=k;
			unchosen_row[t++]=row_mate[l];
		      }
		    // End look at a new zero 22
		  }
	      }
	    else
	      col_inc[l]+=s;
	  // End introduce a new zero into the matrix 21
	}
    breakthru:
      // Begin update the matching 20
      while (1)
	{
	  j=col_mate[k];
	  col_mate[k]=l;
	  row_mate[l]=k;
	  if (j<0)
	    break;
	  k=parent_row[j];
	  l=j;
	}
      // End update the matching 20
      if (--unmatched==0)
	goto done;
      // Begin get ready for another stage 17
      t=0;
      for (l=0;l<n;l++)
	{
	  parent_row[l]= -1;
	  slack[l]=INF;
	}
      for (k=0;k<m;k++)
	if (col_mate[k]<0)
	  {
	    unchosen_row[t++]=k;
	  }
      // End get ready for another stage 17
    }
 done:

  // Begin doublecheck the solution 23
  /* don't bother please
  for (k=0;k<m;k++)
    for (l=0;l<n;l++)
      if (p->cost[k][l]<row_dec[k]-col_inc[l])
	exit(0);
  for (k=0;k<m;k++)
    {
      l=col_mate[k];
      if (l<0 || p->cost[k][l]!=row_dec[k]-col_inc[l])
	exit(0);
    }
  k=0;
  for (l=0;l<n;l++)
    if (col_inc[l])
      k++;
  if (k>m)
    exit(0);
  // End doublecheck the solution 23
  */
  // End Hungarian algorithm 18

  for (i=0;i<m;++i)
    {
      p->assignment[i][col_mate[i]]=HUNGARIAN_ASSIGNED;
      /*TRACE("%d - %d\n", i, col_mate[i]);*/
    }
  for (k=0;k<m;++k)
    {
      for (l=0;l<n;++l)
	{
	  /*TRACE("%d ",p->cost[k][l]-row_dec[k]+col_inc[l]);*/
	  p->cost[k][l]=p->cost[k][l]-row_dec[k]+col_inc[l];
	}
      /*TRACE("\n");*/
    }
  for (i=0;i<m;i++)
    cost+=row_dec[i];
  for (i=0;i<n;i++)
    cost-=col_inc[i];

/* // don't need these for stack-alloc'd arrays
  free(slack);
  free(col_inc);
  free(parent_row);
  free(row_mate);
  free(slack_row);
  free(row_dec);
  free(unchosen_row);
  free(col_mate);
  */
  return cost;
}

 
void refineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph &perg, int B, int* edge_use, bool always_collect_edge_use){
  int i, inf = 4*B, n = g.numVertices();
  //begin block for parallel processing
  {
  int j,k,l,deg_gi,deg_gj,size,u,v;
  bool neighbors;
  bool collect_edge_use = true;
  DenseGraph::Row pergi;
  int biSize = n+B;
  int* Ni; 
  int* Nj;
  int** biGraph;
  int** assignment;
  hungarian_problem_t p;
  int cost;
  biGraph = (int**)malloc(biSize*sizeof(int*));
  assignment = (int**)malloc(biSize*sizeof(int*));
  Ni = (int*)malloc(n*sizeof(int));
  Nj = (int*)malloc(n*sizeof(int));
  for (j = 0; j < biSize; ++j){ 
    biGraph[j] = (int*)malloc(biSize*sizeof(int));
    assignment[j] = (int*)malloc(biSize*sizeof(int));
  } 
  for (i = 0; i < n; i++){
    pergi = perg.row(i);
    for (j = i; (j = pergi.nextVertex(j)) >= 0; ){ //check upper triangle only
      deg_gi = g.deg(i); deg_gj = g.deg(j);
      //set up and solve matching problem
      neighbors = false;
      size = (deg_gi > deg_gj) ? deg_gi : deg_gj ; //return greater of deg_gi and deg_gj
      //printf("B = %d", B);
      size += B;
      neighbors = g.isEdge(i,j);
      if (neighbors) size--; //if i and j are neighbors, we won't add them to the biGraph
      
      buildCostMatrix(i, j, g, g_fixed, perg, size, inf, neighbors, biGraph, Ni, Nj);
      
      for(k = 0; k < size; k++){
	memset(assignment[k], 0, size*sizeof(int));
      } // zero out assignment matrix
	//let's do the initializtion ourselves

      p.num_rows = size; p.num_cols = size;
      p.cost = biGraph; p.assignment = assignment;
	
      cost = hungarian_solve(&p);
	
      // This information needs to be shared between the treads somehow...
      if (cost > 2*B){  //we can kill i,j
	//printf("Cost = %d, killing this permutation\n", cost);
	//printf("Kill permutation (%d, %d) in refineByMatching\n", i, j);
	perg.delEdge(i,j);
	//if we delete an edge from perg we won't be using the edge_use data,
	//so we can save time by not writing it
	if (!always_collect_edge_use) collect_edge_use = false;
      } else if(collect_edge_use) { 
	if (neighbors) {
	  for (k = 0; k < deg_gi-1; ++k){
	    for (l = deg_gj-1; l < size; ++l){
	      if (assignment[k][l] == 1){
	  	u = Ni[k];
		//printf("edge %d, %d got deleted\n", u, i);
		if (i>u) {
		  edge_use[LTIdx(i,u)]++;
		} else { // u>i
		  edge_use[LTIdx(u,i)]++;
		}  
		break; //if we come across 1 in the row, we shouldn't find others
	      }
	    }
	  }
	  for ( ; k < size; ++k){
	    for (l = 0; l < deg_gj-1; ++l){
	      if (assignment[k][l] == 1){
		v = Nj[l];
		//printf("edge %d, %d got deleted\n", v, j);
		if (j>v) {
		  edge_use[LTIdx(j,v)]++;
		} else { //v>j
	 	  edge_use[LTIdx(v,j)]++;
		}
		break; //sim
	      }
	    }
	  }
	} else {
	  for (k = 0; k < deg_gi; ++k){
	    for (l = deg_gj; l < size; ++l){
	      if (assignment[k][l] == 1){
		u = Ni[k];
		//printf("edge %d, %d got deleted\n", u, i);
		if (i>u) {
		  edge_use[LTIdx(i,u)]++;
		} else { // u>i
		  edge_use[LTIdx(u,i)]++;
		}
		break;
	      }	
	    }
	  }
	  for ( ; k < size; ++k){
	    for (l = 0; l < deg_gj; ++l){
	      if (assignment[k][l] == 1){
		v = Nj[l];
		//printf("edge %d, %d got deleted\n", v, j);
		if (j>v) {
		  edge_use[LTIdx(j,v)]++;
		} else { //v>j
	 	  edge_use[LTIdx(v,j)]++;
		}
		break;
	      }
	    }
	  }
	}
	//printf("\n");
      } //end info needing to be shared        
    }
  }
  for (j = 0; j < biSize; ++j){ 
    free(biGraph[j]);
    free(assignment[j]);
  }
  free(Ni); free(Nj);
  free(biGraph);
  free(assignment);
  } //end parallel scope
}

/* TODO: fix logic to handle self-loops in g */
void parRefineByMatching(const DenseGraph &g, const DenseGraph &g_fixed, DenseGraph &perg, int B, int* edge_use, bool always_collect_edge_use ){
  int i, inf = 4*B, n = g.numVertices();
  //begin block for parallel processing
  {
  int j,k,l,deg_gi,deg_gj,size,u,v;
  bool neighbors;
  bool collect_edge_use = true;
  DenseGraph::Row pergi;
  int biSize = n+B;
  int* Ni; 
  int* Nj;
  int** biGraph;
  int** assignment;
  hungarian_problem_t p;
  int cost;
  biGraph = (int**)malloc(biSize*sizeof(int*));
  assignment = (int**)malloc(biSize*sizeof(int*));
  Ni = (int*)malloc(n*sizeof(int));
  Nj = (int*)malloc(n*sizeof(int));
  for (j = 0; j < biSize; ++j){ 
    biGraph[j] = (int*)malloc(biSize*sizeof(int));
    assignment[j] = (int*)malloc(biSize*sizeof(int));
  } 
  //begin figuring out which edges in perg we'll process
  int num_iterations = perg.numEdges(); //- n; //we'll only iterate over perg's 
  int iteration = 0;
  bool * perg_del_edges = new bool[num_iterations](); //itilize to 0
  int iters_per_proc = num_iterations / (uMPI::size);
  int iters_remaining = num_iterations%(uMPI::size);
  int my_start_iter = uMPI::rank*iters_per_proc + ((uMPI::rank < iters_remaining) ? (uMPI::rank) : iters_remaining); //return the smaller of rank and iters remaining
  int my_end_iter = (uMPI::rank+1)*iters_per_proc + (((uMPI::rank+1) < iters_remaining) ? (uMPI::rank+1) : iters_remaining);
  //end 
  for (i = 0; i < n; i++){
    pergi = perg.row(i);
    for (j = i; (j = pergi.nextVertex(j)) >= 0; iteration++){ //check upper triangle only
     if( iteration >= my_start_iter && iteration < my_end_iter) { 
      //    printf("\n(%d, %d) is an edge in the possible permutations graph\n", i, j);
      deg_gi = g.deg(i); deg_gj = g.deg(j);
      //set up and solve matching problem
      neighbors = false;
      size = (deg_gi > deg_gj) ? deg_gi : deg_gj ; //return greater of deg_gi and deg_gj
      //printf("B = %d", B);
      size += B;
      neighbors = g.isEdge(i,j);
      if (neighbors) size--; //if i and j are neighbors, we won't add them to the biGraph
      
      buildCostMatrix(i, j, g, g_fixed, perg, size, inf, neighbors, biGraph, Ni, Nj);
      
      for(k = 0; k < size; k++){
	memset(assignment[k], 0, size*sizeof(int));
      } // zero out assignment matrix
	//let's do the initializtion ourselves

      p.num_rows = size; p.num_cols = size;
      p.cost = biGraph; p.assignment = assignment;
	
      cost = hungarian_solve(&p);
	
      // This information needs to be shared between the treads somehow...
      if (cost > 2*B){  //we can kill i,j
	//printf("Cost = %d, killing this permutation\n", cost);
	//if (uMPI::rank == uMPI::ioProc)  ucout << "Planning to kill permutation (" << i << "," << j << ") in refineByMatching, counter = " << counter << std::endl;
	//perg.delEdge(i,j);
	perg_del_edges[iteration] = true;
	//if we delete an edge from perg we won't be using the edge_use data,
	//so we can save time by not writing it
	if (!always_collect_edge_use)  collect_edge_use = false;
      } else if(collect_edge_use) { 
	if (neighbors) {
	  for (k = 0; k < deg_gi-1; ++k){
	    for (l = deg_gj-1; l < size; ++l){
	      if (assignment[k][l] == 1){
	  	u = Ni[k];
		//printf("edge %d, %d got deleted\n", u, i);
		if (i>u) {
		  edge_use[LTIdx(i,u)]++;
		} else { // u>i
		  edge_use[LTIdx(u,i)]++;
		}  
		break; //if we come across 1 in the row, we shouldn't find others
	      }
	    }
	  }
	  for ( ; k < size; ++k){
	    for (l = 0; l < deg_gj-1; ++l){
	      if (assignment[k][l] == 1){
		v = Nj[l];
		//printf("edge %d, %d got deleted\n", v, j);
		if (j>v) {
		  edge_use[LTIdx(j,v)]++;
		} else { //v>j
	 	  edge_use[LTIdx(v,j)]++;
		}
		break; //sim
	      }
	    }
	  }
	} else {
	  for (k = 0; k < deg_gi; ++k){
	    for (l = deg_gj; l < size; ++l){
	      if (assignment[k][l] == 1){
		u = Ni[k];
		//printf("edge %d, %d got deleted\n", u, i);
		if (i>u) {
		  edge_use[LTIdx(i,u)]++;
		} else { // u>i
		  edge_use[LTIdx(u,i)]++;
		}
		break;
	      }	
	    }
	  }
	  for ( ; k < size; ++k){
	    for (l = 0; l < deg_gj; ++l){
	      if (assignment[k][l] == 1){
		v = Nj[l];
		//printf("edge %d, %d got deleted\n", v, j);
		if (j>v) {
		  edge_use[LTIdx(j,v)]++;
		} else { //v>j
	 	  edge_use[LTIdx(v,j)]++;
		}
		break;
	      }
	    }
	  }
	} 
	//printf("\n");
      } //end info needing to be shared   
     }     
    }
  } 
  //if (num_iterations != iteration) ucout << "You're a dumbass, interation = " << iteration << " and num_iterations = " << num_iterations << std::endl;
  //collect info and act on it
  //set collect_edge_use ... if it's false on some thread it's false
  //everywhere, and we'll need to gather up the deleted edges
  uMPI::reduceCast(MPI_IN_PLACE,&collect_edge_use,1,MPI::BOOL,MPI_LAND);
  if(collect_edge_use){ //in this case, no edges were deleted, so we'll probably be branching
    uMPI::reduceCast(MPI_IN_PLACE,edge_use,LTIdx(n,0),MPI::INT,MPI_SUM);
  } else { //in this case some edges were deleted, so we need to update perg.
    uMPI::reduceCast(MPI_IN_PLACE,perg_del_edges,num_iterations,MPI::BOOL,MPI_LOR);
    iteration = 0;
    for (i = 0; i < n; i++){
      pergi = perg.row(i);
      for (j = i; (j = pergi.nextVertex(j)) >= 0; iteration++){ //check upper triangle only
        if(perg_del_edges[iteration]){
          perg.delEdge(i,j);        
	  //if (uMPI::rank == uMPI::ioProc)  ucout << "Killing permutation (" << i << "," << j << ") in refineByMatching" << std::endl;
        }
      }
    }
  }  
  for (j = 0; j < biSize; ++j){ 
    free(biGraph[j]);
    free(assignment[j]);
  }
  free(Ni); free(Nj);
  free(biGraph);
  free(assignment);

  } //end parallel scope
}


//TODO; This logic of this function should be made to handle the case when g has
//      self-loops
void buildCostMatrix(int i, int j,const DenseGraph & g,const DenseGraph & g_fixed,const DenseGraph & perg, int size, int inf, bool neighbors, int** biGraph, int* Ni, int* Nj){
  int k,l,u,v,deg_gu;
  bool jNeig = true;
  bool * Ni_not_indep; 
  bool * Nj_not_indep;
  DenseGraph::Row gi = g.row(i), gj = g.row(j);
  if (neighbors) {
    //initialize these to false (assume each vertex is independent...we'll exhastively prove this
    // is the case)	  
    Ni_not_indep = new bool[g.deg(i) - 1](); Nj_not_indep = new bool[g.deg(j) - 1]();    
    for (u = -1, k = 0; (u = gi.nextVertex(u)) >= 0; ){ //++k){ //for each neighbor u of i
      if (u != j) { //only add row if neighbor u of i is not j (we can avoid these check if take them out)
	Ni[k] = u; //make this the next element in neighbor list 
	deg_gu = g.deg(u);
	for (v = -1, l = 0; (v = gj.nextVertex(v)) >= 0;){ // ++l ){ //for each neighbor v of j
	   //make next element of neighbor list
	  if (v != i) { //only add column of neighbor v of j is not i (same here)
	    if ( u == v ) {
              biGraph[k][l] = 0;
	    } else { 
	      if ( perg.isEdge(u,v)) { // if v is permutable with u
	        biGraph[k][l] = abs(deg_gu - g.deg(v));
	      } else {
	        biGraph[k][l] = inf; //make this large enough that it's basically +/infty
	      }
	      if (g.isEdge(u,v)) {
	        Ni_not_indep[k] = true;
		Nj_not_indep[l] = true;
	      }
	    }  
	    l++; //only increment l when they're not neighbors
	  }
	}
	for ( ; l < size; ++l) {
	  if (g_fixed.isEdge(i,u)){
	    biGraph[k][l] = inf; //if edge is fixed, it can't be deleted
	  } else {
	    biGraph[k][l] = 2;
	  }
	}
	for (v = u; ((v = gi.nextVertex(v)) >= 0) && !(Ni_not_indep[k]);) { //check each pair only once
		// if Ni_not_indep[k] is already true don't bother to check this and exit when it becomes true
	  if ( (v != j) && g.isEdge(u,v) ) Ni_not_indep[k] = true;
	}
	k++; //only increment k when they're not neighbors
      }
    }
    for ( ; k < size; ++k) {
      for (v = -1, l = 0; (v = gj.nextVertex(v)) >= 0;){// ++l ){ //for each neighbor v of j
	if (v != i) { //only add row of neighbor v of j is not i
	  if (jNeig) Nj[l] = v;
	  if (g_fixed.isEdge(j,v)){
	    biGraph[k][l] = inf; //if edge is fixed, it can't be deleted
	  } else {
	    biGraph[k][l] = 2; //cost of deletion
	  }
	  if (jNeig) { //only do this one pass though the loop
            for (u = v; ((u = gj.nextVertex(u)) >= 0) && !(Nj_not_indep[l]);) {
		    // if Nj_not_indep[l] is already true don't bother to check and exit when it becomes true
	      if( (u != i) && g.isEdge(u,v) ) Nj_not_indep[l] = true;
	    }
          }	    
	  l++; //only increment l when they're not neighbors
	}
      }
      for ( ; l < size; ++l) {
	biGraph[k][l] = 0;
      }
      jNeig = false;
    }
    // we've collected the data, now let's use it!!
    // TODO: This logic needs to be fixed to handle when g has self-loops
    for( k = 0; k < (g.deg(i) - 1); k++ ) {
      if( !(Ni_not_indep[k]) ) { //k is independent from every vertex in biGraph 
        for( l = 0; l < (g.deg(j) - 1); l++) {
          if( g.deg(Ni[k]) > g.deg(Nj[l]) ) { //if deg(Ni[k])>deg(Nj[l]), then we need delete at least the difference of the edges
	    //if the larger degree one is independent, then an edge deletion
	    //on it will only be counted once!
	    biGraph[k][l] *= 2;
	  }
        }
      }	
    }
    for( l = 0; l < (g.deg(j) - 1); l++ ) {
      if( !(Nj_not_indep[l]) ) { //l is independent from every vertex in biGraph 
        for( k = 0; k < (g.deg(i) - 1); k++) {
          if( g.deg(Nj[l]) > g.deg(Ni[k]) ) { //if deg(Nj[l])>deg(Ni[k]), then we need delete at least the difference of the edges
	    //if the larger degree one is independing, then an edge deletion
	    //on it will only be counted once!
	    biGraph[k][l] *= 2;
	  }
        }
      }	
    } 
    delete [] Ni_not_indep; delete [] Nj_not_indep;      
  } else { //they're not neighbors, so we can do things a little faster
    //initialize these to false (assume each vertex is independent...we'll exhastively prove this
    // is the case)	  
    Ni_not_indep = new bool[g.deg(i)](); Nj_not_indep = new bool[g.deg(j)]();    
    for (u = -1, k = 0; (u = gi.nextVertex(u)) >= 0; ++k ){ //for each neighbor u of i
      Ni[k] = u;
      deg_gu = g.deg(u);
      for (v = -1, l = 0; (v = gj.nextVertex(v)) >= 0; ++l ){ //for each neighbor v of j
	if ( u == v ) {
	  biGraph[k][l] = 0;	
	} else {
	  if (perg.isEdge(u,v)) { // if v is permutable with u
	    biGraph[k][l] = abs(deg_gu - g.deg(v));
	  } else {
	    biGraph[k][l] = inf; //make this large enough that it's basically +/infty
	  }
	  if (g.isEdge(u,v)) {
	    Ni_not_indep[k] = true;
	    Nj_not_indep[l] = true;
	  }
        }	  
      }
      for ( ; l < size; ++l) {
	if (g_fixed.isEdge(i,u)){
	  biGraph[k][l] = inf; //if edge is fixed, it can't be deleted
	} else {
	  biGraph[k][l] = 2;
	}
      }
      for(v = u; ((v = gi.nextVertex(v)) >= 0) && !(Ni_not_indep[k]);) { //check each pair only once
	if( g.isEdge(u,v) ) Ni_not_indep[k] = true;
      }	
    }
    for ( ; k < size; ++k) {
      for (v = -1, l = 0; (v = gj.nextVertex(v)) >= 0; ++l ){ //for each neighbor v of j
	if (jNeig) Nj[l] = v; //make next element of neighbor list
	if (g_fixed.isEdge(j,v)){
	  biGraph[k][l] = inf; //if edge is fixed, it can't be deleted
	} else {
	  biGraph[k][l] = 2; //cost of deletion
	}
	if (jNeig) {
	  for (u = v; ((u = gj.nextVertex(u)) >= 0) && !(Nj_not_indep[l]);) {
	    if( g.isEdge(u,v) ) Nj_not_indep[l] = true;
	  }	  
	}
      }
      for ( ; l < size; ++l) {
	biGraph[k][l] = 0;
      }
      jNeig = false;
    }
    // we've collected the data, now let's use it!
    // TODO: This logic need to be fixed to handle when g has self-loops
    for( k = 0; k < g.deg(i); k++ ) {
      if( !(Ni_not_indep[k]) ) { //k is independent from every vertex in biGraph 
        for( l = 0; l < g.deg(j); l++) {
          if( g.deg(Ni[k]) > g.deg(Nj[l]) ) { //if deg(Ni[k])>deg(Nj[l]), then we need delete at least the difference of the edges
	    //if the larger degree one is independent, then an edge deletion
	    //on it will only be counted once!
	    biGraph[k][l] *= 2;
	  }
        }
      }	
    }
    for( l = 0; l < g.deg(j); l++ ) {
      if( !(Nj_not_indep[l]) ) { //l is independent from every vertex in biGraph 
        for( k = 0; k < g.deg(i); k++) {
          if( g.deg(Nj[l]) > g.deg(Ni[k]) ) { //if deg(Nj[l])>deg(Ni[k]), then we need delete at least the difference of the edges
	    //if the larger degree one is independing, then an edge deletion
	    //on it will only be counted once!
	    biGraph[k][l] *= 2;
	  }
        }
      }	
    }  
    delete [] Ni_not_indep; delete [] Nj_not_indep;      
  }
}

void initialRefineByDegrees(const DenseGraph &g, DenseGraph &perg, int B){
  int i,j;
  for(i = 0; i < g.numVertices(); i++){
    for(j = 0; j < i; j++){
      if (abs(g.deg(i) - g.deg(j)) > B) {
	perg.delEdge(i,j);
      }
    }
  }
}


void refineByDegreeDiffFixedEdges(const DenseGraph & g, const DenseGraph & g_fixed, DenseGraph & perg){
  int i,j,n = g.numVertices();
DenseGraph::Row pergi;
  for (i = 0, pergi = perg.row(i) ; i < n; ++i, pergi.nextRow()){
    for (j = i; (j = pergi.nextVertex(j)) >= 0; ){
      if ((g_fixed.deg(i) > g.deg(j)) ||
	  (g_fixed.deg(j) > g.deg(i)))
      {
	//printf("Killing permutation (%d, %d) in refineByDegreeDiffFixedEdges\n", i, j);
	perg.delEdge(i,j);
      }
    }
  }
}

void refineByDegreeDiff(const DenseGraph & g, DenseGraph & perg, int B){
int i,j,n = g.numVertices();
DenseGraph::Row pergi;
  for (i = 0, pergi = perg.row(i) ; i < n; ++i, pergi.nextRow()){
    for (j = i; (j = pergi.nextVertex(j)) >= 0; ){
      if (abs(g.deg(i) - g.deg(j)) > B){
	//printf("Killing permutation (%d, %d) in refineByDegreeDiff\n", i , j);
	  perg.delEdge(i,j);
      }
    }
  }
}

//******************end refinement algorithms******************************//
//******utils******///

int independentSetHeuristic(const DenseGraph & h){
  const int n = h.numVertices();

  //copy the degree list
  std::vector<int> degree(n);
  for(int i=0; i<n; i++) degree[i] = h.deg(i);

  int independenceNum = 0;
  //make a list of uncovered vertices, initialize to everything
  std::vector<int> UncoveredVertices;
  // set asside some space
  UncoveredVertices.reserve(n/4);
  //preprocess these degree 0 vertices
  for(int i=0; i<n; i++) {
    if (degree[i] == 0) independenceNum++;
    else UncoveredVertices.push_back(i);
  }

  while (!UncoveredVertices.empty()) {
    //find an uncovered minimal degree vertex        
    int min_degree = h.numEdges() + 1;
    int cover_vertex = -1; 
    // iterate over uncovered vertices
    for(int v : UncoveredVertices) {
      if (degree[v] < min_degree) {
	min_degree = degree[v];
	cover_vertex = v;
      }
    }
    // add to cover
    {
      std::vector<int>::iterator it;
      for( it = UncoveredVertices.begin(); it!=UncoveredVertices.end(); ++it) {
	if ( cover_vertex == *it ) break;
      }	
      UncoveredVertices.erase(it);
    }
    DenseGraph::Row g_cover = h.row(cover_vertex);
    DenseGraph::Row gj;
    // "remove" the edges in and adjacent to the cover 
    for( int j = -1; (j = g_cover.nextVertex(j)) >= 0; ){   
      // add to cover
      {
        std::vector<int>::iterator it;
        for( it = UncoveredVertices.begin(); it!=UncoveredVertices.end(); ++it) {
          if ( j == *it ) {
	    UncoveredVertices.erase(it);
	    break;
	  }
	  if ( j < *it ) break;
        }	
      }
      // "erase" this vertex's neighbors, including the cover vertex
      gj = h.row(j);
      for (int k = -1; (k = gj.nextVertex(k)) >= 0;){
	degree[k]--;
      }
    }
    // iterate independence number
    independenceNum++;
  }

  return independenceNum;

}

void printArray(int n, int m, int** array){
  int i,j;
  for (i = 0; i < n; ++i){
    printf("[");
    for (j = 0; j < m; ++j){
      printf(" %d ", array[i][j]);
    }
    printf("]\n");
  }
}

