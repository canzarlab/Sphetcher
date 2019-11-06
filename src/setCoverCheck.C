// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include <cstring>
#include "parallel.h"
#include "IO.h"
#include "graph.h"
#include "graphIO.h"
#include "parseCommandLine.h"
using namespace std;
using namespace benchIO;

intT maxElt(graph<intT> G) {
  intT *C = newA(intT,G.n);
  parallel_for (intT i = 0; i < G.n; i++) {
    C[i] = 0;
    for (intT j=0; j < G.V[i].degree; j++) 
      C[i] = max(C[i], G.V[i].Neighbors[j]);
  }
  return sequence::reduce(C, G.n, utils::maxF<intT>()) + 1;
}

int checkSetCover(graph<intT> G, _seq<intT> Sets) {
  intT n = G.n;
  intT m = maxElt(G);
  bool *Elts = newA(bool, m);
  parallel_for (intT i = 0; i < n; i++) Elts[i] = 0;

  parallel_for (intT i = 0; i < Sets.n; i++) {
    intT k = Sets.A[i];
    for (intT j=0; j < G.V[k].degree; j++) 
      Elts[G.V[k].Neighbors[j]] = 1;
  }

  for (intT i=0; i < n; i++)
    for (intT j=0; j < G.V[i].degree; j++) {
      if (!Elts[G.V[i].Neighbors[j]]) {
	  cout << "checkSetCover: uncovered element " 
	       << G.V[i].Neighbors[j] << " from set " << i << endl;
	  return 1;
      }
    }

  return 0;
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  graph<intT> G = readGraphFromFile<intT>(iFile);
  _seq<intT> Out = readIntArrayFromFile<intT>(oFile);
  if (Out.n > G.n) {
    cout << "set Cover: too many sets" << endl;
    return(1);
  }

  return checkSetCover(G, Out);
}
