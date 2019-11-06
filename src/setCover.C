// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2012 Guy Blelloch, Kanat Tangwongsan and the PBBS team
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

// This code is based on the paper:
//    Guy Blelloch, Richard Peng and Kanat Tangwongsan. 
//    Linear-Work Greedy Parallel Approximation Algorithms for Set Covering and Variants. 
//    Proc. ACM Symposium on Parallel Algorithms and Architectures (SPAA), May 2011
// It is a serial version of that code

#include <iostream>
#include <math.h>
#include "sequence.h"
#include "gettime.h"
#include "graph.h"
#include "parallel.h"
#include "blockRadixSort.h"
using namespace std;

timer bucketTime;
timer manisTime;
timer packTime;

typedef int intT;
typedef unsigned int uint;
typedef pair<uint, uint> upair;
typedef graph<intT> Graph;

struct set {
  intT* Elements;
  intT degree;
  intT id;
  set(intT* E, intT d, intT i) : Elements(E), degree(d), id(i) {}
};

struct bucket {
  set *S;
  intT n;
  bucket(set* _S, intT _n) : S(_S), n(_n) {}
};

// maximum element ID in any set + 1
intT maxElt(graph<intT> G) {
  intT *C = newA(intT,G.n);
  for (intT i = 0; i < G.n; i++) {
    C[i] = 0;
    for (intT j=0; j < G.V[i].degree; j++) 
      C[i] = max(C[i], G.V[i].Neighbors[j]);
  }
  return sequence::reduce(C, G.n, utils::maxF<intT>()) + 1;
}

// Puts each vertex into buckets that are multiples of (1+epsilon) based on degree
std::pair<bucket*,intT> putInBuckets(Graph GS, double epsilon) {
  double x = 1.0/log(1.0 + epsilon);
  intT numBuckets = 1 + floor(x * log((double) GS.n));
  cout << "numBuckets = " << numBuckets << endl;
  intT* offsets = newA(intT, numBuckets + 256);

  // determine bucket numbers
  upair *A = newA(upair, GS.n);
  for(intT i=0; i < GS.n; i++) {
    // deal with empty buckets
    intT d = (GS.V[i].degree > 0) ? GS.V[i].degree : 1; 
    A[i] = upair(floor(x * log((double) d)), i);
  }

  // sort based on bucket numbers
  intSort::iSort(A, offsets, (intT) GS.n, numBuckets, utils::firstF<uint, uint>());
  set *S = newA(set, GS.n);
  for(int i=0; i < GS.n; i++) {
    int j = A[i].second;
    S[i] = set(GS.V[j].Neighbors, GS.V[j].degree, j);
  }
  free(A);

  // create each bucket
  bucket *B = newA(bucket, numBuckets);
  for(int i=0; i < numBuckets-1; i++) 
    B[i] = bucket(S+offsets[i],offsets[i+1]-offsets[i]);
  B[numBuckets-1] = bucket(S+offsets[numBuckets-1], GS.n - offsets[numBuckets-1]);
  free(offsets);

  return pair<bucket*,int>(B, numBuckets);
}

void freeBuckets(bucket* B) {
  free(B[0].S);
  free(B);
}

intT processBucket (set* S, intT* elts, intT n, intT threshold) {
  for (intT i = 0; i < n; i++) {
    if (S[i].degree >= threshold) {
      int k = 0;
      for (intT j = 0; j < S[i].degree; j++) {
	intT ngh = S[i].Elements[j];
	if (elts[ngh] == INT_MAX) S[i].Elements[k++] = ngh;
      }
      S[i].degree = k;
      if (k >= threshold) {
	for (intT j = 0; j < k; j++) 
	  elts[S[i].Elements[j]] = -1;
	S[i].degree = -1; // marks that vertex is in set
      }
    }
  }
  return n;
}

_seq<intT> setCover(Graph GS) {
  double epsilon = 0.01;
  intT m = maxElt(GS);
  cout << "m = " << m << endl;

  bucketTime.start();
  pair<bucket*, int> B = putInBuckets(GS, epsilon);
  bucketTime.stop();
  bucket* allBuckets = B.first;
  int numBuckets = B.second;

  set* S = newA(set, GS.n);    // holds sets for current bucket
  set* ST = newA(set, GS.n);   // temporarily S (pack is not inplace)
  int l = 0;                   // size of S
  bool* flag = newA(bool, GS.n);
  intT* inCover = newA(intT, GS.n);
  intT nInCover = 0;
  intT totalWork = 0;
  intT* elts = newA(intT,m);
  intT threshold = GS.n;
  for (int i = 0; i < m; i++) elts[i] = INT_MAX;

  // loop over all buckets, largest degree first
  for (int i = numBuckets-1; i >= 0; i--) {
    bucket currentB = allBuckets[i];

    intT degreeThreshold = ceil(pow(1.0+epsilon,i));
    if (degreeThreshold == threshold && currentB.n == 0) continue;
    else threshold = degreeThreshold;
    packTime.start();

    // pack leftover sets that are below threshold down for the next round
    for (int j = 0; j < l; j++) 
      flag[j] = (S[j].degree > 0 && S[j].degree < threshold);
    intT ln = sequence::pack(S, ST, flag, l);

    // pack leftover sets greater than threshold above for this round
    for (int j = 0; j < l; j++) 
      flag[j] = (S[j].degree >= threshold);
    intT lb = sequence::pack(S, ST+ln, flag, l);

    // copy prebucketed bucket i to end, also for this round
    for (int j = 0; j < currentB.n; j++) 
      ST[j+ln+lb] = currentB.S[j];

    lb = lb + currentB.n;   // total number in this round
    l = ln + lb;            // total number including those for next round
    swap(ST,S);             // since pack is not in place 
    set* SB = S + ln;       // pointer to bottom of sets for this round
    packTime.stop();

    if (lb > 0) { // is there anything to do in this round?

      manisTime.start();
      intT work = processBucket(SB, elts, lb, threshold);
      totalWork += work;
      manisTime.stop();
      packTime.start();

      // check which sets were selected by manis to be in the set cover
      for (int j = 0; j < lb; j++)
	flag[j] = SB[j].degree < 0;

      // add these to inCover and label by their original ID
      int nNew = sequence::packIndex(inCover+nInCover, flag, lb);
      for (int j = nInCover; j < nInCover + nNew; j++) 
	inCover[j] = SB[inCover[j]].id;
      nInCover = nInCover + nNew;
      packTime.stop();
      // cout << "i = " << i << " bc = " << currentB.n << " l = " << l << " lb = " << lb
	  //  << " work = " << work << " new = " << nNew << " threshold = " << threshold << endl;
    }
  }
  cout << "Set cover size = " << nInCover << endl;
  // cout << "Total work = " << totalWork << endl;
  // cout << "Bucket Time = " << bucketTime.total() << endl;
  cout << "Manis Time = " << manisTime.total() << endl;
  // cout << "Pack Time = " << packTime.total() << endl;

  free(elts); free(S); free(ST); free(flag);
  freeBuckets(allBuckets);
  return _seq<intT>(inCover, nInCover);
}
