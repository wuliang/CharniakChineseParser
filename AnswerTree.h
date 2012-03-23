
#ifndef ANSWERTREE_H
#define ANSWERTREE_H

#include <list>
#include "Edge.h"
#include "CntxArray.h"
using namespace std;

class AnswerTree;
class
AnswerTree
{
public:
  AnswerTree(Edge* ed) : e(ed){}
  void extend(AnswerTree* at) { subtrees.push_back(at); }
  void deleteSubTrees();
  Edge* e;
  list<AnswerTree*> subtrees;
};

typedef pair<double,AnswerTree*> AnswerTreePair;
typedef map<CntxArray, AnswerTreePair, less<CntxArray> > AnswerTreeMap;

AnswerTreePair& atpFind(CntxArray& hi, AnswerTreeMap& atm);

#endif
