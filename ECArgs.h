
#ifndef ECARGS_H
#define ECARGS_H

#include <list>
#include "ECString.h"

using namespace std;

class ECArgs
{
 public:
  ECArgs(int argc, char *argv[]);
  int nargs() { return nargs_; }
  bool isset(char c);
  ECString value(char c);
  ECString arg(int n) { return argList[n]; }
 private:
  int nargs_;
  int nopts_;
  ECString argList[32];
  list<ECString> optList;
};
  

#endif /* ! ECARGS_H */
