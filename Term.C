
#include "Term.h"
#include "utils.h"

// wul: 这里的400具有疑惑性。因为句子中的单词数目限制也是400。实际上Term的数目没有这么多。
// wul: 现在Term的总数（包括终结符和非终结符）小于80
Term*  Term::array_[400];
map<ECString, Term*,less<ECString> > Term::termMap_;
int    Term::lastTagInt_ = 0;
int    Term::lastNTInt_ = 0;
Term*  Term::stopTerm;
Term*  Term::startTerm;
Term*  Term::rootTerm;

Term::
Term()
: num_( -1 ), terminal_p_( 0 )
{}

Term::
Term(const ECString s, int terminal, int num )
: name_( s ), num_( num ), terminal_p_( terminal )
{}

Term::
Term( const Term& src )
: name_( src.name() ), 
  //num_(src.toInt()),
  terminal_p_( src.terminal_p() )
{}

ostream&
operator<< ( ostream& os, const Term& t )
{
    os << t.name();
    return os;
}

int
Term::
operator== ( const Term& rhs ) const
{
    if( this == &rhs || name_ == rhs.name() )
	return 1;
    return 0;
}


void
Term::
init(ECString & prefix)
{
  ECString          fileName(prefix);
  // wul: terms.txt中为语法标识符号表
  fileName += "terms.txt";
  ifstream           stream(fileName.c_str(), ios::in);
  if (!stream)
    {
      cerr << "Can't open terms file " << fileName << endl;;
      return;
    }
  
  ECString          termName;
  int ind, n;
  n = 0;
  bool seenNTs = false;
  while (stream >> termName)
    {
      // wul: 第一列(termName)为名称，第二列(ind)为是否为终结符，1为，0不为，其他？
      // 例如 PRP就是终结符，不能继续展开了。
      stream >> ind;
      Term* nextTerm = new Term(termName, ind, n);
      // wul: 名称到Term的映射表
      termMap_[nextTerm->name()] = nextTerm;
      // wul: 为了方便，再记录若干特别的Term
      if(termName == "STOP") Term::stopTerm = nextTerm;
      else if(termName == "G4") Term::startTerm = nextTerm;
      else if(termName == "S1") Term::rootTerm = nextTerm;
      // wul: 按照index记录Term
      array_[n] = nextTerm;
      if(!ind && !seenNTs)
	{
          // wul: 文件的组织格式，为后面都是非终结符号(ind=0)
          // wul: 前面的都是终结符号，可以看作POS（tag?)。记录最后一个tag
	  assert(n > 0);
	  lastTagInt_ = n-1;
	  seenNTs = true;
	}
      n++;
      assert(n < 400);
    }
  assert(!ind);
  // wul: 记录最后一个非终结符号的序号
  lastNTInt_ = n-1;
  //lastNTInt_ = n-4;  //??? hack to ignore G1 and G2 and G3;
  stream.close();
}

Const_Term_p Term::get(const ECString& getname)
{
  TermMap::iterator ti = termMap_.find(getname);
  if( ti == termMap_.end()) return NULL;
  return (*ti).second;
}

