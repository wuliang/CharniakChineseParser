#include "Edge.h"
#include "EdgeHeap.h"

EdgeHeap::
~EdgeHeap()
{
  int i;
  for(i = 0 ; i < unusedPos_ ; i++) delete array[i];
}

EdgeHeap::
EdgeHeap()
{
  int i;
  for (i = 0 ; i < HeapSize ; i++) array[i] = NULL;
  print = false;
  unusedPos_ = 0;
}


void
EdgeHeap::
insert(Edge* edge)
{
  if(print)
    cerr << "heap insertion of " << *edge << " at " << unusedPos_ << endl;
  array[unusedPos_] = edge;
  edge->heapPos() = unusedPos_;
  upheap(unusedPos_);
  unusedPos_++;
  assert(unusedPos_ < HeapSize);
}

bool
EdgeHeap::
upheap(int pos)
{
  // wul: 把merit值高的节点向上移动。这样插入的结果，就是二叉树
  // wul: 按层有次序。即对每一条从叶子到根部的路径，metri值都是递增的。
  if(print) cerr << "in Upheap " << pos << endl;
  if(pos == 0) return false;
  Edge* edge = array[pos];
  assert(edge->heapPos() == pos);
  double merit = edge->merit();
  int   parPos = parent(pos);
  Edge* par = array[parPos];
  assert(par->heapPos() == parPos);
  // wul: add print
  if(print) cerr << "pos " << pos    << ",m=" << merit << ";";
  if(print) cerr << "par " << parPos << ",m=" << par->merit() << endl;

  if(merit > par->merit())
    {
      array[parPos] = edge;
      edge->heapPos() = parPos;
      array[pos] = par;
      par->heapPos() = pos;
      if(print) cerr << "Put " << *edge << " in " << parPos << endl;
      upheap(parPos);
      return true;
    }
  else if(print)
    {
      cerr << "upheap of " << merit << "stopped by "
	<< *par << " " << par->merit() << endl;
    }
  return false;
}


Edge*
EdgeHeap::
pop()
{
  // wul: pop出来的总是index为0的，即堆的根节点。
  if(print)
    cerr << "popping" << endl;
  if(unusedPos_ == 0) return NULL;
  Edge* retVal(array[0]);
  assert(retVal->heapPos() == 0);
  del_(0);
  // wul: 从堆里面弹出的Edge，应该清除pos信息。
  retVal->heapPos() = -1;
  return retVal;
}

void
EdgeHeap::
downHeap(int pos)
{
  if(print) cerr << "downHeap " << pos << endl;
  if(pos >= unusedPos_-1) return;
  Edge* par = array[pos];
  assert(par->heapPos() == pos);
  double merit = par->merit();
  // wul: 计算理论上的左孩子，右孩子节点的序号。
  int lc = left_child(pos);
  int rc = right_child(pos);
  int largec;
  int lcthere = 0;
  Edge* lct;
  if(lc < unusedPos_)
    {
      // 如果实际使用了，就取对应的节点
      lct = array[lc];
      if(lct)
	{ lcthere = 1;
	  assert(lct->heapPos() == lc);
	}
    }
  int rcthere = 0;
  Edge* rct;
  if(rc < unusedPos_)
    {
      rct = array[rc];
      if(rct)
	{
	  rcthere = 1;
	  assert(rct->heapPos() == rc);
	}
    }
  // wul: 至少存在1个子节点
  if(!lcthere && !rcthere) return;
  // wul: 如果存在，左节点必须存在（因为它的序号小）
  assert(lcthere);
  // wul: 选择子节点中merit值较大的。
  if(!rcthere || (lct->merit() > rct->merit()))
    largec = lc;
  else largec = rc;

  Edge* largeEdg = array[largec];
  if(merit >= largeEdg->merit())
    { 
      // wul: 如果本节点的merit更大（比两个子节点都大），则没有必要下移
      if(print) cerr << "downheap of " << merit << " stopped by "
		     << *largeEdg << " " << largeEdg->merit() << endl;
      return;
    }
  // wul: 否则就交换位置。
  array[pos] = largeEdg;
  largeEdg->heapPos() = pos;
  array[largec] = par;
  par->heapPos() = largec;
  // wul: 继续下移
  downHeap(largec);
}

void
EdgeHeap::
del(Edge* edge)
{
  // wul: 在堆中删除该edge的信息。正常而言，在堆中的edge上包含pos信息。
  if(print)
    cerr << "del " << edge << endl;
  int pos = edge->heapPos();
  assert( pos < unusedPos_ && pos >= 0);
  del_( pos );
}

void
EdgeHeap::
del_(int pos)
{
  if(print) cerr << "del_ " << pos << endl;
  assert(unusedPos_);
  if(pos == (unusedPos_ - 1) )
    {
      // wul: 如果删除的正好是堆中最后一个，那么就简单了。除去标志和计数。
      unusedPos_--;
      array[unusedPos_] = NULL;
      return;
    }
  /* move the final edge in heap to empty position */
  // wul: 否则用堆中最后一个节点来填补被删的节点。
  array[pos] = array[unusedPos_ - 1];
  if(!array[pos])
    {
      error("Never get here");
      return;
    }
  array[pos]->heapPos() = pos;
  array[unusedPos_ -1] = NULL;
  unusedPos_--;
  // wul: 因为最后一个节点和被删的节点原来属于的不同的分支，所以节点移动过来之后
  // wul: 为保证原有属性（每个分支有序），需要调整它的位置
  // wul: 尝试向上移动比较容易，而且如果向上移动成功（至少移动1次），也没必要尝试下移了。
  if(upheap(pos)) return;
  downHeap(pos);
}
/*
void
EdgeHeap::
check()
{
  if(size() > 0) array[0]->check();
  for(int i = 1 ; i < unusedPos_ ; i++)
    {
      assert(array[i]);
      array[i]->check();
      if(!(array[parent(i)]->merit() >= array[i]->merit()))
	{
	 cerr << "For i = " << i <<  " parent_i = "
	   << parent(i) << " "
	   << *(array[parent(i)])
	   << " at " << array[parent(i)]->merit() 
	   << " not higher than " << *(array[i])
	   << " at " << array[i]->merit() 
	     << endl;
	 assert(array[parent(i)]->merit() >= array[i]->merit());
       }
    }
}
*/
