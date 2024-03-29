/*
 * Copyright 1999, Brown University, Providence, RI.
 * 
 *                         All Rights Reserved
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose other than its incorporation into a
 * commercial product is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Brown University not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.
 * 
 * BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 * ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <time.h>
#include <sys/resource.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include "ewDciTokStrm.h"
#include "Field.h"
#include "GotIter.h"
#include "Wrd.h"
#include "InputTree.h"
#include "Bchart.h"
#include "ECArgs.h"
#include "MeChart.h"
#include "headFinder.h"

MeChart* curChart;

float endFactor = 1.2;
float midFactor = 1;

class Params 
{
public:
    Params(): 	    file(0),
		    fileString_(),
		    whichSent_( 0 ),
		    ofTotal_( 1 ),
		    field_( 0 ),
                    outputData_(false),
                    stdInput_(false),
		    numString_()
		    {}
    void	    init( ECArgs& args );
    const char *    file;
    const ECString&   fileString()
		    {  return fileString_;  }
    const ECString&   numString()
		    {   return numString_;  }
    const int	    whichSent()
		    {   return whichSent_;   }
    const int	    ofTotal()
		    {   return ofTotal_;   }
    const Field&    field() const
		    {   return *field_;   }
    bool&      stdInput() { return stdInput_; }
    bool&      outputData() { return outputData_; }
private:
    bool       stdInput_;
    bool       outputData_;
    ECString	    fileString_;
    ECString          numString_;
    int             whichSent_;
    int             ofTotal_;
    Field *         field_;
};

InputTree* inputTreeFromAnswerTree(AnswerTree* at, Item* itm);

void
Params::
init( ECArgs& args ) 
{  
  fileString_ = args.arg(0);
  if( args.isset('n') )
    {
      char etemp[16];
      strcpy(etemp,args.value('n').c_str());
      char *	temp = strchr( etemp, '/' );
      if( !temp )
	error( "No terminal '/' found in '-n' argument" );
      *temp = '\0';
      ofTotal_ = atoi( ++temp );
      char *	mask = new char[ ofTotal_ ];
      for( int i = 0; i < ofTotal_; i++ )
	mask[ i ] = 0;
      // fill in mask with valid numbers;
      ECString tmp2 = etemp;
      numString_ = tmp2;		// meaningful id for this process;
      whichSent_ = atoi(tmp2.c_str());
      mask[ whichSent_ ] = 1;
      field_ = new Field( ofTotal_, mask );
    }
    else
    {
	static char mask[1] = { 1 };
	field_ = new Field( 1, mask );
    }
}

// Global objects
//
//		(or object wannabes)
//

Params 		    params;

int
main(int argc, char *argv[])
{
   ECArgs args( argc, argv );
   /* o = basic, but not debugging, output.
      l = length of sentence to be proceeds 0-40 is default
      n = work on each #'th line.
      d = print out debugging info at level #
      W = use wwclasses
      R = use rwclasses
      t = report timings (requires o)
      s = maximum sleep time
      f = f# says multiply ctl2 counts by #
      p = p# use prepFactor #
      P = which types of prob models to use */

   // prevent core file creation;
   struct rlimit 	core_limits;
   core_limits.rlim_cur = 0;
   core_limits.rlim_max = 0;
   setrlimit( RLIMIT_CORE, &core_limits );

   params.init( args );
   if(args.isset('s'))
     {
       int  maxDelay = atoi(args.value('s').c_str());
       srand(params.whichSent());
       int randN = rand();
       int delay = randN%maxDelay;
       sleep(delay);
     }

   if(args.isset('T'))
     {
       int fac = atoi(args.value('T').c_str());
       float ffac = (float)fac;
       ffac /= 10;
       Bchart::timeFactor = ffac;
     }
	 
   int maxSentLen = 70;
   if(args.isset('l'))
     {
       maxSentLen = atoi(args.value('l').c_str());
     }
   int    totEdges = 0;
   int    totPopedEdges = 0;
   double totAccessTime = 0;
   double totParseTime = 0;
   double totSemParseTime = 0;
   clock_t lastTime, currTime;
   double lastTimeSec, currTimeSec, elapsedTime;

   endFactor = 1.2;
   midFactor = (1.0 - (.3684 * endFactor))/(1.0 - .3684);

   if( args.nargs() > 2 || args.nargs() == 0 )	// require path name 
     error( "Need exactly two arg." );
   ECString  path( args.arg( 0 ) );
   readHeadInfo(path);
   Term::init( path );
   InputTree::init();

   ECString testSString( args.arg(1) );

   ewDciTokStrm testSStream(testSString);
   //ifstream testSStream(testSString.c_str());
   if( !testSStream ) error( "No testSstream" );
   int      sentenceCount = 0;  //counts all sentences so we can use 1/50;

   ECString  probSumString( path );
   probSumString += "pSgT.txt";
   ifstream    probSumStream( probSumString.c_str() );
   if( !probSumStream ) error( "Failed to find probSum file" );

   Bchart::readTermProbs(path);

   if( args.isset('d') )
     {
       int lev = atoi(args.value('d').c_str());
       Bchart::printDebug() = lev;
     }
   int totSents = 0;
   int totUnparsed = 0;

   MeChart::init(path);
   Bchart::setPosStarts();
   for( ; !(!testSStream) ; )
     {
       SentRep sr(testSStream, SentRep::SGML); 
       int len = sr.length();
       if(len == 0) continue;
       if(len > maxSentLen) continue;
       if( !params.field().in(sentenceCount) )
	 {
	   sentenceCount++;
	   continue;
	 }
       if(len == 1)
	 {
	   if(sr[0].lexeme() == "</DOC>")
	     {
	       continue;
	     }
	 }
       sentenceCount++;

       //SentRep orgsr( wtList );  // used in precision calc;

       if( args.isset('t') ) lastTime = clock();
       if(args.isset('t') )
	 {
	   currTime = clock();
	   lastTimeSec = (double)lastTime/(double)CLOCKS_PER_SEC;
	   currTimeSec = (double)currTime/(double)CLOCKS_PER_SEC;
	   elapsedTime = currTimeSec - lastTimeSec;
	   if(elapsedTime < 0) elapsedTime += 2147;
	   cerr << "Reading data time = " << elapsedTime << endl;
	   totAccessTime += elapsedTime;
	   lastTime = currTime;
	 }

       MeChart*	chart = new MeChart( sr );
       curChart = chart;
       chart->ruleCountTimeout() = 250000;
       
       totSents++;
       if(args.isset('t') )
	 lastTime = clock();
       double tmpCrossEnt = chart->parse( );
       Item* topS = chart->topS();

       if(!topS)
	 {
	   if(len == 1)
	     {
	       delete chart;
	       continue;
	     }
	   Edge::DemFac = .9;
	   delete chart;
	   chart = new MeChart(sr);
	   chart->ruleCountTimeout() = 350000;
	   curChart = chart;
	   tmpCrossEnt = chart->parse( );
	   topS = chart->topS();
	   Edge::DemFac = .999;
	   if(!topS)
	     {
	       totUnparsed++;
	       cerr << "Parse failed on: " << sr << endl;

	       delete chart;
	       continue;
	     }
	 }
       
       // compute the outside probabilities on the items so that we can
       // skip doing detailed computations on the really bad ones 
       if(args.isset('t') )
	 {
	   currTime = clock();
	   lastTimeSec = (double)lastTime/(double)CLOCKS_PER_SEC;
	   currTimeSec = (double)currTime/(double)CLOCKS_PER_SEC;
	   elapsedTime = currTimeSec - lastTimeSec;
	   if(elapsedTime < 0) elapsedTime += 2147;
	   cerr << "Parsing time = " << elapsedTime
	     << "\tEdges created = " << chart->totEdgeCountAtS()
	       << "\tEdges poped = " << chart->popedEdgeCountAtS() << endl;
	   totParseTime += elapsedTime;
	   //totEdges += chart->totEdgeCountAtS();
	   //totPopedEdges += chart->popedEdgeCountAtS();
	   totEdges += chart->totEdgeCountAtS();
	   totPopedEdges += chart->popedEdgeCountAtS();
	   lastTime = clock();

	 }

       chart->set_Alphas();

       AnswerTree* at = chart->findMapParse();
       if( !at ) 
	 {
	   totUnparsed++;
	   cerr << "MapParse failed on: " << sr << endl;
	   delete chart;
	   continue;
	 }
       InputTree*  mapparse = inputTreeFromAnswerTree(at,topS);
       //at->deleteSubTrees();
       //delete at;
       cout << *mapparse << endl;
       delete mapparse;

       if(args.isset('t') )
	 {
	   currTime = clock();
	   lastTimeSec = (double)lastTime/(double)CLOCKS_PER_SEC;
	   currTimeSec = (double)currTime/(double)CLOCKS_PER_SEC;
	   elapsedTime = currTimeSec - lastTimeSec;
	   if(elapsedTime < 0) elapsedTime += 2147;
	   cerr << "Sem Parsing time = " << elapsedTime << endl;
	   totSemParseTime += elapsedTime;
	 }

       delete chart;
     }
   if( args.isset('t') )
     cout << "Av access time = " << totAccessTime/totSents
       << "\t Av parse time = "
	 << totParseTime/totSents
       << "\t Av stats time = "
	 << totSemParseTime/totSents
       << "\nAv edges created = "
	 << (float)totEdges/totSents
       << "\tAv edges poped = "
	 << (float)totPopedEdges/totSents
	   << endl;

   return 0;
}


InputTree*
inputTreeFromAnswerTree(AnswerTree* at, Item* itm)
{
  int strt = itm->start();
  int fin = itm->finish();
  const Term* trm = itm->term();
  ECString trmString = trm->name();
  ECString wrdString;
  ECString ntString;
  list<InputTree*> subtrs;
  InputTree* ans;
  if(trm->terminal_p())
    {
      wrdString = itm->word()->lexeme();
    }
  else
    {
      assert(at);
      Edge* e = at->e;
      list<AnswerTree*>::iterator ati = at->subtrees.begin();
      Item* subi;
      LeftRightGotIter gi(e);
      for( ; ati != at->subtrees.end() ; ati++)
	{
	  gi.next(subi);
	  //stop terms will now appear in edges.
	  if(subi->term() == Term::stopTerm)
	    {
	      ati--;
	      continue;
	    }
	  AnswerTree* sab = *ati;
	  InputTree* sit = inputTreeFromAnswerTree(sab,subi);
	  subtrs.push_back(sit);
	}
    }
  ans = new InputTree(strt, fin, wrdString, trmString, ntString,
		      subtrs, NULL, NULL);
  return ans;
}
