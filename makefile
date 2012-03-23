.pony: all

all: parseIt parseStdin

.C.o:
	g++ $(CFLAGS) -O3 -c $<


PARSEIT_OBJS = \
	AnswerTree.o \
	Bchart.o \
	BchartSm.o \
	FBinaryArray.o \
	CntxArray.o \
	ChartBase.o \
	ECArgs.o \
	Edge.o \
	EdgeHeap.o \
	Feat.o \
	Feature.o \
	FeatureTree.o \
	Field.o \
	FullHist.o \
	GotIter.o \
	InputTree.o \
	Item.o \
	SentRep.o \
	Term.o \
	ccInd.o \
	ewDciTokStrm.o \
	edgeSubFns.o \
	fhSubFns.o \
	headFinder.o \
	utils.o \
	MeChart.o \
	parseIt.o 


PARSESTDIN_OBJS = \
	AnswerTree.o \
	Bchart.o \
	BchartSm.o \
	FBinaryArray.o \
	CntxArray.o \
	ChartBase.o \
	ECArgs.o \
	Edge.o \
	EdgeHeap.o \
	Feat.o \
	Feature.o \
	FeatureTree.o \
	Field.o \
	FullHist.o \
	GotIter.o \
	InputTree.o \
	Item.o \
	SentRep.o \
	Term.o \
	ccInd.o \
	ewDciTokStrm.o \
	edgeSubFns.o \
	fhSubFns.o \
	headFinder.o \
	utils.o \
	MeChart.o \
	parseStdin.o 


parseIt: $(PARSEIT_OBJS)
	g++ $(CFLAGS) $(PARSEIT_OBJS) -o parseIt


parseStdin: $(PARSESTDIN_OBJS)
	g++ $(CFLAGS) $(PARSESTDIN_OBJS) -o parseStdin


clean:
	@rm -f *.o parseIt parseStdin

