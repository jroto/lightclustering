#include"TROOT.h"
#include <iostream>
//#include"NDK_analyzer.C"
using namespace std;
void macrogrid(int i, int j, int k, int f)
{
  gROOT->ProcessLine(".L Cluster.h+");
  std::cout << "calling " << Form(".x SuperClustering.C(%i,%i,%i,%i)",i,j,k,f)<< std::endl;
  gROOT->ProcessLine(Form(".x SuperClustering.C(%i,%i,%i,%i)",i,j,k,f));
  //NDK_analyzer(imin,imax);
}
