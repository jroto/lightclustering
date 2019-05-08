/*

#include"MyStyle.C"

*/
#include <iostream>
#include <fstream>
#include"TGraphErrors.h"
#include"TPaveText.h"
#include"TTree.h"
#include"TFile.h"
#include"TMath.h"
#include"TAxis.h"
#include"TTimer.h"
#include"TCanvas.h"
#include"TH2D.h"
#include"TLegend.h"
#include"TMultiGraph.h"


using namespace std;
//*
void lets_pause (){
	TTimer * timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	timer->TurnOn();
	timer->Reset();
	std::cout << "q/Q to quit, other to continuee: ";
	char kkey;
	std::cin.get(kkey);
	if( kkey == 'q' || kkey == 'Q') 			std::exit(0); //gSystem->Exit(0);
	timer->TurnOff();
	delete timer;
}//*/
//gInterpreter->GenerateDictionary("Cluster", "Cluster.h");
//gSystem->Load("../lib/libmylibrary.so")

//#include"Cluster.h"
#include"pmtmap.h"
#include"snana2_cluster.h"
#include"snana2_clustering.h"
#include"snana2_clusteringAnalyzer.h"
#include"snana2_radevt.h"
#include"snana2_radrun.h"
#include"snana2_matching.h"

 const string sClusteringD[3]={"3m5","4m5","5m5"};
 const string sClusteringW[3]={"800ns","50ns","100ns"};
 const string sClusteringMW[3]={"1us","1us","100ns"};
 const float fClusteringD[3]={350,450,550};//cm
 const float fClusteringW[3]={0.8,0.05,0.1};//us
 const float fClusteringMW[3]={1.,1.,0.1};//us



 const string fFiles[9]={
                         "/eos/user/j/jsotooto/DUNE/lightana/data/radiologicals/Foils/20190313_v08_12_00_light_Foil_20000.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_Foil_NDK_barnuK+.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/20190313_v08_12_00_light_Foil_SN_Livermore.root",
                          "/eos/user/j/jsotooto/DUNE/lightana/data/radiologicals/HalfFoil/20190313_v08_12_00_light_HalfFoil_20000.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_HalfFoil_NDK_barnuK+.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/20190313_v08_12_00_light_HalfFoil_SN_Livermore.root",
                         "/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root", 
                         "/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/SN_20181210_v07_11_00_Livermore_light_newtree.root"};


string SuperclusterFilename(int f,int i){ return fFiles[f]+"_D"+sClusteringD[i]+"_W"+sClusteringW[i]+"_MW"+sClusteringMW[i]+"_SuperCluster.root";}
string VarName(int i){ return "D"+sClusteringD[i]+"_W"+sClusteringW[i]+"_MW"+sClusteringMW[i];}


double error(double a,double nevents)
{
 return TMath::Sqrt(a*(1.0-a)/nevents);
}

void BGOrigins(string SignalFile, string BGFile, string var, int split, bool genfiles, int DriftDirection)
{

//splitting 0, 1 ,2. 0 - no split (we get the average # of events), 1 split (we do NxNRAD events), 2 (we do NRAD events), 3 (we do N events).

  cout << "Scan on Matchign distance on file " << BGFile << " splitting? " << split << endl;

  snana::matching MyMatch(SignalFile, BGFile, var, genfiles,DriftDirection);

  TFile *ofile;
  ofile = new TFile(Form("%s_BGClusterPerOrigin.root",BGFile.c_str()),"RECREATE");
  ofile->cd();
  double purity, noSIGeff, numberOfBGClusters, nevents;

  TH1F *myhist = MyMatch.OriginOfMatchedEvents(150, 0);

  myhist->Draw("HIST");
  ofile->cd();
  myhist->Write("Origins");
 lets_pause();
}


void BGAnalysis( int i=0,int j=0, int k=4)
{

 int ClusterConf=2;
 BGOrigins(SuperclusterFilename(4,ClusterConf),SuperclusterFilename(3,ClusterConf),VarName(ClusterConf),2,0,1);

}
