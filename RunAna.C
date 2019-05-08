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
#include"snana2_runAnalyzer.h"
#include"snana2_radrun.h"
#include"snana2_matching.h"


void plot(string NoFoilFile, string FoilFile, string HalfFoilFile)
{
 cout << "loading .. " << NoFoilFile << endl;
 cout << "loading .. " << FoilFile << endl;
 cout << "loading .. " << HalfFoilFile << endl;
  snana::runAnalyzer Analyzer0(NoFoilFile, 0, "NDK_NoFoil");
  snana::runAnalyzer Analyzer1(FoilFile, 1, "NDK_Foil");
  snana::runAnalyzer Analyzer2(HalfFoilFile, 1, "NDK_HalfFoil");
 
  TProfile *t0 = Analyzer0.getNPEsVsDrift(); t0->SetLineColor(1);
  TProfile *t1 = Analyzer1.getNPEsVsDrift(); t1->SetLineColor(2);
  TProfile *t2 = Analyzer2.getNPEsVsDrift(); t2->SetLineColor(3);


  TCanvas *c = new TCanvas();

  t0->Draw("HIST");
  t1->Draw("HIST SAME");
  t2->Draw("HIST SAME");
  gPad->BuildLegend();
}

void plot2(string NoFoilFile, string FoilFile, string HalfFoilFile)
{
 cout << "loading .. " << NoFoilFile << endl;
 cout << "loading .. " << FoilFile << endl;
 cout << "loading .. " << HalfFoilFile << endl;
  snana::runAnalyzer Analyzer0(NoFoilFile, 0, "NDK_NoFoil");
  snana::runAnalyzer Analyzer1(FoilFile, 1, "NDK_Foil");
  snana::runAnalyzer Analyzer2(HalfFoilFile, 1, "NDK_HalfFoil");
 
  TH1D *t0 = Analyzer0.getNPEsPerPMT(); t0->SetLineColor(1);
  TH1D *t1 = Analyzer1.getNPEsPerPMT(); t1->SetLineColor(2);
  TH1D *t2 = Analyzer2.getNPEsPerPMT(); t2->SetLineColor(3);


  TCanvas *c = new TCanvas();

  t0->Draw("HIST");
  t1->Draw("HIST SAME");
  t2->Draw("HIST SAME");
  gPad->BuildLegend();
}

void plot3(string NoFoilFile, string FoilFile, string HalfFoilFile)
{
 cout << "loading .. " << NoFoilFile << endl;
 cout << "loading .. " << FoilFile << endl;
 cout << "loading .. " << HalfFoilFile << endl;
  snana::runAnalyzer Analyzer0(NoFoilFile, 0, "NDK_NoFoil");
  snana::runAnalyzer Analyzer1(FoilFile, 1, "NDK_Foil");
  snana::runAnalyzer Analyzer2(HalfFoilFile, 1, "NDK_HalfFoil");
 
  TH1D *t0 = Analyzer0.getNPEs(); t0->SetLineColor(1);
  TH1D *t1 = Analyzer1.getNPEs(); t1->SetLineColor(2);
  TH1D *t2 = Analyzer2.getNPEs(); t2->SetLineColor(3);


  TCanvas *c = new TCanvas();

  t0->Draw("HIST");
  t1->Draw("HIST SAME");
  t2->Draw("HIST SAME");
  gPad->BuildLegend();
}

void RunAna( int i=0,int j=0, int k=4)
{


 string fileRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20181123_v07_11_00_light_newtree.root";
 string fileNDK="/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root";
 string fileLONGRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root"; 
 //string fileSN="/eos/user/a/angalleg/TDR_simulations/data/newtree/SN_20181210_v07_11_00_Livermore_light_newtree.root";


 string fileRADFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_Foil.root";
 string fileRADHalfFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_HalfFoil.root";
 string fileNDKHalfFoil="/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_HalfFoil_NDK_barnuK+.root";
 string fileNDKFoil="/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_Foil_NDK_barnuK+.root";



 string fileRADHalfFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_HalfFoil_20000.root"; 
 string fileRADFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_Foil_20000.root";

// plot(fileNDK,fileNDKFoil,fileNDKHalfFoil);
 if(i==0)plot2(fileLONGRAD,fileRADFoil20k,fileRADHalfFoil20k);
 if(i==1)plot3(fileLONGRAD,fileRADFoil20k,fileRADHalfFoil20k);


//CreateTGraphFile("RAD_SuperCluster_4k_D2m5_W800ns_MW1us.root");
//CreateTGraphFile("RAD_SuperCluster_40k_D2m5_W800ns_MW1us.root");
//MatchDistanceScan("NDK_SuperCluster_D2m5_W800ns_MW1us.root","RAD_SuperCluster_4k_D2m5_W800ns_MW1us.root","D2m5_W800ns_MW1us");

//MatchDistanceScan("NDK_SuperCluster_D2m5_W800ns_MW1us.root","RAD_SuperCluster_40k_D2m5_W800ns_MW1us.root","D2m5_W800ns_MW1us");

///Debugging("/eos/user/a/angalleg/TDR_simulations/data/newtree/SN_20181210_v07_11_00_Livermore_light_newtree.root_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us");


//ploteffvsdistanceNDK("NDK_SuperCluster_D2m5_W800ns_MW1us.root","RAD_SuperCluster_4k_D2m5_W800ns_MW1us.root","D2m5_W800ns_MW1us");
//ploteffvsdistanceNDK("NDK_SuperCluster_D2m5_W800ns_MW1us.root","RAD_SuperCluster_40k_D2m5_W800ns_MW1us.root","D2m5_W800ns_MW1us");

//  for(int ii=0;ii<3;ii++) 
//  for(int jj=0;jj<5;jj++) for(int kk=k;kk<k+1;kk++)
//  genClusters(filename,ii,jj,kk);
//  checkfiles();
//   run_missing(i,j);

//for(int ii=i;ii<i+1;ii++) for(int jj=0;jj<5;jj++)for(int kk=0;kk<9;kk++) barridoBG(ii,jj,kk);
//for(int ii=i;ii<i+1;ii++) for(int jj=0;jj<5;jj++)for(int kk=0;kk<9;kk++) barridoBG_SN(ii,jj,kk);

//  effplotNDK(i);
//  effplotSN();
//    ploteffvsdistanceNDK(1,4,0);

}
