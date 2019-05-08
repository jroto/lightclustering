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


void MatchDistanceScan(string SignalFile, string BGFile, string var, int split, bool genfiles, int DriftDirection, std::vector<double> ddistance)
{

//splitting 0, 1 ,2. 0 - no split (we get the average # of events), 1 split (we do NxNRAD events), 2 (we do NRAD events), 3 (we do N events).

  cout << "Scan on Matchign distance on file " << BGFile << " splitting? " << split << endl;

  snana::matching MyMatch(SignalFile, BGFile, var, genfiles,DriftDirection);

  TFile *ofile;
  switch(split)
  {
    case 0: ofile = new TFile(Form("%s_TGraph_nosplit.root",BGFile.c_str()),"RECREATE");break;
    case 1: ofile = new TFile(Form("%s_TGraph_split.root",BGFile.c_str()),"RECREATE");break;
    case 2: ofile = new TFile(Form("%s_TGraph_split_NRAD.root",BGFile.c_str()),"RECREATE");break;
    case 3: ofile = new TFile(Form("%s_TGraph_split_NEVENTS.root",BGFile.c_str()),"RECREATE"); break;
  }
  ofile->cd();
  TGraphErrors* tg = new TGraphErrors();
  tg->GetXaxis()->SetTitle("Maximum Cluster-MCTruth distance (cm)");
  tg->GetYaxis()->SetTitle("Efficiency");

  TGraphErrors* tg2 = new TGraphErrors();
  tg2->GetXaxis()->SetTitle("Maximum Cluster-MCTruth distance (cm)");
  tg2->GetYaxis()->SetTitle("Purity"); //only signal matching

  TGraphErrors* tg3 = new TGraphErrors();
  tg3->GetXaxis()->SetTitle("Maximum Cluster-MCTruth distance (cm)");
  tg3->GetYaxis()->SetTitle("Efficiency x Purity");//assuming perfect signal matching

  TGraphErrors* tg4 = new TGraphErrors();
  tg4->GetXaxis()->SetTitle("Maximum Cluster-MCTruth distance (cm)");
  tg4->GetYaxis()->SetTitle("Number of BG Clusters in vicinity");//assuming perfect signal matching

  

  for(int i=0; i<ddistance.size(); i++)
  {
    std::cout << "calculating FME at " << ddistance[i] << std::endl;
    double Purity, noSIGeff, numberOfBGClusters;
    double efficiency, nevents;
    switch(split)
    {
      case 0: efficiency = MyMatch.MatchingEfficiencyFast(ddistance[i],Purity,noSIGeff, numberOfBGClusters);break; //using tH2 deprecated
      case 1: efficiency = MyMatch.MatchingEfficiencySplitting(ddistance[i],Purity,noSIGeff, numberOfBGClusters);break; //using th3, too slow
      case 2: efficiency = MyMatch.MatchingEfficiencyRight(ddistance[i],Purity,noSIGeff, numberOfBGClusters, nevents, 0);break; //NRADevents
      case 3: efficiency = MyMatch.MatchingEfficiencyRight(ddistance[i],Purity,noSIGeff, numberOfBGClusters, nevents, 1);break; //N events
    }

    tg->SetPoint(tg->GetN(),ddistance[i],efficiency);
    tg->SetPointError(tg->GetN()-1,0,error(efficiency,nevents));

    tg2->SetPoint(tg2->GetN(),ddistance[i],Purity);
    tg2->SetPointError(tg2->GetN()-1,0,error(Purity,efficiency*nevents));

    tg3->SetPoint(tg3->GetN(),ddistance[i],efficiency*Purity);
    tg3->SetPointError(tg3->GetN()-1,0,error(efficiency*Purity,nevents));

    tg4->SetPoint(tg4->GetN(),ddistance[i],efficiency*numberOfBGClusters);
    tg4->SetPointError(tg4->GetN()-1,0.0,0.0);

    std::cout << ddistance[i] <<" \t" << efficiency  << "\t" << Purity <<  "\t" << efficiency*Purity <<"\t" <<numberOfBGClusters <<std::endl;
  }
   ofile->cd();
   tg->Write("eff");
   tg2->Write("Purity");
   tg3->Write("effNoSIG");
   tg4->Write("numberOfBGClusters");
   ofile->Close();


}


void Debugging(string SignalFile, string var)
{
 cout << "loading .. " << SignalFile << endl;
  snana::clusteringAnalyzer Analyzer(SignalFile, var); //1 to generate densities
 
  std::map<int,TH1D> mymap = Analyzer.NumberOfPEsPerClusterPerPDG(5.0);

   map<int,TH1D>::iterator it;

   for ( it = mymap.begin(); it != mymap.end(); it++ )
   {
     std::cout << it->first  << std::endl ;it->second.Draw();   // string's value 
     lets_pause();
   }
}

void SetAxis(TH2D *h, int nx, int ny)
{

  const char *sClusteringD[5]={"1.5","2.5","3.5","4.5","5.5"};
  const char *sClusteringW[5]={"50","100","200","400","800"};
  const char * sClusteringMW[9]={"1us","2us","4us","8us","16us","500ns","250ns","125ns","63ns"};


   TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,800,800);
   c1->SetGrid();
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);

//   h->SetStats(0);
//   h->GetXaxis()->SetLabelOffset(99);
//   h->GetYaxis()->SetLabelOffset(99);
//   h->Draw("text");
   // draw labels along X
//   Float_t x, y;
//   y = gPad->GetUymin() - 0.2*h->GetYaxis()->GetBinWidth(1);

   for (int i=1;i<=nx;i++) h->GetXaxis()->SetBinLabel(i,sClusteringD[i-1]);
   for (int i=1;i<=ny;i++) h->GetYaxis()->SetBinLabel(i,sClusteringW[i-1]);
   //gPad->SetOptStat(0);
//   h->Draw("colz");lets_pause();

}


void ploteffvsdistanceNDK(string SignalFile, string BGFile, string var, int split, double distance, int DriftDirection, int axis)
{

  std::cout  <<"DOING TProfile of the EFFICIENCIES along X axis! Splitting?=" << split  <<std::endl;
  snana::matching MyMatch(SignalFile, BGFile, var,0,DriftDirection);

  TProfile *tpPurity =new TProfile();
  TProfile *tpEffxPurity=new TProfile();
  TProfile *NumberOfBGClusters = new TProfile();

  TProfile *tpeff;

  switch(split)
  {
    case 0: tpeff = MyMatch.MatchingEfficiencyFastScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters, 20,0); break;
    case 1: tpeff= MyMatch.MatchingEfficiencySplittingScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters, 20,0); break;
    case 2: tpeff = MyMatch.MatchingEfficiencyRightScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters,0, 18,axis,8); break;
    case 3: tpeff = MyMatch.MatchingEfficiencyRightScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters,1, 24,axis,8); break;//MatchingEfficiencyRightScan(double distance, TProfile &Purity, TProfile &EffxPurity, TProfile &NumberOfBGClusters, bool AllowRepeatBackground, int nbins=12, int axis=0,int BGeventsPerSignalEvent=8)
  } 

//  tpnoBGeff->Draw();tpnoSIGeff->SetLineColor(2);
//  tpnoSIGeff->Draw("SAME");
  TFile *ofile;
  switch(split)
  {
    case 0: ofile = new TFile(Form("%s_TProfile_effvsdistance.root",BGFile.c_str()),"RECREATE"); break;
    case 1: ofile = new TFile(Form("%s_TProfile_effvsdistance_split.root",BGFile.c_str()),"RECREATE"); break;
    case 2: ofile = new TFile(Form("%s_TProfile_effvsdistance_NRAD.root",BGFile.c_str()),"RECREATE"); break;
    case 3: ofile = new TFile(Form("%s_TProfile_effvsdistance_NEVENTS.root",BGFile.c_str()),"RECREATE"); break;
  } 
  ofile->cd();
  tpeff->Write("eff");
  tpPurity->Write("Purity"); //no background
  tpEffxPurity->Write("EffxPurity"); //assuming perfect signal matching 
  NumberOfBGClusters->Write("NumberOfBGClusters"); 

  ofile->Close();
  std::cout  <<" \tDONE!"  <<std::endl;

}


void EfficiencyAnalysis( int i=0,int j=0, int k=4)
{
 std::vector<double> ddistance; ddistance.resize(13);
 for (int i=0;i<13;i++) ddistance[i]=50+i*50;
 std::vector<double> ddistance2; ddistance2.resize(15);ddistance2[0]=100;
 for (int i=1;i<15;i++) ddistance2[i]=ddistance2[i-1]+10;

 string fileNDK="/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root";
 string fileLONGRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root"; 


 switch(i)
 {
 case 0:  MatchDistanceScan(SuperclusterFilename(4,0),SuperclusterFilename(3,0),VarName(0),2,0,1,ddistance);break;
 case 1:  MatchDistanceScan(SuperclusterFilename(4,1),SuperclusterFilename(3,1),VarName(1),2,0,1,ddistance);break;
 case 2:  MatchDistanceScan(SuperclusterFilename(4,2),SuperclusterFilename(3,2),VarName(2),2,0,1,ddistance);break;
 case 3:  ploteffvsdistanceNDK(SuperclusterFilename(4,0),SuperclusterFilename(3,0),VarName(0),2,150,1,1);break;
 case 4:  ploteffvsdistanceNDK(SuperclusterFilename(4,1),SuperclusterFilename(3,1),VarName(1),2,150,1,1);break;
 case 5:  ploteffvsdistanceNDK(SuperclusterFilename(4,2),SuperclusterFilename(3,2),VarName(2),2,150,1,1);break;

 case 6:  MatchDistanceScan(fileNDK+"_SuperCluster_D2m5_W1us_MW1us.root",fileLONGRAD+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,0,ddistance);break;
 case 7:  ploteffvsdistanceNDK(fileNDK+"_SuperCluster_D2m5_W1us_MW1us.root",fileLONGRAD+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,250,0,0);break;
}

}
