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


  const string sClusteringD[5]={"1m5","2m5","3m5","4m5","5m5"};
  const string sClusteringW[6]={"50ns","100ns","200ns","400ns","800ns","1us"};
  const string sClusteringMW[9]={"1us","2us","4us","8us","16us","500ns","250ns","125ns","63ns"};
  const float fClusteringD[5]={150,250,350,450,550};//cm
  const float fClusteringW[6]={0.05,0.1,0.2,0.4,0.8,1.0};//us
  const float fClusteringMW[9]={1.,2.,4.,8.,16.,0.5,0.25,0.125,0.063};//us

///:
void genClusters(string filename, int i, int j, int k, string outstr)
{

//  string outfilename = filename+"_clustering_D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k]+".root";
  string outfilename=outstr+"_D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k]+".root";
  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");fout->cd();
  snana::radrun_snana run;
  run=snana::radrun_snana(fout,filename,sClusteringD[i],fClusteringD[i],sClusteringW[j],fClusteringW[j],sClusteringMW[k],fClusteringMW[k],0);
  run.Dump(fout);
  std::cout << "Dumped to " << outfilename << std::endl;
  fout->Close();
}



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


void MatchDistanceScanParallel(string SignalFile, string BGFile, string var, int split, bool genfiles, int dd,int DriftDirection)
{

//splitting 0, 1 ,2. 0 - no split (we get the average # of events), 1 split (we do NxNRAD events), 2 (we do NRAD events), 3 (we do N events).

  cout << "Scan on Matchign distance on file " << BGFile << " splitting? " << split << endl;

  snana::matching MyMatch(SignalFile, BGFile, var, genfiles,DriftDirection);

  ofstream ofile;

  switch(split)
  {
    case 0: ofile.open(Form("%s_TGraph_DEPRECATED_%i.txt",BGFile.c_str(),dd));break;
    case 1: ofile.open(Form("%s_TGraph_SPLIT_%i.txt",BGFile.c_str(),dd));break;
    case 2: ofile.open(Form("%s_TGraph_NRAD_%i.txt",BGFile.c_str(),dd));break;
    case 3: ofile.open(Form("%s_TGraph_NEVENTS_%i.txt",BGFile.c_str(),dd)); std::cout << "printing file " << std::endl; break;
  }

  ofile << "MCTruth-RecoDistance\tPurity\terror\tPurity_noBG\terror\tPurity_onlyBG\terror\tNClusters_in_vicinity" << endl;

    std::cout << "calculating FME at " << dd << std::endl;
    double noBGeff, noSIGeff, numberOfBGClusters;
    double efficiency, nevents;
    switch(split)
    {
      case 0: efficiency = MyMatch.MatchingEfficiencyFast((double)dd,noBGeff,noSIGeff, numberOfBGClusters);break; //using tH2 deprecated
      case 1: efficiency = MyMatch.MatchingEfficiencySplitting((double)dd,noBGeff,noSIGeff, numberOfBGClusters);break; //using th3, too slow
      case 2: efficiency = MyMatch.MatchingEfficiencyRight((double)dd,noBGeff,noSIGeff, numberOfBGClusters, nevents, 0);break; //NRADevents
      case 3: efficiency = MyMatch.MatchingEfficiencyRight((double)dd,noBGeff,noSIGeff, numberOfBGClusters, nevents, 1);break; //N events
    }
    ofile << dd << "\t" << efficiency<< "\t" <<error(efficiency, nevents)<< "\t" <<noBGeff<< "\t" <<error(noBGeff, nevents)<< "\t" <<noSIGeff<< "\t" <<error(noSIGeff, nevents) << "\t" << numberOfBGClusters <<endl;
    std::cout <<"Calculated efficiency: " <<  dd <<" \t" << efficiency <<std::endl;
    
   ofile.close();


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


void ploteffvsdistanceNDK(string SignalFile, string BGFile, string var, int split, double distance, int DriftDirection)
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
    case 2: tpeff = MyMatch.MatchingEfficiencyRightScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters,0, 24,1,8); break;
    case 3: tpeff = MyMatch.MatchingEfficiencyRightScan(distance,*tpPurity,*tpEffxPurity, *NumberOfBGClusters,1, 24,1,8); break;//MatchingEfficiencyRightScan(double distance, TProfile &Purity, TProfile &EffxPurity, TProfile &NumberOfBGClusters, bool AllowRepeatBackground, int nbins=12, int axis=0,int BGeventsPerSignalEvent=8)
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

void plotMaxClusterPerDrift(string f1, string f2, string f3, string endname)
{
   snana::clusteringAnalyzer ana1(f1+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana2(f2+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana3(f3+endname, "D2m5_W1us_MW1us");

  TProfile *p1 = ana1.GetPEvsDrift(0);p1->SetLineColor(1);
  TProfile *p2 = ana2.GetPEvsDrift(1);p2->SetLineColor(2);
  TProfile *p3 = ana3.GetPEvsDrift(1);p3->SetLineColor(3);
  p1->Draw("HIST"); 
  p2->Draw("HIST SAME"); 
  p3->Draw("HIST SAME"); lets_pause();

}

void plotRecoLightRatio(string f1, string f2, string f3, string endname)
{
   snana::clusteringAnalyzer ana1(f1+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana2(f2+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana3(f3+endname, "D2m5_W1us_MW1us");

  TProfile *p1 = ana1.RecoLightRatio(0);p1->SetLineColor(1);
  TProfile *p2 = ana2.RecoLightRatio(1);p2->SetLineColor(2);
  TProfile *p3 = ana3.RecoLightRatio(1);p3->SetLineColor(3);
  p1->Draw("HIST"); 
  p2->Draw("HIST SAME"); 
  p3->Draw("HIST SAME"); lets_pause();

}

void plotBGLightFreq(string f1, string f2, string f3, string endname)
{
   snana::clusteringAnalyzer ana1(f1+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana2(f2+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana3(f3+endname, "D2m5_W1us_MW1us");

  TH1F *p1 = ana1.RecoLightFreq("noFoil");p1->SetLineColor(1);
  TH1F *p2 = ana2.RecoLightFreq("Foil");p2->SetLineColor(2);
  TH1F *p3 = ana3.RecoLightFreq("HalfFoil");p3->SetLineColor(3);
  TCanvas *c =new TCanvas();
  p1->Draw("HIST"); 
  p2->Draw("HIST SAME"); 
  p3->Draw("HIST SAME");
  gPad->BuildLegend();
  lets_pause();

}


void plotRecoLightFreqByBGGen(string f1, string f2, string f3, string endname)
{


  std::cout  <<"Analysing plotRecoLightFreqByBGGen "  <<std::endl;
   snana::clusteringAnalyzer ana1(f1+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana2(f2+endname, "D2m5_W1us_MW1us");
   snana::clusteringAnalyzer ana3(f3+endname, "D2m5_W1us_MW1us");
  int BGGens[6]={0,4,5,6,8,9};
  string sBGGens[6]={"DarkCounts","Ar39","Neutron","Krypton","Radon","Ar42"};

std::cout  <<"Analysing plotRecoLightFreqByBGGen 0"  <<std::endl;
  //TH1F *p1[6]; for(int i=0; i<6; i++) {p1[i] = ana1.RecoLightFreqByBGGen(Form("noFoil_%s",sBGGens[i].c_str() ),BGGens[i]);p1[i]->SetLineColor(2+i);}
std::cout  <<"Analysing plotRecoLightFreqByBGGen 1"  <<std::endl;
  //TH1F *p2[6]; for(int i=0; i<6; i++) {p2[i] = ana2.RecoLightFreqByBGGen(Form("Foil_%s",sBGGens[i].c_str() ),BGGens[i]);p2[i]->SetLineColor(2+i);}
std::cout  <<"Analysing plotRecoLightFreqByBGGen 2"  <<std::endl;
  TH1F *p3[6]; for(int i=0; i<6; i++) {p3[i] = ana3.RecoLightFreqByBGGen(Form("HalfFoil_%s",sBGGens[i].c_str() ),BGGens[i]);p3[i]->SetLineColor(2+i);}

  TCanvas *c =new TCanvas();
  p3[0]->Draw("HIST"); 
  for(int i=1; i<6; i++)p3[i]->Draw("HIST SAME"); 
  gPad->BuildLegend();
  lets_pause();

}


void NDK_analyzer( int i=0,int j=0, int k=4)
{


 //string fileRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20181123_v07_11_00_light_newtree.root";
 string fileNDK="/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root";
 string fileLONGRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root"; 
 //string fileSN="/eos/user/a/angalleg/TDR_simulations/data/SN/SN_20181210_v07_11_00_Livermore_light_newtree.root";


//string fileLONGRADFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/radiologicals/rad_20190313_v08_12_00_light_Foil.root";
 //string fileLONGRADHalfFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_HalfFoil.root";


 string fileRADFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_Foil.root";
 string fileRADHalfFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_HalfFoil.root";
 string fileNDKHalfFoil="/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_HalfFoil_NDK_barnuK+.root";
 string fileNDKFoil="/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_Foil_NDK_barnuK+.root";

// genClusters(fileNDK,1,5,0,fileNDK+"_SuperCluster");
// genClusters(fileSN,1,5,0,fileSN+"_SuperCluster");
// genClusters(fileLONGRAD,1,5,0,fileLONGRAD+"_SuperCluster");
// genClusters(fileRAD,1,5,0,fileRAD+"_SuperCluster");
/*
 switch(i)
 {
   case 0: genClusters(fileRADFoil,1,5,0,fileRADFoil+"_SuperCluster");         break;
   case 1: genClusters(fileRADHalfFoil,1,5,0,fileRADHalfFoil+"_SuperCluster"); break;
   case 2: genClusters(fileNDKFoil,1,5,0,fileNDKFoil+"_SuperCluster");         break;
   case 3: genClusters(fileNDKHalfFoil,1,5,0,fileNDKHalfFoil+"_SuperCluster"); break;
   case 4: genClusters(fileLONGRADFoil,1,5,0,fileNDKFoil+"_SuperCluster");         break;
   case 5: genClusters(fileLONGRADHalfFoil,1,5,0,fileNDKHalfFoil+"_SuperCluster"); break;
 }
*/

 string fileRADHalfFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_HalfFoil_20000.root"; 
 string fileRADFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_Foil_20000.root";
 std::vector<double> ddistance; ddistance.resize(13);
 for (int i=0;i<13;i++) ddistance[i]=50+i*50;

 std::vector<double> ddistance2; ddistance2.resize(15);ddistance2[0]=100;
 for (int i=1;i<15;i++) ddistance2[i]=ddistance2[i-1]+10;


 switch(i)
 {
 case 0:  MatchDistanceScan(fileNDKFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,1,ddistance);break;
 case 1:  MatchDistanceScan(fileNDKHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADHalfFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,1,ddistance);break;
 case 2:  ploteffvsdistanceNDK(fileNDKFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,150,1);break;
 case 3:  ploteffvsdistanceNDK(fileNDKHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADHalfFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,150,1);break;
 case 4:  MatchDistanceScan(fileNDK+"_SuperCluster_D2m5_W1us_MW1us.root",fileLONGRAD+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,1,ddistance);break;
 case 5:  ploteffvsdistanceNDK(fileNDK+"_SuperCluster_D2m5_W1us_MW1us.root",fileLONGRAD+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,250,1);break;
 case 6:  MatchDistanceScan("/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root_D2m5_W1us_MW1us.root",
                            "/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,0,ddistance);break;
 case 7:  MatchDistanceScan(fileNDKFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,1,ddistance2);break;
 case 8:  MatchDistanceScan(fileNDKHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADHalfFoil20k+"_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0,1,ddistance2);break;
}

//MatchDistanceScan(fileNDKHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0);
//MatchDistanceScan(fileNDKFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADFoil+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,0);


//ploteffvsdistanceNDK(fileNDKHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADHalfFoil+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,150);
//ploteffvsdistanceNDK(fileNDKFoil+"_SuperCluster_D2m5_W1us_MW1us.root",fileRADFoil+"_SuperCluster_D2m5_W1us_MW1us.root","D2m5_W1us_MW1us",2,150);

//plotMaxClusterPerDrift(fileNDK,fileNDKFoil,fileNDKHalfFoil,"_SuperCluster_D2m5_W1us_MW1us.root");
//plotRecoLightRatio(fileNDK,fileNDKFoil,fileNDKHalfFoil,"_SuperCluster_D2m5_W1us_MW1us.root");
//plotBGLightFreq(fileLONGRAD,fileRADFoil,fileRADHalfFoil,"_SuperCluster_D2m5_W1us_MW1us.root");
//plotRecoLightFreqByBGGen(fileLONGRAD,fileRADFoil,fileRADHalfFoil,"_SuperCluster_D2m5_W1us_MW1us.root");


//plotBGLightFreq(fileLONGRAD,fileRADFoil20k,fileRADHalfFoil20k,"_D2m5_W1us_MW1us.root");
plotRecoLightFreqByBGGen(fileLONGRAD,fileRADFoil20k,fileRADHalfFoil20k,"_D2m5_W1us_MW1us.root");



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
