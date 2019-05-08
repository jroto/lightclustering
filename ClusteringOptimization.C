#include"MyStyle.C"

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
void lets_pause (){
	TTimer * timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	timer->TurnOn();
	timer->Reset();
	std::cout << "q/Q to quit, other to continuee: ";
	char kkey;
	std::cin.get(kkey);
	if( kkey == 'q' || kkey == 'Q') 			gSystem->Exit(0);
	timer->TurnOff();
	delete timer;
}

#include"pmtmap.h"
#include"snana2_cluster.h"
#include"snana2_clusteringAnalyzerLight.h"
#include"snana2_radevt.h"
#include"snana2_radrun.h"
#include"snana2_clustering.h"
#include"snana2_matching.h"


  const char *cClusteringD[5]={"1.5","2.5","3.5","4.5","5.5"};
  const char *cClusteringW[5]={"50","100","200","400","800"};
  const char * cClusteringMW[9]={"63ns","125ns","250ns","500ns","1us","2us","4us","8us","16us"};


 const string sClusteringD[5]={"1m5","2m5","3m5","4m5","5m5"};
 const string sClusteringW[5]={"50ns","100ns","200ns","400ns","800ns"};
 const string sClusteringMW[9]={"63ns","125ns","250ns","500ns","1us","2us","4us","8us","16us"};
 const float fClusteringD[5]={150,250,350,450,550};//cm
 const float fClusteringW[5]={0.05,0.1,0.2,0.4,0.8};//us
 const float fClusteringMW[9]={0.063,0.125,0.25,0.5,1.,2.,4.,8.,16.};//us

 const string fFiles[9]={
                         "/eos/user/j/jsotooto/DUNE/lightana/data/radiologicals/Foils/20190313_v08_12_00_light_Foil_20000.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_Foil_NDK_barnuK+.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/20190313_v08_12_00_light_Foil_SN_Livermore.root",
                          "/eos/user/j/jsotooto/DUNE/lightana/data/radiologicals/HalfFoil/20190313_v08_12_00_light_HalfFoil_20000.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_HalfFoil_NDK_barnuK+.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/20190313_v08_12_00_light_HalfFoil_SN_Livermore.root",
                         "/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root", 
                         "/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root",
                         "/eos/user/j/jsotooto/DUNE/lightana/data/SN/SN_20181210_v07_11_00_Livermore_light_newtree.root"
};


  const string termination[2]={"_barridoBG_PEs.root", //for NDK
                               "_barridoBG_Hits.root"}; //for SN

void barridoBG(int i, int j, int k, string filerad, bool SN)
{


  string var="D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k];

  filerad=filerad+"_"+var+".root";

  TFile *ofile;
  if (SN) ofile = new TFile(Form("%s%s",filerad.c_str(),termination[1].c_str()),"RECREATE");
  else ofile = new TFile(Form("%s%s",filerad.c_str(),termination[0].c_str()),"RECREATE");

  TGraphErrors *tgBGR;

  snana::clusteringAnalyzerLight clusterconf(filerad,var);

  cout <<  " BGR " << i << " " << j << " " << k<< " " << filerad << endl;
  if (SN) tgBGR=clusterconf.BarridoInversoBGR_Hits({0.05,0.1,0.3,1},1.0);
  else tgBGR=clusterconf.BarridoInversoBGR({0.01,0.1,1,10},8.e-3);

  ofile->cd();
  tgBGR->Write(Form("tgBGR_Invers_%s",var.c_str()));
  ofile->Close(); 

}


void SetAxis(TH2D *h, int nx, int ny, const char *LabelsX[], const char *LabelsY[])
{
// function to make the histogram axis look nice

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

   for (int i=1;i<=nx;i++) h->GetXaxis()->SetBinLabel(i,LabelsX[i-1]);
   for (int i=1;i<=ny;i++) h->GetYaxis()->SetBinLabel(i,LabelsY[i-1]);
   //gPad->SetOptStat(0);
//   h->Draw("colz");lets_pause();

}


TH2D*  OptimizationHist(int bglevel, bool SN, int k, string fileRad, string fileSignal) //SN=true for SN sample, false for NDK
{
// this functions makes a scan of the 1st (i, D) and 2nd (j, W) clustering parameter, at a certain BG level, and at a certain value of the 3rd parameter (k, MW)
cout << " holi " <<endl;
  TFile *ifile[5][5]; TGraphErrors *tg[5][5];
  TMultiGraph *mg = new TMultiGraph();

cout << " holi " <<endl;
  TH2D *Efficiencies = new TH2D("Efficiencies","Efficiencies;Hit distance (m);Hit time distance (ns);Signal efficiency at 95% purity",5,0,5,5,0,5);
    double npe, NumClustersPerWindow;
  Efficiencies->SetStats(0);
  for(int i=0;i<5;i++)for(int j=0;j<5;j++)
  { 
    string var="D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k];
   
    string filesignal_aux=fileSignal+"_"+var+".root";
    string filerad_aux=fileRad+"_"+var+".root";

    if(SN) filerad_aux = fileRad+"_"+var + ".root"+termination[1];
    else filerad_aux = fileRad+"_"+var + ".root"+termination[0];
   
    ifile[i][j] = new TFile(filerad_aux.c_str(),"READ");

    if(!ifile[i][j]){ std::cout << filerad_aux << " - file does not exist "<< std::endl; gSystem->Exit(0);}

    tg[i][j] = (TGraphErrors*)ifile[i][j]->Get(Form("tgBGR_Invers_%s",var.c_str()));
    mg->Add(tg[i][j]);tg[i][j]->Draw(); 

    tg[i][j]->GetPoint(bglevel,NumClustersPerWindow,npe);
//lets_pause();
    snana::clusteringAnalyzerLight clusterconf(filesignal_aux,var);

    double eff;

    if(SN) eff = clusterconf.CalculateDetectionEfficiencySignal( npe);
    else eff = clusterconf.CalculateDetectionEfficiencySignalPECut( npe,1,450,600);
    
    Efficiencies->SetBinContent(Efficiencies->GetBin(i+1,j+1),eff);
    cout << NumClustersPerWindow <<  "Clusters per Window, " << npe << "PE or Hits, eff " << eff<< endl;
//    lets_pause();
    ifile[i][j]->Close();
    clusterconf.Close();

  }
  if(SN)  Efficiencies->SetTitle(Form("Efficiencies - BGR %.2f Hz - MW %s",NumClustersPerWindow,sClusteringMW[k].c_str()));
  else    Efficiencies->SetTitle(Form("%.1f BG Clusters per 8ms - MW %s",NumClustersPerWindow,sClusteringMW[k].c_str()));
  SetAxis(Efficiencies, 5, 5, cClusteringD, cClusteringW);
//  Efficiencies->Draw("colz");lets_pause();
  return Efficiencies; 
}
void effplotNDK(string fileRad, string fileSignal, string outname)
{
  TH2D *hh[4];
  
//this function gets the 2d histogram for all backgrounds and MW paramters, and store it in a pdf.
  TCanvas *c;

  for(int k=0;k<9;k++)
  {
    for(int i=0;i<4;i++) hh[i] = OptimizationHist(i, 0, k, fileRad, fileSignal);   
    c = new TCanvas();

    c->Divide(2,2);
    for(int i=0;i<4;i++)
    {
      c->cd(i+1); hh[i]->Draw("COLZ TEXT");hh[i]->SetMarkerSize(2.5);
      hh[i]->GetXaxis()->SetLabelSize(0.06);
      hh[i]->GetYaxis()->SetLabelSize(0.06);
    }
    c->Update(); c->Modified();
    if(k==0) c->Print(Form("%s(",outname.c_str()),"pdf");
    if(k==8)c->Print(Form("%s)",outname.c_str()),"pdf");
    if(k!=8&&k!=0) c->Print(Form("%s",outname.c_str()),"pdf");
  }
//  lets_pause();

}


void effplotSN( string fileRad, string fileSignal, string outname)
{
  TH2D *hh[4];

  TCanvas *c;

  for(int k=0;k<9;k++)
  {
    for(int i=0;i<4;i++) hh[i] = OptimizationHist(i, 1, k, fileRad, fileSignal);   
    c = new TCanvas();

    c->Divide(2,2);
    for(int i=0;i<4;i++)
    {
      c->cd(i+1); hh[i]->Draw("COLZ TEXT");hh[i]->SetMarkerSize(2.5);
      hh[i]->GetXaxis()->SetLabelSize(0.06);
      hh[i]->GetYaxis()->SetLabelSize(0.06);
    }
    c->Update(); c->Modified();    
    if(k==0) c->Print(Form("%s(",outname.c_str()),"pdf");
    if(k==8)c->Print(Form("%s)",outname.c_str()),"pdf");
    if(k!=8&&k!=0) c->Print(Form("%s",outname.c_str()),"pdf");
  }

  lets_pause();

}


void efficiencyTDR( string fileRad, string fileSignal, string outname, bool SN)
{
  TH2D *hh[4];

  TCanvas *c;
  TFile *ifile; TGraphErrors *tg;

  TH2D *Efficiencies = new TH2D("Efficiencies","Efficiencies; Maximum Time distance (us);Hit distance (m);Signal efficiency at 10BG Clusters per 8ms",9,0,9,5,0,5);


  for(int i=0;i<5;i++)for(int j=1;j<2;j++)for(int k=0;k<9;k++)
  {
     ifile=NULL;
     tg=NULL;
     double npe, NumClustersPerWindow;
    
     string var="D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k];

     cout << " holi " << i  << " " << j << " " << k << "   - " << var << endl;
   
     string filesignal_aux=fileSignal+"_"+var+".root";
     string filerad_aux=fileRad+"_"+var+".root";

     if(SN) filerad_aux = fileRad+"_"+var + ".root"+termination[1];
     else filerad_aux = fileRad+"_"+var + ".root"+termination[0];
   
     ifile = new TFile(filerad_aux.c_str(),"READ");

     if(!ifile){ std::cout << filerad_aux << " - file does not exist "<< std::endl; gSystem->Exit(0);}

     tg = (TGraphErrors*)ifile->Get(Form("tgBGR_Invers_%s",var.c_str()));

     snana::clusteringAnalyzerLight clusterconf(filesignal_aux,var);

     int bglevel=3;
       tg->GetPoint(bglevel,NumClustersPerWindow,npe);
       double eff;
       if(SN) eff = clusterconf.CalculateDetectionEfficiencySignal( npe);
       else eff = clusterconf.CalculateDetectionEfficiencySignalPECut( npe,1,450,600);

       Efficiencies->SetBinContent(Efficiencies->GetBin(k+1,i+1),eff);
       cout << NumClustersPerWindow <<  "Clusters per Window, " << npe << "PE or Hits, eff " << eff<< endl;
     clusterconf.Close();

  }

  SetAxis(Efficiencies, 9, 5, cClusteringMW, cClusteringD);
  Efficiencies->Draw("COLZ");
  lets_pause();
}


void effTree( string fileRad, string fileSignal, string outname, bool SN)
{
  TH2D *hh[4];

  TCanvas *c;
  TFile *ifile; TGraphErrors *tg;
  TFile *ofile = new TFile (Form("%s.root",outname.c_str()),"RECREATE");

  TNtuple *ntp = new TNtuple("EfficienciesOptimization","EfficienciesOptimization","D:W:MW:BG:Eff");

  for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<9;k++)
  {
     ifile=NULL;
     tg=NULL;
     double npe, NumClustersPerWindow;
    
     string var="D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k];

     cout << " holi " << i  << " " << j << " " << k << "   - " << var << endl;
   
     string filesignal_aux=fileSignal+"_"+var+".root";
     string filerad_aux=fileRad+"_"+var+".root";

     if(SN) filerad_aux = fileRad+"_"+var + ".root"+termination[1];
     else filerad_aux = fileRad+"_"+var + ".root"+termination[0];
   
     ifile = new TFile(filerad_aux.c_str(),"READ");

     if(!ifile){ std::cout << filerad_aux << " - file does not exist "<< std::endl; gSystem->Exit(0);}

     tg = (TGraphErrors*)ifile->Get(Form("tgBGR_Invers_%s",var.c_str()));

     snana::clusteringAnalyzerLight clusterconf(filesignal_aux,var);

     for (int bglevel=0; bglevel<4; bglevel++)
     {
       tg->GetPoint(bglevel,NumClustersPerWindow,npe);
       double eff;
       if(SN) eff = clusterconf.CalculateDetectionEfficiencySignal( npe);
       else eff = clusterconf.CalculateDetectionEfficiencySignalPE( npe);
    
       ntp->Fill(fClusteringD[i],fClusteringW[j],fClusteringMW[k],static_cast<float>(NumClustersPerWindow),static_cast<float>(eff));
       cout << NumClustersPerWindow <<  "Clusters per Window, " << npe << "PE or Hits, eff " << eff<< endl;
     }
//    lets_pause();
     ifile->Close();
     clusterconf.Close();

  }

  ofile->cd();
  ntp->Write();
  ofile->Close();
}
void checkfiles()
{
  int counter=0;
  gEnv:gEnv->SetValue("TFile.Recover", 0);
  ofstream ofile("missingfiles2.txt");

  std::map<int,std::string> fentries;

  std::vector<int> filestocheck={0,3,6};
  std::vector<int> counterPerL(6,0);

  int LtoTermination[6]={0,1,0,1,0,1};
  int LtoFile[6]={0,0,3,3,6,6};

    for(int i=0;i<5;i++) 
    for(int j=0;j<5;j++)
    for(int k=0;k<9;k++)
    for(int l=0;l<6;l++)
    {
       string var = "D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k];

       TFile *ff = new TFile(Form("%s_%s.root%s",fFiles[LtoFile[l]].c_str(),var.c_str(),termination[LtoTermination[l]].c_str()),"READ");
       if(!ff || ff->IsZombie()){ ofile << i << " " << j << " " << k<< " " << l << endl; counter++;counterPerL[l]++;cout << i << " " << j << " " << k<< " " << l << " MISSING " << endl;}
       else
       {

            TGraphErrors *tgBGR = (TGraphErrors*) ff->Get(Form("tgBGR_Invers_%s",var.c_str()));

            if (!tgBGR)
            {
               ofile << i << " " << j << " " << k<< " " << l << endl; cout << i << " " << j << " " << k<< " " << l << " CORRUPTED" << endl; counter++;counterPerL[l]++;
            }
            else 
            {
               if (tgBGR->GetN()!=4)
               {
                 ofile << i << " " << j << " " << k<< " " << l << endl; cout << i << " " << j << " " << k<< " " << l << " wrong size! " << tgBGR->GetN() << endl; counter++;counterPerL[l]++;
               }
               else { cout << i << " " << j << " " << k<< " " << l << " WE ARE GOOD!!!"<<endl;}
            } 
       }

       ff->Close();
  }
  ofile.close();
  cout << "DONE!!! Added " << counter << " files " << endl;
  for(int l=0;l<6;l++) cout << "\tAdded " << counterPerL[l] << " files for " << fFiles[LtoFile[l]] << endl;

}
void GenerarBarridos(int i, int j, int k, int l)
{
//Foils
  switch(l)
  {
//Foils
    case 0: barridoBG(i,j,k,fFiles[0],0); break;//NDK
    case 1: barridoBG(i,j,k,fFiles[0],1); break;//SN
//HalfFoil
    case 2: barridoBG(i,j,k,fFiles[3],0);break;
    case 3: barridoBG(i,j,k,fFiles[3],1);break;
//NoFoils
    case 4: barridoBG(i,j,k,fFiles[6],0);break;
    case 5: barridoBG(i,j,k,fFiles[6],1);break;
  }
}

void run_missing(int imin, int imax)
{

  ifstream ifile("missingfiles2.txt");

  int l, i, j, k;
  int counter=0;
  while( !ifile.eof())
  {
    counter++;
    ifile >> i >> j >> k >> l;
    cout << i << " " << j << " " << k << " " << l <<  endl;
     if(counter>=imin && counter<imax)GenerarBarridos(i,j,k,l);
   }  

  ifile.close();

}
void ClusteringOptimization( int i=0,int j=0, int k=0, int l=0)
{

// Functions to analyse the generated cluster files, and optimize the clustering configuration.

//First we established the 4th parameter (cluster size in PEs or hits) for the BG sample:
//This functions will generate new files in the folder curvas (be sure it does exist!), with a scan of cluster size at different background levels.
//Functions are different from SN to NDK because we use hits in SN and PEs in NDK.

//effplotNDK(fFiles[0], fFiles[1],"NDK_Foil_CUT.pdf");
//effplotNDK(fFiles[3], fFiles[4],"NDK_HalfFoil_CUT.pdf");

//effplotNDK(fFiles[6], fFiles[7], "NDK_NoFoil.pdf");


//effplotSN(fFiles[0], fFiles[2],"SNB_Foil.pdf");
//effplotSN(fFiles[3], fFiles[5],"SNB_HalfFoil.pdf");
//effplotSN(fFiles[6], fFiles[8], "SNB_NoFoil.pdf");
//   GenerarBarridos(i,j,k,l);

//effTtree(fFiles[0], fFiles[1],"NDK_Foil");
//effTree(fFiles[3], fFiles[4],"NDK_HalfFoil",0);

efficiencyTDR(fFiles[3], fFiles[4],"NDK_HalfFoil",0);

//checkfiles();

//run_missing(i,j);

//Once we have run the previous scan, we can create our 2d plot of optimization:
// reading the previous created files in folder curvas.

//OptimizationHist(0, 0, 0)->Draw("COLZ"); lets_pause();
//  effplotNDK();
//  effplotSN();


//Next functions gives the detection efficiency vs distance.
//    ploteffvsdistanceNDK(1,4,0);


}
