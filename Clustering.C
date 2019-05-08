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

using namespace std;
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
}

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

void genClusters(string filename, int i, int j, int k)
{

//  string outfilename = filename+"_clustering_D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k]+".root";
  string outfilename=filename+"_D"+sClusteringD[i]+"_W"+sClusteringW[j]+"_MW"+sClusteringMW[k]+".root";
  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");fout->cd();
  snana::radrun_snana run;
  run=snana::radrun_snana(fout,filename,sClusteringD[i],fClusteringD[i],sClusteringW[j],fClusteringW[j],sClusteringMW[k],fClusteringMW[k],0); //1 for light simulation, 0 for heavy one
  run.Dump(fout);
  std::cout << "Dumped to " << outfilename << std::endl;
  fout->Close();
}

void checkfiles(string filename)
{


  gEnv:gEnv->SetValue("TFile.Recover", 0);
  ofstream ofile("missingfiles.txt");
  for(int i=0;i<5;i++) 
  for(int j=0;j<5;j++)
  for(int k=0;k<9;k++)
  {

     std::system(Form("ls %s_D%s_W%s_MW%s.root",filename.c_str(),sClusteringD[i].c_str(),sClusteringW[j].c_str(),sClusteringMW[k].c_str()));

     TFile *ff = new TFile(Form("%s_D%s_W%s_MW%s.root",filename.c_str(),sClusteringD[i].c_str(),sClusteringW[j].c_str(),sClusteringMW[k].c_str()),"READ");
     if(!ff || ff->IsZombie()){ ofile << i << " " << j << " " << k<< endl;}
     else
     {
       cout << "exists - ";

       TTree *tpl = (TTree*)ff->Get("clusteringanatree");
       cout << tpl << endl;
       if (!tpl){ ofile << i << " " << j << " " << k<< endl; cout << "CORRUPTED " << endl;}
       else { cout << "WE ARE GOOD!!!"<<endl;}

     }


   }  

  ofile.close();

}

void run_missing(string filename, int imin, int imax)
{

  ifstream ifile("missingfiles.txt");

  int f, i, j, k;
  int counter=0;
  while( !ifile.eof())
  {
    counter++;
    ifile >> i >> j >> k;
    cout << i << " " << j << " " << k << endl;
     if(counter>=imin && counter<imax)genClusters(filename,i,j,k);
   }  

  ifile.close();

}

void  genstartfile()
{
  ofstream ofile("ScanValues.txt");
  for(int i=0;i<5;i++) 
  for(int j=0;j<5;j++)
  for(int k=0;k<9;k++)
  ofile << i << " " << j << " " << k << endl;
  ofile.close();
}
void Clustering( int i=0,int j=0, int k=4)
{

string fileRAD="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190112_v07_11_00_light_newtree.root";
string fileNDK="/eos/user/a/angalleg/TDR_simulations/data/NDK/NDK_20181130_v07_11_00_light_ana_newtree.root";

 string fileRADFoil="/eos/user/a/angalleg/TDR_simulations/data/radiologicals/rad_20190313_v08_12_00_light_Foil.root";
 string fileNDKFoil="/eos/user/j/jsotooto/DUNE/lightana/data/NDK/20190313_v08_12_00_light_Foil_NDK_barnuK+.root";

// string fileRADFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_HalfFoil_20000.root";
// string fileRADFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_Foil_20000.root";
string fileRADFoil20k="/eos/user/j/jsotooto/DUNE/lightana/data/20190313_v08_12_00_light_Foil_20000.root";
  //we do a scan in (i,j,k) on the clustering paramters, to generate one file per clustering configuration.
  //The function creates one output file called filename_DX_WY_MWZ, in the same directory where the file lives.

   //for (int ii=0;ii<5;ii++)for (int jj=0;jj<5;jj++)
   
   if(i==0)genClusters(fileRAD,1,5,0); 
   else genClusters(fileNDK,1,5,0);
//   genClusters(fileRADFoil,i,j,k); //this function run the Clustering over filename with the ii,jj,kk parameters. 
//   genClusters(fileNDKFoil,i,j,k); //this function run the Clustering over filename with the ii,jj,kk parameters. 

//  genstartfile();
//  checkfiles(filename); // this function checks which cluster configurations are missing, and store it in missingfiles.txt

//   run_missing(filename,i,j); //this function reads from i-th to j-th entry of missingfiles.txt and process it with genClusters.

}
