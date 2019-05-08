namespace snana
{
class radrun_snana
{

 public:

 string filein;
 bool rad;
 bool debug=false;
 TProfile *tp;

 TTree *t;//one entry per event, to store all computed variables.

  string sClusteringD;
  string sClusteringW;
  string sClusteringMW;
  float fClusteringD;//cm
  float fClusteringW;//us
  float fClusteringMW;//us

 radrun_snana(){}
 radrun_snana(TFile *ofile, string ss,string sD, float fD, string sW, float fW, string sMW, float fMW, bool light=false, int levt=0): filein(ss), sClusteringD(sD), fClusteringD (fD), sClusteringW(sW), fClusteringW(fW), sClusteringMW (sMW), fClusteringMW (fMW)
 {
  std::cout << "Processing " << filein  << std::endl;//". split: "<< split << ", level: " << level << std::endl;

  TFile *f =TFile::Open(filein.c_str(),"READ");
  if(debug) f->ls();
  TTree *t2 = (TTree*)f->Get("snanagaushit/SNSimTree");
  if(debug)  t2->Print();


  std::cout << " Ruuning clustering for D" << sClusteringD << " W" << sClusteringW << " MW" << sClusteringMW <<  " light:" << light<< std::endl;
  radevt ev(sClusteringD, fClusteringD , sClusteringW, fClusteringW, sClusteringMW, fClusteringMW, light);
  ofile->cd();
  t = new TTree("clusteringanatree","clusteringanatree");//t->SetDirectory(0);
  ev.SetAnaTreeBranchAdresses(t);
  ev.SetBranchAdresses(t2);if(debug){std::cout << "branches adressed"<<std::endl;lets_pause();}

  int totev = t2->GetEntries();

//  int evini = level*(totev/split-totev%split);
//  int evfin = (level+1)*(totev/split-totev%split);
  int evfin=levt;
  if(levt==0 || totev<levt)evfin=totev;
  
  int evini=0;

  for (int i=evini;i<evfin;i++)
  {
    if(i%20==0){std::cout << i << " out of "<< evfin << " events processed "<< std::endl;}
    t2->GetEntry(i);
    ev.Process();  if(debug){std::cout << "event processed"<<std::endl;lets_pause();} 
    ev.FillAnaTree(t);
    if(debug){std::cout << "end of evt loop"<<std::endl;lets_pause();break;}
   // break;
  }

  if(debug){std::cout << "end of run creation"<<std::endl;lets_pause();}
//  Dump(fout);
  f->Close();

 }


 
 void Dump(TFile *fout)
 {

   std::cout << "dumpling to TFile"<<std::endl; 
//   fout=TFile(ofile.c_str(),"RECREATE");
   fout->cd();
   fout->WriteObject(&filein,"filein");
   std::cout << "writting ttree"<<std::endl; 
   t->Write();
   std::cout << "TTree written"<<std::endl; 

 }
};
}
