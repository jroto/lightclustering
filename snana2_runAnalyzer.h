namespace snana
{
class runAnalyzer
{
  public:

  string filein;
  bool rad;
  bool debug=false;
  radevt ev;
  TTree *t = NULL;
  string title;
  runAnalyzer(){}

  runAnalyzer(string filein, int DriftAxis, string name)
  {
    title=name;
    std::cout << "Processing " << name << ": " << filein << " - drift in " << DriftAxis;
    TFile *f =TFile::Open(filein.c_str(),"READ"); if (!f){ cout << "\tERROR OPENING FILE!" << std::endl; std::exit(0);}
    if(debug) f->ls();
    t = (TTree*)f->Get("snanagaushit/SNSimTree");
    if(debug)  t->Print();
    ev = radevt(DriftAxis);
    ev.SetBranchAdresses(t);
    std::cout << " - " << t->GetEntries() << " entries. " << std::endl;
  }

  TProfile* getNPEsVsDrift()
  {
    TProfile *tp = new TProfile(Form("%s_npe",title.c_str()),Form("%s;Drift Position (cm); # of PEs per event",title.c_str()),25,-600,600,0,1e7);
    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i);
      ev.FillAuxVars(); if(i%1000==0){std::cout << i << " out of "<< t->GetEntries() << " events processed "<< std::endl;}
      tp->Fill(ev.GetDriftPosition(),ev.GetNPEperEvent()); 
    }
    return tp;
  }

  TH1D* getNPEsPerPMT()
  {
    TH1D *tp = new TH1D(Form("%s_npe",title.c_str()),Form("%s;# of PEs per PMT per 1ms readout window; # of PMTs",title.c_str()),5000,0,5000);
    std::vector <double> npe;

    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); ev.GetPEperPMT(npe); if(i%1000==0){std::cout << i << " out of "<< t->GetEntries() << " events processed "<< std::endl;}
      for (int j=0;j<720;j++){ tp->Fill(npe[j]);} 
    }
    tp->Scale(1.0/t->GetEntries());
    return tp;
  }

  TH1D* getNPEs()
  {

    TH1D *tp = new TH1D(Form("%s_npe",title.c_str()),Form("%s;# of PEs per 1ms readout window; Counts",title.c_str()),5000,25000,80000);
    std::vector <double> npe;

    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); if(i%1000==0){std::cout << i << " out of "<< t->GetEntries() << " events processed "<< std::endl;}
      tp->Fill(ev.GetNPEperEvent());
    }
    tp->Scale(1.0/t->GetEntries());
    return tp;
  }

 
};

}
