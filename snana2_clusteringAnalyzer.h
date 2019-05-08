

namespace snana{

class clusteringAnalyzer
{
  public:
  bool debug=false;

  int anaevent;
  float VertX;
  float VertY;
  float VertZ;
  float ENu;
  int _anaevent;


  string filein;
  string branchname;
  TFile *ifile=NULL;
  TTree *t = NULL;
  std::vector<Cluster_t> *ClustersVector = NULL;

  clusteringAnalyzer(){}

  clusteringAnalyzer(string ff, string bb)
  {
    filein=ff;
    branchname=bb;

    ifile = TFile::Open(filein.c_str(),"READ"); if (!ifile) std::cout << "ERROR with file not found!!!! "<< std::endl; std::cout <<  filein << " READED!!!! "<< std::endl;
    t = (TTree*)ifile->Get("clusteringanatree");

    SetBranchesAddresses();
  
  }

  void SetBranchesAddresses()
  {

     t->SetBranchAddress("anaevent"   	         ,&anaevent);
     t->SetBranchAddress("VertX"   	         ,&VertX);
     t->SetBranchAddress("VertY"   	         ,&VertY);
     t->SetBranchAddress("VertZ"   	         ,&VertZ);
     t->SetBranchAddress("E_nu"   	         ,&ENu);

     t->SetBranchAddress(Form("Clusters_%s",branchname.c_str()),&ClustersVector);
  }

  TProfile *GetPEvsDrift(int axis)
  {
    TProfile *tp = new TProfile("tp","Maximum Cluster per Event;Drift direction (cm); # of PEs",30,-600,600,0,30000);
    float nPE=0;
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      float pos[3]={VertX,VertY,VertZ};
      if(ClustersVector->size()>0)
      {
          nPE = ClustersVector->at(0).PEs();
      }
      else
      {
          nPE=0.0;
      }
      tp->Fill(pos[axis],nPE);
    }
   return tp;

  }

  TProfile *RecoLightRatio(int axis)
  {
    TProfile *tp = new TProfile("tp","Maximum Cluster per Event;Drift direction (cm); # of PEs",30,-600,600,0,1);
    float nPE=0;
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      float pos[3]={VertX,VertY,VertZ};
      nPE=0.0;
      if(ClustersVector->size()>0)
      {
          for (int j=0;j<ClustersVector->size();j++) nPE+=ClustersVector->at(j).PEs();
          nPE = ClustersVector->at(0).PEs()/nPE;
      }
      tp->Fill(pos[axis],nPE);
    }
   return tp;

  }

  TH1F *RecoLightFreqByBGGen(string title,int BGGen)
  {
    std::cout << " producing RecoLightFreqByBGGen light freq plot for " << filein <<  " " <<title << " " << BGGen << std::endl;
    TH1F *th = new TH1F(Form("th%s",title.c_str()),Form("%s;# of PEs per cluster; Frequency (Hz)",title.c_str()),3000,0,3000);
    th->Sumw2();std::cout  << t->GetEntries() << "entries"<<  std::endl;
    int nentries=t->GetEntries();
    t->Draw(Form("Clusters_D2m5_W1us_MW1us.nPEs>>th%s",title.c_str()),Form("anaevent<%i&&Clusters_D2m5_W1us_MW1us.Generator==%i",nentries,BGGen));
    double cumulator=0;
    TH1F *th2 = new TH1F(Form("th2%s",title.c_str()),Form("%s;Threshold (# of PEs per Cluster); Frequency (Hz)",title.c_str()),3000,0,3000);
    for(int i=th->GetNbinsX();i>0;i--)
    {
      cumulator+=th->GetBinContent(i);
      th2->SetBinContent(i,cumulator);
    }     
    th2->Scale(1.0/(1e-3*nentries));
    return th2;

  }


  TH1F *RecoLightFreq(string title)
  {
    std::cout << " producing reco light freq plot for " << filein << std::endl;
    TH1F *th = new TH1F(Form("th%s",title.c_str()),Form("%s;# of PEs per cluster; Frequency (Hz)",title.c_str()),3000,0,3000);
    th->Sumw2();std::cout  << t->GetEntries() << "entries"<<  std::endl;
    int nentries=t->GetEntries();
    t->Draw(Form("Clusters_D2m5_W1us_MW1us.nPEs>>th%s",title.c_str()),Form("anaevent<%i",nentries));

/*
    int nClusters;
    for(int i=1;i<=th->GetNbinsX();i++)
    {std::cout << "entry " << i << "out of " << th->GetNbinsX() <<  std::endl;
      nClusters = t->Draw(Form("Clusters_D2m5_W1us_MW1us.nPEs"),Form("anaevent<4000&&Clusters_D2m5_W1us_MW1us.nPEs>%.2f",th->GetBinCenter(i)));
      th->SetBinContent(i,nClusters);
      th->SetBinError(i,TMath::Sqrt(nClusters));
    }
*/
/*
    float nPE=0;
    int entries=4000;
    if(t->GetEntries()<entries) entries=t->GetEntries();
    for (int i=0;i<entries;i++)
    {
      t->GetEntry(i);if(i%500==0)std::cout << "entry " << i << "out of " << t->GetEntries() <<  std::endl;
      nClusters=ClustersVector->size();
      if(ClustersVector->size()>0)
      {
          for (int j=0;j<ClustersVector->size();j++) tp->Fill(ClustersVector->at(j).PEs());
      }
    }
*/   double cumulator=0;
   TH1F *th2 = new TH1F(Form("th2%s",title.c_str()),Form("%s;Threshold (# of PEs per Cluster); Frequency (Hz)",title.c_str()),3000,0,3000);
   for(int i=th->GetNbinsX();i>0;i--)
   {
     cumulator+=th->GetBinContent(i);
     th2->SetBinContent(i,cumulator);
   }     
   th2->Scale(1.0/(1e-3*nentries));
   return th2;

  }

/*
  TH1F ClusterPurity()
  {
   if(debug)  t->Print();
    
   TH1F hist("hist";"hist",100,0,100);

    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      if(ClustersVector->size()>0 && ClustersVector->at(0).Hits()>=hitthreshold)
      {
         CounterPerEventRun++;
         for(unsigned int j=0;j<ClustersVector->size();j++)
         {
           if(ClustersVector->at(j).Hits()>=hitthreshold) {CounterAllRun++;CounterAllPerEventRun++;}
           else {break;}
         }  
      }
      if(CounterAllPerEventRun>1)Counter_multiclusterRun++;
      if(CounterAllPerEventRun>1)CounterAll_multiclusterRun+=CounterAllPerEventRun;
    }
    NumberOfEventsRun=t->GetEntries();

    return hist;
  }
*/

  std::map<int,TH1D> NumberOfPEsPerClusterPerGenerator(int hitthreshold)
  {
   if(debug)  t->Print();
    
   std::map<int,TH1D> mymap;
   for (int i=0; i<10; i++) mymap[i]=TH1D(Form("%i_back",i), ";PE per Cluster;Event", 50000, 0, 50000);

    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      if(ClustersVector->size()>0 && ClustersVector->at(0).Hits()>=hitthreshold)
      {
         for(unsigned int j=0;j<ClustersVector->size();j++)
         {
           if(ClustersVector->at(j).Hits()>=hitthreshold)
           {
              mymap[ ClustersVector->at(j).GetGenType()].Fill(ClustersVector->at(j).PEs());
              if (debug) {cout << " entro "  <<ClustersVector->at(j).GetGenType() << " " <<  ClustersVector->at(j).PEs() << endl; mymap[ ClustersVector->at(j).GetGenType()].Draw("hist");lets_pause();}
           }
           else {break;}
         }  
      }
    }

    return mymap;
  }

  std::map<int,TH1D> NumberOfPEsPerClusterPerPDG(int hitthreshold)
  {
   if(debug)  t->Print();
    
   std::map<int,TH1D> mymap;
   for (int i=0; i<10; i++) mymap[i]=TH1D(Form("%i_back",i), ";PE per Cluster;Event", 50000, 0, 50000);

    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      if(ClustersVector->size()>0 && ClustersVector->at(0).Hits()>=hitthreshold)
      {
         for(unsigned int j=0;j<ClustersVector->size();j++)
         {
           if(ClustersVector->at(j).Hits()>=hitthreshold)
           {
              mymap[ ClustersVector->at(j).GetGenType()].Fill(ClustersVector->at(j).PEs());
              if (debug) {cout << " entro "  <<ClustersVector->at(j).GetGenType() << " " <<  ClustersVector->at(j).PEs() << endl; mymap[ ClustersVector->at(j).GetGenType()].Draw("hist");lets_pause();}
           }
           else {break;}
         }  
      }
    }

    return mymap;
  }

 };

}
