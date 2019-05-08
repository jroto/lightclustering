

namespace snana{

class clustering
{
  public:
  bool debug=false;
  float background;
  float efficiency;

  std::vector<int> nevents;
  int nClusters;
  int anaevent;
  float VertX;
  float VertY;
  float VertZ;
  float ENu;
  float ENu_Lep;
  int _anaevent;


  string filein;
  string branchname;
  TFile *ifile=NULL;
  TTree *t = NULL;
  std::vector<Cluster_t> *ClustersVector = NULL;

  TFile *ClusterFile = NULL;
  TTree *ClusterTree = NULL;
  TH2D *ClusteringHist=NULL;
  TH3D *ClusteringHist3D=NULL;
  int NReadOutWindows;
  int MaxBGCluster=175;
  bool DensitiesLoaded=false;
  bool DriftY=true;

  clustering(){}
  clustering(string ff, string bb, bool dy)
  {
    filein=ff;
    branchname=bb;

    ifile = TFile::Open(filein.c_str(),"READ");
    t = (TTree*)ifile->Get("clusteringanatree");
//    t->SetDirectory(0);

//    cout << "tree loaded " << t << endl;
    SetBranchesAddresses();
//     t->GetEntry(2); cout << ClustersVector->at(0).Hits() <<  endl;
    DriftY = dy;
  }
 
  void SetBranchesAddresses()
  {

//     t->Print();
     t->SetBranchAddress("anaevent"   	         ,&anaevent);
     t->SetBranchAddress("VertX"   	         ,&VertX);
     t->SetBranchAddress("VertY"   	         ,&VertY);
     t->SetBranchAddress("VertZ"   	         ,&VertZ);
     t->SetBranchAddress("E_nu"   	         ,&ENu);
     t->SetBranchAddress("E_nu_Lep"   	         ,&ENu_Lep);

//     gInterpreter->GenerateDictionary("Cluster_t", "Cluster.h");

     t->SetBranchAddress(Form("Clusters_%s",branchname.c_str()),&ClustersVector);
//     cout << "Branches addressed"<<endl;


  }
  double getVertX()
  {return VertX;}
  void GetThisEntry(int i)
  {
    t->GetEntry(i);
  }

  int CounterPerEventRun=0;
  int CounterAllRun=0;
  int Counter_multiclusterRun=0;
  int CounterAll_multiclusterRun=0;
  int CounterAllPerEventRun;
  int NumberOfEventsRun;

  int CounterPerEvent=0;
  int CounterAll=0;
  int Counter_multicluster=0;
  int CounterAll_multicluster=0;
  int CounterAllPerEvent;
  int NumberOfEvents=0;

  double multiplicity=0; //probability of having more than 1 clusters once the event is detected.
  double avmulticluster=0; // average number of clusters in events with more than 1 cluster.

  double CalculateDetectionEfficiencyChain(float hitthreshold)
  {
//    std::cout << Form("ls %s",filein.c_str()) <<std::endl;
//    std::cout << Form("rm out.txt;ls %s>out.txt",filein.c_str()) <<std::endl;
//    std::system(Form("ls %s",filein.c_str()));
//    std::system(Form("rm out.txt;ls %s>out.txt",filein.c_str()));
//    ifstream ifile;
//    ifile.open("out.txt");
    string filename; //getline(ifile,filename);
    CounterPerEvent=0;
    CounterAll=0;
    Counter_multicluster=0; // number of events with more than one cluster found
    CounterAll_multicluster=0; // number of cluster per event in events whith more than one cluster found.
    NumberOfEvents=0;
    int cc=0;
    
//    while( !ifile.eof())
//    {
//      ifile>>filename;//
//      std::cout <<"adding " << filename << "   " << cc << " " << ifile.eof() << " " << hitthreshold << std::endl;//lets_pause();cc++;
      filename=filein;
      double aux = CalculateDetectionEfficiencySignal(hitthreshold);
      CounterPerEvent+=CounterPerEventRun;
      CounterAll+=CounterAllRun;
      Counter_multicluster+=Counter_multiclusterRun; // number of events with more than one cluster found
      CounterAll_multicluster+=CounterAll_multiclusterRun; // number of cluster per event in events whith more than one cluster found.
      NumberOfEvents+=NumberOfEventsRun;

//    }
    //lets_pause();

    return 1.0*CounterPerEvent/NumberOfEvents;
  }

   TPaveText * GetClusterStatistics(int i, float hitthreshold)
  {

    bool debug=false;

    if(debug)  t->Print();

    std::vector<double> clusters;
    double totalPEs=0;
    t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
    for(unsigned int j=0;j<ClustersVector->size();j++)
    {
      totalPEs+=ClustersVector->at(j).PEs();
      if(ClustersVector->at(j).Hits()>=hitthreshold) {clusters.push_back(ClustersVector->at(j).PEs());}
//      else {break;}
    }  

   float width=0.055*(2+clusters.size());
   TPaveText *pt = new TPaveText(.5,0.89-width,.90,.89,"NDC"); pt->SetLineColor(0);


//  pt->SetFillColorAlpha(0,0);
   pt->SetFillStyle(0);
   pt->SetLineStyle(0);
//pt->SetFillAtributes(4000);
   pt->AddText(Form("%.0f PEs in the event",totalPEs));
   pt->AddText(Form("%zu Clusters found above %.0fPE",clusters.size(),hitthreshold));
   for(unsigned int j=0;j<clusters.size();j++) pt->AddText(Form("\tCluster of %.0fPE",clusters[j]));
//   pt->UseCurrentStyle();
   pt->SetFillColor(0);
//    f->Close();

   return pt;
//    return clusters;
  }


  double CalculateDetectionEfficiencySignal(float hitthreshold)
  {
    
    if(debug)  t->Print();
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
//    std::cout <<"Analysing " << t->GetEntries() << std::endl;
    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //std::cout <<"event " << i << std::endl;
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

    std::cout << "Detected " << CounterPerEventRun << " out of " << NumberOfEventsRun <<std::endl;
    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }
  double CalculateDetectionEfficiencySignalPE(double hitthreshold)
  {

    if(debug)  t->Print();
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
//    std::cout <<"Analysing " << t->GetEntries() << std::endl;
    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //std::cout <<"event " << i << std::endl;
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
    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }

  double CalculateDetectionEfficiencySignalPECut(double hitthreshold, float Xmin, float Xmax)
  {
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
    NumberOfEventsRun=0;

    std::cout <<"Analysing " << filein << " between " << Xmin << " " << Xmax <<  std::endl;

    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //std::cout <<"event " << i << std::endl;

      t->GetEntry(i);
//      std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  " VertX" << VertX << std::endl; lets_pause();

      if(VertX<=Xmin || VertX>Xmax) continue;
      else NumberOfEventsRun++;

      if(ClustersVector->size()>0 && ClustersVector->at(0).Hits()>=hitthreshold  )
      {
         CounterPerEventRun++;
      }
      if(CounterAllPerEventRun>1)Counter_multiclusterRun++;
      if(CounterAllPerEventRun>1)CounterAll_multiclusterRun+=CounterAllPerEventRun;
    }
    
    std::cout <<"Detected clusters in " << CounterPerEventRun << " out of  " << NumberOfEventsRun << " above " << hitthreshold << "PEs. "<<   std::endl;

    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }



  double GetMultiplicity(float hitthreshold)
  {
    double efficiency = CalculateDetectionEfficiencySignalPE(hitthreshold);//percentage of events that triggered one cluster.
    std::cout << "multiplicity " << Counter_multiclusterRun  << " " << CounterPerEventRun << std::endl; lets_pause();
    double multip = 1.0*Counter_multiclusterRun / CounterPerEventRun;
    return multip;
  }

  double CalculateBackGroundRate(float hitthreshold)
  {
    double efficiency = CalculateDetectionEfficiencyChain(hitthreshold);//percentage of events that triggered one cluster.

    double backgroundrate = 1.0*CounterAll / (1.e-3*NumberOfEvents ) ;
    return backgroundrate;
  }

  TGraphErrors* BarridoDE(Barrido b, TGraphErrors* tgmultiplicity, TGraphErrors* tgavmulticluster)
  {
     std::cout << "barrido00 DE  " << b.name << std::endl;
     TGraphErrors* tg = new TGraphErrors();
     for(double i=b.min; i<b.max; i+=b.step)
     {
        std::cout << "barrido DE  " << b.name << " " << i << " ";//lets_pause();
        double efficiency = CalculateDetectionEfficiencyChain(i);//percentage of events that triggered one cluster.
        tg->SetPoint(tg->GetN(),i,efficiency);
        tg->SetPointError(tg->GetN()-1,0,efficiency/TMath::Sqrt(NumberOfEvents));

        multiplicity=(double)Counter_multicluster/CounterPerEvent; //probability of having more than 1 clusters once the event is detected.
        avmulticluster=(double)CounterAll_multicluster/Counter_multicluster; // average number of clusters in events with more than 1 cluster.

        tgmultiplicity->SetPoint(tgmultiplicity->GetN(),i,multiplicity);
        tgmultiplicity->SetPointError(tgmultiplicity->GetN()-1,0.0,multiplicity*TMath::Sqrt(1.0/Counter_multicluster+1.0/CounterPerEvent));

        tgavmulticluster->SetPoint(tgavmulticluster->GetN(),i,avmulticluster);
        tgavmulticluster->SetPointError(tgavmulticluster->GetN()-1,0.0,avmulticluster*TMath::Sqrt(1.0/CounterAll_multicluster+1.0/Counter_multicluster));

        std::cout << efficiency << " " << efficiency/TMath::Sqrt(NumberOfEvents) <<  " " << multiplicity <<" " << avmulticluster << std::endl;
     }

     tg->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tg->GetYaxis()->SetTitle("Detection Efficiency");

     tgmultiplicity->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tgmultiplicity->GetYaxis()->SetTitle("Detection Multiplicity");

     tgavmulticluster->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tgavmulticluster->GetYaxis()->SetTitle("Average cluster per MultiClusterEvent");

     return tg;
  }

  TGraphErrors* BarridoDEtest(Barrido b, TGraphErrors* tgmultiplicity, TGraphErrors* tgavmulticluster)
  {
     std::cout << "barrido00 DE  " << b.name << std::endl;
     TGraphErrors* tg = new TGraphErrors();
     for(double i=b.min; i<b.max; i+=b.step)
     {
        std::cout << "barrido DE  " << b.name << " " << i << " ";//lets_pause();
        double efficiency = CalculateDetectionEfficiencySignal(i);//percentage of events that triggered one cluster.
        tg->SetPoint(tg->GetN(),i,efficiency);
        tg->SetPointError(tg->GetN()-1,0,efficiency/TMath::Sqrt(NumberOfEventsRun));

        multiplicity=(double)Counter_multiclusterRun/CounterPerEventRun; //probability of having more than 1 clusters once the event is detected.
        avmulticluster=(double)CounterAll_multiclusterRun/Counter_multiclusterRun; // average number of clusters in events with more than 1 cluster.

        tgmultiplicity->SetPoint(tgmultiplicity->GetN(),i,multiplicity);
        tgmultiplicity->SetPointError(tgmultiplicity->GetN()-1,0.0,multiplicity*TMath::Sqrt(1.0/Counter_multiclusterRun+1.0/CounterPerEventRun));

        tgavmulticluster->SetPoint(tgavmulticluster->GetN(),i,avmulticluster);
        tgavmulticluster->SetPointError(tgavmulticluster->GetN()-1,0.0,avmulticluster*TMath::Sqrt(1.0/CounterAll_multiclusterRun+1.0/Counter_multiclusterRun));

        std::cout << efficiency << " " << efficiency/TMath::Sqrt(NumberOfEventsRun) <<  " " << multiplicity <<" " << avmulticluster << std::endl;
     }

     tg->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tg->GetYaxis()->SetTitle("Detection Efficiency");

     tgmultiplicity->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tgmultiplicity->GetYaxis()->SetTitle("Detection Multiplicity");

     tgavmulticluster->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tgavmulticluster->GetYaxis()->SetTitle("Average cluster per MultiClusterEvent");

     return tg;
  }


  TGraphErrors* BarridoBGR(Barrido b, TGraphErrors* tgBGevtPer2DrifWindow)
  {
     TGraphErrors* tg = new TGraphErrors();
     for(double i=b.min; i<b.max; i+=b.step)
     {
        std::cout << "barrido BGR " << b.name << " " << i << std::endl;//lets_pause();
        double BGR = CalculateBackGroundRate(i);//percentage of events that triggered one cluster.
        tg->SetPoint(tg->GetN(),i,BGR);
        tg->SetPointError(tg->GetN()-1,0,TMath::Sqrt(CounterAll)/(1.e-3*NumberOfEvents));

        tgBGevtPer2DrifWindow->SetPoint(tgBGevtPer2DrifWindow->GetN(),i,CounterAll*15.e-3/(1.e-3*NumberOfEvents));
        tgBGevtPer2DrifWindow->SetPointError(tgBGevtPer2DrifWindow->GetN()-1,0,TMath::Sqrt(CounterAll)*15.e-3/(1.e-3*NumberOfEvents));

        std::cout << " " << BGR << " " << TMath::Sqrt(CounterAll)/(1.e-3*NumberOfEvents) << "Hz  --- " << CounterAll <<" clusters found"<< std::endl;
        if (CounterAll==0)break;
     }
     tg->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tg->GetYaxis()->SetTitle("Background Rate (Hz)");
     tg->SetLineColor(2);

     tgBGevtPer2DrifWindow->GetXaxis()->SetTitle(Form("%s",b.name.c_str()));
     tgBGevtPer2DrifWindow->GetYaxis()->SetTitle("# clusters per 2x Drift Window (Counts/15ms)");
     tgBGevtPer2DrifWindow->SetLineColor(3);
    
     return tg;
  }

  void insert( std::vector<int> &cont, int value ) {
    std::vector<int>::iterator it = std::lower_bound( cont.begin(), cont.end(), value, std::greater<int>() ); // find proper position in descending order
    cont.insert( it, value ); // insert before iterator it
  }

  void insert( std::vector<double> &cont, double value ) {
    std::vector<double>::iterator it = std::lower_bound( cont.begin(), cont.end(), value, std::greater<float>() ); // find proper position in descending order
    cont.insert( it, value ); // insert before iterator it
  }

  double GetMinimumClusterAtLevel(float NumOfEvtsPerWindow, float WindowSize)
  {

    float SampleSize=t->GetEntries()*1e-3; //seconds
    unsigned int NumOfEventsInSample=(unsigned int)(NumOfEvtsPerWindow*SampleSize/WindowSize);

    std::cout << "looking for " << NumOfEventsInSample << " events in the sample of " << SampleSize<< "s."<<   std::endl;
    std::cout << "To get " << NumOfEvtsPerWindow << " events per " << WindowSize<< "s window."<<   std::endl;
    std::vector<double> AllClusters; AllClusters.clear();
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0)->PEs() << " number of clusters: " << ClustersVectorPE->size() <<  std::endl; lets_pause();

      if(ClustersVector->size()>0)
      {
         for(unsigned int j=0;j<ClustersVector->size();j++)
         {
           if(AllClusters.size()<NumOfEventsInSample)
           {
              insert( AllClusters,ClustersVector->at(j).PEs());
           }
           else
           {
             if(ClustersVector->at(j).PEs()<=AllClusters[AllClusters.size()-1]) break;
             if(ClustersVector->at(j).PEs()>AllClusters[AllClusters.size()-1])
             {
               insert(AllClusters,ClustersVector->at(j).PEs()); //std::cout <<"adding" <<std::endl;
               AllClusters.pop_back();
             }
           }
//           for(int k=0;k<AllClusters.size();k++) std::cout << AllClusters[k] << " "; std::cout << std::endl;lets_pause();  
         }
       }
    }

//    for (int i=0;i<AllClusters.size();i++) {std::cout << AllClusters[i] << " ";} std::cout << std::endl ; lets_pause();

    double aux =AllClusters[AllClusters.size()-1];
    return aux;
  }

  int GetMinimumClusterAtLevelHit(float NumOfEvtsPerWindow, float WindowSize)
  {

    float SampleSize=t->GetEntries()*1e-3; //seconds
    unsigned int NumOfEventsInSample=(int)(NumOfEvtsPerWindow*SampleSize/WindowSize);

    std::cout << "looking for " << NumOfEventsInSample << " events in the sample of " << SampleSize<< "s."<<   std::endl;
    std::cout << "To get " << NumOfEvtsPerWindow << " events per " << WindowSize<< "s window."<<   std::endl;
    std::vector<int> AllClusters; AllClusters.clear();
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();

      if(ClustersVector->size()>0)
      {
//         std::cout << "entro" << std::endl;for(int k=0;k<AllClusters.size();k++) std::cout << AllClusters[k] << " "; std::cout << std::endl;
         for(unsigned int j=0;j<ClustersVector->size();j++)
         {
           if(AllClusters.size()<NumOfEventsInSample)
           {
              insert( AllClusters,ClustersVector->at(j).Hits());
           }
           else
           {

//             std::cout << "entro2" << std::endl;for(int k=0;k<AllClusters.size();k++) std::cout << AllClusters[k] << " "; std::cout << std::endl;
//             std::cout << ClustersVector->at(j)->Hits() << " " << AllClusters[AllClusters.size()-1] <<std::endl;lets_pause();
             if(ClustersVector->at(j).Hits()<=AllClusters[AllClusters.size()-1]) break;
             if(ClustersVector->at(j).Hits()>AllClusters[AllClusters.size()-1])
             {
               insert(AllClusters,ClustersVector->at(j).Hits()); //std::cout <<"adding" <<std::endl;
               AllClusters.pop_back();
             }
//std::cout << "entro3" << std::endl;for(int k=0;k<AllClusters.size();k++) std::cout << AllClusters[k] << " "; std::cout << std::endl;
           }
//           for(int k=0;k<AllClusters.size();k++) std::cout << AllClusters[k] << " "; std::cout << std::endl;  
         }
       }
    }

//    for (int i=0;i<AllClusters.size();i++) {std::cout << AllClusters[i] << " ";} std::cout << std::endl ; lets_pause();

    int aux =AllClusters[AllClusters.size()-1];
    return aux;
  }

  TGraphErrors* BarridoInversoBGR(std::vector<float> NumOfEvtsPerWindow, float Window)
  {
     TGraphErrors* tg = new TGraphErrors();
     for(unsigned int i=0; i<NumOfEvtsPerWindow.size(); i++)
     {
        std::cout << "barrido BGR " << branchname << " " << i << std::endl;//lets_pause();
        double NumOfPEs;
        do
        {
          NumOfPEs = GetMinimumClusterAtLevel(NumOfEvtsPerWindow[i], Window);
        }while(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000);
        std::cout << "\t" << NumOfPEs << "PEs per " << NumOfEvtsPerWindow[i] << "evt/"<<Window<<"s" <<std::endl<<std::endl; //if(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000) lets_pause();
        
        tg->SetPoint(tg->GetN(),NumOfEvtsPerWindow[i],NumOfPEs);
        tg->SetPointError(tg->GetN()-1,0,0);
     }
     tg->GetXaxis()->SetTitle(Form("%s",branchname.c_str()));
     tg->GetYaxis()->SetTitle("#PEs)");
     tg->SetLineColor(2);
    
     return tg;
  }

  TGraphErrors* BarridoInversoBGR_Hits(std::vector<float> NumOfEvtsPerWindow, float Window)
  {
     TGraphErrors* tg = new TGraphErrors();
     for(unsigned int i=0; i<NumOfEvtsPerWindow.size(); i++)
     {
        std::cout << "BarridoInversoBGR_Hits " << branchname << " " << i  << " " << NumOfEvtsPerWindow[i] << std::endl;//lets_pause();
        int NumOfHits;
        do
        {
          NumOfHits = GetMinimumClusterAtLevelHit(NumOfEvtsPerWindow[i], Window);
        }while(!NumOfHits);// || NumOfHits>5000 || NumOfHits<-5000);
        std::cout << "\t" << NumOfHits << "Hits per " << NumOfEvtsPerWindow[i] << "evt/"<<Window<<"s" <<std::endl<<std::endl; //if(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000) lets_pause();
        
        tg->SetPoint(tg->GetN(),NumOfEvtsPerWindow[i],(float)NumOfHits);
        tg->SetPointError(tg->GetN()-1,0,0);
     }
     tg->SetTitle(Form("%s",branchname.c_str()));
     tg->GetXaxis()->SetTitle(Form("# Events in %.2f s",Window));
     tg->GetYaxis()->SetTitle("#Hits");
     tg->SetLineColor(2);
    
     return tg;
  }

  double GetNumberOfClustersAt(float TruePos[3],double distance,double PEThres)
  {
    double counter=0; double dist2 = distance*distance;
    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      for(unsigned int j=0;j<ClustersVector->size();j++)
      {
        if(ClustersVector->at(j).PEs()<PEThres) break;
        //double dx=TruePos[0]-ClustersVector->at(j).pos[0];
        double dy=TruePos[1]-ClustersVector->at(j).pos[1];
        double dz=TruePos[2]-ClustersVector->at(j).pos[2];
        //cout << TruePos[0] << " " << ClustersVector->at(j).pos[0]<< endl;
        //cout << dy*dy+dz*dz << " vs " << dist2 << " && " << ClustersVector->at(j).PEs() << " vs " << PEThres << endl;
        if(dy*dy+dz*dz<dist2 && ClustersVector->at(j).PEs()>PEThres) counter++;

      }
    }
    return counter;
  }
  double GetNumberOfClustersAtBetween(float TruePos[3],double distance,double PEThres, int i_start, int i_end)
  {
    double counter=0; double dist2 = distance*distance;
    for (int i=i_start;i<i_end;i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << " number of clusters: " << ClustersVector->size() << ", first cluster: "  <<   std::endl; ClustersVector->at(0).Print();lets_pause();
      for(unsigned int j=0;j<ClustersVector->size();j++)
      {
        if(ClustersVector->at(j).PEs()<PEThres) break;
        //double dx=TruePos[0]-ClustersVector->at(j).pos[0];
        double dy=TruePos[1]-ClustersVector->at(j).pos[1];
        if(DriftY) dy=TruePos[0]-ClustersVector->at(j).pos[0];
        double dz=TruePos[2]-ClustersVector->at(j).pos[2];
        //cout << TruePos[0] << " " << ClustersVector->at(j).pos[0]<< endl;
        //cout << dy*dy+dz*dz << " vs " << dist2 << " && " << ClustersVector->at(j).PEs() << " vs " << PEThres << endl;
        if(dy*dy+dz*dz<dist2 && ClustersVector->at(j).PEs()>=PEThres) counter++;
      }
    }
    return counter;
  }

  int GetOriginOfLargerClusterAtBetween(float TruePos[3],double distance,double PEThres, int i_start, int i_end)
  {
    std::vector<int> ClusterOrigin, ClusterSize;
    double counter=0; double dist2 = distance*distance;
    for (int i=i_start;i<i_end;i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << " number of clusters: " << ClustersVector->size() << ", first cluster: "  <<   std::endl; ClustersVector->at(0).Print();lets_pause();
      for(unsigned int j=0;j<ClustersVector->size();j++)
      {
        if(ClustersVector->at(j).PEs()<PEThres) break;
        //double dx=TruePos[0]-ClustersVector->at(j).pos[0];
        double dy=TruePos[1]-ClustersVector->at(j).pos[1];
        if(DriftY) dy=TruePos[0]-ClustersVector->at(j).pos[0];
        double dz=TruePos[2]-ClustersVector->at(j).pos[2];
        //cout << TruePos[0] << " " << ClustersVector->at(j).pos[0]<< endl;
        //cout << dy*dy+dz*dz << " vs " << dist2 << " && " << ClustersVector->at(j).PEs() << " vs " << PEThres << endl;
        if(dy*dy+dz*dz<dist2 && ClustersVector->at(j).PEs()>=PEThres) { ClusterSize.push_back( ClustersVector->at(j).PEs()); ClusterOrigin.push_back( ClustersVector->at(j).GetGenType()); break;}
      }
    }
    int j=0;
    for (int i=0; i<ClusterOrigin.size(); i++)
    {
      if(ClusterSize[i]>ClusterSize[j]) j=i;
    }
    return ClusterOrigin[j];
  }

  void GenDensityOfClusters()
  {


    std::cout << "creating density histograms for " << filein << endl;

    TFile *fout = new TFile(Form("%s_densities.root",filein.c_str()),"RECREATE");
    TH2D *DensityOfClusters[400];
    double PEThres[400];
    for(int h=0;h<400;h++)
    {
      DensityOfClusters[h] = new TH2D(Form("dth2_%i",h),Form("dth2_%i",h),120,-600,600,600,0,6000);
    }
    double counter=0; 

    for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i); 	//std::cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  std::endl; lets_pause();
      for(unsigned int j=0;j<ClustersVector->size();j++)
      {
        int hh=0;
        while(ClustersVector->at(j).PEs()>hh)
        {
          DensityOfClusters[hh]->Fill(ClustersVector->at(j).pos[1],ClustersVector->at(j).pos[2]);
          hh++;
        }
      }
    }

    for(int h=0;h<400;h++)
    {
      DensityOfClusters[h]->Write(Form("dth2_%i",h));
    }

    fout->Close();
    std::cout << "density histograms created for " << filein << endl;
  }


   void LoadClusterDensities()
  {
      cout <<"LOADING density histogram files!!!!. "  <<endl;

    ClusterFile = TFile::Open(Form("%s_densities.root",filein.c_str()),"READ");// ifile->ls();
    if(!ClusterFile)
    {
      cout <<"Problem with the density histogram files!!!!. EXITING... "  <<endl;
//      char option; std::cin >> option;
//      if(option=='Y'||option=='y'){ GenDensityOfClusters();ClusterFile = TFile::Open(Form("%s_densities.root",filein.c_str()),"READ"); }
//      ifile = TFile::Open(Form("%s_densities.root",filein.c_str()),"READ");
//      else
 std::exit(0);
    }

    ClusterTree = (TTree*)ClusterFile->Get("densityTree");//ClusterTree->Print();

    ClusterTree->SetBranchAddress("bth2",&ClusteringHist);
    ClusterTree->SetBranchAddress("bth3",&ClusteringHist3D);
    ClusterTree->GetEntry(0);NReadOutWindows=ClusteringHist3D->GetZaxis()->GetNbins();
    DensitiesLoaded=true;
      cout <<"density histogram files LOADED!!!!. "  <<endl;

  }
  int GetNReadOutWindows(){return NReadOutWindows;}
  double GetNumberOfClustersAtFast(float TruePos[3],double distance,double PEThres)
  {
/*
    bool debug = true;

    if(PEThres>=2*MaxBGCluster)PEThres=2*(MaxBGCluster-1);

    int pt = (int)(PEThres/2);

    double counter=0; //double dist2 = distance*distance;
    int ifirst= ClusteringHist[pt]->GetXaxis()->FindBin(TruePos[1]-distance);
    int ilast = ClusteringHist[pt]->GetXaxis()->FindBin(TruePos[1]+distance);

    int jfirst= ClusteringHist[pt]->GetYaxis()->FindBin(TruePos[2]-distance);
    int jlast = ClusteringHist[pt]->GetYaxis()->FindBin(TruePos[2]+distance);
    
    for (int i=ifirst;i<=ilast;i++)for (int j=jfirst;j<=jlast;j++)
    {
        double dy = TruePos[1]-ClusteringHist[pt]->GetXaxis()->GetBinCenter(i);
        double dz = TruePos[2]-ClusteringHist[pt]->GetYaxis()->GetBinCenter(j);
        if(dy*dy+dz*dz<distance) counter+=ClusteringHist[pt]->GetBinContent(ClusteringHist[pt]->GetBin(i,j));

    }
    if(counter==0)counter++;
    return counter;
*/
   return 0;
  }
  
  void GetEntry(double PEThres)
  {

    if(PEThres>=2*MaxBGCluster)PEThres=2*(MaxBGCluster-1);

    int pt = (int)(PEThres/2);
    ClusterTree->GetEntry(pt);


  }
  std::vector<double> GetNumberOfClustersAtFast3D(float TruePos[3],double distance,double PEThres)
  {

    bool debug = true;

    if(PEThres>=2*MaxBGCluster)PEThres=2*(MaxBGCluster-1);

    int pt = (int)(PEThres/2);
    ClusterTree->GetEntry(pt);


    std::vector<double> counter(NReadOutWindows,0.0); //double dist2 = distance*distance;
    int ifirst= ClusteringHist3D->GetXaxis()->FindBin(TruePos[1]-distance);
    int ilast = ClusteringHist3D->GetXaxis()->FindBin(TruePos[1]+distance);

    int jfirst= ClusteringHist3D->GetYaxis()->FindBin(TruePos[2]-distance);
    int jlast = ClusteringHist3D->GetYaxis()->FindBin(TruePos[2]+distance);
    
    for (int i=ifirst;i<=ilast;i++)for (int j=jfirst;j<=jlast;j++)
    {
        double dy = TruePos[1]-ClusteringHist3D->GetXaxis()->GetBinCenter(i);
        double dz = TruePos[2]-ClusteringHist3D->GetYaxis()->GetBinCenter(j);
        if(dy*dy+dz*dz<distance) for(int k=0;k<NReadOutWindows;k++) counter[k]+=ClusteringHist3D->GetBinContent(ClusteringHist3D->GetBin(i,j,k+1));

    }
//   if(counter==0)counter++;
    return counter;


//   return 0;
  }


 };

}
