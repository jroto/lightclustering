namespace snana{

class clusteringAnalyzerLight
{
  public:
  float background;
  float efficiency;

  string filein;
  string branchname;
  TFile *ifile=NULL;
  TTree *t = NULL;

  std::vector<int> nevents;
  int nClusters;
  int anaevent;
  float VertX;
  float VertY;
  float VertZ;
  float ENu;
  float ENu_Lep;
  int _anaevent;

  std::vector<int> *ClustersVector = NULL;
  std::vector<double> *ClustersVectorPE = NULL;

  
  clusteringAnalyzerLight(){}
  clusteringAnalyzerLight(string ff, string bb)
  {
    filein=ff;
    branchname=bb;

    ifile = TFile::Open(filein.c_str(),"READ");
    t = (TTree*)ifile->Get("clusteringanatree");
//    t->SetDirectory(0);

//    cout << "tree loaded " << t << endl;
    SetBranchesAddresses();
//     t->GetEntry(2); cout << ClustersVector->at(0).Hits() <<  endl;
  }
  void Close(){ifile->Close();}

 
  void SetBranchesAddresses()
  {
     t->SetBranchAddress("anaevent"   	         ,&anaevent);
     t->SetBranchAddress("VertX"   	         ,&VertX);
     t->SetBranchAddress("VertY"   	         ,&VertY);
     t->SetBranchAddress("VertZ"   	         ,&VertZ);
     t->SetBranchAddress("E_nu"   	         ,&ENu);

     //t->Print();
     t->SetBranchAddress(Form("ClustersHits_%s",branchname.c_str()),&ClustersVector);
     t->SetBranchAddress(Form("ClustersPE_%s",branchname.c_str()),&ClustersVectorPE);
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


   TPaveText * GetClusterStatistics(int i, float hitthreshold)
  {

    bool debug=false;

    std::vector<double> clusters;
    double totalPEs=0;
    t->GetEntry(i); 	//cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  endl; lets_pause();
    for(int j=0;j<ClustersVectorPE->size();j++)
    {
      totalPEs+=ClustersVectorPE->at(j);
      if(ClustersVectorPE->at(j)>=hitthreshold) {clusters.push_back(ClustersVectorPE->at(j));}
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
   for(int j=0;j<clusters.size();j++) pt->AddText(Form("\tCluster of %.0fPE",clusters[j]));
//   pt->UseCurrentStyle();
   pt->SetFillColor(0);

   return pt;
//    return clusters;
  }


  double CalculateDetectionEfficiencySignal(float hitthreshold)
  {

    bool debug=false;

    if(debug)  t->Print();
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
//    cout <<"Analysing " << t->GetEntries() << endl;
    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //cout <<"event " << i << endl;
      t->GetEntry(i); 	//cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  endl; lets_pause();
      if(ClustersVector->size()>0 && ClustersVector->at(0)>=hitthreshold)
      {
         CounterPerEventRun++;
         for(int j=0;j<ClustersVector->size();j++)
         {
           if(ClustersVector->at(j)>=hitthreshold) {CounterAllRun++;CounterAllPerEventRun++;}
           else {break;}
         }  
      }
      if(CounterAllPerEventRun>1)Counter_multiclusterRun++;
      if(CounterAllPerEventRun>1)CounterAll_multiclusterRun+=CounterAllPerEventRun;
    }
    NumberOfEventsRun=t->GetEntries();

    cout << "Detected " << CounterPerEventRun << " out of " << NumberOfEventsRun <<endl;
    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }
  double CalculateDetectionEfficiencySignalPE(double hitthreshold)
  {

    bool debug=false;
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
//    cout <<"Analysing " << t->GetEntries() << endl;
    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //cout <<"event " << i << endl;
      t->GetEntry(i); 	//cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  endl; lets_pause();
      if(ClustersVectorPE->size()>0 && ClustersVectorPE->at(0)>=hitthreshold)
      {
         CounterPerEventRun++;
         for(int j=0;j<ClustersVectorPE->size();j++)
         {
           if(ClustersVectorPE->at(j)>=hitthreshold) {CounterAllRun++;CounterAllPerEventRun++;}
           else {break;}
         }  
      }
      if(CounterAllPerEventRun>1)Counter_multiclusterRun++;
      if(CounterAllPerEventRun>1)CounterAll_multiclusterRun+=CounterAllPerEventRun;
    }
    NumberOfEventsRun=t->GetEntries();
    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }

  double CalculateDetectionEfficiencySignalPECut(double hitthreshold, int axis, float Xmin, float Xmax)
  {


    bool debug=false;
    CounterPerEventRun=0;
    CounterAllRun=0;
    Counter_multiclusterRun=0; // number of events with more than one cluster found
    CounterAll_multiclusterRun=0; // number of cluster per event in events whith more than one cluster found.
    NumberOfEventsRun=0;

    cout <<"Analysing " << filein << " between " << Xmin << " " << Xmax <<  endl;

    for (int i=0;i<t->GetEntries();i++)
    {

      CounterAllPerEventRun=0; //cout <<"event " << i << endl;

      t->GetEntry(i);
//      cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  " VertX" << VertX << endl; lets_pause();

      float position[3]={VertX,VertY,VertZ};
      if(position[axis]<=Xmin || position[axis]>Xmax) continue;
      else NumberOfEventsRun++;

      if(ClustersVector->size()>0 && ClustersVectorPE->at(0)>=hitthreshold  )
      {
         CounterPerEventRun++;
      }
      if(CounterAllPerEventRun>1)Counter_multiclusterRun++;
      if(CounterAllPerEventRun>1)CounterAll_multiclusterRun+=CounterAllPerEventRun;
    }
    
    cout <<"Detected clusters in " << CounterPerEventRun << " out of  " << NumberOfEventsRun << " above " << hitthreshold << "PEs. "<<   endl;

    return 1.0*CounterPerEventRun/NumberOfEventsRun;
  }



  double GetMultiplicity(float hitthreshold)
  {
    double efficiency = CalculateDetectionEfficiencySignalPE(hitthreshold);//percentage of events that triggered one cluster.
    cout << "multiplicity " << Counter_multiclusterRun  << " " << CounterPerEventRun << endl; lets_pause();
    double multip = 1.0*Counter_multiclusterRun / CounterPerEventRun;
    return multip;
  }

  double CalculateBackGroundRate(float hitthreshold)
  {
    double efficiency = CalculateDetectionEfficiencySignal(hitthreshold);//percentage of events that triggered one cluster.

    double backgroundrate = 1.0*CounterAll / (1.e-3*NumberOfEvents ) ;
    return backgroundrate;
  }

  TGraphErrors* BarridoDE(string filein, Barrido b, TGraphErrors* tgmultiplicity, TGraphErrors* tgavmulticluster)
  {
     cout << "barrido00 DE  " << b.name << endl;
     TGraphErrors* tg = new TGraphErrors();
     for(double i=b.min; i<b.max; i+=b.step)
     {
        cout << "barrido DE  " << b.name << " " << i << " ";//lets_pause();
        double efficiency = CalculateDetectionEfficiencySignal(i);//percentage of events that triggered one cluster.
        tg->SetPoint(tg->GetN(),i,efficiency);
        tg->SetPointError(tg->GetN()-1,0,efficiency/TMath::Sqrt(NumberOfEvents));

        multiplicity=(double)Counter_multicluster/CounterPerEvent; //probability of having more than 1 clusters once the event is detected.
        avmulticluster=(double)CounterAll_multicluster/Counter_multicluster; // average number of clusters in events with more than 1 cluster.

        tgmultiplicity->SetPoint(tgmultiplicity->GetN(),i,multiplicity);
        tgmultiplicity->SetPointError(tgmultiplicity->GetN()-1,0.0,multiplicity*TMath::Sqrt(1.0/Counter_multicluster+1.0/CounterPerEvent));

        tgavmulticluster->SetPoint(tgavmulticluster->GetN(),i,avmulticluster);
        tgavmulticluster->SetPointError(tgavmulticluster->GetN()-1,0.0,avmulticluster*TMath::Sqrt(1.0/CounterAll_multicluster+1.0/Counter_multicluster));

        cout << efficiency << " " << efficiency/TMath::Sqrt(NumberOfEvents) <<  " " << multiplicity <<" " << avmulticluster << endl;
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
        cout << "barrido BGR " << b.name << " " << i << endl;//lets_pause();
        double BGR = CalculateBackGroundRate(i);//percentage of events that triggered one cluster.
        tg->SetPoint(tg->GetN(),i,BGR);
        tg->SetPointError(tg->GetN()-1,0,TMath::Sqrt(CounterAll)/(1.e-3*NumberOfEvents));

        tgBGevtPer2DrifWindow->SetPoint(tgBGevtPer2DrifWindow->GetN(),i,CounterAll*15.e-3/(1.e-3*NumberOfEvents));
        tgBGevtPer2DrifWindow->SetPointError(tgBGevtPer2DrifWindow->GetN()-1,0,TMath::Sqrt(CounterAll)*15.e-3/(1.e-3*NumberOfEvents));

        cout << " " << BGR << " " << TMath::Sqrt(CounterAll)/(1.e-3*NumberOfEvents) << "Hz  --- " << CounterAll <<" clusters found"<< endl;
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

   bool debug=false;


    float SampleSize=t->GetEntries()*1e-3; //seconds
    int NumOfEventsInSample=(int)(NumOfEvtsPerWindow*SampleSize/WindowSize);

    cout << "looking for " << NumOfEventsInSample << " events in the sample of " << SampleSize<< "s."<<   endl;
    cout << "To get " << NumOfEvtsPerWindow << " events per " << WindowSize<< "s window."<<   endl;
    std::vector<double> AllClusters; AllClusters.clear();
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//cout << "entry " << i << ", first cluster: " << ClustersVectorPE->at(0) << " number of clusters: " << ClustersVectorPE->size() <<  endl; lets_pause();

      if(ClustersVectorPE->size()>0)
      {
         for(int j=0;j<ClustersVectorPE->size();j++)
         {
           if(AllClusters.size()<NumOfEventsInSample)
           {
              insert( AllClusters,ClustersVectorPE->at(j));
           }
           else
           {
             if(ClustersVectorPE->at(j)<=AllClusters[AllClusters.size()-1]) break;
             if(ClustersVectorPE->at(j)>AllClusters[AllClusters.size()-1])
             {
               insert(AllClusters,ClustersVectorPE->at(j)); //cout <<"adding" <<endl;
               AllClusters.pop_back();
             }
           }
//           for(int k=0;k<AllClusters.size();k++) cout << AllClusters[k] << " "; cout << endl;lets_pause();  
         }
       }
    }

//    for (int i=0;i<AllClusters.size();i++) {cout << AllClusters[i] << " ";} cout << endl ; lets_pause();

    double aux =AllClusters[AllClusters.size()-1];
    return aux;
  }

  int GetMinimumClusterAtLevelHit(float NumOfEvtsPerWindow, float WindowSize)
  {

    float SampleSize=t->GetEntries()*1e-3; //seconds
    int NumOfEventsInSample=(int)(NumOfEvtsPerWindow*SampleSize/WindowSize);

    cout << "looking for " << NumOfEventsInSample << " events in the sample of " << SampleSize<< "s."<<   endl;
    cout << "To get " << NumOfEvtsPerWindow << " events per " << WindowSize<< "s window."<<   endl;
    std::vector<int> AllClusters; AllClusters.clear();
    for (int i=0;i<t->GetEntries();i++)
    {

      t->GetEntry(i); 	//cout << "entry " << i << ", first cluster: " << ClustersVector->at(0) << " number of clusters: " << ClustersVector->size() <<  endl; lets_pause();

      if(ClustersVector->size()>0)
      {
//         cout << "entro" << endl;for(int k=0;k<AllClusters.size();k++) cout << AllClusters[k] << " "; cout << endl;
         for(int j=0;j<ClustersVector->size();j++)
         {
           if(AllClusters.size()<NumOfEventsInSample)
           {
              insert( AllClusters,ClustersVector->at(j));
           }
           else
           {

//             cout << "entro2" << endl;for(int k=0;k<AllClusters.size();k++) cout << AllClusters[k] << " "; cout << endl;
//             cout << ClustersVector->at(j) << " " << AllClusters[AllClusters.size()-1] <<endl;lets_pause();
             if(ClustersVector->at(j)<=AllClusters[AllClusters.size()-1]) break;
             if(ClustersVector->at(j)>AllClusters[AllClusters.size()-1])
             {
               insert(AllClusters,ClustersVector->at(j)); //cout <<"adding" <<endl;
               AllClusters.pop_back();
             }
//cout << "entro3" << endl;for(int k=0;k<AllClusters.size();k++) cout << AllClusters[k] << " "; cout << endl;
           }
//           for(int k=0;k<AllClusters.size();k++) cout << AllClusters[k] << " "; cout << endl;  
         }
       }
    }

//    for (int i=0;i<AllClusters.size();i++) {cout << AllClusters[i] << " ";} cout << endl ; lets_pause();
   int aux =AllClusters[AllClusters.size()-1];
    return aux;
  }

  TGraphErrors* BarridoInversoBGR(std::vector<float> NumOfEvtsPerWindow, float Window)
  {
     TGraphErrors* tg = new TGraphErrors();
     for(int i=0; i<NumOfEvtsPerWindow.size(); i++)
     {
        cout << "barrido BGR " << branchname << " " << i << endl;//lets_pause();
        double NumOfPEs;
        do
        {
          NumOfPEs = GetMinimumClusterAtLevel(NumOfEvtsPerWindow[i], Window);
        }while(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000);
        cout << "\t" << NumOfPEs << "PEs per " << NumOfEvtsPerWindow[i] << "evt/"<<Window<<"s" <<endl<<endl; //if(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000) lets_pause();
        
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
     for(int i=0; i<NumOfEvtsPerWindow.size(); i++)
     {
        cout << "BarridoInversoBGR_Hits " << branchname << " " << i  << " " << NumOfEvtsPerWindow[i] << endl;//lets_pause();
        int NumOfHits;
        do
        {
          NumOfHits = GetMinimumClusterAtLevelHit(NumOfEvtsPerWindow[i], Window);
        }while(!NumOfHits);// || NumOfHits>5000 || NumOfHits<-5000);
        cout << "\t" << NumOfHits << "Hits per " << NumOfEvtsPerWindow[i] << "evt/"<<Window<<"s" <<endl<<endl; //if(!NumOfPEs || NumOfPEs>5000 || NumOfPEs<-5000) lets_pause();
        
        tg->SetPoint(tg->GetN(),NumOfEvtsPerWindow[i],(float)NumOfHits);
        tg->SetPointError(tg->GetN()-1,0,0);
     }
     tg->SetTitle(Form("%s",branchname.c_str()));
     tg->GetXaxis()->SetTitle(Form("# Events in %.2f s",Window));
     tg->GetYaxis()->SetTitle("#Hits");
     tg->SetLineColor(2);
    
     return tg;
  }

  
 };

}
