

namespace snana{

class radevt
{

 bool light;
 bool debug=false;
 public:
 string filein;
 int DriftAxis=0;

 float initime=0; //startime to split the waveform in non radiological sample
 int evtnum;
 int Run;
 int SubRun;
 int Event;
 std::vector<int> *PDS_OpHit_OpChannel= NULL;
 std::vector<double> *PDS_OpHit_X= NULL;
 std::vector<double> *PDS_OpHit_Y= NULL;
 std::vector<double> *PDS_OpHit_Z= NULL;
 std::vector<double> *PDS_OpHit_PeakTime= NULL;
 std::vector<double> * PDS_OpHit_PeakTimeAbs= NULL;
 std::vector<unsigned short> * PDS_OpHit_Frame= NULL;
 std::vector<double> * PDS_OpHit_Width= NULL;
 std::vector<double> * PDS_OpHit_Area= NULL;
 std::vector<double> * PDS_OpHit_Amplitude= NULL;
 std::vector<double> * PDS_OpHit_PE= NULL;
 std::vector<double> * PDS_OpHit_FastToTotal= NULL;
 std::vector<int> * PDS_OpHit_True_GenType= NULL;
 std::vector<double> * PDS_OpHit_True_Energy= NULL;
 std::vector<int> * PDS_OpHit_True_TrackID= NULL;
 std::vector<int> * PDS_OpHit_True_GenTypeAll= NULL;
 std::vector<double> * PDS_OpHit_EnergyAll= NULL;
 std::vector<int> * PDS_OpHit_TrackIDAll= NULL;
 std::vector<int> * PDS_OpHit_IndexAll= NULL;
 std::vector<int> * True_Bck_Mode= NULL;
 std::vector<double> * True_Bck_VertX= NULL;
 std::vector<double> * True_Bck_VertY= NULL;
 std::vector<double> * True_Bck_VertZ= NULL;
 std::vector<double> * True_Bck_Time= NULL;
 std::vector<double> * True_Bck_Energy= NULL;
 std::vector<int> * True_Bck_PDG= NULL;
 std::vector<int> * True_Bck_ID= NULL;
 int TotGen_Marl;
 int TotGen_Ar39;
 int TotGen_Neut;
 int TotGen_Kryp;
 int TotGen_Plon;
 int TotGen_Rdon;
 int TotGen_Ar42;
 std::vector<int> sorted; 
 std::vector<int> sortedBySize; 

 // pmtmap mypmts = pmtmap(720);

  int _anaevent;
  float _VertX;
  float _VertY;
  float _VertZ;
  float _E_nu;
  float _E_nu_Lep;

  std::vector<Cluster_t> Clusters;
  std::vector<double> LightClustersPE;
  std::vector<int> LightClustersHits;

  string sClusteringD;
  string sClusteringW;
  string sClusteringMW;
  float fClusteringD;//cm
  float fClusteringW;//us
  float fClusteringMW;//us
 radevt() { }
 radevt(int da) { DriftAxis=da; }
 radevt(string sD, float fD, string sW, float fW, string sMW, float fMW, bool lt):  sClusteringD(sD), fClusteringD (fD), sClusteringW(sW), fClusteringW(fW), sClusteringMW (sMW), fClusteringMW (fMW), light(lt)
 {
    std::cout << "ligth boolean " << light << endl;
 }
 double GetDriftPosition()
 {
   if(DriftAxis==0) return _VertX;
   if(DriftAxis==1) return _VertY;
   else return _VertZ;
 }

 void Clear()
 {
  Clusters.clear();
  LightClustersPE.clear();
  LightClustersHits.clear();
  sorted.clear();
  sortedBySize.clear();
  Map_TrackIdToPDG.clear();
  Map_TrackIdToTrueX.clear();
  Map_TrackIdToTrueY.clear();
  Map_TrackIdToTrueZ.clear();
  multiplicity.clear();
  True_TrackPDG.clear();
  True_Multiplicity.clear();
 }
 std::map<int, int> Map_TrackIdToPDG;
 std::map<int, double> Map_TrackIdToTrueX;
 std::map<int, double> Map_TrackIdToTrueY;
 std::map<int, double> Map_TrackIdToTrueZ;

 std::map<int, int> multiplicity;
 std::vector<int>  True_TrackPDG;
 std::vector<int>  True_Multiplicity;

 void FillAuxVars()
 {
  
   int nupdg[] = {-18,-16,-14,-12,12,14,16,18};
   std::vector<double> nupos[3], nuE;
   for(int i=0;i<True_Bck_ID->size();i++)
   {
     if(debug) std::cout << "Particle " << i << " - TrackID: " << True_Bck_ID->at(i) << " - PDG:" << True_Bck_PDG->at(i) <<std::endl;
     Map_TrackIdToPDG[True_Bck_ID->at(i)]=True_Bck_PDG->at(i);
     Map_TrackIdToTrueX[True_Bck_ID->at(i)]=True_Bck_VertX->at(i);
     Map_TrackIdToTrueY[True_Bck_ID->at(i)]=True_Bck_VertY->at(i);
     Map_TrackIdToTrueZ[True_Bck_ID->at(i)]=True_Bck_VertZ->at(i);
     if(std::find(std::begin(nupdg), std::end(nupdg), True_Bck_PDG->at(i)) != std::end(nupdg))
     {
       nupos[0].push_back(True_Bck_VertX->at(i));
       nupos[1].push_back(True_Bck_VertY->at(i));
       nupos[2].push_back(True_Bck_VertZ->at(i));
       nuE.push_back(True_Bck_Energy->at(i));
     }
   }
   switch(nuE.size())
   {
     case 0:
       if(debug) std::cout <<"No neutrino found!, are we running without signal>??? " << std::endl;
       break;
     case 1:
       _VertX=nupos[0][0];
       _VertY=nupos[1][0];
       _VertZ=nupos[2][0];
       break;
     default:
       if(debug) std::cout <<"More than one neutrino found!, we will take the first one!!" << std::endl;
       _VertX=nupos[0][nuE.size()-1];
       _VertY=nupos[1][nuE.size()-1];
       _VertZ=nupos[2][nuE.size()-1];
       break;
   }
 }

 std::vector<Cluster_t> PureAdjacentClustering(float distance, float timedistance, float maxtimedistance)
 {
   std::vector<bool> processed;
   processed.resize(PDS_OpHit_X->size(),false);

//   float counterPE=0;for(int i=0;i<PDS_OpHit_PE->size();i++)counterPE+=PDS_OpHit_PE->at(i);

//   std::cout << PDS_OpHit_X->size() << " total hits - " << counterPE << "PEs. " << std::endl;

   std::vector<Cluster_t> clusterlist; clusterlist.clear();//nhits per cluster

   for (unsigned int h=0;h<PDS_OpHit_X->size();h++)
   { 
//     std::cout << " Analysing hit " << h <<std::endl;
     if(!processed[h])
     {
       Cluster_t cl; //create an empty cluster
       processed[h]=true;
       double pos[3]={PDS_OpHit_X->at(sorted[h]),PDS_OpHit_Y->at(sorted[h]),PDS_OpHit_Z->at(sorted[h])};
       cl.FillCluster(PDS_OpHit_PeakTime->at(sorted[h]),PDS_OpHit_OpChannel->at(sorted[h]),pos,PDS_OpHit_PE->at(sorted[h]),PDS_OpHit_True_TrackID->at(sorted[h]),PDS_OpHit_True_GenType->at(sorted[h]),PDS_OpHit_True_GenType->at(sorted[h]));
       clusterlist.push_back(cl);

       std::vector<int> N; N.clear(); //vector of neighbours

       std::vector<int> nb =getNeighbors(h,distance, timedistance, maxtimedistance, processed, PDS_OpHit_PeakTime->at(sorted[h])); 
       N.insert( N.end(), nb.begin(), nb.end() );//std::cout << "added " << nb.size() << "pmts "<< std::endl;

       for (unsigned int i=0;i<N.size();i++)
       {
           double pos[3]={  PDS_OpHit_X->at(sorted[N[i]]),
                                                             PDS_OpHit_Y->at(sorted[N[i]]),
                                                             PDS_OpHit_Z->at(sorted[N[i]])   };

           clusterlist[clusterlist.size()-1].FillCluster( PDS_OpHit_PeakTime->at(sorted[N[i]]),
                                                          PDS_OpHit_OpChannel->at(sorted[N[i]]),
                                                          pos,
                                                          PDS_OpHit_PE->at(sorted[N[i]]),PDS_OpHit_True_TrackID->at(sorted[h]),PDS_OpHit_True_GenType->at(sorted[h]),PDS_OpHit_True_GenType->at(sorted[h])        );
           std::vector<int> nb2 =getNeighbors(N[i], distance, timedistance, maxtimedistance, processed, PDS_OpHit_PeakTime->at(sorted[h])); 
           N.insert( N.end(), nb2.begin(), nb2.end() );// std::cout << "added " << nb2.size() << "pmts "<< std::endl;
       }
     } 
   }
  std::sort(  std::begin(clusterlist), std::end(clusterlist), [&](Cluster_t c1, Cluster_t c2) { return c1.nPEs > c2.nPEs; } );
  return clusterlist;
 }

 void Process()
 {

  Clear();


  int ntothits=PDS_OpHit_OpChannel->size();

  if(debug) std::cout << "Processing event with " << ntothits << " hits." <<std::endl;


  std::size_t n(0);
  sorted.resize(ntothits);
  std::generate(std::begin(sorted), std::end(sorted), [&]{ return n++; });
  std::sort(  std::begin(sorted), std::end(sorted), [&](int i1, int i2) { return PDS_OpHit_PeakTime->at(i1) < PDS_OpHit_PeakTime->at(i2); } );

  FillAuxVars();

//  std::size_t n2(0);
//  sortedBySize.resize(ntothits);
//  std::generate(std::begin(sortedBySize), std::end(sortedBySize), [&]{ return n22++; });
//  std::sort(  std::begin(sortedBySize), std::end(sortedBySize), [&](int i1, int i2) { return PDS_OpHit_PE->at(i1) > PDS_OpHit_PE->at(i2); } );


//  std::cout << "holi" << std::endl;lets_pause();
  //distance, time, maxtime
  Clusters = PureAdjacentClustering(fClusteringD,fClusteringW,fClusteringMW);

  if(light)
  {
    for(int i =0; i<Clusters.size();i++) LightClustersPE.push_back(Clusters[i].PEs());
    for(int i =0; i<Clusters.size();i++) LightClustersHits.push_back(Clusters[i].Hits());
    std::sort(  std::begin(LightClustersHits), std::end(LightClustersHits), [&](int i1, int i2) { return i1 > i2; } );
    std::sort(  std::begin(LightClustersPE),std::end(LightClustersPE), [&](int i1, int i2) { return i1 > i2; } );

  }
  else
  {
    for(int i =0; i<Clusters.size();i++) Clusters[i].FillPurity();
    for(int i =0; i<Clusters.size();i++) {Clusters[i].SetPDG(Map_TrackIdToPDG[Clusters[i].GetTrackID()]);multiplicity[Clusters[i].GetTrackID()]=0;}
    for(int i =0; i<Clusters.size();i++) multiplicity[Clusters[i].GetTrackID()]++;
    for(int i =0; i<Clusters.size();i++) Clusters[i].SetPosition(Map_TrackIdToTrueX[Clusters[i].GetTrackID()],Map_TrackIdToTrueY[Clusters[i].GetTrackID()],Map_TrackIdToTrueZ[Clusters[i].GetTrackID()]);
    for(int i =0; i<Clusters.size();i++) Clusters[i].Clear();

    map<int,int>::iterator it;
    for ( it = multiplicity.begin(); it != multiplicity.end(); it++ )
    {
      True_TrackPDG.push_back( Map_TrackIdToPDG[it->first]);
      True_Multiplicity.push_back( it->second);
    }
  }
 }

//, distance, timedistance, maxtimedistance, processed); 
 std::vector<int> getNeighbors(int hitnumber, float distthres, float timedistance, float maxtimedistance, std::vector<bool> &processed, float initimecluster)
//std::vector<int> hits, std::map<int,bool> &processed)
 {

  int itTmin=0;
  int itTmax=PDS_OpHit_PeakTime->size();
  for (int h=hitnumber;h>0;h--) if(PDS_OpHit_PeakTime->at(sorted[h]) < PDS_OpHit_PeakTime->at(sorted[hitnumber])-timedistance)
  {
    itTmin=h; break;
  }
  for (unsigned int h=hitnumber;h<PDS_OpHit_PeakTime->size();h++) if(PDS_OpHit_PeakTime->at(sorted[h]) > PDS_OpHit_PeakTime->at(sorted[hitnumber])+timedistance)
  {
    itTmax=h; break;
  }

   std::vector<int> neighbors;
//   std::cout << "hitnumber " << hitnumber <<" min "<< itTmin << " - max " << itTmax << std::endl; lets_pause();
   float dx, dy, dz, dt, dtmax;
   for (int h=itTmin;h<itTmax;h++)
//   for (int h=0;h<PDS_OpHit_PeakTime->size();h++)
   {
     if(!processed[h])
// && TMath::Abs(initimecluster - PDS_OpHit_PeakTime->at(sorted[h]))<maxtimedistance)
     {
       dx  = PDS_OpHit_X->at(sorted[hitnumber])        - PDS_OpHit_X->at(sorted[h]);
       dy  = PDS_OpHit_Y->at(sorted[hitnumber])        - PDS_OpHit_Y->at(sorted[h]);
       dz  = PDS_OpHit_Z->at(sorted[hitnumber])        - PDS_OpHit_Z->at(sorted[h]);
       dt  = TMath::Abs(PDS_OpHit_PeakTime->at(sorted[hitnumber]) - PDS_OpHit_PeakTime->at(sorted[h]) );
       dtmax  = TMath::Abs(initimecluster - PDS_OpHit_PeakTime->at(sorted[h]) );
       float distance = dx*dx + dy*dy + dz*dz;
//       std::cout << distance << " " << dt << " " << dtmax << std::endl; lets_pause();
       if(distance<distthres*distthres && dt<timedistance && dtmax<maxtimedistance){ neighbors.push_back(h);processed[h]=true;}
     }
   }

   return neighbors;
 }


 void printmap(std::map<int,bool>mymap)
 {
   map<int,bool>::iterator it;

   for ( it = mymap.begin(); it != mymap.end(); it++ )
   {
     std::cout << it->first  // string (key)
              << ':'
              << it->second   // string's value 
              << std::endl ;
   }
 }

 void FillAnaTree(TTree *t)
 {
     _anaevent=t->GetEntries();
     t->Fill();//lets_pause();

 }
 void SetAnaTreeBranchAdresses(TTree *t)
 {
   t->Branch("anaevent"       ,&_anaevent);
   t->Branch("VertX"   	      ,&_VertX);
   t->Branch("VertY"   	      ,&_VertY);
   t->Branch("VertZ"   	      ,&_VertZ);
   t->Branch("E_nu"   	      ,&_E_nu);

//   gInterpreter->GenerateDictionary("Cluster_t", "Cluster.h");
   if (light)
   {
     t->Branch(Form("ClustersHits_D%s_W%s_MW%s",sClusteringD.c_str(),sClusteringW.c_str(),sClusteringMW.c_str()), &LightClustersHits   );
     t->Branch(Form("ClustersPE_D%s_W%s_MW%s",sClusteringD.c_str(),sClusteringW.c_str(),sClusteringMW.c_str()), &LightClustersPE   );
   }
   else
   {
     t->Branch(Form("Clusters_D%s_W%s_MW%s",sClusteringD.c_str(),sClusteringW.c_str(),sClusteringMW.c_str()), &Clusters   );
     t->Branch(Form("True_TrackPDG"), &True_TrackPDG   );
     t->Branch(Form("True_Multiplicity"), &True_Multiplicity   );

   }

 }
 void SetBranchAdresses(TTree *t)
 {
  t->SetBranchAddress("Run"       ,  &Run );
  t->SetBranchAddress("SubRun"       ,  &SubRun );
  t->SetBranchAddress("Event"       ,  &Event );

  t->SetBranchAddress("PDS_OpHit_OpChannel", &PDS_OpHit_OpChannel );
  t->SetBranchAddress("PDS_OpHit_X"       ,  &PDS_OpHit_X );
  t->SetBranchAddress("PDS_OpHit_Y"       ,  &PDS_OpHit_Y );
  t->SetBranchAddress("PDS_OpHit_Z"       ,  &PDS_OpHit_Z );
  t->SetBranchAddress("PDS_OpHit_PeakTime"       ,  &PDS_OpHit_PeakTime );
  t->SetBranchAddress("PDS_OpHit_PeakTimeAbs"       ,  &PDS_OpHit_PeakTimeAbs );
  t->SetBranchAddress("PDS_OpHit_Frame"       ,  &PDS_OpHit_Frame );
  t->SetBranchAddress("PDS_OpHit_Width"       ,  &PDS_OpHit_Width );
  t->SetBranchAddress("PDS_OpHit_Area"       ,  &PDS_OpHit_Area );
  t->SetBranchAddress("PDS_OpHit_Amplitude"       ,  &PDS_OpHit_Amplitude );
  t->SetBranchAddress("PDS_OpHit_PE"       ,  &PDS_OpHit_PE );
  t->SetBranchAddress("PDS_OpHit_FastToTotal"       ,  &PDS_OpHit_FastToTotal );
  t->SetBranchAddress("PDS_OpHit_True_GenType"       ,  &PDS_OpHit_True_GenType );
  t->SetBranchAddress("PDS_OpHit_True_Energy"       ,  &PDS_OpHit_True_Energy );
  t->SetBranchAddress("PDS_OpHit_True_TrackID"       ,  &PDS_OpHit_True_TrackID );

  t->SetBranchAddress("True_Bck_VertX"    ,  &True_Bck_VertX );
  t->SetBranchAddress("True_Bck_VertY"    ,  &True_Bck_VertY );
  t->SetBranchAddress("True_Bck_VertZ"    ,  &True_Bck_VertZ );
  t->SetBranchAddress("True_Bck_Time"     ,  &True_Bck_Time );
  t->SetBranchAddress("True_Bck_Energy"   ,  &True_Bck_Energy );
  t->SetBranchAddress("True_Bck_PDG"      ,  &True_Bck_PDG );
  t->SetBranchAddress("True_Bck_ID"       ,  &True_Bck_ID );
  t->SetBranchAddress("True_Bck_PDG"      ,  &True_Bck_PDG );

  t->SetBranchAddress("TotGen_Marl"       ,  &TotGen_Marl );
  t->SetBranchAddress("TotGen_Ar39"       ,  &TotGen_Ar39 );
  t->SetBranchAddress("TotGen_Neut"       ,  &TotGen_Neut );
  t->SetBranchAddress("TotGen_Kryp"       ,  &TotGen_Kryp );
  t->SetBranchAddress("TotGen_Plon"       ,  &TotGen_Plon );
  t->SetBranchAddress("TotGen_Rdon"       ,  &TotGen_Rdon );
  t->SetBranchAddress("TotGen_Ar42"       ,  &TotGen_Ar42 );

 }

  double GetNPEperEvent()
  {
    double nPEs=0;
    for(int i=0;i<PDS_OpHit_PE->size();i++)
    {
      nPEs+=PDS_OpHit_PE->at(i);
    }
    return nPEs;
  }

  void GetPEperPMT(std::vector<double> &nPEs)
  {
    nPEs.clear();
    nPEs.resize(720,0.0);
//    std::vector<double> nPEs(nchannels,0.0);

    for(int i=0;i<PDS_OpHit_PE->size();i++)
    {
      nPEs[PDS_OpHit_OpChannel->at(i)]+=PDS_OpHit_PE->at(i);
    }
    //return nPEs;    
  }

 };

}
 
