namespace snana{


struct ClusterFinderPack
{
  double TimeWindow; // in nanosecons
  float clusterjump; //number of times to the timewindow to jump after a cluster is found;

   int nHitsPerEvent;
   double avHitsPerPMTperEvent;
   double nrecoPEperEvent;
   double avAmpPerEvent;
   int nPMTperEvent;
   int nPMTperEvent_more45ADC;
   int nPMTperEvent_more75ADC;
   int nPMTperEvent_more105ADC;
   int nPMTperEvent_more1hit;
   int nPMTperEvent_more2hit;
   int nPMTperEvent_more3hit;
   double avrecoPEperPMTperEvent;
   float evtsizeY;
   float evtsizeZ;
   float evtsizeY_more1hit;
   float evtsizeZ_more1hit; 
   float evtsizeY_more2hit;
   float evtsizeZ_more2hit; 
   float evtsizeY_more3hit;
   float evtsizeZ_more3hit;
   int nHitsperCluster_adjacents_1m5;
   int nHitsperCluster_adjacents_2m;
   int nHitsperCluster_adjacents_2m5;
   int nHitsperCluster_adjacents_3m;
   int nHitsperCluster_adjacents_3m5;
   int nHitsperCluster_adjacents_4m;
   float xpos;
};

struct Barrido
{
  string name;
  float min;
  float step;
  float max;
};

}


