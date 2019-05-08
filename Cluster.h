#include<iostream>
#include<map>
#include<algorithm>

using namespace std;
class Hit_t
{
 public:
 
  double time;
  int pmt;
  double pos[3];
  float nPE;
  int trackID; //particle that originated the hit
  int GenType; //physical origin of the hit
  Hit_t() { }
  Hit_t(double t, int np, double p[3], double pe, int ti, int gt)
  {
    time=t;
    pmt=np;
    pos[0]=p[0];
    pos[1]=p[1];
    pos[2]=p[2];
    nPE=pe;
    trackID=ti;
    GenType=gt;
  }
 
  ClassDef(Hit_t,1)

};

class Cluster_t
{

  public:
  double peaktime=0; //in us
  int nhits=0;
  int npmts=0;
  float pos[3]={0,0,0};
  double nPEs=0;
  double starttime=20e3; // in us
  double endtime=0; // in us
  std::vector<Hit_t> hits; 

  std::map<int,double> trackID;
  std::map<int,double> GenType;
  
  double purity;
  int PDG;
  int Generator;
  int track;

  double TrueX;
  double TrueZ;
  double TrueY;

 
  Cluster_t(){}
  void FillCluster(double t, int pm, double p[3], double pe, int tI, int gt, int pp)
  {
    hits.push_back(Hit_t(t,pm,p,pe,tI, gt));
    peaktime = (peaktime*nPEs + pe*t)/(pe+nPEs);
    pos[0]   = (pos[0]*nPEs + pe*p[0])/(pe+nPEs);
    pos[1]   = (pos[1]*nPEs + pe*p[1])/(pe+nPEs);
    pos[2]   = (pos[2]*nPEs + pe*p[2])/(pe+nPEs);
    trackID[tI]+=pe; //std::cout << " adding " << pe << "PEs to trackID "<< tI << " with post/value " << trackID[tI] << std::endl; 
    GenType[gt]+=pe;
    if(starttime>t)starttime=t;
    if(endtime<t)endtime=t;

    peaktime = (peaktime*nPEs + pe*t)/(pe+nPEs);
    nhits++;
    nPEs+=pe;
//      std::cout << "filling cluster " << std::endl; 
    if(hits.size()==1)npmts=1;
    else
    {
      for(unsigned int i=0;i<hits.size()-1;i++)
      {
//      std::cout << i << " " << pm << " " << hits[i].pmt << std::endl;

        if(hits[i].pmt==pm) break;
        if(i==hits.size()-2) npmts++;
      }
    }
  }
  void FillPurity()
  {
    std::map<int, double>::iterator x = std::max_element(trackID.begin(), trackID.end(), [](const pair<int, double>& p1, const pair<int, double>& p2) { return p1.second < p2.second; });
    std::map<int, double>::iterator y = std::max_element(GenType.begin(), GenType.end(), [](const pair<int, double>& p1, const pair<int, double>& p2) { return p1.second < p2.second; });

    purity=x->second/nPEs;
    track=x->first;
    Generator=y->first;
  }
  double PEs(){return nPEs;}
  int Hits() {return nhits;}
  int GetTrackID() {return track;}
  int GetGenType() {return Generator;}
  void SetPDG(int p) {PDG=p;}
  void SetPosition(double x, double y , double z) {TrueX=x; TrueY=y; TrueZ=z;}
  void Clear(){trackID.clear(); GenType.clear();hits.clear();}
  void Print()
  {
    std::cout << "Cluster: " << peaktime << "us, "<<nhits << "hits, "<< nPEs <<"PEs, ("<< pos[0]<< ","<< pos[1]<< ","<< pos[2]<< ")"<<std::endl;
  }
  ClassDef(Cluster_t,1)
};
class ClusterCollection
{
  public:
  std::vector<Cluster_t> cl;

  ClassDef(ClusterCollection,1)

};

