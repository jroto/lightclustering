class pmtmap
{

public:
 std::vector<float> PMTPosition[3];
 std::vector<int> PMTIndexZ;
 std::vector<int> PMTIndexY; int NPMTs;
 pmtmap(int n, int drift): NPMTs(n)
{
 PMTIndexY.resize(NPMTs);
 PMTIndexZ.resize(NPMTs);
      PMTPosition[0].resize(NPMTs,0);
      PMTPosition[1].resize(NPMTs,0);
      PMTPosition[2].resize(NPMTs,0);
 if(drift==0)
 {
   for(int pm=0;pm<NPMTs;pm++)PMTIndexZ[NPMTs-pm-1]=(pm-pm%12)/12+1;
   for(int pm=0;pm<NPMTs;pm++)PMTIndexY[NPMTs-pm-1]=pm%12+1;
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[0][pm]=-7.8;//in m
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[1][pm]=1.0*PMTIndexY[pm]-6.5;//in m
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[2][pm]=1.0*PMTIndexZ[pm]-0.5;// in m
 }
 else
 {
   for(int pm=0;pm<NPMTs;pm++)PMTIndexZ[NPMTs-pm-1]=(pm-pm%12)/12+1;
   for(int pm=0;pm<NPMTs;pm++)PMTIndexY[NPMTs-pm-1]=pm%12+1;
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[1][pm]=-7.8;//in m
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[0][pm]=1.0*PMTIndexY[pm]-6.5;//in m
   for(int pm=0;pm<NPMTs;pm++)PMTPosition[2][pm]=1.0*PMTIndexZ[pm]-0.5;// in m
 }
}

 float pos(int pm, int i) {return PMTPosition[i][pm];}
};
