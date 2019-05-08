

namespace snana{

  class matching
  {
    public:
    bool debug=false;

    clustering *BG=NULL;
    clustering *Signal=NULL;
    
    matching(string sg, string bg, string var, bool GenDensities, bool dy)
    {
      BG=new clustering(bg,var,dy);
      Signal=new clustering(sg, var,dy);
      if(GenDensities)BG->GenDensityOfClusters();
    }

    double MatchingEfficiency(double distance, string mode, double &noBGeff, double &noSIGeff, double &numberOfBGClusters)
    {

      std::cout << Signal->t->GetEntries()<< "events of signal "<< endl;
      std::cout << BG->t->GetEntries()<< "events of BG "<< endl;
      noBGeff=0;
      noSIGeff=0;
      double MatchEff=0;
      numberOfBGClusters=0;

      for(int i=0;i<Signal->t->GetEntries();i++)
      {
         Signal->t->GetEntry(i);
         if(Signal->ClustersVector->size()>0)
         {
           if(debug){ std::cout <<"Entry " <<i << " "; Signal->ClustersVector->at(0).Print();}
        
           float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};
           double NumberOfClustersAt;
	   if(mode=="full") NumberOfClustersAt= BG->GetNumberOfClustersAt(TruePos,distance,0.8*Signal->ClustersVector->at(0).PEs()); //amount of clusters that matches with this event
           else NumberOfClustersAt = BG->GetNumberOfClustersAtFast(TruePos,distance,0.5*Signal->ClustersVector->at(0).PEs()); //amount of clusters that matches with this event
          
           double TotalBGWindow = BG->t->GetEntries()*1e-3; //seconds of BG sample
           double ReadoutWindow = 16.5e-3; //seconds
           double AverageNumberOfBackgroundClustersPerReadoutWindow = NumberOfClustersAt*ReadoutWindow/TotalBGWindow;
           double dy =Signal->VertY - Signal->ClustersVector->at(0).pos[1];
           double dz =Signal->VertZ - Signal->ClustersVector->at(0).pos[2];
           double ProbabilityToBeMatchedWell;

           if(dy*dy+dz*dz > distance*distance) {ProbabilityToBeMatchedWell = 0;}
           else {ProbabilityToBeMatchedWell = 1.0/(1.0+AverageNumberOfBackgroundClustersPerReadoutWindow);noBGeff++;}

           noSIGeff+=1.0/(1.0+AverageNumberOfBackgroundClustersPerReadoutWindow);
           numberOfBGClusters+=AverageNumberOfBackgroundClustersPerReadoutWindow;
           if(i%10000==0) std::cout << "Event " << i << " out of " << Signal->t->GetEntries() << endl; //<<NumberOfClustersAt << "BG Clusters, "<< AverageNumberOfBackgroundClustersPerReadoutWindow << "BG Clusters per RW, Prob=" << ProbabilityToBeMatchedWell << std::endl; 
           MatchEff+=ProbabilityToBeMatchedWell;
           if(debug){ std::cout <<"AverageBGClusters " <<AverageNumberOfBackgroundClustersPerReadoutWindow << " ("<<NumberOfClustersAt<< " at " << BG->t->GetEntries()*1.e-3 << "s" <<endl;}
         }
         else {if(debug)std::cout << "Event without clusters. " << std::endl;}
         //break;
      }
      noBGeff/=(double)Signal->t->GetEntries();
      noSIGeff/=(double)Signal->t->GetEntries();
      MatchEff/=(double)Signal->t->GetEntries();
      numberOfBGClusters/=(double)Signal->t->GetEntries();
      
      return MatchEff;
    }
    double MatchingEfficiencyRight(double distance, double &Purity, double &noSIGeff, double &numberOfBGClusters, double &nevents, bool AllowRepeatBackground, int BGeventsPerSignalEvent=8)
    {

      int nSig=Signal->t->GetEntries();
      int nBG=(BG->t->GetEntries()-(BG->t->GetEntries()%BGeventsPerSignalEvent))/BGeventsPerSignalEvent;

      std::cout << nSig<< "events of signal "<< endl;
      std::cout << BG->t->GetEntries()<< "events of BG, corresponding to " << nBG  << " full readout windows." << endl;

      int ntot;
      if (AllowRepeatBackground)
      { 
        if(nSig>nBG)ntot=nSig;
      }
      else
      {
        if(nSig>nBG)ntot=nBG;
      }
      nevents=0;
      double Efficiency=0;
      Purity=0;
      noSIGeff=0;
      numberOfBGClusters=0;
      double numberOfDetectedEvents=0;
      int i=0;
      int entrySIG=0;
      int entryBG=0;

      do
      {
         Signal->t->GetEntry(entrySIG);
         float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};

         if(Signal->ClustersVector->size()>0)
         {
           if(debug){ std::cout <<"Entry " <<i << " "; Signal->ClustersVector->at(0).Print(); lets_pause();}
           int j=0;
           for(j=0;j<Signal->ClustersVector->size();j++)
           {
             double dy =Signal->VertY - Signal->ClustersVector->at(j).pos[1];
             if (Signal->DriftY)dy = Signal->VertX - Signal->ClustersVector->at(j).pos[0];
             double dz =Signal->VertZ - Signal->ClustersVector->at(j).pos[2];
             if(dy*dy+dz*dz < distance*distance)
             {
     
               double NumberOfClustersAt = BG->GetNumberOfClustersAtBetween(TruePos,distance,0.9*Signal->ClustersVector->at(0).PEs(),entryBG*BGeventsPerSignalEvent,(entryBG+1)*(BGeventsPerSignalEvent)); //amount of clusters that matches with this event
               double ProbabilityToBeMatchedWell = 1.0/(1.0+NumberOfClustersAt);
               Purity+=ProbabilityToBeMatchedWell;
               Efficiency++;
          
               noSIGeff+=1.0/(1.0+NumberOfClustersAt);
               numberOfBGClusters+=NumberOfClustersAt;

               if(i%200==0||debug) std::cout << "Event " << i << " out of " << ntot << ", " <<NumberOfClustersAt << "BG Clusters, Prob=" << ProbabilityToBeMatchedWell << std::endl; 
               break;
             }
           }
         }
         else {if(debug)std::cout << "Event without clusters. " << std::endl;}
         //break;

        nevents++;
        i++; entrySIG++; entryBG++;
        if(entrySIG==nSig)entrySIG=0;
        if(entryBG==nBG)entryBG=0;
      } while(i<ntot);

      Purity/=(double)Efficiency;
      noSIGeff/=(double)Efficiency;
      numberOfBGClusters/=(double)Efficiency;
      Efficiency/=ntot;
      
      return Efficiency;
    }



    TH1F* OriginOfMatchedEvents(double distance, bool AllowRepeatBackground, int BGeventsPerSignalEvent=8)
    {

      TH1F *Origins= new TH1F("Origin","Origin;Origin;Counts",12,0,12);
      int nSig=Signal->t->GetEntries();
      int nBG=(BG->t->GetEntries()-(BG->t->GetEntries()%BGeventsPerSignalEvent))/BGeventsPerSignalEvent;

      std::cout << nSig<< "events of signal "<< endl;
      std::cout << BG->t->GetEntries()<< "events of BG, corresponding to " << nBG  << " full readout windows." << endl;

      int ntot;
      if (AllowRepeatBackground)
      { 
        if(nSig>nBG)ntot=nSig;
      }
      else
      {
        if(nSig>nBG)ntot=nBG;
      }
      int nevents=0;
      int i=0;
      int entrySIG=0;
      int entryBG=0;
      do
      {
       Signal->t->GetEntry(entrySIG);
         float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};

         if(Signal->ClustersVector->size()>0)
         {
           if(debug){ std::cout <<"Entry " <<i << " "; Signal->ClustersVector->at(0).Print();}
           int j=0;
           for(j=0;j<Signal->ClustersVector->size();j++)
           {
             double dy =Signal->VertY - Signal->ClustersVector->at(j).pos[1];
             if (Signal->DriftY)dy = Signal->VertX - Signal->ClustersVector->at(j).pos[0];
             double dz =Signal->VertZ - Signal->ClustersVector->at(j).pos[2];
             if(dy*dy+dz*dz < distance*distance)
             {
     
               double NumberOfClustersAt = BG->GetNumberOfClustersAtBetween(TruePos,distance,Signal->ClustersVector->at(0).PEs(),entryBG*BGeventsPerSignalEvent,(entryBG+1)*(BGeventsPerSignalEvent)); //amount of clusters that matches with this event
               int ClusterOrigin=0;
               
               if (NumberOfClustersAt>0) ClusterOrigin = BG->GetOriginOfLargerClusterAtBetween(TruePos,distance,Signal->ClustersVector->at(0).PEs(),entryBG*BGeventsPerSignalEvent,(entryBG+1)*(BGeventsPerSignalEvent));
               if(i%200==0||debug) std::cout << "NumberOfClustersAt= "<< NumberOfClustersAt << ", Cluster Origin " << ClusterOrigin << std::endl;
               if(debug && NumberOfClustersAt!=0) lets_pause();
               Origins->Fill(ClusterOrigin);
               break;
             }
           }
         }
         else {if(debug)std::cout << "Event without clusters. " << std::endl;}
         //break;

        nevents++;
        i++; entrySIG++; entryBG++;
        if(entrySIG==nSig)entrySIG=0;
        if(entryBG==nBG)entryBG=0;
      } while(i<ntot);
      
      return Origins;
    }

    TProfile* MatchingEfficiencyRightScan(double distance, TProfile &Purity, TProfile &EffxPurity, TProfile &NumberOfBGClusters, bool AllowRepeatBackground, int nbins=12, int axis=0,int BGeventsPerSignalEvent=8)
    {

      string axisname; double min, max;
      switch(axis)
      {
         case 0: axisname="X Drift (cm)"; min=-600.; max=600.; break;
         case 1: axisname="Y (cm)"; min=-600.; max=600.; break;
         case 2: axisname="Z (cm)"; min=0.; max=6000.;   break;
         default: axis=0; axisname="X (cm)"; min=-600.; max=600.;
      }

      TProfile *eff = new TProfile("Efficiency",Form(" Efficiency Scan;%s;Effiency",axisname.c_str()),nbins, min, max);
      Purity= TProfile("Purity",Form("Purity Scan - No Background;%s;Purity",axisname.c_str()),nbins, min, max);
      EffxPurity= TProfile("EffxPurity",Form("EffxPurity Scan - Assuming perfect signal detection;%s;EffxPurity effiency",axisname.c_str()),nbins, min, max);
      NumberOfBGClusters= TProfile("NumberOfBGClusters",Form("NumberOfBGClusters Scan;%s;Number Of BGClusters in vecinity",axisname.c_str()),nbins, min, max);

      bool debug=false;
      int nSig=Signal->t->GetEntries();
      int nBG=(BG->t->GetEntries()-(BG->t->GetEntries()%BGeventsPerSignalEvent))/BGeventsPerSignalEvent;

      std::cout << nSig<< "events of signal "<< endl;
      std::cout << BG->t->GetEntries()<< "events of BG, corresponding to " << nBG  << " full readout windows." << endl;

      int ntot;
      if (AllowRepeatBackground)
      { 
        if(nSig>nBG)ntot=nSig;
      }
      else
      {
        if(nSig>nBG)ntot=nBG;
      }
      int i=0;
      int entrySIG=0;
      int entryBG=0;
      do
      {
         Signal->t->GetEntry(entrySIG);
         float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};
         if(Signal->ClustersVector->size()>0)
         {
           if(debug){ std::cout <<"Entry " <<i << " "; Signal->ClustersVector->at(0).Print(); lets_pause();}
           int j=0;
           for(j=0;j<Signal->ClustersVector->size();j++)
           {
             double dy =Signal->VertY - Signal->ClustersVector->at(j).pos[1];
             if (Signal->DriftY)dy = Signal->VertX - Signal->ClustersVector->at(j).pos[0];
             double dz =Signal->VertZ - Signal->ClustersVector->at(j).pos[2];
             if(dy*dy+dz*dz < distance*distance)
             {
     
               double NumberOfClustersAt;
	       NumberOfClustersAt= BG->GetNumberOfClustersAtBetween(TruePos,distance,0.9*Signal->ClustersVector->at(0).PEs(),entryBG*BGeventsPerSignalEvent,(entryBG+1)*(BGeventsPerSignalEvent)); //amount of clusters that matches with this event
               double ProbabilityToBeMatchedWell = 1.0/(1.0+NumberOfClustersAt);
               
               eff->Fill(TruePos[axis],1.0);
               Purity.Fill(TruePos[axis],ProbabilityToBeMatchedWell);

               if(i%200==0||debug) std::cout << "Event " << i << " out of " << ntot << ", " <<NumberOfClustersAt << "BG Clusters, Prob=" << ProbabilityToBeMatchedWell << std::endl; 

               NumberOfBGClusters.Fill(TruePos[axis],NumberOfClustersAt);
               EffxPurity.Fill(TruePos[axis],ProbabilityToBeMatchedWell);
               break;
             }
           }
           if(j==Signal->ClustersVector->size())
           {
             EffxPurity.Fill(TruePos[axis],0.0);
             eff->Fill(TruePos[axis],0.0);
           }
         
        }
        else
        {
          EffxPurity.Fill(TruePos[axis],0.0);
          eff->Fill(TruePos[axis],0.0);
        }

        i++; entrySIG++; entryBG++;
        if(entrySIG==nSig)entrySIG=0;
        if(entryBG==nBG)entryBG=0;

      } while(i<ntot);
      std::cout << "leaving"<<std::endl;
     return eff;
    }

    double MatchingEfficiencyFast(double distance, double &noBGeff, double &noSIGeff, double &noBGClusters)
    {

      if(!BG->DensitiesLoaded)BG->LoadClusterDensities();
     return MatchingEfficiency(distance,"fast", noBGeff,noSIGeff, noBGClusters);
    }
    double MatchingEfficiencySplitting(double distance, double &noBGeff, double &noSIGeff, double &noBGClusters)
    {

      if(!BG->DensitiesLoaded)BG->LoadClusterDensities();



      std::cout << Signal->t->GetEntries()<< "events of signal "<< endl;
      std::cout << BG->t->GetEntries()<< "events of BG "<< endl;

      double efficiency=0;
      noBGeff=0;
      noSIGeff=0;
      noBGClusters=0;
      int NReadOutWindows = BG->GetNReadOutWindows();
      std::vector<double> NumberOfClustersVector;
      std::cout << "Looping over " << NReadOutWindows << " readout windows." << endl;
      double NumberOfClustersAt;
      double dz, dy, ProbabilityToBeMatchedWell;
      for(int i=0;i<Signal->t->GetEntries();i++)
      {
         Signal->t->GetEntry(i);
         if(Signal->ClustersVector->size()>0)
         {
           if(debug){ std::cout <<"Entry " <<i << " "; Signal->ClustersVector->at(0).Print();}
        
           float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};

           NumberOfClustersVector = BG->GetNumberOfClustersAtFast3D(TruePos,distance,0.8*Signal->ClustersVector->at(0).PEs());

           for( int b=0; b<NReadOutWindows; b++)
           {
              NumberOfClustersAt = NumberOfClustersVector[b]; //amount of clusters that matches with this event
 
              dy =Signal->VertY - Signal->ClustersVector->at(0).pos[1];
              dz =Signal->VertZ - Signal->ClustersVector->at(0).pos[2];

              if(dy*dy+dz*dz > distance*distance) {ProbabilityToBeMatchedWell = 0;}
              else {ProbabilityToBeMatchedWell = 1.0/(1.0+NumberOfClustersAt);noBGeff++;}

              efficiency+=ProbabilityToBeMatchedWell;
              noSIGeff+=1.0/(1.0+NumberOfClustersAt);
              noBGClusters+=NumberOfClustersAt;

              if(debug){ std::cout <<"NumberOfClustersAt " <<NumberOfClustersAt << " "<<b <<endl;}
           }
           
           if(i%500==0) std::cout << "Event " << i << " out of " << Signal->t->GetEntries() << endl; //<<NumberOfClustersAt << "BG Clusters, "<< AverageNumberOfBackgroundClustersPerReadoutWindow << "BG Clusters per RW, Prob=" << ProbabilityToBeMatchedWell << std::endl; 

         }
         else {if(debug)std::cout << "Event without clusters. " << std::endl;}
         //break;
      }

      efficiency/=(double)(Signal->t->GetEntries()*NReadOutWindows);
      noBGeff/=(double)(Signal->t->GetEntries()*NReadOutWindows);
      noSIGeff/=(double)(Signal->t->GetEntries()*NReadOutWindows);
      noBGClusters/=(double)(Signal->t->GetEntries()*NReadOutWindows);
      
      return efficiency;

    }


    TProfile* MatchingEfficiencyFastScan(double distance, TProfile &noBGeff, TProfile &noSIGeff, TProfile &NumberOfBGClusters, int nbins=12, int axis=0)
    {

      BG->LoadClusterDensities();

      string mode="fast";
      string axisname; double min, max;
      switch(axis)
      {
         case 0: axisname="X Drift (cm)"; min=-600.; max=600.; break;
         case 1: axisname="Y (cm)"; min=-600.; max=600.; break;
         case 2: axisname="Z (cm)"; min=0.; max=6000.;   break;
         default: axis=0; axisname="X (cm)"; min=-600.; max=600.;
      }
      TProfile *eff = new TProfile("MatchEffScan",Form("Matching Efficiency Scan;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      noBGeff= TProfile("MatchEffScanNoBG",Form("Matching Efficiency Scan - No Background;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      noSIGeff= TProfile("MatchEffScanNoSIG",Form("Matching Efficiency Scan - Assuming perfect signal detection;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      NumberOfBGClusters= TProfile("NumberOfBGClusters",Form("NumberOfBGClusters Scan;%s;Number Of BGClusters in vecinity",axisname.c_str()),nbins, min, max);

      for(int i=0;i<Signal->t->GetEntries();i++)
      {
         Signal->t->GetEntry(i);
         float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};
         if(Signal->ClustersVector->size()>0)
         {
        
           double NumberOfClustersAt;
	   if(mode=="full") NumberOfClustersAt= BG->GetNumberOfClustersAt(TruePos,distance,0.8*Signal->ClustersVector->at(0).PEs()); //amount of clusters that matches with this event
           else NumberOfClustersAt = BG->GetNumberOfClustersAtFast(TruePos,distance,0.5*Signal->ClustersVector->at(0).PEs()); //amount of clusters that matches with this event
           double TotalBGWindow = BG->t->GetEntries()*1e-3; //seconds of BG sample
           double ReadoutWindow = 16.5e-3; //seconds
           double AverageNumberOfBackgroundClustersPerReadoutWindow = NumberOfClustersAt*ReadoutWindow/TotalBGWindow;
           double dy =Signal->VertY - Signal->ClustersVector->at(0).pos[1];
           double dz =Signal->VertZ - Signal->ClustersVector->at(0).pos[2];
           double ProbabilityToBeMatchedWell;

//           cout << AverageNumberOfBackgroundClustersPerReadoutWindow <<endl; lets_pause();

           if(dy*dy+dz*dz > distance*distance) {ProbabilityToBeMatchedWell = 0;noBGeff.Fill(TruePos[axis],0.0);}
           else {ProbabilityToBeMatchedWell = 1.0/(1.0+AverageNumberOfBackgroundClustersPerReadoutWindow);noBGeff.Fill(TruePos[axis],1.0);}

           if(i%10000==0) std::cout << "Event " << i << " out of " << Signal->t->GetEntries() << endl;
           NumberOfBGClusters.Fill(TruePos[axis],AverageNumberOfBackgroundClustersPerReadoutWindow);
           eff->Fill(TruePos[axis],ProbabilityToBeMatchedWell);
           noSIGeff.Fill(TruePos[axis],1.0/(1.0+AverageNumberOfBackgroundClustersPerReadoutWindow));

         }
         else
         {
           eff->Fill(TruePos[axis],0.0);
           noBGeff.Fill(TruePos[axis],0.0);
           noSIGeff.Fill(TruePos[axis],0.0);
         }
         //if (i==10)break;
      }
      cout << "leaving"<<endl;
     return eff;
    }
    TProfile* MatchingEfficiencySplittingScan(double distance, TProfile &noBGeff, TProfile &noSIGeff, TProfile &NumberOfBGClusters, int nbins=12, int axis=0)
    {

      BG->LoadClusterDensities();
      string axisname; double min, max;
      switch(axis)
      {
         case 0: axisname="X Drift (cm)"; min=-600.; max=600.; break;
         case 1: axisname="Y (cm)"; min=-600.; max=600.; break;
         case 2: axisname="Z (cm)"; min=0.; max=6000.;   break;
         default: axis=0; axisname="X (cm)"; min=-600.; max=600.;
      }
      TProfile *eff = new TProfile("MatchEffScan",Form("Matching Efficiency Scan;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      noBGeff= TProfile("MatchEffScanNoBG",Form("Matching Efficiency Scan - No Background;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      noSIGeff= TProfile("MatchEffScanNoSIG",Form("Matching Efficiency Scan - Assuming perfect signal detection;%s;Matching effiency",axisname.c_str()),nbins, min, max);
      NumberOfBGClusters= TProfile("NumberOfBGClusters",Form("NumberOfBGClusters Scan;%s;Number Of BGClusters in vecinity",axisname.c_str()),nbins, min, max);

      int NReadOutWindows = BG->GetNReadOutWindows();
   
      std::cout << "Looping over " << NReadOutWindows << " readout windows." << endl;
      std::vector<double> NumberOfClustersVector;
      for(int i=0;i<Signal->t->GetEntries();i++)
      {
         Signal->t->GetEntry(i);

         float TruePos[3]={Signal->VertX,Signal->VertY,Signal->VertZ};

         if(Signal->ClustersVector->size()>0)
         {
           NumberOfClustersVector = BG->GetNumberOfClustersAtFast3D(TruePos,distance,0.8*Signal->ClustersVector->at(0).PEs());
           if(i%10000==0) std::cout << "Event " << i << " out of " << Signal->t->GetEntries() << endl;
           for( int b=0; b<NReadOutWindows; b++)
           {
              double NumberOfClustersAt = NumberOfClustersVector[b]; //amount of clusters that matches with this event
 
              double dy =Signal->VertY - Signal->ClustersVector->at(0).pos[1];
              double dz =Signal->VertZ - Signal->ClustersVector->at(0).pos[2];
              double ProbabilityToBeMatchedWell;

              if(dy*dy+dz*dz > distance*distance) {ProbabilityToBeMatchedWell = 0; noBGeff.Fill(TruePos[axis],0.0);}
              else {ProbabilityToBeMatchedWell = 1.0/(1.0+NumberOfClustersAt); noBGeff.Fill(TruePos[axis],1.0);}

              

              eff->Fill(TruePos[axis],ProbabilityToBeMatchedWell);
              noSIGeff.Fill(TruePos[axis],1.0/(1.0+NumberOfClustersAt));
              NumberOfBGClusters.Fill(TruePos[axis],NumberOfClustersAt);
  

              if(debug){ std::cout <<"NumberOfClustersAt " <<NumberOfClustersAt << " "<<b <<endl;}
           }

         }
         else
         {
           for( int b=0; b<NReadOutWindows; b++)
           {
             eff->Fill(TruePos[axis],0.0);
             noBGeff.Fill(TruePos[axis],0.0);
             noSIGeff.Fill(TruePos[axis],0.0);
           }
         }
         //if (i==10)break;
      }
      cout << "leaving"<<endl;
     return eff;
    }

  
  };

}
