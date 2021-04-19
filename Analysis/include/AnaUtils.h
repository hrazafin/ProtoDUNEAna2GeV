#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_

#include "AnaIO.h"

class AnaUtils
{
  public:
    // Constructor & Destructor
    AnaUtils() {};
    ~AnaUtils() {};
    // Tell you the particle type from pdg (only primary particle type)
    int GetParticleType(const int pdg);
    // Set if event is a signal based on phase space cuts on protons 
    void SetFullSignal();
    // Get a vector of true Final-State(FS) particles 4 vector with Pi0 + leading + subleading proton 
    vector<TLorentzVector> GetFSParticlesTruth(); 
    // Check event topology for signal
    bool IsSignal(const int nProton, const int nPiZero, const int nPiPlus, const int nParticleBkg);
    // Check if event is Signal/EvtBkg/BmBkg(beamBkg) using Truth info
    int GetFillEventType();

    void FillBeamKinematics(const int kMC);
    TVector3 GetTruthBeamFull();
    TVector3 GetRecBeamFull();

    vector<double> GetdEdxVector(const vector<double> &arraydEdx, const bool kForward);
    double GetTruncatedMean(const vector<double> & dEdxVec, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac);
    double GetFSParticleTME(const unsigned int ii, const bool kForward);
    double GetChi2NDF(const int ii);
 
    TVector3 GetTruthMatchedTrackVectLab(const int ii);
    TVector3 GetRecTrackVectLab(const int ii, const bool kProton);
    TLorentzVector GetMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton);
    void FillFSParticleKinematics(const int recIndex, const int truthParticleType, const int recParticleType);      
    
    TVector3 GetRecShowerDistVector(const int ii);
    TVector3 GetTruthMatchedShowerVectLab(const int ii);
    TVector3 GetRecShowerVectLab(const int ii);
    TLorentzVector GetRecShowerLTVectLab(const int ii);
    TLorentzVector GetShowerMomentumRefBeam(const bool isTruth, const int recIndex);
    
    TLorentzVector GetPiZero();
    // Define particle types
    enum parType{

      //1-4
      gkProton=1,
      gkPiPlus,
      gkPiMinus,
      gkGamma,
      
      //5-10
      gkSecondaryProton,
      gkSecondaryPiPlus,
      gkSecondaryPiMinus,
      gkSecondaryGamma,
      gkSecondaryEplusEminus,
      gkSecondaryMuon,
      
      //11
      gkOthers,
      
      //12-17
      gkNeutron,
      gkPiZero,
      gkElectron,
      gkPositron,
      gkMuPlus,
      gkMuMinus,
      
      //18
      gkKaon,
      
      //19
      gkNeutrino,
      
      //20
      gkHyperon,
      
      //21
      gkNucleus,
      
      //22
      gkNoTruth
    }; // End of parType

    // Define event types (signal/eventBk/BeamBk)
    enum evtType{
      gkSignal = 0,
      gkEvtBkg,
      gkBmBkg
     };

    void SavePiZeroShower(TLorentzVector shower, double showerE, TVector3 pos, int showerType)
    {
      showerArray.push_back(shower);
      showerEarr.push_back(showerE);
      showerTypeArray.push_back(showerType);
      showerPos.push_back(pos);
    }
    void CleanShowerArray(){showerArray.clear();showerEarr.clear();showerPos.clear();showerTypeArray.clear();}
  private:
    PlotUtils plotUtils;
    int nProton;
    int nNeutron;
    int nPiPlus;
    int nPiZero;
    int nGamma;
    int nParticleBkg;

    vector<TLorentzVector> showerArray;
    vector<int> showerTypeArray;
    vector<double> showerEarr;
    vector<TVector3> showerPos;
};
#endif
