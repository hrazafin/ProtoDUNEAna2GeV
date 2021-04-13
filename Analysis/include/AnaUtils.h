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

    bool GetgkSignal(){return gkSignal;}
  private:
    PlotUtils plotUtils;
    bool Signal;
    int nProton;
    int nNeutron;
    int nPiPlus;
    int nPiZero;
    int nGamma;
    int nParticleBkg;
};
#endif
