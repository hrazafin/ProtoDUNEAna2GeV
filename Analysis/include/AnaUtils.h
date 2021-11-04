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
    // Fill reco beam theta and momentum and truth-matched resolution
    void FillBeamKinematics(const int kMC);
    // Get the truth beam momentum directly at the end point 
    TVector3 GetTruthBeamFull();
    // Get the reco beam direction and KE at the end point and then find the beam momentum (end point)
    TVector3 GetRecBeamFull();
    // Get meaningful dEdx vector
    vector<double> GetdEdxVector(const vector<double> &arraydEdx, const bool kForward);
    // Calculate the truncated mean dEdx
    double GetTruncatedMean(const vector<double> & dEdxVec, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac);
    // Get the TME of a FS reco particle
    double GetFSParticleTME(const unsigned int ii, const bool kForward);
    // Get Chi2/NDF
    double GetChi2NDF(const int ii);
    // Get the truth-matched track vector in the lab frame
    TVector3 GetTruthMatchedTrackVectLab(const int ii);
    // Get the reco track vector in the lab frame
    TVector3 GetRecTrackVectLab(const int ii, const bool kProton);
    // Return a LT vector relative to beam direction (only consider theta)
    TLorentzVector GetMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton);
    // Return a Pi0 LT vector relative to beam direction (only consider theta)
    TLorentzVector GetPi0MomentumRefBeam(const TLorentzVector dummyPi0);

    void TruthMatchingTKI(TLorentzVector dummypi0, TLorentzVector dummyproton);

    // Fill different reco FS particle theta and momentum and truth-matched resolution
    void FillFSParticleKinematics(const int recIndex, const int truthParticleType, const int recParticleType);      
    // Get reco shower distance vector from vertex to start point 
    TVector3 GetRecShowerDistVector(const int ii);
    // Get LT vector of truth-matched shower
    TLorentzVector GetTruthMatchedShowerLTVectLab(const int ii);
    // Get LT vector of reco shower
    TLorentzVector GetRecShowerLTVectLab(const int ii, bool DoCorrection = true);
    // Combine two showers to reconstruct pi0
    TLorentzVector GetPiZero();
    // Get the info for Fitting
    void GetPi0Showers();
    // Truth TKI calculation
    void DoTruthTKICalculation();

    //void Chi2FCN(int &npars, double *grad, double &value, double *par, int flag);
    //void KinematicFitting(double openAngle, double E1, double E2, double sigmaE1, double sigmaE2);
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
    // Save good shower candidates info for pi0 reco
    void SavePiZeroShower(TLorentzVector shower, TLorentzVector showerRaw, TLorentzVector showerTruth, double showerE, double showerTruthE, TVector3 pos, int showerType)
    {
      showerArray.push_back(shower);
      showerArrayRaw.push_back(showerRaw);
      showerTruthArray.push_back(showerTruth);
      showerEarr.push_back(showerE);
      showerTruthEarr.push_back(showerTruthE);
      showerTypeArray.push_back(showerType);
      showerPos.push_back(pos);
    }
    // Clean vectors
    void CleanShowerArray(){
      showerArray.clear(); showerArrayRaw.clear(); showerTruthArray.clear();showerEarr.clear();showerTruthEarr.clear();showerPos.clear();showerTypeArray.clear();
    }
    // Get bufferType
    vector<int> GetBufferType(){
      return bufferType;
    }
    // Get FS particle number 
    vector<double> GetNParticles(){
      vector<double> NParList;
      NParList.push_back(nProton);
      NParList.push_back(nNeutron);
      NParList.push_back(nPiZero);
      return NParList;
    }

    static vector<double> LdShowerEnergyTruth;
    static vector<double> SlShowerEnergyTruth; 
    static vector<double> OpenAngleTruth; 
    static vector<double> LdShowerEnergyRaw; 
    static vector<double> SlShowerEnergyRaw; 
    static vector<double> OpenAngle;
    // TKI analysis
    static TLorentzVector RecPi0LTVet;
    static TLorentzVector RecProtonLTVet;

  private:
    PlotUtils plotUtils;
    int nProton;
    int nNeutron;
    int nPiPlus;
    int nPiZero;
    int nGamma;
    int nParticleBkg;
    
    vector<int> bufferType;

    vector<TLorentzVector> showerArray;
    vector<TLorentzVector> showerArrayRaw;
    vector<TLorentzVector> showerTruthArray;
    vector<int> showerTypeArray;
    vector<double> showerEarr;
    vector<double> showerTruthEarr;
    vector<TVector3> showerPos;
};

vector<double> AnaUtils::LdShowerEnergyTruth;
vector<double> AnaUtils::SlShowerEnergyTruth;
vector<double> AnaUtils::OpenAngleTruth;
vector<double> AnaUtils::LdShowerEnergyRaw;
vector<double> AnaUtils::SlShowerEnergyRaw;
vector<double> AnaUtils::OpenAngle;

TLorentzVector AnaUtils::RecPi0LTVet;
TLorentzVector AnaUtils::RecProtonLTVet;

#endif
