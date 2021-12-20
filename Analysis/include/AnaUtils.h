#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_

#include "AnaIO.h"

class AnaUtils
{
  public:
    // Constructor & Destructor
    AnaUtils() {};
    ~AnaUtils() {};

    int GetBeamParticleType(const int pdg);
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
    // Return the transverse momentum relative to the beam direction
    double GetTransverseMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton);
    // Return a Pi0 LT vector relative to beam direction (only consider theta)
    TLorentzVector GetPi0MomentumRefBeam(const TLorentzVector RecPi0, const TLorentzVector TruthPi0, const bool isTruth);

    void TruthMatchingTKI(TLorentzVector dummypi0, TLorentzVector dummyproton, TLorentzVector dummypi0Truth, TLorentzVector dummyprotonTruth, const bool kMC);

    // Fill different reco FS particle theta and momentum and truth-matched resolution
    void FillFSParticleKinematics(const int recIndex, const int truthParticleType, const int recParticleType);      
    // Get reco shower distance vector from vertex to start point 
    TVector3 GetRecShowerDistVector(const int ii);

    TVector3 GetTruthMatchedShower3VectLab(const int ii);
    // Get LT vector of truth-matched shower in lab frame
    TLorentzVector GetTruthMatchedShowerLTVectLab(const int ii);

    TVector3 GetRecShower3VectLab(const int ii, bool DoCorrection = true);
    // Get LT vector of reco shower in lab frame
    TLorentzVector GetRecShowerLTVectLab(const int ii, bool DoCorrection = true);
    // Get LT vector of reco shower relative to the beam direction
    TLorentzVector GetRecShowerRefBeam(const bool isTruth, const int ii, bool DoCorrection = true);

    // Combine two showers to reconstruct pi0
    TLorentzVector GetRecPiZeroFromShowers();

    vector<TLorentzVector> GetTwoPi0Showers();

    vector<TLorentzVector> GetTwoTruthMatchedPi0Showers(int &truthPi0Type);

    // Get the info for Fitting
    void SavePi0ShowersForKF();
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

    enum BeamparType{

      //1-5
      gkBeamProton=1,
      gkBeamPiPlus,
      gkBeamMuPlus,
      gkBeamElectronGamma,
      gkBeamPiMinus,

      //6
      gkBeamOthers,
      
    }; // End of parType


    // Define event types (signal/eventBk/BeamBk)
    enum evtType{
      gkSignal = 0,
      gkEvtBkg,
      gkBmBkg
    };

    enum truthPi0Type{
      gkTwoGammas = 0,
      gkOneGamma,
      gkNoGammas
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

    static Double_t CauchyDens(Double_t *x, Double_t *par)
    {
      Double_t pi   = TMath::Pi();
      Double_t mean = par[0];
      Double_t fwhm = par[1];

      Double_t arg = x[0]-mean;
      Double_t top = fwhm;
      Double_t bot = pi*(arg*arg+top*top);

      Double_t func = top/bot;
      return func;
    }

    static Double_t CauchyPeak(Double_t *x, Double_t *par)
    {
      Double_t height = par[2];
      Double_t func = height*CauchyDens(x,par);
      return func;
    }

    double GetProtonCorrectedMom(double rawMom){
      TF1 *fpCor = new TF1("fpCor","pol5",0,2);
      fpCor->SetParameter(0,3.63507);
      fpCor->SetParameter(1,-28.0287);
      fpCor->SetParameter(2,81.5051);
      fpCor->SetParameter(3,-113.576);
      fpCor->SetParameter(4,76.2824);
      fpCor->SetParameter(5,-19.8619);
      const double factor = fpCor->Eval(rawMom);
      return rawMom/(1+factor);
    }

    static Double_t ShowerEenergyFCN(Double_t *x, Double_t *par)
    {
      return par[3] + (par[0] - par[3])/(1 + pow((x[0]/par[2]),par[1]));
    }

    double GetShowerCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",ShowerEenergyFCN,0,1,4);
      fpCor->SetParameter(0,-0.5921); 
      fpCor->SetParameter(1,2.236); 
      fpCor->SetParameter(2,0.1076);
      fpCor->SetParameter(3,-0.1334);

      //fpCor->SetParameter(0,-65.13); 
      //fpCor->SetParameter(1,2.722); 
      //fpCor->SetParameter(2,0.01652);
      //fpCor->SetParameter(3,-0.1459);

      const double factor = fpCor->Eval(rawE);
      return rawE/(1+factor);
    }

    static vector<double> LdShowerEnergyTruth;
    static vector<double> SlShowerEnergyTruth; 
    static vector<double> OpenAngleTruth; 
    static vector<double> LdShowerEnergyRaw; 
    static vector<double> SlShowerEnergyRaw; 
    static vector<double> OpenAngle;
    static vector<TVector3> LdShowerDirTruth;
    static vector<TVector3> SlShowerDirTruth;
    static vector<TVector3> LdShowerDir;
    static vector<TVector3> SlShowerDir;

    // Energy Correction
    static vector<double> ProtonMomTruth;
    static vector<double> ProtonMomRaw;

    // TKI analysis
    static TLorentzVector RecPi0LTVet;
    static TLorentzVector TruthPi0LTVet;
    static TLorentzVector RecProtonLTVet;
    static TLorentzVector TruthProtonLTVet;
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

vector<TVector3> AnaUtils::LdShowerDirTruth;
vector<TVector3> AnaUtils::SlShowerDirTruth;

vector<TVector3> AnaUtils::LdShowerDir;
vector<TVector3> AnaUtils::SlShowerDir;

vector<double> AnaUtils::ProtonMomTruth;
vector<double> AnaUtils::ProtonMomRaw;

TLorentzVector AnaUtils::RecPi0LTVet;
TLorentzVector AnaUtils::RecProtonLTVet;
TLorentzVector AnaUtils::TruthPi0LTVet;
TLorentzVector AnaUtils::TruthProtonLTVet;


#endif
