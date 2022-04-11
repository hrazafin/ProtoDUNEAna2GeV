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
    // Get a vector of true pi0 showers 4 vector 
    vector<TLorentzVector> GetFSPiZeroDecayDaughterTruth(); 
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
    TVector3 GetRecTrackVectLab(const int ii, const bool kProton, bool DoCorrection = true);
    // Return a LT vector relative to beam direction (only consider theta)
    TLorentzVector GetMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton, bool DoProtonMomCorrection = true);
    // Return the transverse momentum relative to the beam direction
    double GetTransverseMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton, bool DoProtonMomCorrection = true);
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
    TLorentzVector GetRecPiZeroFromShowers(double &OA, bool kMC, bool kFill, bool kBefore, int &truthPi0Type);

    vector<TLorentzVector> GetTwoPi0Showers(double &separation, bool kMC, bool kFill, bool kBefore);

    vector<TLorentzVector> GetTwoTruthMatchedPi0Showers(int &truthPi0Type, const bool &kFill);

    // Get the info for Fitting
    void SavePi0ShowersForKF();
    // Truth TKI calculation
    void DoTruthTKICalculation();

    // Do Kinematic Fitting
    void DoKinematicFitting();

    void KF(const TLorentzVector &ldShower, const TLorentzVector &slShower, vector<double> &FittedVars);

    // Set the value of CVM
    void SetCVM();

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
      gkMesons,
      
      //23
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
      gkTwoGammasSamePi0 = 0,
      gkTwoGammasDiffPi0,
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

    

    static Double_t CorrectionFCN(Double_t *x, Double_t *par)
    {
      return par[3] + (par[0] - par[3])/(1 + pow((x[0]/par[2]),par[1]));
    }

    double GetProtonCorrectedMom(double rawMom){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,2,4);
      //TF1 *fpCor = new TF1("fpCor","pol5",0,2);
      fpCor->SetParameter(0,-3.815e4);
      fpCor->SetParameter(1,7.181);
      fpCor->SetParameter(2,0.07183);
      fpCor->SetParameter(3,-0.005464);

      const double factor = fpCor->Eval(rawMom);
      return rawMom/(1+factor);
    }

    double GetShowerCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,1,4);
      /*fpCor->SetParameter(0,-0.9257); // Old
      fpCor->SetParameter(1,1.591); 
      fpCor->SetParameter(2,0.06717);
      fpCor->SetParameter(3,-0.1116);

      fpCor->SetParameter(0,-0.8409); // Pandora Type
      fpCor->SetParameter(1,1.442); 
      fpCor->SetParameter(2,0.06301);
      fpCor->SetParameter(3,-0.1096);
*/
      fpCor->SetParameter(0,-0.8995); 
      fpCor->SetParameter(1,1.442); 
      fpCor->SetParameter(2,0.06068);
      fpCor->SetParameter(3,-0.1096);

      const double factor = fpCor->Eval(rawE);
      return rawE/(1+factor);
    }

    double GetShowerCorrectedTheta(double rawTheta){
      TF1 *fpCor = new TF1("fpCor","gaus",0,180);
      /*fpCor->SetParameter(0,5.163);
      fpCor->SetParameter(1,66.37);
      fpCor->SetParameter(2,26.14); 

      fpCor->SetParameter(0,5.169);
      fpCor->SetParameter(1,66.9);
      fpCor->SetParameter(2,28.14);  
*/
      fpCor->SetParameter(0,5.102);
      fpCor->SetParameter(1,65.83);
      fpCor->SetParameter(2,26.4); 

      const double factor = fpCor->Eval(rawTheta);
      return rawTheta - factor;
    }

    double GetShowerCorrectedPhi(double rawPhi){
      TF1 *fpCor = new TF1("fpCor","pol7",0,180);
      /*fpCor->SetParameter(0,-0.0635279);
      fpCor->SetParameter(1,0.0941388);
      fpCor->SetParameter(2,-4.78192e-05);
      fpCor->SetParameter(3,-2.0647e-05);
      fpCor->SetParameter(4,7.65803e-09);
      fpCor->SetParameter(5,1.12047e-09);
      fpCor->SetParameter(6,-1.93702e-13); 
      fpCor->SetParameter(7,-1.76401e-14);

      fpCor->SetParameter(0,-0.298625);
      fpCor->SetParameter(1,0.100689);
      fpCor->SetParameter(2,6.19076e-05);
      fpCor->SetParameter(3,-2.20382e-05);
      fpCor->SetParameter(4,-2.50278e-09);
      fpCor->SetParameter(5,1.20338e-09);
      fpCor->SetParameter(6,5.23649e-14); 
      fpCor->SetParameter(7,-1.9246e-14);
*/
      fpCor->SetParameter(0,-0.0942805);
      fpCor->SetParameter(1,0.0915092);
      fpCor->SetParameter(2,-3.00198e-05);
      fpCor->SetParameter(3,-2.01208e-05);
      fpCor->SetParameter(4,5.16471e-09);
      fpCor->SetParameter(5,1.0919e-09);
      fpCor->SetParameter(6,-1.19572e-13); 
      fpCor->SetParameter(7,-1.72078e-14);

      const double factor = fpCor->Eval(rawPhi);
      return rawPhi - factor;
    }

    double GetProtonCorrectedTheta(double rawTheta){
      TF1 *fpCor = new TF1("fpCor","gaus",0,180);
      fpCor->SetParameter(0,4.601);
      fpCor->SetParameter(1,50.44);
      fpCor->SetParameter(2,24.42);  
      const double factor = fpCor->Eval(rawTheta);
      return rawTheta - factor;
    }

    double GetProtonCorrectedPhi(double rawPhi){
      TF1 *fpCor = new TF1("fpCor","pol5",0,180);
      fpCor->SetParameter(0,-1.16553);
      fpCor->SetParameter(1,-0.00766996);
      fpCor->SetParameter(2,0.00019561);
      fpCor->SetParameter(3,1.012e-07);
      fpCor->SetParameter(4,-4.78252e-09);
      fpCor->SetParameter(5,6.90576e-12);   
      const double factor = fpCor->Eval(rawPhi);
      return rawPhi - factor;
    }

    double GetLDShowerCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,1,4);
      fpCor->SetParameter(0,-0.412); 
      fpCor->SetParameter(1,2.002); 
      fpCor->SetParameter(2,0.1559);
      fpCor->SetParameter(3,-0.1185);

      const double factor = fpCor->Eval(rawE);
      return rawE/(1+factor);
    }

    double GetSLShowerCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,1,4);
      fpCor->SetParameter(0,-0.4897); 
      fpCor->SetParameter(1,3.457); 
      fpCor->SetParameter(2,0.1216);
      fpCor->SetParameter(3,-0.1419);

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

    // KF CVM 
    static std::map<std::pair <int,int>,vector<double>> CVM;
    static int selected;
    static int total;

  private:
    PlotUtils plotUtils;
    int nProton;
    int nNeutron;
    int nPiPlus;
    int nPiZero;
    int nGamma;
    int nParticleBkg;
    
    vector<int> bufferType;

    int nPiZeroGamma;
    int nPiZeroElectron;
    int nPiZeroPositron;

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

std::map<std::pair <int,int>,vector<double>> AnaUtils::CVM;

int AnaUtils::selected = 0;
int AnaUtils::total = 0;

#endif
