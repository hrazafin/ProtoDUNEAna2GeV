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
    int GetChannelType(const TString beamEndprocess, const int beamPDG, const int npi0Daughters, const int npiplusDaughters, const int npiminusDaughters);

    int GetFillXSEventType();

    // Tell you the particle type from pdg (only primary particle type)
    int GetParticleType(const int pdg);
    // Set if event is a signal based on phase space cuts on protons 
    void SetFullSignal();
    // Get a vector of true Final-State(FS) particles 4 vector with Pi0 + leading + subleading proton 
    vector<TLorentzVector> GetFSParticlesTruth(bool kFill = true); 
    // Get a vector of true pi0 showers 4 vector 
    vector<TLorentzVector> GetFSPiZeroDecayDaughterTruth(bool kFill = true); 
    // Get truth interacting histogram for xsec calculation
    double MakeTrueIncidentEnergies(vector<double> *true_beam_traj_Z, vector<double> *true_beam_traj_KE, vector<double> *true_beam_new_incidentEnergies);
    double MakeRecoIncidentEnergies(vector<double> *reco_beam_traj_Z, vector<double> *reco_beam_traj_KE, vector<double> *reco_beam_new_incidentEnergies);
    
    // Get reco track length from the space points
    double GetRecoTrackLength();
    // Get true track length from caloXYZ
    double GetTrueTrackLength();

    int GetInitialSliceID();
    int GetInteractionSliceID();

    bool IsFakeData(){
      // Add fake data type for xsec study
      const int event = AnaIO::event;
      if(event%2) return false;
      else return true;
    }
    

    // Check event topology for signal
    //bool IsSignal(const int nProton, const int nPiZero, const int nPiPlus, const int nPiMinus, const int nParticleBkg);
    bool IsSignal(const int nProton, const int nPiZero, const int nPiPlus, const int nPiMinus, const int nParticleBkg, const int nPionAboveThreshold);

    // Check if event is Signal/EvtBkg/BmBkg(beamBkg) using Truth info
    int GetFillEventType();
    // Decompose events into NpMn etc. for TKI analysis
    int GetFillTKIEventType();
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

    void TruthMatchingTKI(TLorentzVector dummypi0, TLorentzVector dummyproton, TLorentzVector dummypi0Truth, TLorentzVector dummyprotonTruth, const bool kMC,  const bool GoodTruthMatch);

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

    double CalWeight(const bool & kMC);
    double CalBckWeight(const bool & kMC);
    double CalXSEvtWeight(const bool & kMC, const double & intE, const int & evtXStype);
    double CalBeamIniWeight(const double & iniE);
    double CalBeamIntWeight(const double & intE);
    double CalCEXIntWeight(const double & intE);
    double CalCEXPi0KEWeight(const double & intE);

    double GetUpStreamEnergyLoss(const bool & kMC, const double & InstE);

    // Get truth particle PDG from it's ID
    int GetTruthPDGFromID(const int inID, const vector<int> * idarray, const vector<int> * pdgarray);
    string GetTruthProcessFromID(const int inID, const vector<int> * idarray, const vector<string> * processarray);
    string GetTruthEndProcessFromID(const int inID, const vector<int> * idarray, const vector<string> * Endprocessarray);
    // Get truth-matched particle info using reco particle index
    int GetTruthParticleInfoFromRec(const int recidx);

    // Used in AnaCut beam energy related
    // Fill beam quality cuts histograms (use fitted value for beam quality cut latter)
    void FillBeamQualityHist();
    // Set beam energy for upstream energy loss study
    void SetBeamInstKEandFrontFaceKE(double &beam_inst_KE, double &true_ffKE, bool kFill = false);
    void FillUpStreamEnergyLossHistBeforeCut(double beam_inst_KE, double true_ffKE);
    void FillUpStreamEnergyLossHistAfterCut(double beam_inst_KE, double true_ffKE);
    void FillBeamVariablesAfterAllCuts(int parType, int channelType);

    double GetLoss(double InstE){
      TF1 *fpCor = new TF1("fpCor","pol2",-400,1400);
      fpCor->SetParameter(0,170.777);
      fpCor->SetParameter(1,-0.595308);
      fpCor->SetParameter(2,0.000456932);
      const double factor = fpCor->Eval(InstE);
      return factor;
    }


    void FillXSTrueHistograms();
    void FillXSRecoHistograms();

    //void FillEsliceHistograms(double KE_init, double KE_end, double KE_int, double weight, const vector<double> binning, int N_bin, bool fill_int);
    void FillEsliceHistograms(TH1D* hinit, TH1D* hend, TH1D* hinc, TH1D* hint, double KE_init, double KE_end, double KE_int, double weight, const vector<double> binning, int N_bin, bool fill_int);

    //vector<vector<double>> ComputeTrueIncidentHist(const double & initialE, const double & interactingE, const double & ffe, const double & intE);
    vector<vector<double>> ComputeTrueIncidentHist(const double & initialE, const double & interactingE, const double & ffe, const double & intE, vector<double> *true_beam_incidentEnergies);


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
      gkBeamPiPlus=0,
      gkBeamMuon,
      gkMisIDProton,
      gkMisIDPion,
      gkMisIDMuon,
      gkMisIDElectronGamma,
      gkCosmicMuon,
      gkBeamOthers,
      
    }; // End of BeamparType


    enum ChannelType{

      //1-4
      gkChargeExchange=0,
      gkInelastic,
      gkAbsorption,
      gkPionDecays,
      gkOtherChannels,
      gkBeamMuons,
      gkBackground,
      
    }; // End of ChannelType


    // Define event types (signal/eventBk/BeamBk)
    enum evtType{
      gkSignal = 0,
      gkEvtBkg,
      gkBmBkg
    };

    enum evtXSType{
      gkXSSignal = 0,
      gkXSEvtBkgAbs,
      gkXSEvtBkgInel,
      gkXSEvtBkgSinglePi0,
      gkXSEvtBkgMultiPi0,
      //gkXSEvtBkgOthers,
      gkXSBmBkg
    };

    enum truthPi0Type{
      gkTwoGammasSamePi0 = 0,
      gkTwoGammasDiffPi0,
      gkOneGamma,
      gkNoGammas
    };

    enum TKIevtType{
      gk1p0n = 0,
      gk1pMn,
      //gkNp0n,
      gkNpMn,
      gkEvtTKIBkg,
      gkBmTKIBkg
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

    vector<double> GetRecoTrackLengthAccumVect(){
      return reco_trklen_accum;
    }

    vector<double> GetTrueTrackLengthAccumVect(){
      return true_trklen_accum;
    }

    double GetTrueFrontFaceEnergy(){
      return true_ffKE;
    }

    double GetTrueIntEnergy(){
      return int_energy_true;
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
      fpCor->SetParameter(0,-8.35387e+04);
      fpCor->SetParameter(1,7.79370e+00);
      fpCor->SetParameter(2,7.67495e-02);
      fpCor->SetParameter(3,-5.68081e-03);

      const double factor = fpCor->Eval(rawMom);
      return rawMom/(1+factor);
    }

    double GetShowerCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,1,4);
      fpCor->SetParameter(0,-0.8323); 
      fpCor->SetParameter(1,1.665); 
      fpCor->SetParameter(2,0.06923);
      fpCor->SetParameter(3,-0.1219);

      const double factor = fpCor->Eval(rawE);
      return rawE/(1+factor);
    }

    double GetShowerCorrectedTheta(double rawTheta){
      TF1 *fpCor = new TF1("fpCor","gaus",0,180);
      fpCor->SetParameter(0,5.033);
      fpCor->SetParameter(1,66.61);
      fpCor->SetParameter(2,26.44); 

      const double factor = fpCor->Eval(rawTheta);
      return rawTheta - factor;
    }

    double GetShowerCorrectedPhi(double rawPhi){
      TF1 *fpCor = new TF1("fpCor","pol7",0,180);
      fpCor->SetParameter(0,-0.00806448);
      fpCor->SetParameter(1,0.0978646);
      fpCor->SetParameter(2,-0.000136381);
      fpCor->SetParameter(3,-2.12771e-05);
      fpCor->SetParameter(4,1.61586e-08);
      fpCor->SetParameter(5,1.14406e-09);
      fpCor->SetParameter(6,-3.82649e-13); 
      fpCor->SetParameter(7,-1.78852e-14);

      const double factor = fpCor->Eval(rawPhi);
      return rawPhi - factor;
    }

    double GetProtonCorrectedTheta(double rawTheta){
      TF1 *fpCor = new TF1("fpCor","gaus",0,180);
      fpCor->SetParameter(0,4.329);
      fpCor->SetParameter(1,50.98);
      fpCor->SetParameter(2,25.21);  
      const double factor = fpCor->Eval(rawTheta);
      return rawTheta - factor;
    }

    double GetProtonCorrectedPhi(double rawPhi){
      TF1 *fpCor = new TF1("fpCor","pol5",0,180);
      fpCor->SetParameter(0,-1.12529);
      fpCor->SetParameter(1,-0.00592945);
      fpCor->SetParameter(2,0.000186931);
      fpCor->SetParameter(3,-1.4181e-07);
      fpCor->SetParameter(4,-4.49753e-09);
      fpCor->SetParameter(5,1.38988e-11);   
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

    double GetBeamCorrectedPhi(double rawPhi){
      TF1 *fpCor = new TF1("fpCor","pol3",-180,-50);
      fpCor->SetParameter(0,494.8);
      fpCor->SetParameter(1,10.18);
      fpCor->SetParameter(2,0.07259);
      fpCor->SetParameter(3,0.0001799);   
      const double factor = fpCor->Eval(rawPhi);
      return rawPhi - factor;
    }

    double GetBeamCorrectedTheta(double rawTheta){
      TF1 *fpCor = new TF1("fpCor","pol3",0,80);
      fpCor->SetParameter(0,1.865);
      fpCor->SetParameter(1,-0.6969);
      fpCor->SetParameter(2,0.04245);
      fpCor->SetParameter(3,-0.0003842);  
      const double factor = fpCor->Eval(rawTheta);
      return rawTheta - factor;
    }

    double GetBeamCorrectedE(double rawE){
      TF1 *fpCor = new TF1("fpCor",CorrectionFCN,0,1,4);
      fpCor->SetParameter(0,5.443e6); 
      fpCor->SetParameter(1,10.64); 
      fpCor->SetParameter(2,0.009267);
      fpCor->SetParameter(3,0.02225);

      const double factor = fpCor->Eval(rawE);
      return rawE/(1+factor);
    }

    double GetSliceWidth(){
      return slice_width;
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

    static bool GoodTruthMatch;
    static bool GoodFit;

    // KF CVM 
    static std::map<std::pair <int,int>,vector<double>> CVM;
    static int selected;
    static int total;

  private:
    PlotUtils plotUtils;
    int nProton;
    int nNeutron;
    int nPiPlus;
    int nPiMinus;
    int nPiZero;
    int nGamma;
    int nParticleBkg;
    // for xs analysis
    int nPionAboveThreshold;
    
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

    vector<double> reco_trklen_accum;
    vector<double> true_trklen_accum;

    double true_ffKE;
    double int_energy_true;

    // E-slice parameters
    const double Eslicewidth = 50; //MeV
    const double Eklim = 1000; //MeV
    const int nthinslices = 20;

    const double slice_width = 1.0;

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

bool AnaUtils::GoodTruthMatch = false;
bool AnaUtils::GoodFit = false;

#endif
