#ifndef _ANACUT_H_
#define _ANACUT_H_

#include "AnaUtils.h"
class AnaCut
{
  public:
    //Constructor & Destructor
    AnaCut() {};
    ~AnaCut() {};
    // All beam cut for both MC and data including 1.beam ID cut; 2.primary beam type cut; 3.beam position cut; 4.APA3 cut
    bool CutBeamAllInOne(const bool kMC, const bool kStopMuon = false, const bool kFill = false);
    bool CutBeamLongTrack(const bool kMC, const bool kFill = false);

    // Cut on beam PDG code
    bool CutBeamPDG(const bool kMC);
    bool CutBeamPDG_StopMuon(const bool kMC);

    // Cut on Pandora slice
    bool CutPandoraSlice();
    // Cut on Calo size
    bool CutCaloSize();
    // Cut on beam quality
    bool CutBeamQuality(const bool kMC, bool DoAngleCut = true, const bool kFill = false);
    bool CutBeamInstQuality(bool kMC);
    // Cut on end Z APA3
    bool CutAPA3EndZ(const bool kMC, const bool kFill = false, const double weight = 1.0);
    // Cut on Michel Score (remove muons)
    bool CutMichelScore(const bool kMC, const bool kFill = false, const double weight = 1.0);
    // Cut on Michel Score (stop muons)
    bool CutMichelScore_StopMuon(const bool kMC, const bool kFill = false, const double weight = 1.0);
    // Cut on median dE/dx (remove protons)
    bool CutMediandEdx(const bool kMC);
    // Cut on median dE/dx (remove protons)
    bool CutProtonChi2DOF(const bool kMC, bool kFill = false, const double weight = 1.0);
    // Cut on beam scrapers
    bool CutBeamScraper(const bool kMC);

    // --------------------------------- Legacy Cuts (not used anymore replaced by beam quality cuts) ----------------------- //
    // Cut on beam position for MC
    bool CutMCBeamPos(); 
    // Cut on beam position for data
    bool CutDataBeamPos();
    // Beam position for MC 
    bool Manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, 
                           const double beam_dirX, const double beam_dirY,   const double beam_dirZ, 
                           const double true_dirX,   const double true_dirY, const double true_dirZ,   
                           const double true_startX, const double true_startY, const double true_startZ);
    // Beam position for data
    bool Manual_beamPos_data(const int event,            const double data_startX,
                             const double data_startY,   const double data_startZ,
                             const double data_dirX,     const double data_dirY,
                             const double data_dirZ,     const double beam_inst_X,
                             const double beam_inst_Y,     const double beam_inst_dirX,
                             const double beam_inst_dirY,  const double beam_inst_dirZ,
                             const int beam_inst_nMomenta, const int beam_inst_nTracks);

    // --------------------------------- Legacy Cuts (not used anymore replaced by beam quality cuts) ----------------------- //

    // Event and final state particle selections
    bool CutTopology(const bool kMC, const bool kFill = false);
    void CountPFP(const bool kMC, const bool kFill);
    bool IsProton(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);
    bool IsTrack(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);  
    bool IsPionTrack(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);  
    bool PassProtonSubPID(const int ii, const double weight = 1.0, const bool kFill = false);
    bool PassPionSubPID(const int ii, const double weight = 1.0, const bool kFill = false);
    bool IsPiplus(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);
    bool IsShower(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);
    bool IsMichel(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false);
    bool IsPiZeroShower(const int ii, const bool kMC, const double weight = 1.0, const bool kFill = false); 
    bool IsPizero(const bool kMC, const bool kFill = false, const double weight = 1.0);

    static double pi0KineticE;
    static double pi0costheta;

  private:

    AnaUtils anaUtils;
    PlotUtils plotUtils;
    int nproton;
    int npiplus;
    int nshower;
    int npi0shower;
    int nmichel;
    int npi0;
    int truthParticleType;
    int recParticleType;

    double OA;
    int truthPi0Type;
    TLorentzVector PiZeroVec;
    

    // Cut values
    const double cut_EndZ_APA3 = 220.0; // in cm
    const double cut_MichelScore = 0.55;
    // Beam quality cuts
    const double MCmeanStartX = -30.85; // in cm
    const double MCmeanStartY = 422.4;  
    const double MCmeanStartZ = 0.1146; 
    const double MCmeanThetaX = 101.6 * TMath::DegToRad(); // in rad.
    const double MCmeanThetaY = 101.2 * TMath::DegToRad();
    const double MCmeanThetaZ = 16.61 * TMath::DegToRad();
    const double MCsigmaStartX = 5.043;
    const double MCsigmaStartY = 4.535;
    const double MCsigmaStartZ = 0.2166;

    const double DATAmeanStartX = -28.35; // in cm
    const double DATAmeanStartY = 424.6;
    const double DATAmeanStartZ = 3.175;
    const double DATAmeanThetaX = 100.6 * TMath::DegToRad(); // in rad.
    const double DATAmeanThetaY = 103.6 * TMath::DegToRad();
    const double DATAmeanThetaZ = 17.79 * TMath::DegToRad();
    const double DATAsigmaStartX = 4.636;
    const double DATAsigmaStartY = 5.217;
    const double DATAsigmaStartZ = 1.326;

    const double MCmeanInstX = -29.28;
    const double MCsigmaInstX = 4.127;
    const double DATAmeanInstX = -30.66;
    const double DATAsigmaInstX = 4.163;
    const double MCmeanInstY = 421.7;
    const double MCsigmaInstY = 3.983;
    const double DATAmeanInstY = 422.3;
    const double DATAsigmaInstY = 3.867;


    
};

double AnaCut::pi0KineticE = -999;
double AnaCut::pi0costheta = -999;

#endif
