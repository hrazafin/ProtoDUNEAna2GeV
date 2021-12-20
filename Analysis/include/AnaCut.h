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
    bool CutBeamAllInOne(const bool kMC);
    // Fill beam quality cuts histograms (use fitted value for beam quality cut latter)
    void FillBeamQualityHist();
    // Cut on beam PDG code
    bool CutBeamPDG(const bool kMC);
    // Cut on Pandora slice
    bool CutPandoraSlice();
    // Cut on Calo size
    bool CutCaloSize();
    // Cut on beam quality
    bool CutBeamQuality(const bool kMC, bool DoAngleCut = true);
    // Cut on end Z APA3
    bool CutAPA3EndZ();
    // Cut on Michel Score (remove muons)
    bool CutMichelScore(const bool kMC);
    // Cut on median dE/dx (remove protons)
    bool CutMediandEdx(const bool kMC);

    
    // --------------------------------- Legacy Cuts (not used anymore) ----------------------- //
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

    // --------------------------------- Legacy Cuts (not used anymore) ----------------------- //

    // Get truth particle PDG from it's ID
    int GetTruthPDGFromID(const int inID, const vector<int> * idarray, const vector<int> * pdgarray);
    // Get truth-matched particle info using reco particle index
    int GetTruthParticleInfoFromRec(const int recidx);
    bool CutTopology(const bool kMC);
    void CountPFP(const bool kMC, const bool kFill);
    bool IsProton(const int ii, const bool kMC);
    bool IsTrack(const int ii, const bool kMC);  
    bool PassProtonSubPID(const int ii);
    bool IsPiplus(const int ii, const bool kMC);
    bool IsShower(const int ii, const bool kMC);
    bool IsMichel(const int ii, const bool kMC);
    bool IsPiZeroShower(const int ii, const bool kMC); 


  private:

    AnaUtils anaUtils;
    PlotUtils plotUtils;
    int nproton;
    int npiplus;
    int nshower;
    int npi0shower;
    int nmichel;
    int truthParticleType;
    int recParticleType;

    // Cut values
    const double cut_EndZ_APA3 = 220.0; // in cm
    const double cut_MichelScore = 0.55;

    const double MCmeanStartX = -30.85; // in cm
    const double MCmeanStartY = 422.4;  
    const double MCmeanStartZ = 0.1145; 
    const double MCmeanThetaX = 101.6 * TMath::DegToRad(); // in rad.
    const double MCmeanThetaY = 101.2 * TMath::DegToRad();
    const double MCmeanThetaZ = 16.61 * TMath::DegToRad();
    const double MCsigmaStartX = 5.047;
    const double MCsigmaStartY = 4.528;
    const double MCsigmaStartZ = 0.2164;

    const double DATAmeanStartX = -28.44; // in cm
    const double DATAmeanStartY = 424.6;
    const double DATAmeanStartZ = 2.959;
    const double DATAmeanThetaX = 100.6 * TMath::DegToRad(); // in rad.
    const double DATAmeanThetaY = 103.4 * TMath::DegToRad();
    const double DATAmeanThetaZ = 17.65 * TMath::DegToRad();
    const double DATAsigmaStartX = 4.802;
    const double DATAsigmaStartY = 5.355;
    const double DATAsigmaStartZ = 1.332;

    
};


#endif
