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
    // Cut on data beam ID only
    bool CutBeamID(const std::vector<int> & pidCandidates);
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
    int GetTruthPDGFromID(const int inID, const vector<int> * idarray, const vector<int> * pdgarray);
    int GetTruthParticleInfoFromRec(const int recidx);
    bool CutTopology(const bool kMC);
    void CountPFP(const bool kMC);
    bool IsProton(const int ii);
    bool IsTrack(const int ii);  
    bool PassProtonSubPID(const int ii);
    bool IsPiplus(const int ii);
    bool IsShower(const int ii);
    bool IsMichel(const int ii);
    bool IsPiZeroShower(const int ii); 
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
    
};


#endif
