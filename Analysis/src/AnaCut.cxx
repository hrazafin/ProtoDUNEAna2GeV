#include "../include/AnaCut.h"

bool AnaCut::CutBeamAllInOne(const bool kMC)
{
  // Standard procedure from pion analyses
  
  // 1. Beam ID cut (pion beam/other)
  if(kMC){ // MC
    // In data the beam instrumentation is able to filter for these events but is not inside the MC
    const bool mc_beamID =(AnaIO::true_beam_PDG == 211 || AnaIO::true_beam_PDG == -13);
    AnaIO::hCutBeamIDPass->Fill(mc_beamID);
    if(!mc_beamID){
      return false;
    }
  }
  else{ // Data
    const bool data_beamID = CutBeamID((*AnaIO::beam_inst_PDG_candidates));
    AnaIO::hCutBeamIDPass->Fill(data_beamID);
    if(!data_beamID){
      return false;
    }
  }
  
  //2. Primary beam type cut (track like/shower like)
  const int beam_recoType = AnaIO::reco_beam_type;
  AnaIO::hCutBeamType->Fill(beam_recoType);
  if(beam_recoType!=13){//13: Pandora "track like"
    return false;
  }

  //3. beam position cut
  const bool kBeamPosPass = kMC ? CutMCBeamPos() : CutDataBeamPos();
  AnaIO::hCutBeamPosPass->Fill(kBeamPosPass);
  if(!kBeamPosPass){
    return false;
  }

  //4. APA3 cut
  const double endzcut = 226;
  AnaIO::hCutBeamEndZ->Fill(AnaIO::reco_beam_endZ);
  AnaIO::hCutBeamEndZPass->Fill(!(AnaIO::reco_beam_endZ>=endzcut));
  if(AnaIO::reco_beam_endZ >= endzcut){
    return false;
  }
  
  return true;
}


bool AnaCut::CutBeamID(const std::vector<int> & pidCandidates)
{
  return ( (std::find(pidCandidates.begin(), pidCandidates.end(), 211)) != pidCandidates.end());
}

bool AnaCut::CutMCBeamPos()
{
  return Manual_beamPos_mc(AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                           AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                           AnaIO::true_beam_startDirX, AnaIO::true_beam_startDirY, AnaIO::true_beam_startDirZ,
                           AnaIO::true_beam_startX, AnaIO::true_beam_startY, AnaIO::true_beam_startZ);
}

bool AnaCut::CutDataBeamPos()
{
  return Manual_beamPos_data(AnaIO::event, AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                             AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                             AnaIO::beam_inst_X, AnaIO::beam_inst_Y, AnaIO::beam_inst_dirX, AnaIO::beam_inst_dirY,
                             AnaIO::beam_inst_dirZ, AnaIO::beam_inst_nMomenta, AnaIO::beam_inst_nTracks);
}

bool AnaCut::Manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, const double beam_dirX, const double beam_dirY,   const double beam_dirZ, const double true_dirX,   const double true_dirY, const double true_dirZ,   const double true_startX, const double true_startY, const double true_startZ)
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.h

  //For MC from Owen Goodwins studies
  const double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  const double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  const double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  const double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  const double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;

  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;

  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;

  if ( cos < coslow)
    return false;
    
  return true;
}

bool AnaCut::Manual_beamPos_data(const int event,            const double data_startX,
                         const double data_startY,   const double data_startZ,
                         const double data_dirX,     const double data_dirY,
                         const double data_dirZ,     const double beam_inst_X,
                         const double beam_inst_Y,     const double beam_inst_dirX,
                         const double beam_inst_dirY,  const double beam_inst_dirZ,
                         const int beam_inst_nMomenta, const int beam_inst_nTracks)
{

  //For Data from Owen Goodwin
  const double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
  const double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

  const double deltaX = data_startX - beam_inst_X;
  const double deltaY = data_startY - beam_inst_Y;
  const double cos = beam_inst_dirX*data_dirX + beam_inst_dirY*data_dirY +
               beam_inst_dirZ*data_dirZ;

  if(beam_inst_nMomenta != 1 || beam_inst_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

} 
