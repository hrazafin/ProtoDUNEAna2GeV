#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

namespace pi{
  const unsigned int nCuts = 7;

  // can change order of cuts
  const char cutName[nCuts][100] = {
    "Nocut",
    "PandoraSlice",
    "CaloSize",
    "BeamQuality",
    "ProtonCut",
    "MichelScore",
    "APA3",
  };
  
  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality,
    kProtonCut,
    kMichelScore,
    kAPA3,
  };
}

namespace p{
  const unsigned int nCuts = 4;

  const char cutName[nCuts][100] = {"Nocut",
                                    "PandoraSlice",
                                    "CaloSize",
                                    "BeamQuality"};

  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality
  };
}

const double beam_startX_data = -28.3483;
const double beam_startY_data = 424.553;
const double beam_startZ_data = 3.19841;
const double beam_startX_rms_data = 4.63594;
const double beam_startY_rms_data = 5.21649;
const double beam_startZ_rms_data = 1.2887;

const double beam_startX_mc = -30.7834;
const double beam_startY_mc = 422.422;
const double beam_startZ_mc = 0.113008;
const double beam_startX_rms_mc = 4.97391;
const double beam_startY_rms_mc = 4.47824;
const double beam_startZ_rms_mc = 0.214533;

const double beam_angleX_data = 100.464;
const double beam_angleY_data = 103.442;
const double beam_angleZ_data = 17.6633;

const double beam_angleX_mc = 101.579;
const double beam_angleY_mc = 101.212;
const double beam_angleZ_mc = 16.5822;

#endif
