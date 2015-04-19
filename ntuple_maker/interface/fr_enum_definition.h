#ifndef FR_ENUM_DEFINITION_H
#define FR_ENUM_DEFINITION_H

enum Flags {
  GenElectronInTheEvent  = 1UL<<1,
  GenMuonInTheEvent  = 1UL<<2,
  LepTightSelectionV1  = 1UL<<3,
  LepTightSelectionV2  = 1UL<<4,
  LepTightSelectionV3  = 1UL<<5,
  LepTightSelectionV4  = 1UL<<6,
  LepTightSelectionV5  = 1UL<<7,
  LepLooseSelectionV1  = 1UL<<8,
  LepLooseSelectionV2  = 1UL<<9,
  LepLooseSelectionV3  = 1UL<<10,
  LepLooseSelectionV4  = 1UL<<11,
  LepLooseSelectionV5  = 1UL<<12,
};

#endif
