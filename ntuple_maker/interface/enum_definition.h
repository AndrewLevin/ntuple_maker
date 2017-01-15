#ifndef NTUPLEMAKER_ENUM_DEFINITION_H
#define NTUPLEMAKER_ENUM_DEFINITION_H

enum LeptonSelections {
  Lep1TightSelectionV1  = 1UL<<1,
  Lep1TightSelectionV2  = 1UL<<2,
  Lep1TightSelectionV3  = 1UL<<3,
  Lep1TightSelectionV4  = 1UL<<4,
  Lep1TightSelectionV5  = 1UL<<5,
  Lep2TightSelectionV1  = 1UL<<6,
  Lep2TightSelectionV2  = 1UL<<7,
  Lep2TightSelectionV3  = 1UL<<8,
  Lep2TightSelectionV4  = 1UL<<9,
  Lep2TightSelectionV5  = 1UL<<10,
  Lep1LooseSelectionV1  = 1UL<<11,
  Lep1LooseSelectionV2  = 1UL<<12,
  Lep1LooseSelectionV3  = 1UL<<13,
  Lep1LooseSelectionV4  = 1UL<<14,
  Lep1LooseSelectionV5  = 1UL<<15,
  Lep2LooseSelectionV1  = 1UL<<16,
  Lep2LooseSelectionV2  = 1UL<<17,
  Lep2LooseSelectionV3  = 1UL<<18,
  Lep2LooseSelectionV4  = 1UL<<19,
  Lep2LooseSelectionV5  = 1UL<<20,

};

enum Cuts {

  WLLJJVetoV1 = 1UL << 1,
  WLLJJVetoV2 = 1UL << 2,
  WLLJJVetoV3 = 1UL << 3,
  WLLJJVetoV4 = 1UL << 4,
  WLLJJVetoV5 = 1UL << 5,
  WLLJJVetoV6 = 1UL << 6,
  WLLJJVetoV7 = 1UL << 7,
  WLLJJVetoV8 = 1UL << 8,
  WLLJJVetoV9 = 1UL << 9,
  WLLJJVetoV10 = 1UL << 10,
  WLLJJVetoV11 = 1UL << 11,
  PassTriggerV1 = 1UL << 12,
  PassTriggerV2 = 1UL << 13,
  PassTriggerV3 = 1UL << 14,
  PassTriggerV4 = 1UL << 15,
  PassTriggerV5 = 1UL << 16,
  PassTriggerV6 = 1UL << 17,
};

#endif
