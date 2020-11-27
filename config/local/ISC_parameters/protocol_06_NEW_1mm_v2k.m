% 6)NEW  1mm protocol - seq version v2k:
% PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
global hmri_def
hmri_def.MPMacq_set.name = 'v2k protocol';
hmri_def.MPMacq_set.tag  = 'v2k';
hmri_def.MPMacq_set.val  = [24.5 24.5 6 21];
hmri_def.imperfectSpoilCorr.v2k.tag = 'v2k protocol';
hmri_def.imperfectSpoilCorr.v2k.P2_a = [71.2817617982844,-92.2992876164017,45.8278193851731];
hmri_def.imperfectSpoilCorr.v2k.P2_b = [-0.137859046784839,0.122423212397157,0.957642744668469];

