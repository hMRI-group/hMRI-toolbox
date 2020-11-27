% 7) 800um protocol - seq version v3* released used by MEG group:
% TR = 25ms for all volumes; flipAngles = [6, 21 deg] for PDw and T1w
% Correction parameters below were determined via Bloch-Torrey 
% simulations but end result agrees well with EPG-derived correction 
% for this RF spoiling increment of 137 degrees.
% See: Callaghan et al. ISMRM, 2015, #1694
global hmri_def
hmri_def.MPMacq_set.names{7} = 'v3star protocol';
hmri_def.MPMacq_set.tags{7}  = 'v3star';
hmri_def.MPMacq_set.vals{7}  = [25 25 6 21];
hmri_def.imperfectSpoilCorr.v3star.tag = 'v3star protocol';
hmri_def.imperfectSpoilCorr.v3star.P2_a = [57.427573706259864,-79.300742898810441,39.218584751863879];
hmri_def.imperfectSpoilCorr.v3star.P2_b = [-0.121114060111119,0.121684347499374,0.955987357483519];
