% 4) High-res (0.8mm) FIL protocol:
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
global hmri_def
hmri_def.MPMacq_set.name = 'High-res FIL protocol';
hmri_def.MPMacq_set.tag  = 'HResFIL';
hmri_def.MPMacq_set.val  = [23.7 23.7 6 28];
hmri_def.imperfectSpoilCorr.HResFIL.tag = 'High-res FIL protocol';
hmri_def.imperfectSpoilCorr.HResFIL.P2_a = [1.317257319014170e+02,-1.699833074433892e+02,73.372595677371650];
hmri_def.imperfectSpoilCorr.HResFIL.P2_b = [-0.218804328507184,0.178745853134922,0.939514554747592];
