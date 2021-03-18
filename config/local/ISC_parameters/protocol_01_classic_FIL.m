% 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
% PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
global hmri_def
hmri_def.MPMacq_set.name = 'Classic FIL protocol';
hmri_def.MPMacq_set.tag  = 'ClassicFIL';
hmri_def.MPMacq_set.val  = [23.7 18.7 6 20];
hmri_def.imperfectSpoilCorr.ClassicFIL.tag = 'Classic FIL protocol';
hmri_def.imperfectSpoilCorr.ClassicFIL.P2_a = [78.9228195006542,-101.113338489192,47.8783287525126];
hmri_def.imperfectSpoilCorr.ClassicFIL.P2_b = [-0.147476233142129,0.126487385091045,0.956824374979504];
