% 2) new FIL/Helms protocol
% PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
global hmri_def
hmri_def.MPMacq_set.name = 'New FIL/Helms protocol';
hmri_def.MPMacq_set.tag  = 'NewFILHelms';
hmri_def.MPMacq_set.val  = [24.5 24.5 5 29];
hmri_def.imperfectSpoilCorr.NewFILHelms.tag = 'New FIL/Helms protocol';
hmri_def.imperfectSpoilCorr.NewFILHelms.P2_a = [93.455034845930480,-120.5752858196904,55.911077913369060];
hmri_def.imperfectSpoilCorr.NewFILHelms.P2_b = [-0.167301931434861,0.113507432776106,0.961765216743606];
