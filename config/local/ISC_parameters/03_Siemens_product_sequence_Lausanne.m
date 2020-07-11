% 3) Siemens product sequence protocol used in Lausanne (G Krueger)
% PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
global hmri_def
hmri_def.MPMacq_set.name = 'Siemens product Lausanne (GK) protocol';
hmri_def.MPMacq_set.tag  = 'SiemPrLausGK';
hmri_def.MPMacq_set.val  = [24.0 19.0 6 20];
hmri_def.imperfectSpoilCorr.SiemPrLausGK.tag = 'Siemens product Lausanne (GK) protocol';
hmri_def.imperfectSpoilCorr.SiemPrLausGK.P2_a = [67.023102027100880,-86.834117103841540,43.815818592349870];
hmri_def.imperfectSpoilCorr.SiemPrLausGK.P2_b = [-0.130876849571103,0.117721807209409,0.959180058389875];
