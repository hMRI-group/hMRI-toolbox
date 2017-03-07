% Fixing the MPM-based tissue probability maps
%=============================================
%
% Following some discussions with John Ashburner about TPMs definition for
% the US it follows that (copy-paste from John's email):
% 1) Make sure the TPMs have values between (just above) zero and one, and
%    that they sum to one.
% 2) Registration is driven by the gradients of the logarithm of the tissue
%    probabilities, so check that all those values close to zero are 
%    reasonably smooth when you look at their logs.
% 
% In practice all tissue classes have a minimal value of 10^-6 anywhere in
% the TPMs
% 
% PROBLEM:
% When checking the TPMs named nwTPM_sl2.nii, I noticed that in some places
% the tissue probability was exactly equal to 0.
% The places where the TPMs are =0 are 
% - outside the head volume for the GM/WM/CSF tissue classes,
% - inside the brain for the skull tissue class
% - inside the brain for the 'other' tissue class
% - moreover images are stored in int16
%
% FIX:
% Set the voxels of GM/WM/CSF/skull/other that have a TPM<10^-6 to 10^-6 
% and adjust the 'scalp' tissue class to keep the sum of TPMs = 1.It's
% working because we're talking about really tiny values ~[10^-5 10^-6]
% And save everything as float.
%
% One could also explicitly smooth some log-TPMs to make it smoother in
% areas of low values, where there is a bit of 'noise'.
%
%_______________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% For sanity check, visualize all TPMs together and display intensities
spm_check_registration(...
    fullfile(spm('dir'),'tpm','TPM.nii,1'),...
    fullfile(spm('dir'),'tpm','TPM.nii,2'),...
    fullfile(spm('dir'),'tpm','TPM.nii,3'),...
    fullfile(spm('dir'),'tpm','TPM.nii,4'),...
    fullfile(spm('dir'),'tpm','TPM.nii,5'),...
    fullfile(spm('dir'),'tpm','TPM.nii,6'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,1'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,2'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,3'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,4'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,5'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,6')...
    );
% Then display the intensities and click in a 'corner'.

%% Fixing the TPms
thr = 1e-6;
Ptpm = fullfile(spm('dir'),'tpm','nwTPM_sl2.nii');
Vtpm = spm_vol(Ptpm);
tpm_orig = spm_read_vols(Vtpm);
tpm_updt = tpm_orig;

% Deal with skull & Other -> set to thr (inside head/brain)
for ii = [4 6]
    tpm_ii = squeeze(tpm_orig(:,:,:,ii));
    tpm_ii(tpm_ii<thr) = thr;
    tpm_updt(:,:,:,ii) = tpm_ii;
end
% % Smooth a wee bit (logs of) 'other' inside the brain
% tpm_ot = tpm_updt(:,:,:,6) ;
% stpm_ot = 10.^(smooth3(log10(tpm_ot),'gaussian'));
% tpm_ot(tpm_ot<thr*100) = stpm_ot(tpm_ot<thr*100);
% tpm_updt(:,:,:,6) = tpm_ot;

% Deal with GM, WM, CSF -> set to thr (outside head)
for ii = 1:3
    tpm_ii = squeeze(tpm_orig(:,:,:,ii));
    tpm_ii(tpm_ii<thr) = thr;
    tpm_updt(:,:,:,ii) = tpm_ii;
end
% -> adjust 'scalp' outside the brain
tpm_updt(:,:,:,5) = 1 - sum(tpm_updt(:,:,:,[1:4 6]),4);

% Check intensities
SZ = size(tpm_updt);
vtpm_updt = reshape(tpm_updt,[prod(SZ(1:3)) SZ(4)])';
fplot(sum(vtpm_updt))

% Save updated TPMs
Vtpm_u = spm_vol( fullfile(spm('dir'),'tpm','TPM.nii'));
fn_TPM_upd = fullfile(spm('dir'),'tpm','unwTPM_sl2.nii');
for ii=1:6
    Vtpm_u(ii).fname = fn_TPM_upd;
    Vtpm_u(ii) = spm_create_vol(Vtpm_u(ii));
    Vtpm_u(ii) = spm_write_vol(Vtpm_u(ii),tpm_updt(:,:,:,ii));
end


%% Visual check
spm_check_registration(...
    fullfile(spm('dir'),'tpm','TPM.nii,1'),...
    fullfile(spm('dir'),'tpm','TPM.nii,2'),...
    fullfile(spm('dir'),'tpm','TPM.nii,3'),...
    fullfile(spm('dir'),'tpm','TPM.nii,4'),...
    fullfile(spm('dir'),'tpm','TPM.nii,5'),...
    fullfile(spm('dir'),'tpm','TPM.nii,6'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,1'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,2'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,3'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,4'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,5'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,6'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,1'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,2'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,3'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,4'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,5'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,6')...
    );

