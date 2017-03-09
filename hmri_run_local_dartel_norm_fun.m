function out = hmri_run_local_dartel_norm_fun(job)
% function out = hmri_dartel_norm_fun_local(job)
%  derived from spm_dartel_norm_fun_local

% TO DO:
%  Build the output 'out' structure with created file names! 

out = struct();

% MFC - Setting up feds structure which has a copy of the (reordered) subject info:
job = dartel_perimage_to_persubject(job);

feds.template = job.template;
feds.vox      = job.vox;
feds.bb       = job.bb;
feds.fwhm     = [0 0 0];
for m = 1:length(job.subjd)
    feds.data.subj(m).flowfield = job.subjd(m).flowfield;
    feds.data.subj(m).images    = job.subjd(m).tc_vols;
end

% MFC - Jacobian modulation correction to preserve total signal intensity:
feds.preserve = 1;

% MFC - Produces mwc* images, i.e. modulated, spatially normalised images.
% This produces w = |Dphi|t(phi), the product of the Jacobian determinants
% of deformation phi and the tissue class image warped by phi, as per
% Draganski 2011, NeuroImage.
spm_dartel_norm_fun(feds);

% MFC - Now take the MPMs and do a regular normalisation but don't apply
% Jacobian modulation. Produces ws* images.
for mm=1:length(job.subjd)
    feds.data.subj(mm).images = job.subjd(mm).mp_vols;
end
feds.preserve = 0;
spm_dartel_norm_fun(feds);

% Specify the out structure with created file names

end

% ========================================================================
%% SUBFUNCTIONS
% ========================================================================

function job = dartel_perimage_to_persubject(job)
% Rearrange data per subject for dartel processing.
for i=1:numel(job.multsdata.vols_field)
    job.subjd(i).tc_vols = {};
    for j=1:numel(job.multsdata.vols_tc)
        job.subjd(i).tc_vols{j} = job.multsdata.vols_tc{j}{i};
    end
    job.subjd(i).flowfield = {};
    job.subjd(i).flowfield{1} = job.multsdata.vols_field{i};
    job.subjd(i).mp_vols = {};
    for j=1:numel(job.multsdata.vols_pm)
        job.subjd(i).mp_vols{j} = job.multsdata.vols_pm{j}{i};
    end
end
end
%_______________________________________________________________________


% % Smoothing part
% for nm=1:length(job.subjd)
%     for i=1:length(job.subjd(nm).mp_vols)
%         chk = check_entry(job.subjd(nm));
%         if ~isempty(chk)
%             error(chk)
%         end
%         c1 = insert_pref(job.subjd(nm).images{1},'mw');  % removed s
%         c2 = insert_pref(job.subjd(nm).images{2},'mw');  % removed s
%         p  = spm_str_manip(job.subjd(nm).mp_vols{1},'h'); %#ok<*NASGU>
%         f  = insert_pref(job.subjd(nm).mp_vols{i},'w');  % removed s
%         c  = spm_imcalc(char(char(c1),char(c2)),insert_pref(f,'bb_'),'(i1+i2)');
%         c  = c.fname;
%         m_c1 = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',1'];
%         m_c2 = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',2'];
%         m_c  = [spm_select('FPList',fullfile(spm('Dir'),'tpm'),'^TPM.nii') ',6'];
%         p1 = spm_imcalc(char(char(c1),char(f),m_c1),insert_pref(f,'p1_'), ...
%             '(i1.*i2).*(i3>0.05)');
%         p1 = p1.fname;
%         p2 = spm_imcalc(char(char(c2),char(f),m_c2),insert_pref(f,'p2_'), ...
%             '(i1.*i2).*(i3>0.05)');
%         p2 = p2.fname;
%         pp = spm_imcalc(char(char(c),char(f),m_c),insert_pref(f,'p_'), ...
%             '(i1.*i2).*((1-i3)>0.05)');
%         pp = pp.fname;
%         m1 = insert_pref(c1,'s'); spm_smooth(c1,m1,job.fwhm);
%         m2 = insert_pref(c2,'s'); spm_smooth(c2,m2,job.fwhm);
%         m  = insert_pref(c, 's');  spm_smooth(c,m,job.fwhm);
%         n1 = insert_pref(p1,'s'); spm_smooth(p1,n1,job.fwhm);
%         n2 = insert_pref(p2,'s'); spm_smooth(p2,n2,job.fwhm);
%         n  = insert_pref(pp,'s'); spm_smooth(pp,n,job.fwhm);
%         q1 = spm_imcalc(char(n1,m1,m1),insert_pref(p1,'fin_dart_'), ...
%             '(i1./i2).*(i3>0.05)');
%         q2 = spm_imcalc(char(n2,m2,m2),insert_pref(p2,'fin_dart_'), ...
%             '(i1./i2).*(i3>0.05)');
%         q = spm_imcalc(char(n,m,m),insert_pref(pp,'fin_dart_bb_'), ...
%             '(i1./i2).*((i3)>0.05)');
%         q = q.fname;
%         delfiles=strrep({c,p1,p2,m1,m2,n1,n2,pp,m,n,q},'.nii,1','.nii');
%         for ii=1:numel(delfiles)
%             delete(delfiles{ii});
%         end
%     end
% end


% function fout = insert_pref(f,p)
% fout = strcat(spm_str_manip(f,'h'),filesep,p,spm_str_manip(f,'t'));
% end
% %_______________________________________________________________________
% 
% function chk = check_entry(job)
% n1 = numel(job.images);
% chk = '';
% n2 = sum(~cellfun('isempty',regexp(spm_str_manip(job.images,'t'),'(^c1|^c2).*.nii')));
% if n1 ~= 2,
%     chk = 'Wrong input - should be c1 and c2';
% end
% if n2 ~= 2,
%     chk = 'Wrong input - should be c1 and c2';
% end
% end