function outfiles = hmri_autoreorient(ref, template, other)

% FORMAT outfiles = hmri_autoreorient(ref, template, other)
%
% Function to automatically (but approximately) rigid-body reorient a T1
% image (or any other usual image modality) in the MNI space, i.e. mainly
% set the AC location and correct for head rotation, in order to further
% proceed with the segmentation/normalisation of the image.
% A set of other images can be reoriented along the 1st one. They should be
% specified as "Other Images".
% This is useful as the "unified segmentation" process is rather sensitive
% to the starting orientation of the image.
%
% IN:
% - ref      : filename of the reference image to reorient,
% - template : a template image, already in the MNI space, to which the
%              image is reoriented,
% - other    : filenames of other images to be reoriented along with the
%              reference image.
%
% OUT:
% - outfiles : the list (cell array of strings) of reoriented images 
%              listed in the same order as the input (ref, then other).
%__________________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Code originally written by Carlton Chu, FIL, UCL, London, UK
% Modified and extended by Christophe Phillips & Evelyne Balteau, CRC, ULg,
% Liege, Belgium

% If no input provided, possibility to select images here (otherwise, all
% parameters should be provided)
if nargin<1 || isempty(ref)
    ref = spm_select(inf,'image','Select the reference image to be reoriented');
    other = spm_select(inf,'image','Select the other images be reoriented along with the reference image');
    template = fullfile(spm('dir'),'canonical','avg152T1.nii');
end

if iscell(ref), ref = char(ref); end
if iscell(other), other = char(other); end

Vtempl = spm_vol(template);
flags.regtype = 'mni';
Nother = size(other,1);
method = 1;

switch method
    case 1
        % get image and smooth
        ref = strtrim(ref);
        spm_smooth(ref,'temp.nii',[8 8 8]);
        vf = spm_vol('temp.nii');
        % estimate reorientation
        [M, scal] = spm_affreg(Vtempl,vf,flags);
%         for citer = 2:8
%             [M,scal] = spm_affreg(Vtempl,vf,flags,M,scal);
%         end

        M3 = M(1:3,1:3);
        [u s v] = svd(M3);
        M3 = u*v';
        M(1:3,1:3) = M3;
        % apply it on image to reorient
        N = nifti(ref);
        N.mat = M*N.mat;
        create(N);
        % apply it on other images
        for cother = 1:Nother
            fo = strtrim(other(cother,:));
            if ~isempty(fo) && ~strcmp(ref,fo)
                % allow case where :
                % - noname was passed
                % - name is same as the image used for the reorient
                % => skip
                N = nifti(fo);
                N.mat = M*N.mat;
                create(N);
            end
        end
        % clean up
        spm_get_space(deblank('temp.nii'),M*spm_get_space('temp.nii'));
        %delete('temp.nii');
        
    case 2
        Niter = 8;
        sep = 8./[1 2 4*ones(1,Niter)];
        flags = struct('WG'      ,[]    ,...
            'WF'      ,[]    ,...
            'sep'     ,8     ,...
            'regtype' ,'mni' ,...
            'globnorm',0);
        
        VRef = spm_smoothto8bit(spm_vol(ref),8);
        VTmp = spm_smoothto8bit(spm_vol(template),0);
        VRef.pinfo(1:2,:) = VRef.pinfo(1:2,:)/spm_global(VRef);
        VTmp.pinfo(1:2,:) = VTmp.pinfo(1:2,:)/spm_global(VTmp);
        
        [M,scal] = spm_affreg(VTmp,VRef,flags,eye(4));
%         for citer = 2:Niter
%             flags.sep = sep(citer);
%             [M,scal] = spm_affreg(VTmp,VRef,flags,M,scal);
%         end
        [A,B,C] = svd(M(1:3,1:3)); R = A*C'; %#ok<ASGLU>
        R(:,4) = R*(M(1:3,1:3)\M(1:3,4)); R(4,4) = 1;
        
        spm_progress_bar('Init',Nother,'Auto-Reorient to MNI space','volumes completed');
        % apply reorientation to ref image
        spm_get_space(deblank(ref),R*spm_get_space(deblank(ref)));
        % apply reorientation to all the other images
        for cother = 1:Nother
            spm_get_space(deblank(other(cother,:)),...
                R*spm_get_space(deblank(other(cother,:))));
            spm_progress_bar('Set',cother);
        end
        spm_progress_bar('Clear');
        
end

outfiles = cellstr(char(ref,other));

end