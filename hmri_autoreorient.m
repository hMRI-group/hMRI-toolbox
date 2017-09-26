function out = hmri_autoreorient(ref, template, other)

% FORMAT outfiles = hmri_autoreorient(ref, template, other)
%
% Function to automatically (but approximately) rigid-body reorient a T1
% image (or any other usual image modality) in the MNI space, i.e. mainly
% set the AC location and correct for head rotation, in order to further
% proceed with the segmentation/normalisation of the image.
% A set of other images can be reoriented along the 1st (reference) one.
% They should be specified as "Other Images".
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

flags.regtype = 'mni';
Nother = size(other,1);

% get image and smooth (8mm for reference, no smoothing for template)
ref = strtrim(ref);
template = strtrim(template);
Vref = spm_smoothto8bit(spm_vol(ref),8);
Vtempl = spm_smoothto8bit(spm_vol(template),0);

% estimate reorientation
[M, scal] = spm_affreg(Vtempl,Vref,flags); %#ok<*ASGLU>
% extract rigid body part of the transform 
[U,S,V] = svd(M(1:3,1:3));
R = U*V'; % pure rotation without scaling
T = R*(M(1:3,1:3)\M(1:3,4)); % = R(1:3,1:3)*M(1:3,1:3)^-1*M(1:3,4) translation corrected for scaling
% replace in M
M(1:3,1:3) = R; 
M(1:3,4) = T;

spm_progress_bar('Init',Nother+1,'Auto-Reorient to MNI space','volumes completed');
        
% apply to reference image
spm_get_space(ref,M*spm_get_space(ref));
spm_progress_bar('Set',1);

% apply to other images
for cother = 1:Nother
    fo = strtrim(other(cother,:));
    if ~isempty(fo) && ~strcmp(ref,fo)
        % allow case where :
        % - noname was passed
        % - name is same as the image used for the reorient
        % => skip
        spm_get_space(fo,M*spm_get_space(fo));
    end
    spm_progress_bar('Set',cother+1);
end
spm_progress_bar('Clear');        

out.files = cellstr(char(ref,other));
out.M = M;

end