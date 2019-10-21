% ========================================================================
% function hMRI_read_vols(V,VG,res,p)
%
% read image volume
% 
% Input:
%   VG     - structure containing image volume information of ith image
%   V      - structure containing image volume information of target image
%   res    - resampling function / order of sinc (-7 recommended)
%   p      - z position
%   dim    - defines the axis along which to slice
%
% ========================================================================
% S.Mohammadi 18/10/2019
function Atmp = hMRI_read_vols(V,VG,res,p,dim)
    if ~exist('dim','var')
        dim = 3;
    end
    if ~exist('p','var')
        p = [];
    end
    
    dm      = VG.dim;
    if ~isempty(p)
        switch dim
            case 1
                M = VG.mat(:,[2 3 1 4])*spm_matrix([0 0 p]);
                Atmp = spm_slice_vol(V,V.mat\M,dm(2:3),res);
            case 2
                M = VG.mat(:,[1 3 2 4])*spm_matrix([0 0 p]);
                Atmp = spm_slice_vol(V,V.mat\M,dm([1 3]),res);
            case 3
                M = VG.mat*spm_matrix([0 0 p]);
                Atmp = spm_slice_vol(V,V.mat\M,dm(1:2),res);
        end
    else
        Atmp    = zeros(dm);
        for p=1:dm(3)
            M = VG.mat*spm_matrix([0 0 p]);
            Atmp(:,:,p) = spm_slice_vol(V,V.mat\M,dm(1:2),res);
        end
        switch dim
            case 1
                Atmp = permute(Atmp,[2 3 1]);
            case 2
                Atmp = permute(Atmp,[1 3 2]);
        end
    end
end


