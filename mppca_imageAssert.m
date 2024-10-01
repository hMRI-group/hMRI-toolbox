function [image,mask] = mppca_imageAssert(image,mask)
%{
REFERENCES:
--Veraart et al., NeuroImage (2016) 142, p 394-406 (https://doi.org/10.1016/j.neuroimage.2016.08.016)  
--Does, MD al. Evaluation of principal component analysis image denoising on multi‐exponential MRI relaxometry. Magn Reson Med. 2019; 81: 3503– 3514. https://doi.org/10.1002/mrm.27658
--https://github.com/Neurophysics-CFIN/MP-PCA-Denoising/tree/master
%}
% Want first image indices to discriminate between pixels and last
% dimension to hold data for each pixel. This function puts the image data
% on the form: row x col x slice x rest.

dims = size(image);
assert(length(dims)<=4,'image data array must not have more than 4 dimensions')

% construct mask if not given
if numel(mask)==0
    if sum(dims~=1)==1
        mask = true;
    else
        mask = true([dims(1:end-1) 1]);
    end
end
maskDims = size(mask);

% voxel count must match
assert(numel(mask)==prod(dims(1:end-1)),'mask dimensions does not match image dimensions')

% if image only contains data from single pixel
if (all(dims(1:end-1)==1) && length(dims)<4) || (dims(1)>1 && dims(2)==1 && length(dims)<3) % (1xn || 1x1xn) || (nx1)
    assert(numel(mask)==1,'mask dimensions does not match image dimensions')
    dummy(1,1,1,:) = image;
    image = dummy;
end

% if image is on form rxcxn
if length(dims)==3
    assert(all(maskDims==dims(1:end-1)),'mask dimensions does not match image dimensions')
    dummy(:,:,1,:) = image;
    image = dummy;
end

% if image is on form rxn
if length(dims)==2 && dims(2)~=1
    assert(all(maskDims(1)==dims(1)),'mask dimensions does not match image dimensions')
    dummy(:,1,1,:) = image;
    image = dummy;
end