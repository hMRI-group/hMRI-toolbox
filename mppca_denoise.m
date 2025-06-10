%{
REFERENCES:
--Veraart et al., NeuroImage (2016) 142, p 394-406 (https://doi.org/10.1016/j.neuroimage.2016.08.016)  
--Does, MD al. Evaluation of principal component analysis image denoising on multi‐exponential MRI relaxometry. Magn Reson Med. 2019; 81: 3503– 3514. https://doi.org/10.1002/mrm.27658
--https://github.com/Neurophysics-CFIN/MP-PCA-Denoising/tree/master
%}
function [denoised,S2,P] = mppca_denoise(image,window,mask)
% MP-PCA denoising
%
% input:
% image:  images to be denoised. Must have 3 or 4 indices with MRI images 
%         along the last index and voxels in the first 2 or 3.
% window: sliding window
% mask:   is true for all voxels per default but can be manually set to
%         mask out regions.
%
% output:
% denoised: denoised images
% S2:       map of estimated noise variance
% P:        number of detected signal principal components


%% adjust image dimensions and assert
dimsOld = size(image);
if ~exist('mask','var')
    mask = [];
end
[image,mask] = mppca_imageAssert(image,mask);

dims = size(image);
assert(length(window)>1 && length(window)<4,'window must have 2 or 3 dimensions')
assert(all(window>0),'window values must be strictly positive')
assert(all(window<=dims(1:length(window))),'window values must not exceed image dimensions')
if length(window)==2
    window(3) = 1;
end


%% denoise image
denoised = zeros(size(image));
P = zeros(dims(1:3));
S2 = zeros(dims(1:3));
M = dims(1)-window(1)+1;
N = dims(2)-window(2)+1;
O = dims(3)-window(3)+1;
count = zeros(dims(1:3));
for index = 0:M*N*O-1
    k = 1 + floor(index/M/N);
    j = 1 + floor(mod(index,M*N)/M);
    i = 1 + mod(mod(index,M*N),M);
    rows = i:i-1+window(1);
    cols = j:j-1+window(2);
    slis = k:k-1+window(3);

    % Create X data matrix
    X = reshape(image(rows,cols,slis,:),[],dims(4))';

    % remove masked out voxels
    maskX = reshape(mask(rows,cols,slis),[],1)';
    if nnz(maskX)==0 || nnz(maskX)==1
        continue
    end

    % denoise X
    [X(:,maskX),s2,p] = mppca_denoiseMatrix(X(:,maskX));

    % assign
    X(:,~maskX) = 0;
    denoised(rows,cols,slis,:) = denoised(rows,cols,slis,:) + reshape(X',[window dims(4)]);
    P(rows,cols,slis) = P(rows,cols,slis) + p;
    S2(rows,cols,slis) = S2(rows,cols,slis) + s2;
    count(rows,cols,slis) = count(rows,cols,slis) + 1;
end
skipped = count==0 | ~mask;
denoised = denoised + image.*skipped; % Assign original data to denoisedImage outside of mask and at skipped voxels
count(skipped) = 1;
denoised = denoised./count;
P = P./count;
S2 = S2./count;
P(~mask) = nan;
S2(~mask) = nan;


%% adjust output to match input dimensions
denoised = reshape(denoised,dimsOld);
P = reshape(P,dimsOld(1:end-1));
S2 = reshape(S2,dimsOld(1:end-1));


end


% Want first image indices to discriminate between pixels and last
% dimension to hold data for each pixel. This function puts the image data
% on the form: row x col x slice x rest.
function [image,mask] = mppca_imageAssert(image,mask)

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
end


% X: denoised matrix, s2: original noise variance, p: number of signal components, s2_after: noise variance after denoising
function [X,s2,p,s2_after] = mppca_denoiseMatrix(X) 
M = size(X,1);
N = size(X,2);
if M<N
    [U,lambda] = eig(X*X','vector');
else
    [U,lambda] = eig(X'*X,'vector');
end
[lambda,order] = sort(lambda,'descend');
U = U(:,order);
csum = cumsum(lambda,'reverse');
p = (0:length(lambda)-1)';
p = -1 + find((lambda-lambda(end)).*(M-p).*(N-p) < 4*csum*sqrt(M*N),1);
if p==0
    X = zeros(size(X));
elseif M<N
    X = U(:,1:p)*U(:,1:p)'*X;
else
    X = X*U(:,1:p)*U(:,1:p)';
end
s2 = csum(p+1)/((M-p)*(N-p));
s2_after = s2 - csum(p+1)/(M*N);
end

