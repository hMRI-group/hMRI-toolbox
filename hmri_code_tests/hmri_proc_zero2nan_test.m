% Test file for hmri_proc_zero2nan.m
% Run with
%   results = runtests('hmri_proc_zero2nan_test.m')

%% Main function to generate tests
function tests = hmri_proc_zero2nan_test
tests = functiontests(localfunctions);
end

%% Test functions
% Check all data format available in SPM: 
%   uint8, int16, int32, float32, float64, int8, uint16, uint32

% uint8
function testFunction_uint8(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_uint8.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% int16
function testFunction_int16(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_int16.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% int32
function testFunction_int32(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_int32.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% float32
function testFunction_float32(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_float32.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% float64
function testFunction_float64(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_float64.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% int8
function testFunction_int8(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_int8.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% uint16
function testFunction_uint16(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_uint16.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end

% uint32
function testFunction_uint32(testCase)
% Temporaty folder with the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');

% Pick the 2 files, original & expected results
fn_Orig = fullfile(pth_Dat,'Dat_uint32.nii');
fn_ExpR = spm_file(fn_Orig,'suffix','_ExpRes');

% Apply zero-to-NaN conversion function
hmri_proc_zero2nan(fn_Orig);

% Load image values
actVal = spm_read_vols(spm_vol(fn_Orig));
expVal = spm_read_vols(spm_vol(fn_ExpR));
verifyEqual(testCase,actVal,expVal)
end


%% Those will be run once the test file is loaded/unloaded
function setupOnce(testCase)  %#ok<*INUSD> % do not change function name
% Generate the random data files used for the tests
% Put these files in Matlab's (predefined) temporary folder

% Get SPM's data types for images
Dtypes = spm_type;

% Temporaty folder to save the temporary data files
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');
if ~exist(pth_Dat,'dir'), mkdir(pth_Dat), end % CP: Need this check ?

% Then create some 3D synthetic images, z-axis is of size 3 such that
% - z=1 made of random non-zero values, positive & negative
% - z=2 made of zeros
% - z=3 made of NaNs
Img_sz = [2 4 3]; % 3D image contains 3 slices of size 2x4
Img_val = zeros(Img_sz);
Img_val(:,:,1) = 10.^(randn(Img_sz(1:2))).*sign(randn(Img_sz(1:2)));
Img_val(:,:,3) = NaN;
% Expected results for integer/float images after applying the
% zero-to-nan conversion:
% - integer: NaN's were converted to zero's at image creation & zeros
%            remained untouched at conversion
% - unsigned integer: like integer but negative values turned into zero at
%                     image creation
% - float: zero's were turned into NaN's by conversiotn
% -> integer formatted data
Img_val_expInteger = Img_val;
Img_val_expInteger(:,:,3) = 0;
% -> unsigned integer formatted data
Img_val_expUInteger = Img_val_expInteger;
Img_val_expUInteger(Img_val_expUInteger(:)<0) = 0;
% -> float formatted data
Img_val_expFloat = Img_val;
Img_val_expFloat(:,:,2) = NaN;

% Save the same data in all formats with SPM functions:
% use spm_vol, spm_create_vol, and spm_write_vol
for ii=1:numel(Dtypes)
    % Create original data image
    fn_ii = fullfile(pth_Dat,sprintf('Dat_%s.nii',spm_type(Dtypes(ii))));
    V_ii = struct( ...
        'fname', fn_ii, ...
        'dim',   Img_sz, ...
        'dt',    [Dtypes(ii) 0], ...
        'mat',   eye(4) , ...
        'descrip', sprintf( 'Test %s data',spm_type(Dtypes(ii)) ));
    V_ii = spm_create_vol(V_ii);
    spm_write_vol(V_ii,Img_val); %#ok<*NASGU,*AGROW>
    
    % Create expected resulting data image
    fn_exp_ii = spm_file(fn_ii,'suffix','_ExpRes');
    V_exp_ii = struct( ...
        'fname', fn_exp_ii, ...
        'dim',   Img_sz, ...
        'dt',    [Dtypes(ii) 0], ...
        'mat',   eye(4) , ...
        'descrip', sprintf( 'ExpRes %s data',spm_type(Dtypes(ii)) ));
    V_exp_ii = spm_create_vol(V_exp_ii);
    if spm_type(Dtypes(ii),'nanrep') % with NaN rep -> float
        spm_write_vol(V_exp_ii,Img_val_expFloat);
    else % No NaN rep -> signed/unsigned Integer
        if spm_type(Dtypes(ii),'minval')<0 % -> Integer
            spm_write_vol(V_exp_ii,Img_val_expInteger);
        else % unsigned Integer
            spm_write_vol(V_exp_ii,Img_val_expUInteger);
        end
    end
end

end

function teardownOnce(testCase)  % do not change function name
% Delete the random data files used for the tests, as well as any other
% temporary files

% Temporaty folder where the temporary data files are saved
pth_Dat = fullfile(tempdir,'hMRI_zero2nan_test');
if exist(pth_Dat,'dir')
    fn_2del = spm_select('FPlist',pth_Dat,'.*');
    for ii=1:size(fn_2del,1)
        delete(deblank(fn_2del(ii,:)));
    end
    rmdir(pth_Dat,'s'); 
end

end

%% Those will be run at the beginning/end of each test function
function setup(testCase)  % do not change function name

end
function teardown(testCase)  % do not change function name

end
