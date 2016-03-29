function filename = GetImgFromMosaic( sOutDir, ENCLOSING_path )

Header_textFile = fullfile(sOutDir, 'GetImgFromMosaic_header_dump.txt');

filename = '';
DICOM_filename = getFirstDicomInfolder( ENCLOSING_path );
NEchoes = 1;

[~, ~, ext] = fileparts(DICOM_filename);
if isempty(DICOM_filename) || strcmp(ext, '.nii')
    % Ignore if no DICOM file or if Nifti is already there
    return
end

% Extract Dicom Header
disp(DICOM_filename);
info_first = spm_dicom_headers(DICOM_filename);
info_first = info_first{1};

if ~isfield(info_first, 'CSAImageHeaderInfo')
    % error('CSAImageHeaderInfo tag is missing. Are you sure your scans are in MOSAIC format?');
    return
end

num_imag_in_mosaic = get_numaris4_val(info_first.CSAImageHeaderInfo, 'NumberOfImagesInMosaic');
if isempty(num_imag_in_mosaic)
    % Not a real mosaic, it should be fine to leave it alone
    return
end

if ~isfield(info_first, 'CSASeriesHeaderInfo')
    error('CSASeriesHeaderInfo tag is missing. Are you sure your scans are in MOSAIC format?');
end

phoenix = get_numaris4_val(info_first.CSASeriesHeaderInfo, 'MrPhoenixProtocol');
if isempty(phoenix)
    error('NUMARIS4 value ''MrPhoenixProtocol''  is missing. Are you sure your scans are in MOSAIC format?');
end

if ~exist('dicomread', 'file')
    error('You need Image Processing Toolbox to process MOSAIC scans.')
end

fid = fopen(Header_textFile, 'w');
fwrite(fid, phoenix);
fclose(fid);

% Open Dicom header of the first Dicom image
%% ======================
fid = fopen(Header_textFile);
k=[];
TE=[];

tline = fgetl(fid);
p = 0;
inc = 0;


% Find number of slices
%% ======================
Found = 0;

while ischar(tline) && ~Found
    tline = fgetl(fid);
    
    k = strfind(tline, 'sKSpace.lImagesPerSlab');
    if ~isempty(k)
          NSlices = sscanf(tline, 'sKSpace.lImagesPerSlab = %i' );
          Found = 1;
    end
end

fseek(fid, 0, -1);
tline = fgetl(fid);


% Find number of Echoes
%% ======================
while ischar(tline) && inc < NEchoes
    tline = fgetl(fid);
    
    k = strfind(tline, 'lContrasts');
    if ~isempty(k)
          NEchoes = sscanf(tline, 'lContrasts = %i' );
    end
end

fseek(fid, 0, -1);
tline = fgetl(fid);

% Find Echo time for each Echo
%% ======================
while ischar(tline) && inc < NEchoes
  tline = fgetl(fid);
  if NEchoes > 0
        for h=inc:NEchoes
            name = sprintf('alTE[%i]', h);
            j = strfind(tline, name);
            if ~isempty(j)
                te = sscanf(tline, [ name ' = %i'] );
                TE(h+1) = te;
                inc = inc + 1;
            end
        end
    end
end

% Compute Number of Echoes per Row of the mosaic image
%% ======================
Num = sqrt(double(NEchoes));
EchoPerRow = round(Num+.49);
H = info_first.Rows/EchoPerRow; W = info_first.Columns/EchoPerRow;

D = dir( ENCLOSING_path );
Path = ENCLOSING_path; % fullfile(ENCLOSING_path, 'Echoes');
% mkdir(Path);

it = 1; % iterator to give consecutive numbers to each echo image

% Write Dicom image for each Echo
%% ======================
for j = 1:NEchoes
    mkdir(Path, sprintf('Echo%d', j));
    for i = 3:length(D)
        
        if isdir(    fullfile(ENCLOSING_path,D(i).name) ), continue, end
        if ~isdicom( fullfile(ENCLOSING_path,D(i).name) ), continue, end
        
        % Modify info from Dicom header for each Echo 
        filename = fullfile(ENCLOSING_path,D(i).name);
        I = dicomread(filename);
        info = dicominfo(filename);
        info.InstanceNumber =  info.InstanceNumber + (j-1)*NSlices;                                  %it;
        info.EchoTime = TE(j)*1e-3;
        info.width = W; info.Height = H;
        info.Rows = H; info.Columns = W;
        info.WindowCenter = 0;
        info.WindowWidth = 0;
        
        % Extract one echo image from th emosaic
        X = mod((j-1)*W, EchoPerRow*W)+1; Y = floor((j-1)/EchoPerRow)*H+1;
        I2 = imcrop(I,[X Y W-1 H-1]);
        
        % Write Dicom image corresponding to each echo
        EchoPath = fullfile(Path,sprintf('Echo%d', j));
        dicomwrite(I2, ...
            fullfile(EchoPath, D(i).name), ...
            info ...
            );
        
        it = it + 1;
        
    end
end

for i = 3:length(D)
    filename = fullfile(ENCLOSING_path,D(i).name);
    delete(filename);
end

function filename = getFirstDicomInfolder(path)
        x = dir(path);
        filename='';
        for i=3:numel(x)
            if x(i).isdir
                continue
            end
            
            [~, ~, ext] = fileparts(x(i).name);
            
            if strcmp(ext, '.nii') || strcmp(ext, '.txt')
                continue
            end
            
            filename = fullfile(path, x(i).name);
        end

function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str),
    if strcmp(deblank(str(i).name),name),
        for j=1:str(i).nitems,
            if  str(i).item(j).xx(1),
                val = {val{:} str(i).item(j).val};
            end;
        end;
        break;
    end;
end;
val = strvcat(val{:});
