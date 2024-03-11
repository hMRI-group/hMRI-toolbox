function Ni = hmri_create_nifti(Pout, VG, dt, txtdescrip)
%hmri_create_nifti Create a nifti file and populate the header
% S. Mohammadi 06/09/2019
%
% In:
%   Pout       - file path and name
%   VG         - target structure with orientation information
%   dt         - file type
%   txtdescrip - description of file
%
% Out:
%   Ni         - nifti file

% Initialise nifti object
Ni = nifti;

% Get orientation information from nifti object in VG if possible
if isfield(VG,'private')
    NG             = VG.private;
    Ni.mat         = NG.mat;
    Ni.mat_intent  = NG.mat_intent;
    Ni.mat0        = NG.mat0;
    Ni.mat0_intent = NG.mat0_intent;
else
    Ni.mat         = VG.mat;
    Ni.mat0        = VG.mat;
end

% Fill other fields from inputs
Ni.descrip = txtdescrip;
Ni.dat     = file_array(Pout,VG.dim,dt,0,1,0);

% Save nifti object to disk
create(Ni);

end