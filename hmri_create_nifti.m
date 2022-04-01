function Ni = hmri_create_nifti(Pout, VG, dt, txtdescrip)
% This function creates a nifti file.
% S. Mohammadi 06/09/2019
%
% In:
% Pout          - file path and name
% VG            - target structure of
% dt            - file type
% txtdescrip    - description of file
%
% Out:
% Ni            - nifti file

Ni          = nifti;
Ni.mat      = VG.mat;
Ni.mat0     = VG.mat;
Ni.descrip  = txtdescrip;
Ni.dat      = file_array(Pout,VG.dim,dt,0,1,0);
create(Ni);

end
