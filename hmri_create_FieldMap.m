function varargout = hmri_create_FieldMap(varargin)
% This is a fraction of the FieldMap script (case 'createfieldmap')
% rewritten for the hMRI toolbox in order to make use of the original
% SPM12's FieldMap script wherever no modification is required. The
% modification implies the use of new segmentation tools to create the
% brain mask for unwrapping.
%
% Implementation following FieldMap version:
% $Id: FieldMap.m 6994 2017-01-26 16:19:14Z guillaume $
%
% FORMAT IP.fm = hmri_create_FieldMap('CreateFieldMap',IP);
%=======================================================================

persistent IP            % Input and results

if nargin == 0
    warning('Can only be called for map creation.\nFORMAT: FM = hmri_create_FieldMap(''CreateFieldMap'',IP);');
    return;
else
    Action = varargin{1};
end

switch lower(Action)
    %=======================================================================
    %
    % Create unwrapped phase-map - NO gui
    %
    %=======================================================================

    case 'createfieldmap'

        IP=varargin{2};

        % First check that images are in same space.

        if size([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],2)==4
            ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim' IP.P{3}.dim' IP.P{4}.dim']');
            %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]');
            ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]'));
        else
            ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim']');
            %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:)]');
            ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:)]'));
        end

        if any(any(diff(ip_dim,1,1),1)&[1,1,1])
            errordlg({'Images don''t all have same dimensions'});
            drawnow;
            varargout{1}=[];
        elseif any(any(abs(diff(ip_mat,1,1))>1e-4))
            errordlg({'Images don''t all have same orientation & voxel size'});
            drawnow;
            varargout{1}=[];
        else
            % Update flags for unwarping (in case TEs have been adjusted
            IP.uflags.etd = IP.et{2}-IP.et{1};

            % Clear any brain mask
            IP.uflags.bmask = [];

            % SPM5 Update
            % If flag selected to mask brain and the field map is not based on EPI
            if IP.maskbrain==1
                IP.fmagP = FieldMap('Magnitude',IP);
                % IP.uflags.bmask = pm_brain_mask(IP.fmagP,IP.mflags);
                IP.uflags.bmask = hmri_create_pm_brain_mask(IP.fmagP,IP.mflags);
                varargout{2} = IP.fmagP;
            end

            IP.fm = pm_make_fieldmap([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],IP.uflags);
            varargout{1} = IP.fm;
        end

        %=======================================================================
        %
        % Create unwarped epi - NO gui
        %
        %=======================================================================

    case 'unwarpepi'

        %
        % Update unwarped EPI
        %
        IP=varargin{2};
        IP.uepiP = struct('fname',   'Image in memory',...
            'dim',     IP.epiP.dim,...
            'dt',[64 spm_platform('bigend')],...
            'pinfo',   IP.epiP.pinfo(1:2),...
            'mat',     IP.epiP.mat);

        % Need to sample EPI and voxel shift map in space of EPI...
        [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
        xyz = [x(:) y(:) z(:)];

        % Space of EPI is IP.epiP{1}.mat and space of
        % voxel shift map is IP.vdmP{1}.mat
        tM = inv(IP.epiP.mat\IP.vdmP.mat);

        x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
        y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
        z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
        xyz2 = [x2(:) y2(:) z2(:)];

        %
        % Make mask since it is only meaningful to calculate undistorted
        % image in areas where we have information about distortions.
        %
        msk = reshape(double(xyz2(:,1)>=1 & xyz2(:,1)<=IP.vdmP.dim(1) &...
            xyz2(:,2)>=1 & xyz2(:,2)<=IP.vdmP.dim(2) &...
            xyz2(:,3)>=1 & xyz2(:,3)<=IP.vdmP.dim(3)),IP.epiP.dim(1:3));

        % Read in voxel displacement map in correct space
        tvdm = reshape(spm_sample_vol(spm_vol(IP.vdmP.fname),xyz2(:,1),...
            xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

        % Voxel shift map must be added to the y-coordinates.
        uepi = reshape(spm_sample_vol(IP.epiP,xyz(:,1),...
            xyz(:,2)+tvdm(:),xyz(:,3),1),IP.epiP.dim(1:3));% TEMP CHANGE

        % Sample Jacobian in correct space and apply if required
        if IP.ajm==1
            if IP.epifm==1 % If EPI, use inverted jacobian

                IP.jim = reshape(spm_sample_vol(IP.vdm.ijac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            else
                IP.jim = reshape(spm_sample_vol(IP.vdm.jac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            end
            uepi = uepi.*(1+IP.jim);
        end

        IP.uepiP.dat=uepi.*msk;
        varargout{1}=IP.uepiP;

        %=======================================================================
        %
        % Create unwarped epi - NO gui ***FOR XY***
        %
        %=======================================================================

    case 'unwarpepixy'
        %
        % Update unwarped EPI
        %
        IP=varargin{2};
        IP.uepiP = struct('fname',   'Image in memory',...
            'dim',     IP.epiP.dim,...
            'dt',[64 spm_platform('bigend')],...
            'pinfo',   IP.epiP.pinfo(1:2),...
            'mat',     IP.epiP.mat);

        % Need to sample EPI and voxel shift map in space of EPI...
        [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
        xyz = [x(:) y(:) z(:)];

        % Space of EPI is IP.epiP{1}.mat and space of
        % voxel shift map is IP.vdmP{1}.mat
        tM = inv(IP.epiP.mat\IP.vdmP.mat);

        x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
        y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
        z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
        xyz2 = [x2(:) y2(:) z2(:)];

        %
        % Make mask since it is only meaningful to calculate undistorted
        % image in areas where we have information about distortions.
        %
        msk = reshape(double(xyz2(:,1)>=1 & xyz2(:,1)<=IP.vdmP.dim(1) &...
            xyz2(:,2)>=1 & xyz2(:,2)<=IP.vdmP.dim(2) &...
            xyz2(:,3)>=1 & xyz2(:,3)<=IP.vdmP.dim(3)),IP.epiP.dim(1:3));

        % Read in voxel displacement map in correct space
        tvdm = reshape(spm_sample_vol(spm_vol(IP.vdmP.fname),xyz2(:,1),...
            xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

        % Voxel shift map must be added to the x-coordinates.
        uepi = reshape(spm_sample_vol(IP.epiP,xyz(:,1)+tvdm(:),...
            xyz(:,2),xyz(:,3),1),IP.epiP.dim(1:3));% TEMP CHANGE

        % Sample Jacobian in correct space and apply if required
        if IP.ajm==1
            if IP.epifm==1 % If EPI, use inverted jacobian

                IP.jim = reshape(spm_sample_vol(IP.vdm.ijac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            else
                IP.jim = reshape(spm_sample_vol(IP.vdm.jac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            end
            uepi = uepi.*(1+IP.jim);
        end

        IP.uepiP.dat=uepi.*msk;
        varargout{1}=IP.uepiP;

        %=======================================================================
        %
        % Create unwarped epi - YZ not implemented in original FieldMap toolbox
        %
        %=======================================================================

    case 'unwarpepiyz'

        %
        % Update unwarped EPI
        %
        IP=varargin{2};
        IP.uepiP = struct('fname',   'Image in memory',...
            'dim',     IP.epiP.dim,...
            'dt',[64 spm_platform('bigend')],...
            'pinfo',   IP.epiP.pinfo(1:2),...
            'mat',     IP.epiP.mat);

        % Need to sample EPI and voxel shift map in space of EPI...
        [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
        xyz = [x(:) y(:) z(:)];

        % Space of EPI is IP.epiP{1}.mat and space of
        % voxel shift map is IP.vdmP{1}.mat
        tM = inv(IP.epiP.mat\IP.vdmP.mat);

        x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
        y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
        z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
        xyz2 = [x2(:) y2(:) z2(:)];

        %
        % Make mask since it is only meaningful to calculate undistorted
        % image in areas where we have information about distortions.
        %
        msk = reshape(double(xyz2(:,1)>=1 & xyz2(:,1)<=IP.vdmP.dim(1) &...
            xyz2(:,2)>=1 & xyz2(:,2)<=IP.vdmP.dim(2) &...
            xyz2(:,3)>=1 & xyz2(:,3)<=IP.vdmP.dim(3)),IP.epiP.dim(1:3));

        % Read in voxel displacement map in correct space
        tvdm = reshape(spm_sample_vol(spm_vol(IP.vdmP.fname),xyz2(:,1),...
            xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

        % Voxel shift map must be added to the y-coordinates.
        uepi = reshape(spm_sample_vol(IP.epiP,xyz(:,1),...
            xyz(:,2),xyz(:,3)+tvdm(:),1),IP.epiP.dim(1:3));

        % Sample Jacobian in correct space and apply if required
        if IP.ajm==1
            if IP.epifm==1 % If EPI, use inverted jacobian

                IP.jim = reshape(spm_sample_vol(IP.vdm.ijac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            else
                IP.jim = reshape(spm_sample_vol(IP.vdm.jac,xyz2(:,1),...
                    xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
            end
            uepi = uepi.*(1+IP.jim);
        end

        IP.uepiP.dat=uepi.*msk;
        varargout{1}=IP.uepiP;


        %=======================================================================
        %
        % Coregister fieldmap magnitude image to EPI to unwarp
        %
        %=======================================================================

    case 'matchvdm'

        IP=varargin{2};

        %
        % Need a fieldmap magnitude image
        %

        if isempty(IP.pP) && ~isempty(IP.P{1})

            IP.fmagP = struct(...
                'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
                'dim',   IP.P{1}.dim,...
                'dt',    IP.P{1}.dt,...
                'pinfo', IP.P{1}.pinfo,...
                'mat',   IP.P{1}.mat);

            % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
            % If using phase and magnitude, use magnitude image.
            if strcmp(IP.uflags.iformat,'RI')
                IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
            else
                IP.fmagP = IP.P{2};
            end
        elseif ~isempty(IP.pP) && ~isempty(IP.fmagP)
            fprintf('Using %s for matching\n',IP.fmagP.fname);
        else
            IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
        end

        % Now we have field map magnitude image, we want to coregister it to the
        % EPI to be unwarped.
        % If using an EPI field map:
        % 1) Coregister magnitude image to EPI.
        % 2) Apply resulting transformation matrix to voxel shift map
        % If using a non-EPI field map:
        % 1) Forward warp magnitude image
        % 2) Coregister warped magnitude image to EPI.
        % 3) Apply resulting transformation matrix to voxel shift map

        if IP.epifm==1
            [~,M] = FieldMap('Coregister',IP.epiP,IP.fmagP);
            MM = IP.fmagP.mat;
        else
            % Need to sample magnitude image in space of EPI to be unwarped...
            [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
            xyz = [x(:) y(:) z(:)];

            % Space of EPI is IP.epiP{1}.mat and space of fmagP is IP.fmagP.mat
            tM = inv(IP.epiP.mat\IP.fmagP.mat);
            x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
            y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
            z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
            xyz2 = [x2(:) y2(:) z2(:)];
            wfmag = reshape(spm_sample_vol(IP.fmagP,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

            % Need to sample voxel shift map in space of EPI to be unwarped
            tvdm = reshape(spm_sample_vol(IP.vdm.vdm,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),0),IP.epiP.dim(1:3));

            % Now apply warps to resampled forward warped magnitude image...
            wfmag = reshape(spm_sample_vol(wfmag,xyz(:,1),xyz(:,2)-tvdm(:),...
                xyz(:,3),1),IP.epiP.dim(1:3));

            % Write out forward warped magnitude image
            IP.wfmagP = struct('dim',  IP.epiP.dim,...
                'dt',[64 spm_platform('bigend')],...
                'pinfo',   IP.epiP.pinfo,...
                'mat',     IP.epiP.mat);
            IP.wfmagP = FieldMap('Write',IP.epiP,wfmag,'wfmag_',4,'Voxel shift map');

            % Now coregister warped magnitude field map to EPI
            [~,M] = FieldMap('Coregister',IP.epiP,IP.wfmagP);

            % Update the .mat file of the forward warped mag image
            spm_get_space(deblank(IP.wfmagP.fname),M*IP.wfmagP.mat);

            % Get the original space of the fmap magnitude
            MM = IP.fmagP.mat;
        end

        % Update .mat file for voxel displacement map
        IP.vdmP.mat=M*MM;
        spm_get_space(deblank(IP.vdmP.fname),M*MM);

        varargout{1} = IP.vdmP;

        %=======================================================================
        %
        % Coregister fieldmap magnitude image to EPI to do unwarpxy
        %
        %=======================================================================

    case 'matchvdmxy'

        IP=varargin{2};

        %
        % Need a fieldmap magnitude image
        %

        if isempty(IP.pP) && ~isempty(IP.P{1})

            IP.fmagP=struct(...
                'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
                'dim',   IP.P{1}.dim,...
                'dt',    IP.P{1}.dt,...
                'pinfo', IP.P{1}.pinfo,...
                'mat',   IP.P{1}.mat);

            % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
            % If using phase and magnitude, use magnitude image.
            if strcmp(IP.uflags.iformat,'RI')
                IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
            else
                IP.fmagP = IP.P{2};
            end
        elseif ~isempty(IP.pP) && ~isempty(IP.fmagP)
            fprintf('Using %s for matching\n',IP.fmagP.fname);
        else
            %IP.fmagP = spm_vol(spm_get(1,'*.img','Select field map magnitude image'));
            % SPM5 Update
            IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
        end

        % Now we have field map magnitude image, we want to coregister it to the
        % EPI to be unwarped.
        % If using an EPI field map:
        % 1) Coregister magnitude image to EPI.
        % 2) Apply resulting transformation matrix to voxel shift map
        % If using a non-EPI field map:
        % 1) Forward warp magnitude image
        % 2) Coregister warped magnitude image to EPI.
        % 3) Apply resulting transformation matrix to voxel shift map

        if IP.epifm==1
            [~,M] = FieldMap('Coregister',IP.epiP,IP.fmagP);
            MM = IP.fmagP.mat;
        else
            % Need to sample magnitude image in space of EPI to be unwarped...
            [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
            xyz = [x(:) y(:) z(:)];

            % Space of EPI is IP.epiP{1}.mat and space of fmagP is IP.fmagP.mat
            tM = inv(IP.epiP.mat\IP.fmagP.mat);
            x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
            y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
            z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
            xyz2 = [x2(:) y2(:) z2(:)];
            wfmag = reshape(spm_sample_vol(IP.fmagP,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

            % Need to sample voxel shift map in space of EPI to be unwarped
            tvdm = reshape(spm_sample_vol(IP.vdm.vdm,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),0),IP.epiP.dim(1:3));

            % Now apply warps to resampled forward warped magnitude image...
            wfmag = reshape(spm_sample_vol(wfmag,xyz(:,1)-tvdm(:),xyz(:,2),...
                xyz(:,3),1),IP.epiP.dim(1:3));

            % Write out forward warped magnitude image
            IP.wfmagP = struct('dim',  IP.epiP.dim,...
                'dt',[64 spm_platform('bigend')],...
                'pinfo',   IP.epiP.pinfo,...
                'mat',     IP.epiP.mat);
            IP.wfmagP = FieldMap('Write',IP.epiP,wfmag,'wfmag_',4,'Voxel shift map');

            % Now coregister warped magnitude field map to EPI
            [~,M] = FieldMap('Coregister',IP.epiP,IP.wfmagP);

            % Update the .mat file of the forward warped mag image
            spm_get_space(deblank(IP.wfmagP.fname),M*IP.wfmagP.mat);

            % Get the original space of the fmap magnitude
            MM = IP.fmagP.mat;
        end

        % Update .mat file for voxel displacement map
        IP.vdmP.mat=M*MM;
        spm_get_space(deblank(IP.vdmP.fname),M*MM);

        varargout{1} = IP.vdmP;

        %=======================================================================
        %
        % Coregister fieldmap magnitude image to EPI to do unwarpyz
        %
        %=======================================================================

    case 'matchvdmyz'

        IP=varargin{2};

        %
        % Need a fieldmap magnitude image
        %

        if isempty(IP.pP) && ~isempty(IP.P{1})

            IP.fmagP=struct(...
                'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
                'dim',   IP.P{1}.dim,...
                'dt',    IP.P{1}.dt,...
                'pinfo', IP.P{1}.pinfo,...
                'mat',   IP.P{1}.mat);

            % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
            % If using phase and magnitude, use magnitude image.
            if strcmp(IP.uflags.iformat,'RI')
                IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
            else
                IP.fmagP = IP.P{2};
            end
        elseif ~isempty(IP.pP) && ~isempty(IP.fmagP)
            fprintf('Using %s for matching\n',IP.fmagP.fname);
        else
            IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
        end

        % Now we have field map magnitude image, we want to coregister it to the
        % EPI to be unwarped.
        % If using an EPI field map:
        % 1) Coregister magnitude image to EPI.
        % 2) Apply resulting transformation matrix to voxel shift map
        % If using a non-EPI field map:
        % 1) Forward warp magnitude image
        % 2) Coregister warped magnitude image to EPI.
        % 3) Apply resulting transformation matrix to voxel shift map

        if IP.epifm==1
            [~,M] = FieldMap('Coregister',IP.epiP,IP.fmagP);
            MM = IP.fmagP.mat;
        else
            % Need to sample magnitude image in space of EPI to be unwarped...
            [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
            xyz = [x(:) y(:) z(:)];

            % Space of EPI is IP.epiP{1}.mat and space of fmagP is IP.fmagP.mat
            tM = inv(IP.epiP.mat\IP.fmagP.mat);
            x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
            y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
            z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
            xyz2 = [x2(:) y2(:) z2(:)];
            wfmag = reshape(spm_sample_vol(IP.fmagP,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

            % Need to sample voxel shift map in space of EPI to be unwarped
            tvdm = reshape(spm_sample_vol(IP.vdm.vdm,xyz2(:,1),...
                xyz2(:,2),xyz2(:,3),0),IP.epiP.dim(1:3));

            % Now apply warps to resampled forward warped magnitude image...
            wfmag = reshape(spm_sample_vol(wfmag,xyz(:,1),xyz(:,2),...
                xyz(:,3)-tvdm(:),1),IP.epiP.dim(1:3));

            % Write out forward warped magnitude image
            IP.wfmagP = struct('dim',  IP.epiP.dim,...
                'dt',[64 spm_platform('bigend')],...
                'pinfo',   IP.epiP.pinfo,...
                'mat',     IP.epiP.mat);
            IP.wfmagP = FieldMap('Write',IP.epiP,wfmag,'wfmag_',4,'Voxel shift map');

            % Now coregister warped magnitude field map to EPI
            [~,M] = FieldMap('Coregister',IP.epiP,IP.wfmagP);

            % Update the .mat file of the forward warped mag image
            spm_get_space(deblank(IP.wfmagP.fname),M*IP.wfmagP.mat);

            % Get the original space of the fmap magnitude
            MM = IP.fmagP.mat;
        end

        % Update .mat file for voxel displacement map
        IP.vdmP.mat=M*MM;
        spm_get_space(deblank(IP.vdmP.fname),M*MM);

        varargout{1} = IP.vdmP;
end