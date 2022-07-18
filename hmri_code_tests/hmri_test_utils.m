classdef hmri_test_utils
    
    % Class with general utility methods for unit testing the hMRI toolbox.
    
    methods(Static)
        
        function S=ernst(alpha,TR,R1)
            % Ernst equation
            
            S=sin(alpha).*(1-exp(-TR.*R1))./(1-cos(alpha).*exp(-TR.*R1));
            
        end
        
        % TODO: warn about activating legacy RNG
        function seedRandomNumberGenerator
            % set randomness method and seed for reproducibility
            %
            % it would be best to use a separate random stream for each test
            % for simplicity we set the global stream
            % this will seed the RNG separately on every worker
            rng(319, 'twister');
        end
        
        function w_TEs=createDecaySignal(w_TE0,TEs,R2s)
            
            dims=size(w_TE0);
            
            % Account for 1D case
            if (length(dims)==2)&&(dims(2)==1), dims=dims(1); end
            
            TEs=reshape(TEs,[ones(1,length(dims)),length(TEs)]);
            
            w_TEs=w_TE0.*exp(-R2s.*TEs);
            
        end
        
        
        function [P,V]=makePhantom(output,transform)
            % Synthetic (in silico) phantom for testing tools requiring nifti volumes
            % in the hMRI toolbox.
            % Roughly approximates the shape of a human head in 3D by replicating the
            % Modified Shepp-Logan Phantom in the "slice" direction.
            %
            % Params:
            % output - filename, must have .nii, e.g. "silico_phantom.nii"
            %          Without a path the phantom is written to the pwd.
            % transform - 4x4 affine matrix to be applied to the standard transform, optional.
            %
            % P - phantom as an array
            % V - SPM volume structure.
            
            % Dimensions - arbitrarily chosen; final size (2*Wxy)^3
            Wxy=30; % Wx=Wy=Wxy
            Wz=25;   % different size in z
            
            P=zeros(2*Wxy,2*Wxy,Wxy);
            xy=Wxy;
            for z=1:Wz+1
                % Modified Shepp-Logan Phantom size in final matrix changes across z
                p=rot90(phantom(2*xy));
                
                % Pad so that each slice has a constant size
                P(:,:,z)=padarray(p,[1,1]*(Wxy-xy),0,'both');
                
                % invert equation of ovoid:
                %   (x/Wx)^2+(y/Wy)^2+(z/Wz)^2 = 1
                % to find the in-plane length, x=y=xy, given the position z, Wz, and
                % Wx=Wy=Wxy
                xy=round(real(Wxy*sqrt(1-(z/Wz)^2)));
            end
            
            % Symmetrise result in z to make phantom "round"
            P=cat(3,P(:,:,end:-1:1),P);
            
            % Scale values to sensible range for scanner reconned data and convert to
            % integer as some tools may break with input values < 1...
            P=int32((5000/max(P(:)))*P);
            
            %% Compute transform
            voxelsize=3; % mm; chosen so that size of volume approximates brain
            
            dim=size(P);
            rotscale=voxelsize*diag([-1,1,1]);
            translation=0.5*voxelsize*[dim(1);-dim(2);-dim(3)];
            
            transform0=[[rotscale,translation];[0,0,0,1]];
            
            if ~exist('transform','var')
                % Identity transform
                transform=transform0;
            else
                transform=transform*transform0;
            end
            
            %% Output file using SPM's nifti facilities
            % Note we need to check endianness
            % https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/17569/versions/1/previews/ByteConversion.html
            byteSeq = typecast(uint32(1), 'uint8'); isBigEndian = (byteSeq(end) == 1);
            V.dt=[spm_type('int32'),isBigEndian];
            
            V.dim=size(P);
            V.pinfo=[1;0;352]; % 352 means there is no extended header
            
            V.mat=transform;
            V.descrip='Synthetic phantom for testing components of the hMRI toolbox';
            
            if exist('output','var')
                V.fname=output;
                spm_write_vol(V,P);
            end
            
        end
        
        function [B1map,V]=makePhantom_B1(dims,B1lim,output)
            % Synthetic (in silico) phantom for testing tools requiring B1
            % maps in the hMRI toolbox. The "Doming" effect associated
            % with scanning brains at high and ultra-high field is 
            % represented by Gaussian decay of B1 away from the centre of 
            % the volume.
            %
            % Params:
            %  dims   - size of the phantom. Either a scalar (isotropic) or
            %           3-element array (anisotropic).
            %  B1lim  - B1 values at the edge (lower value) and centre 
            %           (higher value) of the volume in p.u. The minimum B1
            %           will be reached 10 mm away from the edge.
            %  output - filename, must have .nii, e.g. "silico_phantom.nii"
            %           Without a path the phantom is written to pwd. If
            %           not provided, no file is written out.
            %
            % Outputs:
            %  B1map - phantom B1 map in p.u. as an array
            %  V     - SPM volume structure.
            
            if isscalar(dims)
                dims=dims*[1,1,1];
            end
            
            assert(numel(B1lim)==2,'B1lim must contain two values: a minimum B1 and a maximum B1')
            
            minB1=min(B1lim);
            maxB1=max(B1lim);
            
            % Choose voxel size so that the size of the volume approximates that used
            % for brain imaging
            voxelsize=round(240/max(dims)); % mm
            
            % dim(1) and dim(2) swapped as this is what Matlab expects.
            [X,Y,Z]=meshgrid(linspace(-1,1,dims(2)),linspace(-1,1,dims(1)),linspace(-1,1,dims(3)));
            
            % Set Gaussian parameter k so that minB1 is reached at a sensible point;
            % here B1 = minB1 10 mm away from the volume edge
            voxelsaway=10/voxelsize;
            XYZat10mm=(1-voxelsaway*(2/min(dims)));
            k=log(maxB1/minB1)/(XYZat10mm.^2);
            
            % Gaussian decay normalised to maxB1
            B1map=maxB1*exp(-k*(X.^2+Y.^2+Z.^2));
            
            %% nifti header
            V.descrip='Synthetic B1 phantom for testing components of the hMRI toolbox';
            V.dim=dims;
            
            % Sets centre of coordinate system at centre of volume.
            % Signs of terms chosen to match a real nifti volume.
            rotscale=voxelsize*diag([-1,1,1]);
            translation=0.5*voxelsize*[dims(1);-dims(2);-dims(3)];
            V.mat=[[rotscale,translation];[0,0,0,1]];
            
            % Set magic numbers
            byteSeq = typecast(uint32(1), 'uint8'); isBigEndian = (byteSeq(end) == 1);
            V.dt=[spm_type('float32'),isBigEndian];
            V.pinfo=[1;0;352]; % 352 means there is no extended header
            
            %% Save output in p.u.
            if exist('output','var')&&~isempty(output)
                V.fname=output;
                spm_write_vol(V,B1map);
            end
            
        end
    end
end
