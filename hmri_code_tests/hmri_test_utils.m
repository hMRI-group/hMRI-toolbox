classdef hmri_test_utils
    
    % Class with general utility methods for unit testing the hMRI toolbox.
    
    methods(Static)
        
        function w_TEs=createDecaySignal(w_TE0,TEs,R2s)
            
            dims=size(w_TE0);
            
            % Account for 1D case
            if (length(dims)==2)&&(dims(2)==1), dims=dims(1); end
            
            TEs=reshape(TEs,[ones(1,length(dims)),length(TEs)]);
            
            w_TEs=w_TE0.*exp(-R2s.*TEs);
            
        end
        
        % Synthetic (in silico) phantom for testing tools requiring nifti volumes
        % in the hMRI toolbox.
        % Roughly approximates the shape of a human head in 3D by replicating the
        % Modified Shepp-Logan Phantom in the "slice" direction.
        %
        % Params:
        % output - filename, must have .nii, e.g. "silico_phantom.nii"
        %          Without a path the phantom is written to the pwd.
        % transform - 4x4 affine matrix, optional.
        %
        % P - phantom as an array
        % V - SPM structure.
        %
        % ledwards@cbs.mpg.de
        
        function [P,V]=makePhantom(output,transform)
            
            if ~exist('transform','var')
                % Identity transform
                transform=eye(4);
            end
            
            % Dimensions - arbitrarily chosen; final size (2*Wxy)^3
            Wxy=130; % Wx=Wy=Wxy
            Wz=90;   % different size in z
            
            P=zeros(2*Wxy,2*Wxy,Wxy);
            xy=Wxy;
            for z=1:Wz+1
                % Modified Shepp-Logan Phantom size in final matrix changes across z
                p=phantom(2*xy);
                
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
            
            % Add text to break rotational symmetry. Mean converts RGB to white.
            % Text was rotated and shifted so that it looks nice in
            % spm_check_registration.
            hMRI=mean(insertText(zeros(2*Wxy),[50,Wxy-35],"hmri.info",'BoxColor','black','TextColor','white','FontSize',32),3);
            P(:,:,10:35)=P(:,:,10:35)+fliplr(hMRI');
            hMRI=mean(insertText(zeros(2*Wxy),[20,Wxy-35],"hMRI toolbox",'BoxColor','black','TextColor','white','FontSize',32),3);
            P(:,:,end+1-(10:35))=P(:,:,end+1-(10:35))+fliplr(hMRI');
            
            % Scale values to sensible range for scanner reconned data and convert to
            % integer as some tools may break with input values < 1...
            P=int32(5000*P);
            
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
    end
end