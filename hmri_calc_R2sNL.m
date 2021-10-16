function [R2s,extrapolated]=hmri_calc_R2sNL(weighted_data,method,opt)
% R2* estimation using raw signals.
%
% ledwards@cbs.mpg.de
%
% Based on the ESTATICS implementation in the hMRI toolbox:
%   https://github.com/hMRI-group/hMRI-toolbox
%
% Input:
%   array of structures (one per contrast) in the form:
%     weighted(contrast).data (NvoxelsX x NvoxelsY x ... x Nechoes)
%     weighted(contrast).TE  (1 x Nechoes)
%   -Voxels must correspond between the weightings (i.e. the images should
%    have been resliced to the same space), but the sampled TEs may be
%    different.
%   -Nechoes must be at least 2 for each weighting.
%
%   method:
%
% Outputs:
%   R2s (NvoxelsX x NvoxelsY x ...): the voxelwise-estimated
%       common R2* of the weightings.
%   extrapolated: cell array containing data extrapolated to TE=0 in the
%       same order as the input (e.g. matching contrast order).
%
% Examples:
%   
%
% References:
%   Weiskopf et al. Front. Neurosci. (2014), "Estimating the apparent
%     transverse relaxation time (R2*) from images with different contrasts
%     (ESTATICS) reduces motion artifacts",
%     https://doi.org/10.3389/fnins.2014.00278

assert(isstruct(weighted_data),'hmri:structError',['inputs must be structs; see help ' mfilename])

dims=size(weighted_data(1).data);
Nvoxels=prod(dims(1:end-1));
Nweighted=numel(weighted_data);

%% Build regression arrays
% Design matrix
D=[];
for w=1:Nweighted
    d=zeros(length(weighted_data(w).TE),Nweighted+1);
    d(:,1)=-weighted_data(w).TE;
    d(:,w+1)=1;
    D=[D;d]; %#ok<AGROW>
end

% Response variable vector
y=[];
for w=1:Nweighted
    
    nTE=length(weighted_data(w).TE);
    assert(nTE>1,'each weighting must have more than one TE')
    
    localDims=size(weighted_data(w).data);
    assert(localDims(end)==nTE,'echoes must be in the final dimension')
    assert(prod(localDims(1:end-1))==Nvoxels,'all input data must have the same number of voxels');
    
    rData=reshape(weighted_data(w).data,Nvoxels,nTE);
    
    y=[y;rData.']; %#ok<AGROW>
end

%% R2* estimation

switch lower(method)
    case 'arlo'
        % Pei, M., et al. (2015), Algorithm for fast monoexponential
        % fitting based on Auto-Regression on Linear Operations (ARLO) of
        % data. Magn. Reson. Med., 73: 843-850.
        % https://doi.org/10.1002/mrm.25137
        % (2019), Erratum to: Algorithm for fast monoexponential fitting
        % based on Auto-Regression on Linear Operations (ARLO) of data
        % (Magn Reson Med 2015;73:843-850). Magn Reson Med, 82: 1576-1576.
        % https://doi.org/10.1002/mrm.27807
        
        nTE=sum(D(:,2:end),1);
        
        %TODO: check all delta TE values are the same!
        TE=D(:,1).*D(:,2:end);
        deltaTE=-(TE(2)-TE(1));
        
        deltay=zeros(sum(nTE-2),Nvoxels);
        inty=zeros(sum(nTE-2),Nvoxels);
        N=0;
        for w=1:Nweighted
            V=y(logical(D(:,w+1)),:);
            nTElocal=sum(D(:,w+1));
            inty((1:nTElocal-2)+N,:)=(deltaTE/3)*(V(1:end-2,:)+4*V(2:end-1,:)+V(3:end,:));
            deltay((1:nTElocal-2)+N,:)=V(1:end-2,:)-V(3:end,:);
            N=N+nTE(w)-2;
        end
        
        T2s=(sum(inty.^2,1)+(deltaTE/3)*sum(inty.*deltay,1))...
            ./((deltaTE/3)*sum(deltay.^2,1)+sum(inty.*deltay,1));
        
        beta=zeros(Nweighted+1,Nvoxels);
        beta(1,:)=1./T2s;
        
        for w=1:Nweighted
            mask=logical(D(:,w+1));
            beta(w+1,:)=median(y(mask,:)./exp(beta(1,:).*D(mask,1)),1,'omitnan');
        end
    case 'darlo'
        % Pei, M., et al. (2015), Algorithm for fast monoexponential
        % fitting based on Auto-Regression on Linear Operations (ARLO) of
        % data. Magn. Reson. Med., 73: 843-850.
        % https://doi.org/10.1002/mrm.25137
        % (2019), Erratum to: Algorithm for fast monoexponential fitting
        % based on Auto-Regression on Linear Operations (ARLO) of data
        % (Magn Reson Med 2015;73:843-850). Magn Reson Med, 82: 1576-1576.
        % https://doi.org/10.1002/mrm.27807
        
        nTE=sum(D(:,2:end),1);
        
        %TODO: check all delta TE values are the same!
        TE=D(:,1).*D(:,2:end);
        deltaTE=-(TE(2)-TE(1));
        
        deltay=zeros(sum(nTE-2),Nvoxels);
        suby=zeros(sum(nTE-2),Nvoxels);
        N=0;
        for w=1:Nweighted
            V=y(logical(D(:,w+1)),:);
            nTElocal=sum(D(:,w+1));
            suby((1:nTElocal-2)+N,:)=V(2:end-1,:);
            deltay((1:nTElocal-2)+N,:)=V(1:end-2,:)-V(3:end,:);
            N=N+nTE(w)-2;
        end
        
        beta=zeros(Nweighted+1,Nvoxels);
        beta(1,:)=(0.5/deltaTE)*sum(suby.*deltay,1)./sum(suby.^2,1);
        
        for w=1:Nweighted
            mask=logical(D(:,w+1));
            beta(w+1,:)=median(y(mask,:)./exp(beta(1,:).*D(mask,1)),1,'omitnan');
        end
        
    otherwise
        % Check for NLLS case, where specification of the log-linear 
        % initialisation method is in the method string following a hyphen
        r=regexp(lower(method),'^nlls-(.*)$','tokens');
        if ~isempty(r)
            initmethod=r{1}{1};
            
            % Initial estimate using log-linear method
            [R2s0,A0]=hmri_calc_R2sLL(weighted_data,initmethod);
            
            if ~exist('opt','var')
                opt=optimset('lsqcurvefit');
            end
            opt=optimset(opt,'Display','off');
            
            beta0=zeros(Nweighted+1,Nvoxels);
            beta0(1,:)=R2s0(:);
            for w=1:Nweighted
                beta0(w+1,:)=A0{w}(:);
            end
            
            % Loop over voxels
            parfor n=1:size(y,2)
                beta(:,n)=lsqcurvefit(@(x,D)expDecay(x,D),beta0(:,n),D,y(:,n),[],[],opt);
            end
        else
            error(['method ' method ' not recognised'])
        end
end

%% Output
% extra unity in reshape argument avoids problems if size(dims)==2.
R2s=reshape(beta(1,:),[dims(1:end-1),1]);

% Extrapolate weightings to TE=0
if nargout>1
    extrapolated=cell(size(weighted_data)); % cell element per contrast
    for w=1:Nweighted
        extrapolatedData=beta(w+1,:);
        extrapolated{w}=reshape(extrapolatedData(:),[dims(1:end-1),1]);
    end
end

end

function S=expDecay(x,D)
S=(D(:,2:end)*x(2:end)).*exp(x(1)*D(:,1));
end