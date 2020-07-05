function [bDL, bDT] = EPG_diffusion_weights(G,tau,D,Nvals)
%   function [bDL bDT] = EPG_diffusion_weights(G,tau,D,Nvals)
%       
%   Computes array of diffusion weighting values exp(-bD)
%   Inputs:
%       - G = gradient amplitude(s) mT/m - can be one or more
%       - tau = duration of gradients ms - same length as G
%       - diffusion coefficient - m^2 / s (i.e. order of 10^-9 expected)
%       - Nvals = array of integers for EPG orders that need to be computed
%
%   Outputs:
%       - bDL/bDT are arrays of exp(-bD) for each value of N specified, or
%       longitudinal/transverse states respectively
%
%   Shaihan Malik 2017-07-19
    
%%% Compute the total dephasing between two EPG states
gmT = 42.58e6 * 1e-3 * 2*pi; % rad s^-1 mT^-1
tau = tau(:)*1e-3; % convert to sec
dur = sum(tau); % total duration of dephasing period to consider
G = G(:)';
dk = gmT*G*tau; %<--- dk is total dephasing between two EPG states

%%% Define helper functions
bLong = @(n)((n*dk).^2*dur);
bTrans = @(n)(bfactors(n));
    
%%% Compute arrays
NN = length(Nvals);
bDL = zeros([NN 1]);
bDT = zeros([NN 1]);

for ii=1:NN
    bDL(ii)=exp(-bLong(Nvals(ii))*D);
    bDT(ii)=exp(-bTrans(Nvals(ii))*D);
end

    %%% Helper function to compute b-value for transverse states as function 
    % of N (order of EPG state).This is a nested function: it sees arguments 
    % from above. See Weigel review paper 2015 (JMRI).
    function bT = bfactors(N)
        
        k0 = N.*dk; %<--- initial k-space is the order of the state * dk
        %%% for each segment we need the initial k-value ki and final value kf
        %%% initialize with ki=k0.
        ki = k0;
        bT=0;
        % The overall b value depends on the k-space value, hence we
        % account for each gradient lobe in order
        for jj=1:length(G)
            kf = ki+gmT*G(jj)*tau(jj);
            bT = bT + (tau(jj)/3)*(ki.^2+kf.^2+ki.*kf);
            
            % now update ki for next gradient lobe
            ki=kf;
        end
    end

end