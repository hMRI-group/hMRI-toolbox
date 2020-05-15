function [F0,Fn,Zn,F] = EPG_GRE(theta,phi,TR,T1,T2,varargin)
%   [F0,Fn,Zn,F] = EPG_GRE(theta,phi,TR,T1,T2,varargin)
%
%   Single pool EPG (classic version) for gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition time, ms
%               T1:         T1, ms
%               T2:         T2, ms
%
%   optional arguments (use string then value as next argument)
%
%               kmax:       maximum EPG order to include. Can be used to
%                           accelerate calculation. 
%                           Setting kmax=inf ensures ALL pathways are
%                           computed
%              diff:        structure with fields:
%                           G    - Gradient amplitude(s)
%                           tau  - Gradient durations(s)
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
%               prep:       can be used to simulate prep pulse prior to
%                           gradient echo. Assume all transverse
%                           magnetization after prep is spoiled.
%                           structure with fields:
%                           flip    -   flip angle, rad
%                           t_delay -   time delay, ms
%
%   Outputs:                
%               F0:         signal (F0 state) directly after each
%                           excitation
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0 F0* Z0 F1 F-1* Z1 F2 F-2* Z2 ...] etc
%
%
%   Shaihan Malik 2017-07-20


%% Extra variables

for ii=1:length(varargin)
    
    % Kmax = this is the maximum EPG 'order' to consider
    % If this is infinity then don't do any pruning
    if strcmpi(varargin{ii},'kmax')
        kmax = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    
    % Prep pulse - struct contains flip (rad), t_delay
    if strcmpi(varargin{ii},'prep')
        prep = varargin{ii+1};
    end
end

%%% The maximum order varies through the sequence. This can be used to speed up the calculation    
np = length(theta);
% if not defined, assume want max
if ~exist('kmax','var')
    kmax = np - 1;
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = np-1; 
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = (0:kmax) + 1; %<-- +1 because (0:kmax) is correct after each RF pulse, but we must increase order by one to also deal with subsequent shift
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; %<-- don't exceed kmax, we break after last RF pulse
else
    kmax_per_pulse = [1:ceil(np/2) (floor(np/2)):-1:1];
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=3*(kmax+1);

%%% Build Shift matrix, S
S = EPG_shift_matrices(kmax);
S = sparse(S);

%% Set up matrices for Relaxation

E1 = exp(-TR/T1);
E2 = exp(-TR/T2);
E = diag([E2 E2 E1]);

%%% regrowth
b = zeros([N 1]);
b(3) = 1-E1;%<--- just applies to Z0

%%% Add in diffusion at this point 
if exist('diff','var')
    E = E_diff(E,diff,kmax,N);
else
    % If no diffusion, E is the same for all EPG orders
    E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N);
end
    

%%% Composite relax-shift
ES=E*S;
ES=sparse(ES);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 3x3 corner, this helps build_T
i1 = [];
for ii=1:3
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end


%% F matrix (many elements zero, not efficient)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(3)=1;   % M0 - could be variable


%% Prep pulse - execute here
if exist('prep','var')
    %%% Assume the prep pulse leaves NO transverse magnetization, or that
    %%% this is spoiled so that it cannot be refocused. Only consider
    %%% z-terms
    
    % RF rotation just cos(theta) on  Mz term
    FF(3)=cos(prep.flip)*FF(3);
    
    % Now apply time evolution during delay period
    E1p = exp(-prep.t_delay/T1); 
    FF(3) = E1p * FF(3) + (1-E1p);
    
end

%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot(theta(jj),phi(jj));
   
    %%% Variable order of EPG, speed up calculation
    kmax_current = kmax_per_pulse(jj);
    kidx = 1:3*(kmax_current+1); %+1 because states start at zero
    
    %%% Replicate A to make large transition matrix
    build_T(A);
    
    %%% Apply flip and store this: splitting these large matrix
    %%% multiplications into smaller ones might help
    F(kidx,jj)=T(kidx,kidx)*FF(kidx);
    
    if jj==np
        break
    end
    
    %%% Now deal with evolution
    FF(kidx) = ES(kidx,kidx)*F(kidx,jj)+b(kidx);
    
    % Deal with complex conjugate after shift
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate 
end


%%% Return signal
F0=F(1,:);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi(:)) *1i;


%%% Construct Fn and Zn
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; 
kvals = -kmax:kmax;

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
Zn = F(3:3:end,:);



    %%% NORMAL EPG transition matrix as per Weigel et al JMR 2010 276-285 
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
