function E = E_diff(E,diff,kmax,N)
% E = E_diff(E,diff,kmax,N)
%
%    Function to build E operator with diffusion effects for standard EPG
%    3 states per k-value). 
%   
%       E = relaxation matrix (diag(E2 E2 E1))
%       diff = structure with fields:
%              G    - Gradient amplitude(s)
%              tau  - Gradient durations(s)
%              D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
% Shaihan Malik July 2017

[bDL, bDT] = EPG_diffusion_weights(diff.G,diff.tau,diff.D,0:kmax);

% E is a simple diagonal matrix - just need to compute this diagonal
Ed = diag(E);
EdT = Ed(1:2);
EdL = Ed(3);

EdT = EdT*bDT';
EdL = EdL*bDL';

% Combine them
Ed = cat(1,EdT,EdL);
Ed=Ed(:);

%%% Now use sparse diagonal function to define overall matrix
E = spdiags(Ed,0,N,N);

end