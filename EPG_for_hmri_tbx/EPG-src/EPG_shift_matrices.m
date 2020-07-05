function S = EPG_shift_matrices(Nmax)
% Generates shift matrices up to order Nmax for 3 standard EPG
% Layout of state vector is [F0 F0* Z0 F1 F-1 Z1 F2 F-2 Z2 ... etc]
% i.e. 3 states per n-value
%
%   Shaihan Malik 2017-07-20

N = (Nmax+1) * 3;
S = zeros([N N]);

%%% F(k>=1) 
kidx = 4:3:N; 
sidx = kidx-3;%<-- this is the state that we shift FROM
idx = kidx + N*(sidx-1);% linear indices of S(kidx,sidx)
S(idx)=1;


%%% F(k<1) <--- start at F-1 
kidx = 5:3:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+3;%<- Negative states come from more negative states
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states don't shift
kidx = 3:3:N ;
ix = kidx + N*(kidx-1);
S(ix)=1;


%%% finally F0+ - relates to F-1
S([1 2],5)=1; % also F0* - this isn't used in practice


end