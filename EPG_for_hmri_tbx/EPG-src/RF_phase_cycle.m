function [ phi ] = RF_phase_cycle( N, arg2 )
%   RF_phase_cycle(N, args) Supplies phase cycling pattern to EPG code
%
%   SPGR: phi = RF_phase_cycle(N,Phi0) 
%               N    = number of RF pulses
%               Phi0 = (quadratic) spoiling phase increment DEGREES
%
%   balanced: phi = RF_phase_cycle(N,'balanced')
%             produces [0 pi] cycling
%
%   Shaihan Malik July 2017

if ~ischar(arg2)
    spgr = true;
    Phi0 = arg2*pi/180; %<-- degree to radian
else
    spgr = false;
end

if spgr
    % Quadratic 
    p = 1:N;
    phi = p.*(p-1)/2 * Phi0;
else
    % balanced case
    if mod(N,2)==0
        phi = repmat([0 pi],[1 N/2]);
    else
        phi = [repmat([0 pi],[1 floor(N/2)]) 0];
    end
end

end

