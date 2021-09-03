% Calculate A (PD) map from PDw and T1w data.
% 
% Inputs of form:
%   {PDw,T1w}.data (array)
%   {PDw,T1w}.fa   (rad)
%   {PDw,T1w}.TR   (s or ms)
%   {PDw,T1w}.B1   (actual fa / nominal fa)
%
% Output:
%   A (in arbitrary units)

function A=hmri_calc_A(PDw,T1w,small_angle_approx)

if isempty(PDw.B1), PDw.B1=1; end
if isempty(T1w.B1), T1w.B1=1; end

if ~small_angle_approx
    PDw.t=2*tan(PDw.B1.*PDw.fa/2);
    T1w.t=2*tan(T1w.B1.*T1w.fa/2);
else
    PDw.t=PDw.B1.*PDw.fa;
    T1w.t=T1w.B1.*T1w.fa;
end

A=T1w.data.*PDw.data.*( T1w.TR.*PDw.t./T1w.t - PDw.TR.*T1w.t./PDw.t )...
    ./ ( PDw.data.*T1w.TR.*PDw.t - T1w.data.*PDw.TR.*T1w.t );
    
% Zero data points where we cannot estimate result
zeromask=(T1w.data==0)|(PDw.data==0)|(T1w.t==0)|(PDw.t==0);
A(zeromask)=0;

end