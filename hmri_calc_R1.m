% Calculate R1 map from PDw and T1w data.
% 
% Inputs of form:
%   {PDw,T1w}.data (array)
%   {PDw,T1w}.fa   (rad)
%   {PDw,T1w}.TR   (s)
%   {PDw,T1w}.B1   (actual fa / nominal fa)

function R1=hmri_calc_R1(PDw,T1w,small_angle_approx)

if isempty(PDw.B1), PDw.B1=1; end
if isempty(T1w.B1), T1w.B1=1; end

if ~small_angle_approx
    PDw.t=2*tan(PDw.B1.*PDw.fa/2);
    T1w.t=2*tan(T1w.B1.*T1w.fa/2);
    
    R1=0.5*( PDw.data.*PDw.t/PDw.TR - T1w.data.*T1w.t/T1w.TR )...
        ./( T1w.data./T1w.t - PDw.data./PDw.t + 0.25*(PDw.data.*PDw.t - T1w.data.*T1w.t) );
else    
    PDw.t=PDw.B1.*PDw.fa;
    T1w.t=T1w.B1.*T1w.fa;
    R1=0.5*( PDw.data.*PDw.t/PDw.TR - T1w.data.*T1w.t/T1w.TR )...
        ./( T1w.data./T1w.t - PDw.data./PDw.t );
end

% Zero data points where we cannot estimate result
zeromask=(T1w.data==0)|(PDw.data==0)|(T1w.t==0)|(PDw.t==0);
R1(zeromask)=0;

end