% Define reasonable parameter set
dims=[20,50,30];

R1min=0.5; % / s
R1max=2; % / s
R1=R1min+(R1max-R1min)*rand(prod(dims),1); % / s

PDscale=1000;
PD=PDscale*rand(prod(dims),1);

B1min=0.4;
B1max=1.6;
B1map=B1min+(B1max-B1min)*rand([dims,1]);

reltol=1e-3;

% Ernst equation
ernst=@(alpha,TR,R1) sin(alpha).*(1-exp(-TR.*R1))./(1-cos(alpha).*exp(-TR.*R1));

small_angle_approx=false;

%% Typical 7T protocol
PDw.TR=28.5e-3; % s
PDw.fa=deg2rad(5); % rad
T1w.TR=28.5e-3; % s
T1w.fa=deg2rad(26); % rad

PDw.data=reshape(bsxfun(@times,PD,ernst(PDw.fa*B1map(:),PDw.TR,R1)),dims);
T1w.data=reshape(bsxfun(@times,PD,ernst(T1w.fa*B1map(:),T1w.TR,R1)),dims);

PDw.B1=B1map;
T1w.B1=B1map;

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

assert(all(abs(PDest(:)-PD)./PD<reltol),'Estimated PD has large error!')

%% Typical 3T protocol
PDw.TR=25e-3; % s
PDw.fa=deg2rad(6); % rad
T1w.TR=25e-3; % s
T1w.fa=deg2rad(21); % rad

PDw.data=reshape(bsxfun(@times,PD,ernst(PDw.fa*B1map(:),PDw.TR,R1)),dims);
T1w.data=reshape(bsxfun(@times,PD,ernst(T1w.fa*B1map(:),T1w.TR,R1)),dims);

PDw.B1=B1map;
T1w.B1=B1map;

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

assert(all(abs(PDest(:)-PD)./PD<reltol),'Estimated PD has large error!')

%% Compare with old small angle approximation implementation

% Protocol from Weiskopf, et al. (2013).
PDw.TR=23.7e-3; % s
PDw.fa=deg2rad(6); % rad
T1w.TR=18.7e-3; % s
T1w.fa=deg2rad(20); % rad

PDw.data=reshape(bsxfun(@times,PD,ernst(PDw.fa*B1map(:),PDw.TR,R1)),dims);
T1w.data=reshape(bsxfun(@times,PD,ernst(T1w.fa*B1map(:),T1w.TR,R1)),dims);

PDw.B1=B1map;
T1w.B1=B1map;

PDest=hmri_calc_A(PDw,T1w,small_angle_approx);

assert(all(abs(PDest(:)-PD)./PD<reltol),'Estimated PD has large error!')

% Old implementation of PD calculation in hMRI toolbox
PDsa=zeros(size(PDest));
for p = 1:dims(3)
    
    PDwSlice = PDw.data(:,:,p);
    
    f_T = B1map(:,:,p); 
    
    T1wSlice = T1w.data(:,:,p);
    
    % Transmit bias corrected quantitative T1 values
    % correct T1 for transmit bias f_T with fa_true = f_T * fa_nom
    % T1corr = T1 / f_T / f_T
    T1 = ((((PDwSlice / PDw.fa) - (T1wSlice / T1w.fa)+eps) ./ ...
        max((T1wSlice * T1w.fa / 2 / T1w.TR) - (PDwSlice * PDw.fa / 2 / PDw.TR),eps))./f_T.^2);
    
    % A values proportional to PD
    PDsa(:,:,p) = T1 .* (T1wSlice .*(T1w.fa*f_T) / 2 / T1w.TR) + (T1wSlice ./ (T1w.fa*f_T));
    
end

errinerr = abs(PDest(:)-PD)./PD - abs(PDsa(:)-PD)./PD;

assert(all(errinerr(errinerr>0)<reltol),'Small angle approximation gives substantially smaller residuals in some cases!')
