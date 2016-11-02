% Code for processing of the multi-parameter QA data
% Datasets required: two runs of the MTw, PDw and T1w each accompanied by
% their own B1 and B0 mapping data. Here the B1 mapping data is assumed to
% be acquired with the 3D EPI SE/STE method (Lutti A. et al PLoS ONE 2012)
% with the 'v2b' settings standard for 3T acquisitions
% The code runs as follows:
% 1. Loads all images. For each image contrast load all acquired echoes in
% increasing order
% 2. Calculation of the 2 B1 maps
% 3. Calculation of the MT, R1 and R2s maps.
% 4. Quantitative analysis of the stability of the MT, R1, R2s and B1
% values across the two runs.
% Notes: a. Specific 'MTProtQA.m' version of vbq_MTProt.m with settings
% optimized for processing of phantom data: i. No coregistration between the
% runs ii. Dynamic range of the R1 map set between -1e4 and 1e4 to account for
% the shorter R1 values of the phantom iii. Calculation of the proton
% density maps is disabled.
%       b. This code was developed for processing of data acquired on a
%       Siemens doped water bottle. The ROI for the quantitative analysis
%       was designed accordingly, accounting for the strong motion artefacts
%       visible in the data due to vibration of the water during the scan.
%       c. All images (.nii and .fig) are saved in a separate 'MPMQA_analysis'
%       folder located in the folder containing the MTw images of the 1st run.
% Output of the analysis:
% 1. DMT.fig: changes in MT across the two repetitions. Sagittal, Coronal
% and Axial views. Also contains a histogram of the MT changes with the
% mean and standard deviation
% 2. DR1.fig: same as 1. for the R1 values.
% 3. DR2s.fig: same as 1. for the R2s values.
% 4. DMT, DR1 and DR2s are also available in .nii format for visualization
% Antoine Lutti, Wellcome Trust Centre for Neuroimaging at UCL, London
% 03/2012

function vbq_MPMQA(ImagedObject,ImageLoading,QAMode)

% Check inputs
while (~strcmp(ImagedObject,'Fbirn')&&~strcmp(ImagedObject,'WaterBottle'))
    ImagedObject = input('Unknown phantom. Enter phantom type (Fbirn or WaterBottle)  ','s');
end

while (~strcmp('Manual',ImageLoading)&&~strcmp('Auto',ImageLoading))
    ImageLoading = input('Unknown image loading. Enter image loading (AutoLoading or ManualLoading)  ','s');
end

while (~strcmp('Offline',QAMode)&&~strcmp('Online',QAMode))
    QAMode = input('Unknown QA mode. Enter QA mode (Offline or Online)  ','s');
end
[~,host_name]=unix('hostname');
if ((strcmp('mr2',host_name)||strcmp('mr2b',host_name))&&strcmp('Offline',QAMode))
    warning('Offline mode not allowed on this machine; QA Mode changed to Online')
    QAMode='Online';
end
if strcmp('Online',QAMode)&&strcmp('Manual',ImageLoading)
    warning('Manual loading; Loading mode changed to Auto')
    ImageLoading='Auto';
end    
    
if strcmp('Fbirn',ImagedObject)
    T1=570;%in ms
elseif strcmp('WaterBottle',ImagedObject)
    T1=130;%in ms
end

% Set up of the directories, dicom conversion for online mode
if strcmp('Offline',QAMode)
    scanner = host_name;
    images_dir=pwd;
    target_dir=images_dir;
    qafilename=[target_dir filesep sprintf('mpm_qa_%s.txt',deblank(scanner))];
elseif strcmp('Online',QAMode)
    if strcmp(lower(deblank(host_name)),'mr2')
        scanner = 'quattro';
        qa_folder =[filesep 'data' filesep scanner filesep 'MPMqa'];%Output folder containing QA results.
        images_dir='/dicom_from_scanner';
        target_dir='/home/dataman/physics/qa_tmp';
        qafilename=[qa_folder filesep sprintf('mpm_qa_%s.txt',deblank(scanner))];
    elseif strcmp(lower(deblank(host_name)),'mr2b')
        scanner = 'trio';
        qa_folder =[filesep 'data' filesep scanner filesep 'MPMqa'];
        images_dir='/dicom_from_scanner';
        target_dir='/home/dataman/physics/qa_tmp';
        qafilename=[qa_folder filesep sprintf('mpm_qa_%s.txt',deblank(scanner))];
    else%for debugging of 'online' mode. No difference with offline mode
        scanner = host_name;
        images_dir=pwd;
        target_dir=images_dir;
%         target_dir=fullfile(images_dir,'niftis');%fullfile('..',qa_folder);
%         if(~exist(target_dir,'dir'))
%             mkdir(target_dir);
%         end
        qafilename=[target_dir filesep sprintf('mpm_qa_%s.txt',deblank(scanner))];
    end
    P=spm_select(Inf,'^.*\.ima$','Select ima files',[],images_dir);
    if ~strcmp(lower(deblank(host_name)),'mr2')&&~strcmp(lower(deblank(host_name)),'mr2b')
        if ~strcmp(images_dir,spm_str_manip(P,'H'))% when host name is not mr2 or mr2b, images_dir and target_dir were set to pwd. This is changed here to ensure that the niftis are saved somewhere appropriate
            images_dir=spm_str_manip(P(1,:),'H');target_dir=spm_str_manip(P(1,:),'H');
        end
    end
    hdrs = spm_dicom_headers(P);
    current_dir=pwd;
    eval(['cd ' target_dir]);
    spm_dicom_convert(hdrs);
    eval(['cd ' current_dir]);
end

% Load images
if strcmp('Manual',ImageLoading)
    [MTw1,PDw1,T1w1,B11,B01,MTw2,PDw2,T1w2,B12,B02]=ManualLoading(images_dir);
elseif strcmp('Auto',ImageLoading)
%     if strcmp('Offline',QAMode)
%         ImageFolder=spm_str_manip(spm_select(Inf,'dir','Image folder',[],images_dir),'H');
%     else
%         ImageFolder=target_dir;
%     end
    Niftis=spm_select('FPList',images_dir,'^s.*.(img|nii)$');
    while (size(Niftis,1)~=90)
        images_dir=spm_str_manip(spm_select(Inf,'dir','!!Wrong number of files. Image folder'),'H');
        Niftis=spm_select('FPList',images_dir,'^s.*.(img|nii)$');
    end
    B11 = Niftis(1:22,:);B01 = Niftis(23:25,:);
    MTw1 = Niftis(26:31,:);PDw1 = Niftis(32:39,:);T1w1 = Niftis(40:45,:);
    B12 = Niftis(1+45:22+45,:);B02 = Niftis(23+45:25+45,:);
    MTw2 = Niftis(26+45:31+45,:);PDw2 = Niftis(32+45:39+45,:);T1w2 = Niftis(40+45:45+45,:);
end


disp('----- Calculation of B1 maps -----');
vbq_B1map_v2(B11,B01,T1);
if (~isequal(B11,B12))
    vbq_B1map_v2(B12,B02,T1);
end
disp('----- Calculation of quantitative maps -----');
P_trans=spm_select('FPList',spm_str_manip(B11(1,:),'h'),'^(uSumOfSq|smu).*\.(img|nii)$');
P_trans(1:end,:)=P_trans(end:-1:1,:);%uSumOfSq has to appear first
MT_analysis_QA(MTw1,PDw1,T1w1,P_trans,'')

P_trans=spm_select('FPList',spm_str_manip(B12(1,:),'h'),'^(uSumOfSq|smu).*\.(img|nii)$');
P_trans(1:end,:)=P_trans(end:-1:1,:);
MT_analysis_QA(MTw2,PDw2,T1w2,P_trans,'');%uSumOfSq has to appear first

disp('----- Analysis - quantitative maps -----');

MT1=spm_select('FPList',spm_str_manip(MTw1(1,:),'h'),'.*\_MT.(img|nii)$');
MT2=spm_select('FPList',spm_str_manip(MTw2(1,:),'h'),'.*\_MT.(img|nii)$');

R11=spm_select('FPList',spm_str_manip(MTw1(1,:),'h'),'.*\_R1.(img|nii)$');
R12=spm_select('FPList',spm_str_manip(MTw2(1,:),'h'),'.*\_R1.(img|nii)$');

R2s1=spm_select('FPList',spm_str_manip(MTw1(1,:),'h'),'.*\_R2s.(img|nii)$');
R2s2=spm_select('FPList',spm_str_manip(MTw2(1,:),'h'),'.*\_R2s.(img|nii)$');

% Below is useful when quantitative maps from both runs are saved at the
% same level
MT1=squeeze(MT1(1,:));R11=squeeze(R11(1,:));R2s1=squeeze(R2s1(1,:));
if (size(R2s2,1)==2)
    R2s2=squeeze(R2s2(2,:));
end
if (size(R12,1)==2)
R12=squeeze(R12(2,:));
end
if (size(MT2,1)==2)
MT2=squeeze(MT2(2,:));
end

[meanMT,DMT,meanR1,DR1,meanR2s,DR2s,avgMT,SDMT,stabMT,avgR1,SDR1,stabR1,avgR2s,SDR2s,stabR2s]=MPMQA_calc(MT1,MT2,R11,R12,R2s1,R2s2,ImagedObject);

disp('----- Analysis - B1 maps -----');
[meanB1,DB1,avgB1,SDB1,stabB1]=B1anal(B11(1,:),B12(1,:),ImagedObject);

% if (~isequal(B11,B12))
%     disp('----- Analysis - B1 maps -----');
%     [meanB1,DB1,avgB1,SDB1,stabB1]=B1anal(B11(1,:),B12(1,:),ImagedObject);
%     % MaskSSQ(spm_select('FPList',spm_str_manip(B12(1,:),'h'),'^(uSumOfSq).*\.(img|nii)$'));
% else
% %     meanB1=zeros(size(spm_read_vols(spm_vol(B11(1,:))),1),size(spm_read_vols(spm_vol(B11(1,:))),2),size(spm_read_vols(spm_vol(B11(1,:))),3));
%     meanB1=spm_read_vols(spm_vol(spm_select('FPList',spm_str_manip(B11(1,:),'h'),'^(smu).*\.(img|nii)$')));
%     DB1=zeros(size(spm_read_vols(spm_vol(B11(1,:))),1),size(spm_read_vols(spm_vol(B11(1,:))),2),size(spm_read_vols(spm_vol(B11(1,:))),3));
%     avgB1=zeros(size(spm_read_vols(spm_vol(B11(1,:))),1),size(spm_read_vols(spm_vol(B11(1,:))),2),size(spm_read_vols(spm_vol(B11(1,:))),3));
%     SDB1=zeros(size(spm_read_vols(spm_vol(B11(1,:))),1),size(spm_read_vols(spm_vol(B11(1,:))),2),size(spm_read_vols(spm_vol(B11(1,:))),3));
%     stabB1=zeros(size(spm_read_vols(spm_vol(B11(1,:))),1),size(spm_read_vols(spm_vol(B11(1,:))),2),size(spm_read_vols(spm_vol(B11(1,:))),3));
% end
% Write out QA results to to_network/qa/qa_stability.txt
fid = fopen(qafilename,'a');
fprintf(fid,'%s\t mean MT: %1.2f\t SD MT: %1.2f\t stab MT: %1.2f\t mean R1: %4.0f\t SD R1: %3.0f\t stab R1: %3.0f\t mean R2s (s-1): %1.4f\t SD R2s (s-1): %1.4f\t stab R2s (s-1): %1.4f\t mean B1: %1.4f\t SD B1: %1.4f\t stab B1: %1.4f\n',...
    datestr(now),avgMT,SDMT,stabMT,avgR1,SDR1,stabR1,avgR2s*1000,SDR2s*1000,stabR2s*1000,avgB1,SDB1,stabB1);
fclose(fid);

processingtime=clock;
% DrawImage(DMT,DR1,DR2s,DB1,target_dir,['MPMstability-',[datestr(now,'yyyymmdd'),'-',sprintf('%2.0f',processingtime(4)),'h',sprintf('%2.0f',processingtime(5))]])
% DrawImage(meanMT,meanR1,meanR2s,meanB1,target_dir,['MPMaccuracy-',[datestr(now,'yyyymmdd'),'-',sprintf('%2.0f',processingtime(4)),'h',sprintf('%2.0f',processingtime(5))]])
DrawImage(DMT,[-0.15 0.15],DR1,[-100 100],DR2s,[-3 3],DB1,[-0.8 0.8],target_dir,['MPMstability-',[datestr(now,'yyyymmdd'),'-',sprintf('%2.0f',processingtime(4)),'h',sprintf('%2.0f',processingtime(5))]])
DrawImage(meanMT,[0 0.3],meanR1,[1500 2000],meanR2s,[15 30],meanB1,[75 125],target_dir,['MPMaccuracy-',[datestr(now,'yyyymmdd'),'-',sprintf('%2.0f',processingtime(4)),'h',sprintf('%2.0f',processingtime(5))]])

% % Removes all created niftis from the 'MPMQAresults' folder
% Niftis=spm_select('FPList',images_dir,'.(hdr|img|nii)$');
% for counter=1:size(Niftis,1)
%     delete(deblank(Niftis(counter,:)));
% end

end

function DrawImage(Mat1,lim1,Mat2,lim2,Mat3,lim3,Mat4,lim4,target_dir,ImageName)


FS=6;
L1=[0.0 0.25 0.52 0.75];
H1=[0.04 0.29 0.54 0.78];
L2=[0.2 0.2 0.2 0.4];
H2=[0.2 0.2];
rect11=[L1(1) H1(4) L2(1) H2(1)];
rect12=[L1(1) H1(3) L2(1) H2(1)];
rect21=[L1(2) H1(4) L2(1) H2(1)];
rect22=[L1(2) H1(3) L2(1) H2(1)];
rect13=[L1(1) H1(2) L2(1) H2(1)];
rect14=[L1(1) H1(1) L2(1) H2(1)];
rect23=[L1(2) H1(2) L2(3) H2(1)];
rect24=[L1(2) H1(1) L2(3) H2(1)];
rect31=[L1(3) H1(4) L2(2) H2(1)];
rect32=[L1(3) H1(3) L2(2) H2(1)];
rect41=[L1(4) H1(4) L2(2) H2(1)];
rect42=[L1(4) H1(3) L2(2) H2(1)];
rect33=[L1(3) H1(2) L2(2) H2(2)];
rect34=[L1(3) H1(1) L2(2) H2(2)];
rect43=[L1(4) H1(2) L2(2) H2(2)];
rect44=[L1(4) H1(1) L2(2) H2(2)];

figure
% Top left corner
axes('Position',rect11);
imagesc(squeeze(Mat1(size(Mat1,1)/2,end:-1:1,:)),lim1);
% title([inputname(1) '-Coronal'],'FontSize',8,'Position',[rect11(1)+80 rect11(2)+50])
title([inputname(1) '-Coronal'],'FontSize',8,'Position',[rect11(1)+80 rect11(2)+50],'FontWeight','bold','Color','red')
% title('\color[rgb]{1 0 0} Coronal','FontSize',8,'FontWeight','bold')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'ZTick',-1:1:1,'FontSize',FS);
colorbar('East','FontSize',8)
axes('Position',rect21);
imagesc(squeeze(Mat1(:,end:-1:1,size(Mat1,3)/2)'),lim1);
title([inputname(1) '-Sagittal'],'FontSize',8,'Position',[rect21(1)+80 rect21(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect12);
imagesc(squeeze(Mat1(:,size(Mat1,2)/2,:)),lim1);
title([inputname(1) '-Transverse'],'FontSize',8,'Position',[rect12(1)+80 rect12(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect22);
test=Mat1(Mat1~=0);
hist(test,100,[-1 1])
title(sprintf('%1.2f\t +- %1.2f\t', mean(test), std(test,[],1)),'FontSize',8,'FontWeight','bold','Color','red')%,'Position',[rect22(1) rect22(2)+22000]
% set(gca,'XTick',-1:1:1);
% xlim([-1 1])
set(gca,'FontSize',FS);
    
% Top right corner
axes('Position',rect31);
imagesc(squeeze(Mat2(size(Mat2,1)/2,end:-1:1,:)),lim2);
title([inputname(3) '-Coronal'],'FontSize',8,'Position',[rect31(1)+80 rect31(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'ZTick',-1:1:1,'FontSize',FS);
colorbar('East','FontSize',8)
axes('Position',rect41);
imagesc(squeeze(Mat2(:,end:-1:1,size(Mat2,3)/2)'),lim2);
title([inputname(3) '-Sagittal'],'FontSize',8,'Position',[rect41(1)+80 rect41(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect32);
imagesc(squeeze(Mat2(:,size(Mat2,2)/2,:)),lim2);
title([inputname(3) '-Transverse'],'FontSize',8,'Position',[rect32(1)+80 rect32(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect42);
test=Mat2(Mat2~=0);
hist(test,100,[-1 1])
title(sprintf('%1.2f\t +- %1.2f\t', mean(test), std(test,[],1)),'FontSize',8,'FontWeight','bold','Color','red')%,'Position',[rect42(1)+80 rect42(2)+50]
% set(gca,'XTick',-1:1:1);
% xlim([-1 1])
set(gca,'FontSize',FS);
    
Mat3=Mat3*1000;%For convenience. Mat3 is R2s or DR2s
% Bottom left corner
axes('Position',rect13);
imagesc(squeeze(Mat3(size(Mat3,1)/2,end:-1:1,:)),lim3);
title([inputname(5) '-Coronal'],'FontSize',8,'Position',[rect13(1)+80 rect13(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'ZTick',-1:1:1,'FontSize',FS);
colorbar('East','FontSize',8)
axes('Position',rect23);
imagesc(squeeze(Mat3(:,end:-1:1,size(Mat3,3)/2)'),lim3);
title([inputname(5) '-Sagittal'],'FontSize',8,'Position',[rect23(1)+80 rect23(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect14);
imagesc(squeeze(Mat3(:,size(Mat3,2)/2,:)),lim3);
title([inputname(5) '-Transverse'],'FontSize',8,'Position',[rect14(1)+80 rect14(2)+50],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect24);
test=Mat3(Mat3~=0);
hist(test,100,[-1 1])
title(sprintf('%1.2f\t +- %1.2f\t (s-1)', mean(test), std(test,[],1)),'FontSize',8,'FontWeight','bold','Color','red')%,[rect24(1)+80 rect24(2)+50]
% set(gca,'XTick',-1:1:1);
% xlim([-1 1])
set(gca,'FontSize',FS);

% Bottom right corner
axes('Position',rect33);
imagesc(permute(squeeze(Mat4(end:-1:1,size(Mat4,2)/2,end:-1:1)),[2 1]),lim4);
title([inputname(7) '-Coronal'],'FontSize',8,'Position',[rect33(1)+25 rect33(2)+10],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'ZTick',-1:1:1,'FontSize',FS);
colorbar('East','FontSize',8)
axes('Position',rect43);
imagesc(permute(squeeze(Mat4(size(Mat4,1)/2,end:-1:1,end:-1:1)),[2 1]),lim4);
title([inputname(7) '-Sagittal'],'FontSize',8,'Position',[rect43(1)+25 rect43(2)+10],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS)
axes('Position',rect34);
imagesc(permute(squeeze(Mat4(end:-1:1,end:-1:1,size(Mat4,3)/2)),[2 1]),lim4);
title([inputname(7) '-Transverse'],'FontSize',8,'Position',[rect34(1)+25 rect34(2)+10],'FontWeight','bold','Color','red')
axis equal
set(gca,'XTick',5000:1000:6000,'YTick',5000:1000:6000,'FontSize',FS,'FontWeight','bold','Color','red')
axes('Position',rect44);
test=Mat4(Mat4~=0);
hist(test,100,[-1 1])
title(sprintf('%1.2f\t +- %1.2f\t', mean(test), std(test,[],1)),'FontSize',8,'FontWeight','bold','Color','red')%,'Position',[rect44(1)+80 rect44(2)+50]
% set(gca,'XTick',-1:1:1);
% xlim([-1 1])
set(gca,'FontSize',FS);
print(gcf,'-dtiff',fullfile(target_dir,[ImageName,'.tiff']))

end

function [MTw1,PDw1,T1w1,B11,B01,MTw2,PDw2,T1w2,B12,B02]=ManualLoading(images_dir)

    MTw1 = spm_select(Inf,'image','Select MTw images - run1',[],images_dir);
    while (size(MTw1,1)~=6)
        MTw1 = spm_select(Inf,'image','!!Wrong number of files!! Select MTw images - run1',[],images_dir);
    end
    PDw1 = spm_select(Inf,'image','Select PDw images - run1',[],images_dir);
    while (size(PDw1,1)~=8)
        PDw1 = spm_select(Inf,'image','!!Wrong number of files!! Select PDw images - run1',[],images_dir);
    end
    T1w1 = spm_select(Inf,'image','Select T1w images - run1',[],images_dir);
    while (size(T1w1,1)~=6)
        T1w1 = spm_select(Inf,'image','!!Wrong number of files!! Select T1w images - run1',[],images_dir);
    end
    B11 = spm_select(Inf,'image','Select B1 mapping data - run1',[],images_dir);
    while (size(B11,1)~=22)
        B11 = spm_select(Inf,'image','!!Wrong number of files!! Select B1 mapping data - run1',[],images_dir);
    end
    B01 = spm_select(Inf,'image','Select B0 mapping data - run1',[],images_dir);
    while (size(B01,1)~=3)
        B01 = spm_select(Inf,'image','!!Wrong number of files!! Select B0 mapping data - run1',[],images_dir);
    end
    
    MTw2 = spm_select(Inf,'image','Select MTw images - run2',[],images_dir);
    while (size(MTw2,1)~=6)
        MTw2 = spm_select(Inf,'image','!!Wrong number of files!! Select MTw images - run2',[],images_dir);
    end
    PDw2 = spm_select(Inf,'image','Select PDw images - run2',[],images_dir);
    while (size(PDw2,1)~=8)
        PDw2 = spm_select(Inf,'image','!!Wrong number of files!! Select PDw images - run2',[],images_dir);
    end
    T1w2 = spm_select(Inf,'image','Select T1w images - run2',[],images_dir);
    while (size(T1w2,1)~=6)
        T1w2 = spm_select(Inf,'image','!!Wrong number of files!! Select T1w images - run2',[],images_dir);
    end
    B12 = spm_select(Inf,'image','Select B1 mapping data - run2',[],images_dir);
    while (size(B12,1)~=22)
        B12 = spm_select(Inf,'image','!!Wrong number of files!! Select B1 mapping data - run2',[],images_dir);
    end
    B02 = spm_select(Inf,'image','Select B0 mapping data - run2',[],images_dir);
    while (size(B02,1)~=3)
        B02 = spm_select(Inf,'image','!!Wrong number of files!! Select B0 mapping data - run2',[],images_dir);
    end

end
function [meanMT,DMT,meanR1,DR1,meanR2s,DR2s,avgMT,SDMT,stabMT,avgR1,SDR1,stabR1,avgR2s,SDR2s,stabR2s]=MPMQA_calc(MT1,MT2,R11,R12,R2s1,R2s2,ImagedObject)
% OutputDir=fullfile(spm_str_manip(MT1(1,:),'h'),'MPMQA_analysis');
% if(~exist(OutputDir,'dir'))
%     mkdir(OutputDir);
% end

% Creates mask for analysis of the quantitative maps
Mask=zeros(size(spm_read_vols(spm_vol(R11)),1),size(spm_read_vols(spm_vol(R11)),2),size(spm_read_vols(spm_vol(R11)),3));
if strcmp('WaterBottle',ImagedObject)
    Mask(90:150,80:200,60:120)=1;
else
    SegmentData(MT1(1,:));
    V=spm_vol(spm_select('FPList',spm_str_manip(MT1(1,:),'h'),'^(c1|c2|c3).*\.(img|nii)$'));
    Tissue=sum(spm_read_vols(V),4);
    Mask(Tissue>0.8)=1;
    temp=spm_select('FPList',spm_str_manip(MT1,'h'),'^c.*\.(img|nii)$');
    for counter=1:size(temp,1)
        delete(deblank(temp(counter,:)));
    end
    delete(spm_select('FPList',spm_str_manip(MT1,'h'),'^.*_seg8'));

end

DMT=(spm_read_vols(spm_vol(MT1))-spm_read_vols(spm_vol(MT2))).*Mask;
meanMT=(spm_read_vols(spm_vol(MT1))+spm_read_vols(spm_vol(MT2)))/2.*Mask;
DR1=(spm_read_vols(spm_vol(R11))-spm_read_vols(spm_vol(R12))).*Mask;
meanR1=(spm_read_vols(spm_vol(R11))+spm_read_vols(spm_vol(R12)))/2.*Mask;
DR2s=(spm_read_vols(spm_vol(R2s1))-spm_read_vols(spm_vol(R2s2))).*Mask;
meanR2s=(spm_read_vols(spm_vol(R2s1))+spm_read_vols(spm_vol(R2s2)))/2.*Mask;

avgMT=mean(meanMT(meanMT~=0));
SDMT=std(meanMT(meanMT~=0),[],1);%spatial homogeneity of MT
stabMT=mean(DMT(DMT~=0));%stability of MT across repetitions
avgR1=mean(meanR1(meanR1~=0));SDR1=std(meanR1(meanR1~=0),[],1);stabR1=mean(DR1(DR1~=0));
avgR2s=mean(meanR2s(meanR2s~=0));SDR2s=std(meanR2s(meanR2s~=0),[],1);stabR2s=mean(DR2s(DR2s~=0));

end
function MT_analysis_QA(P_mtw,P_pdw,P_t1w,P_trans,P_receiv)

% $Id$

if nargin==0,
    P_mtw = spm_select(Inf,'nifti','MT-weighted');
    P_pdw = spm_select(Inf,'nifti','PD-weighted');
    P_t1w = spm_select(Inf,'nifti','T1-weighted');
    P_trans = spm_select([0 2],'nifti','B1 map: T1w+map');
    P_receiv = spm_select([0 2],'nifti','Sensitivity map: T1w+map');
end

p = hinfo(P_mtw);
TE_mtw = cat(1,p.te);
TR_mtw = p(1).tr;
fa_mtw = p(1).fa;

p = hinfo(P_pdw);
TE_pdw = cat(1,p.te);
TR_pdw = p(1).tr;
fa_pdw = p(1).fa;

p = hinfo(P_t1w);
TE_t1w = cat(1,p.te);
TR_t1w = p(1).tr;
fa_t1w = p(1).fa;

vbq_MTProtQA(P_mtw, P_pdw, P_t1w, TE_mtw, TE_pdw, TE_t1w, TR_mtw, TR_pdw, TR_t1w, fa_mtw, fa_pdw, fa_t1w, P_trans, P_receiv);
end

function [meanB1,DB1,avgB1,SDB1,stabB1]=B1anal(B11files,B12files,ImagedObject)
SSQ=spm_select('FPList',spm_str_manip(B11files,'h'),'^(uSumOfSq).*\.(img|nii)$');
SegmentData(SSQ)

V=spm_vol(spm_select('FPList',spm_str_manip(SSQ,'h'),'^(c1|c2|c3).*\.(img|nii)$'));
Tissue=sum(spm_read_vols(V),4);

Mask=zeros(size(Tissue,1),size(Tissue,2),size(Tissue,3));
Mask(Tissue>0.8)=1;
if strcmp('WaterBottle',ImagedObject)
    Mask(:,:,37:end)=0;%zeroes the end of the bottle where the B1 is wrong
end
Vsave=spm_vol(spm_select('FPList',spm_str_manip(B11files(1,:),'h'),'^(smu).*\.(img|nii)$'));
Vsave=Vsave(1);
Vsave=rmfield(Vsave,'pinfo');
% Vsave.fname=fullfile(spm_str_manip(B11files(1,:),'h'),'mask.img');
% spm_write_vol(Vsave,Mask);

temp=spm_select('FPList',spm_str_manip(SSQ,'h'),'^c.*\.(img|nii)$');
for counter=1:size(temp,1)
    delete(deblank(temp(counter,:)));
end
delete(spm_select('FPList',spm_str_manip(SSQ,'h'),'^.*_seg8'));

if (~isequal(B11files,B12files))
    if strcmp(spm_str_manip(B11files(1,:),'h'),spm_str_manip(B12files(1,:),'h'))%i.e. if the two B1 maps are saved at the same level
        B11=spm_read_vols(spm_vol(spm_select('FPList',spm_str_manip(B11files(1,:),'h'),'^(smu).*\.(img|nii)$')));
        B12=squeeze(B11(:,:,:,2));B11=squeeze(B11(:,:,:,1));
    else
        B11=spm_read_vols(spm_vol(spm_select('FPList',spm_str_manip(B11files(1,:),'h'),'^(smu).*\.(img|nii)$')));
        B12=spm_read_vols(spm_vol(spm_select('FPList',spm_str_manip(B12files(1,:),'h'),'^(smu).*\.(img|nii)$')));
    end
else%i.e. only one B1 map was acquired
    B11=spm_read_vols(spm_vol(spm_select('FPList',spm_str_manip(B11files(1,:),'h'),'^(smu).*\.(img|nii)$')));
    B12=B11;
end
DB1=(B11-B12).*Mask;meanB1=(B11+B12)/2.*Mask;
avgB1=mean(meanB1(meanB1~=0));SDB1=std(meanB1(meanB1~=0),[],1);stabB1=mean(DB1(DB1~=0));

% % Vsave=spm_vol(spm_select('FPList',spm_str_manip(B11files(1,:),'h'),'^(smu).*\.(img|nii)$'));
% Vsave.fname=fullfile(SaveDir,'DB1.nii');
% spm_write_vol(Vsave,DB1);
% 
% % Creates and saves B1-change figure
% figure
% title('Change in B1')
% subplot(2,2,1)
% imagesc(permute(squeeze(DB1(end:-1:1,size(DB1,2)/2,end:-1:1)),[2 1]));
% title('Coronal')
% axis equal
% subplot(2,2,2)
% imagesc(permute(squeeze(DB1(size(DB1,1)/2,end:-1:1,end:-1:1)),[2 1]));
% title('Sagittal')
% axis equal
% subplot(2,2,3)
% % imagesc(squeeze(DB1(:,:,size(DB1,3)/2)));
% imagesc(permute(squeeze(DB1(end:-1:1,end:-1:1,size(DB1,3)/2)),[2 1]));
% title('Transverse')
% axis equal
% subplot(2,2,4)
% test=DB1(DB1~=0);
% hist(test,100)
% title(['Mean: ' num2str(mean(test)) '; SD: ' num2str(std(test,[],1))])
% saveas(gcf,fullfile(SaveDir,'DB1.fig'))
% close all
end

function SegmentData(Path)
clear matlabbatch
matlabbatch{1}.spm.tools.preproc8.channel.vols = {Path};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,1'};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,2'};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,3'};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,4'};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,5'};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {'D:\home\Software\SPM8\toolbox\Seg\TPM.nii,6'};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end

function RenameHeader(MTw2,PDw2,T1w2)
V=spm_vol(spm_select('FPList',spm_str_manip(MTw2(1,:),'h'),'^(s|r).*\.(img|nii)$'));
for Runctr=1:size(V,1)/2
    V(Runctr).descrip=V(Runctr+size(V,1)/2).descrip;
    spm_write_vol(V(Runctr),spm_read_vols(V(Runctr)));
end
V=spm_vol(spm_select('FPList',spm_str_manip(PDw2(1,:),'h'),'^(s|r).*\.(img|nii)$'));
for Runctr=1:size(V,1)/2
    V(Runctr).descrip=V(Runctr+size(V,1)/2).descrip;
    spm_write_vol(V(Runctr),spm_read_vols(V(Runctr)));
end
V=spm_vol(spm_select('FPList',spm_str_manip(T1w2(1,:),'h'),'^(s|r).*\.(img|nii)$'));
for Runctr=1:size(V,1)/2
    V(Runctr).descrip=V(Runctr+size(V,1)/2).descrip;
    spm_write_vol(V(Runctr),spm_read_vols(V(Runctr)));
end
end
function CoregRuns(MTw1,PDw1,T1w1,MTw2,PDw2,T1w2)

for CtrstCtr=1:nargin/2
    if(CtrstCtr==1)
        % MTw
        RefImage=cellstr(MTw1(1,:));
        SourceImage=cellstr(MTw2(1,:));
        OtherImage=cellstr(MTw2(2:end,:));
    elseif(CtrstCtr==2)
        % PDw
        RefImage=cellstr(PDw1(1,:));
        SourceImage=cellstr(PDw2(1,:));
        OtherImage=cellstr(PDw2(2:end,:));
    elseif(CtrstCtr==2)
        % T1w
        RefImage=cellstr(T1w1(1,:));
        SourceImage=cellstr(T1w2(1,:));
        OtherImage=cellstr(T1w2(2:end,:));
    end
    
    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {RefImage{1}};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {SourceImage{1}};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {OtherImage{1:end}};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);
    
end
end

