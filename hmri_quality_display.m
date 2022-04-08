function h = hmri_quality_display(disp_list)
%==========================================================================
% PURPOSE: to display a series of volumes using spm_check_registration.
% USAGE:  h = hmri_quality_display(disp_list)
% where
% - disp_list is an array of structures with fields:
%    - fnam: name of the image file name
%    - title: title for the image
%    - range: intensity range [min max]
% - h is the figure handle
%==========================================================================
% Written by Evelyne Balteau
% Cyclotron Research Centre, University of Liege
%==========================================================================

global st; 
% st is a global variable storing the display parameters.
% st is a structure with fields:
%            n: 0
%         vols: {24x1 cell}
%           bb: [2x3 double]
%        Space: [4x4 double]
%       centre: [-2.4096 30.6205 -6.1265]
%     callback: ';'
%       xhairs: 1
%          hld: 1
%          fig: [1x1 Figure]
%         mode: 1
%      plugins: {'browser'  'contour'  'display'  'goto_max'  'mesh'  'movie'  'reorient'  'rgb'  'roi'  'save'  'dwprofile'  'pli'  'quiver'  'quiver3d'  'tracepath'}
%         snap: []
% st.vols is a cell array whose elements (st.vols{i}) are structures with fields:
%       fname: 'D:\home\data\truc.nii'
%         dim: [240 256 176]
%          dt: [512 0]
%       pinfo: [3x1 double]
%         mat: [4x4 double]
%           n: [1 1]
%     descrip: 'NIfTI description field'
%     private: [1x1 nifti]
%          ax: {3x1 cell}
%      premul: [4x4 double]
%      window: 'auto'
%     mapping: 'linear'
%        area: [0.0100 0.5100 0.4800 0.4800]



% gather all file names
fnamlist = char(disp_list.fnam);
% display all images
spm_check_registration(fnamlist);
% to get the section shifted in x-direction 
% (to avoid displaying the plane between left and right hemispheres)
% NB: depending on the history of manipulating st.centre, the latter can
% either be a column or a line vector. So...
st.centre(1) = st.centre(1) + 10;
st.centre(2) = st.centre(2) + 0;
st.centre(3) = st.centre(3) + 0;
% st.centre = st.centre + [10 0 0];

% to zoom in
%st.Space = [0.25 0 0 0;0 0.25 0 0;0 0 0.25 0;0 0 0 1];

% to give a title to each orthviews and adjust intensity range:
for cim=1:length(disp_list)
    txt = disp_list(cim).title;
    % rearrange title if too long 
    % (split over several lines of maximum maxchar characters length)  
    maxchar = 20;
    if length(txt)>maxchar
        try
            txtsplit = strsplit(txt); % problem when using Matlab < R2013a (when strsplit and strjoin were added)
        catch %#ok<*CTCH>
            txtsplit = regexp(txt,regexptranslate('escape',' '),'split');
        end
        clear linsplit;
        if length(txtsplit)> 1 % text naturally splittable
            linsplit{1} = txtsplit{1};
            cw = 2;
            cl = 1;
            while cw<length(txtsplit)+1
                tmp = [linsplit{cl} ' ' txtsplit{cw}];
                if length(tmp)<maxchar
                    linsplit{cl} = tmp; %#ok<*AGROW>
                else
                    cl = cl+1;
                    linsplit{cl} = txtsplit{cw};
                end
                cw = cw+1;
            end
        else % force split
            cl = 1;
            while length(txt)>maxchar
                linsplit{cl} = txt(1:maxchar);
                txt = txt(maxchar+1:end);
                cl = cl+1;
            end
            linsplit{cl} = txt;
        end
            
        try 
            txt = strjoin(linsplit,'\n'); % problem when using Matlab < R2013a (when strsplit and strjoin were added)
        catch
            txt = linsplit{1};
            for cc = 2:length(linsplit)
                txt = sprintf('%s\n%s', txt, linsplit{cc});
            end
        end
    end
    
    %htit = get(st.vols{cim}.ax{3}.ax,'Title');
    hxlab = get(st.vols{cim}.ax{3}.ax,'Xlabel');
    if any( disp_list(cim).range - round( disp_list(cim).range))
        set(hxlab, 'String',sprintf('%s\n\nIntensity range\n[%.1f %.1f]', txt, disp_list(cim).range),'FontSize',10);
        % set(st.vols{cim}.ax{3}.ax, 'Xlabel',text('String',sprintf('Intensity range\n[%.1f %.1f]', disp_list(cim).range)));
    else
        set(hxlab, 'String',sprintf('%s\n\nIntensity range\n[%d %d]', txt, disp_list(cim).range),'FontSize',10);
        % set(st.vols{cim}.ax{3}.ax, 'Xlabel',text('String',sprintf('Intensity range\n[%d %d]', disp_list(cim).range)));
    end
    %set(htit, 'String',sprintf('%s',disp_list(cim).title));
    % set(st.vols{cim}.ax{3}.ax, 'Title',text('String',sprintf('%s',disp_list(cim).title)));
    st.vols{cim}.window = disp_list(cim).range;
end

spm_orthviews('Redraw');

h = st.fig;

end