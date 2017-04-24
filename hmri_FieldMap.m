function varargout = hmri_FieldMap(varargin)
% This is a fraction of the FieldMap script (case 'createfieldmap')
% rewritten for the hMRI toolbox in order to make use of the original
% SPM12's FieldMap script wherever no modification is required.
%
% Implementation following FieldMap version:
% $Id: FieldMap.m 6994 2017-01-26 16:19:14Z guillaume $
% 
% FORMAT IP.fm = hmri_FieldMap('CreateFieldMap',IP);
%=======================================================================

persistent IP            % Input and results

if nargin == 0
   warning('Can only be called for map creation.\nFORMAT: FM = hmri_FieldMap(''CreateFieldMap'',IP);');
   return;
else
   Action = varargin{1};
end

switch lower(Action)
%=======================================================================
%
% Create unwrapped phase-map - NO gui
%
%=======================================================================

   case 'createfieldmap'

      IP=varargin{2};

      % First check that images are in same space. 

      if size([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],2)==4
         ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim' IP.P{3}.dim' IP.P{4}.dim']');
         %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]');
         ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]'));
      else
         ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim']');
         %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:)]');
         ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:)]'));
      end

      if any(any(diff(ip_dim,1,1),1)&[1,1,1])
        errordlg({'Images don''t all have same dimensions'});
        drawnow;
        varargout{1}=[];
      elseif any(any(abs(diff(ip_mat,1,1))>1e-4))
        errordlg({'Images don''t all have same orientation & voxel size'});
        drawnow;
        varargout{1}=[];
      else
         % Update flags for unwarping (in case TEs have been adjusted
         IP.uflags.etd = IP.et{2}-IP.et{1};      

         % Clear any brain mask 
         IP.uflags.bmask = [];

         % SPM5 Update
         % If flag selected to mask brain and the field map is not based on EPI
         if IP.maskbrain==1 
            IP.fmagP = FieldMap('Magnitude',IP);
            % IP.uflags.bmask = pm_brain_mask(IP.fmagP,IP.mflags);
            IP.uflags.bmask = hmri_pm_brain_mask(IP.fmagP,IP.mflags);
            varargout{2} = IP.fmagP;
         end

         IP.fm = pm_make_fieldmap([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],IP.uflags);
         varargout{1} = IP.fm;
      end
      
end