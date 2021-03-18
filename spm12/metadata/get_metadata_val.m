function [parValue, parLocation] = get_metadata_val(varargin)
% This is get_metadata_val, part of the metadata library
%
% PURPOSE
% To retrieve parameter values (mainly acquisition parameters, but this can
% be developed further in the future) from nifti extended header (extended
% nifti format) or simple nifti files with associated JSON metadata. The
% first argument can either be a file name (extended nifti image, normal
% nifti image accompanied by a JSON metadata file, or a JSON metadata file)
% or the matlab structure previously extracted from such files. The second
% argument must be either a predefined parameter name (see list below), or
% a potential match for a given field in the metadata structure. The script
% can therefore be used to efficiently search any (metadata) structure.
%
% IMPORTANT NOTE FOR FUTURE DEVELOPMENTS
% Many cases below are "Siemens-specific", i.e. implemented for Siemens
% data and only valid with Siemens data. Future development will aim at
% making these cases valid for all vendors.
%
% USAGE AND EXAMPLES
% [parValue, parLocation] = get_metadata_val(mstruc, inParName)
% [parValue, parLocation] = get_metadata_val(filenam, inParName)
% parValue = get_metadata_val(filenam, 'RepetitionTime')
%
% ARGUMENTS
% - mstruc is the matlab structure containing the metadata
% - filenam is the name of a file containing or associated with JSON
%   metadata (see above).
% - inParName can be an arbitrary string that will be searched in the
%   metadata structure or one of the following predefined parameter name.
%   Names listed with an '=' are valid for all vendors, while '-' indicate
%   Siemens-only:
%    = 'RepetitionTime' TR [ms]
%    = 'RepetitionTimes' TR [ms]
%    = 'EchoTime' TE(s) [ms]
%    = 'EchoTimes' TE(s) [ms]
%    = 'FlipAngle' [deg]
%    = 'ProtocolName'
%    = 'SequenceName'
%    - 'MT' (1/0 = ON/OFF)
%    = 'FieldStrength' [T]
%    = 'Frequency' [Hz]
%    = 'ScanningSequence'
%    = 'BandwidthPerPixelRO' [Hz/Px]
%    - 'BandwidthPerPixelPE' [Hz/Px]
%    - 'PATparameters' [struct]
%    - 'AccelFactorPE'
%    - 'AccelFactor3D'
%    - 'PELines'
%    - 'PELinesPF'
%    - 'PELinesPAT'
%    - 'PhaseEncodingDirectionPositive' A>>P & R>>L = 1; P>>A & L>>R = 0.
%    = 'PhaseEncodingDirection' A>>P/P>>A = 'COL' & R>>L/L>>R = 'ROW'
%    - 'NumberOfMeasurements'
%    - 'NumberOfSlices'
%    - 'epiReadoutDuration' [ms]
%    - 'WipParameters' structure containing fields alFree & adFree
%    - 'B1mapNominalFAValues' [deg]
%    - 'RFSpoilingPhaseIncrement' [deg]
%    - 'B1mapMixingTime' [ms]
%    - 'AllDiffusionDirections' list of diffusion directions
%    - 'AllBValues' list of b-value
%    - 'DiffusionDirection' diffusion direction of individual DW image
%    - 'BValue' b-value of individual DW image
%    - 'MultiBandFactor' for CMRR multiband and Siemens SMS sequences
%
% - parValue is the value of the parameter. The type of the output can vary
%   according to the request (single real/char/complex value, an array, a
%   cellarray, ...). If inParName is not leading to a unique value, a cell
%   array of values is returned, each element corresponding to a different
%   location in the header structure, as specified by the returned
%   parLocation.
% - parLocation is a string (or a cellarray of strings if the solution is
%   not unique) giving the location of the retrieved value(s) in the
%   metadata structure. For example, if 'RepetitionTime' value is found
%   in the mstruc as mstruc.acqpar.RepetitionTime, parLocation =
%   'acqpar.RepetitionTime'.
%
%==========================================================================
% Written by Evelyne Balteau - June 2016 - Cyclotron Research Centre
%
% December 2016: modified to deal with JSON metadata contained either in
% extended header or separate JSON files.
%==========================================================================

if nargin~=2
    error('Wrong number of arguments. Type: help get_metadata_val for usage.');
end

mstruc = varargin{1};
inParName = varargin{2};

% filename as argument, first retrieve the matlab structure:
if ischar(mstruc)
  mstruc = get_metadata(mstruc);
  mstruc = mstruc{1};
end

if ~isstruct(mstruc)
  % there seems to be no metadata available (mstruc is not a structure
  % and is likely to be empty). Still try a few tricks in case either
  % TR/TE/FA is the target field and available in the description field:
  [parValue, parLocation] = get_metadata_val_header(varargin{:});
  return;
end

if isfield(mstruc, 'acqpar')
  [parValue, parLocation] = get_metadata_val_classic(varargin{:});
else
  [parValue, parLocation] = get_metadata_val_bids(varargin{:});
end

end
