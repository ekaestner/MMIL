
function [fieldnames,offsets,types] = getoffsetsX(ver)
%function [fieldnames,offsets,types] = getoffsetsX(ver)
%
% 	Returns names, types and offsets for certain variables within
%	a p-file header, given different versions.
%
%	INPUT:
%		ver	= GE software version number
%
%	OUTPUT:
%		fieldnames = list of names of each field.
%		offsets = array of byte offsets in p-file, for given version.
%		types = data types of each field.
%
%	B.Hargreaves -- June 2009.
%       Atsushi Takahashi added field 'ileaves' -- Aug 2012
%
%       now handles 12.X, 14.X, and 20.X Pfiles
%
%

% CVS Version:  (Do NOT Edit Line Below!)
% $Id: getoffsetsX.m,v 1.3 2009-09-09 01:22:48 brian Exp $
%

%
% ===========================================================
% Notes on Version.  It seems like 7 or less is before ESE 9.1
%  ver          Signa release name
%  ---          ------------------ 
%  5.0          5.X 
%  7.0          LX 
%  8.0          "New" LX 
%  9.0          EXCITE (11.0) 
%  11.0         12.0 
%  14.2         14.0 
%  14.3         14.0 M3 and above 
%  20.005       20.X 
%  20.006       20.X    (latest release as of 6/2008)

%
%	This determines the sizes of parts of the header, 
%	so it's nice to do something about it.
% ===========================================================
%
if (ver < 11) error('Versions prior to 11 not supported'); end;


% -- Each parameter in the header has a name (that will be
% -- used in the structure) a byte offset (depending on version) and a type.  
% -- We generate a "table" of fieldnames, offsets and types to pass to
% -- programs that then read/write these values from/to a header.

fieldnames = {'version','run','npasses','frsize','nframes','nslices','nechoes','rawhdrsize','hnover'};
offsets = [0,4,64,80,74,68,70,1468,78];
types = {'float','uint32','uint16','uint16','uint16','uint16','uint16','uint32','uint16'};

% -- Just adding more fields to read.
fieldnames = {fieldnames{:},'ptsize','nex','startrcvr','endrcvr','rhimsize','rhrecon','rhtype','rhdayres'};
offsets = [offsets,82,72,200,202,106,60,56,104];
types = {types{:},'uint16','uint16','uint16','uint16','uint16','uint16','uint16','uint16'};

% -- Add TEs
	fieldnames = {fieldnames{:},'te','te2'};
	offsets = [offsets,1212,1216];
	types = {types{:},'uint32','uint32'};
%
% version-specific header values
%
if (ver >= 20) 
	% -- Even more fields to read.
	fieldnames = {fieldnames{:},'rawsize','exam','series','image'};
	offsets = [offsets,1660,148712,148724,148726];
	types = {types{:},'uint64','uint16','uint16','uint16'};
    
	% -- Atsushi added interleaves....
	fieldnames = {fieldnames{:},'ileaves'};
	offsets = [offsets,914];
	types = {types{:},'uint16'};
        
        % -- Kangrong added slquant in image header.
	fieldnames = {fieldnames{:},'image_slquant'};
	offsets = [offsets,148782];
	types = {types{:},'int16'};

elseif (ver >= 14.3)
	% -- Even more fields to read.
	fieldnames = {fieldnames{:},'rawsize','exam','series','image'};
	offsets = [offsets,116,144884,144896,144898];
	types = {types{:},'uint32','uint16','uint16','uint16'};


elseif (ver >= 14.0)
%
% Note that Matlab reads 14.2 as 14.19999980926514, which is *not* greater than 14.2!
%
	% -- Even more fields to read.
	fieldnames = {fieldnames{:},'rawsize','exam','series','image'};
	offsets = [offsets,116,143384,143396,143398];
	types = {types{:},'uint32','uint16','uint16','uint16'};


else	%%% ver = 11.0
	% -- Even more fields to read.
	fieldnames = {fieldnames{:},'rawsize','exam','series','image'};
	offsets = [offsets,116,65200,65212,65214];
	types = {types{:},'uint32','uint16','uint16','uint16'};
        
        % -- Kangrong added interleaves for ver = 11.0.
	fieldnames = {fieldnames{:},'ileaves'};
	offsets = [offsets,914];
	types = {types{:},'uint16'};
end;





