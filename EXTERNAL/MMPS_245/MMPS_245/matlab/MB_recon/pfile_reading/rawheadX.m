function [header,rhuser,byteheader] = rawheadX(fname,displayinfo);
%function [header,rhuser,byteheader] = rawheadX(fname,displayinfo);
%
%	Function returns header information from the given P-file.
%	the header information is primarily returned in the structure
%	header, although the rhuser variables are returned as an array
%	in rhuser.
%	See 'help geX' for more information.
%
%	INPUT:
%		fname   = path and file name.
%		displayinfo = 1 to display information as it is read.
%
%	OUTPUT:
%		header = header structure with certain fields.
%		rhuser = array of rhuser variables
%		byteheader = full header in bytes.
%
%       SEE ALSO:  rawloadX, writepfileX, reconX
%
%	B.Hargreaves -- April 2008.
%       minor mods 7/21/08  Tom Brosnan
%
%       now handles 12.X, 14.X, and 20.X Pfiles
%
%

% CVS Version:  (Do NOT Edit Line Below!)
% $Id: rawheadX.m,v 1.11 2009-06-12 14:49:42 brian Exp $
%

if (nargin < 2) displayinfo=0; end;

% -- Check if file exists, and error if not.
if (~exist(fname))
	tt = sprintf('%s:  file not found.',fname);
	error(tt);
end;

%
%

% -- Open file.
if exist('octave_config_info', 'builtin')
    % Octave allows us to use gzip-mode, in case the file is compressed.
    fip = fopen(fname,'rb','l');
else
    fip = fopen(fname,'r','l');
end
if fip == -1
  tt = sprintf('File %s not found\n',fname);
  error(tt);
end

header = struct('filename',fname);

%
% ===========================================================
% Extract Version.  It seems like 7 or less is before ESE 9.1
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
ver = fread(fip,1,'float');
%if (ver ~= 12) error('Only version 12 currently supported'); end;
if (ver < 11) error('Versions prior to 11 not supported'); end;


% -- Each parameter in the header has a name (that will be
% -- used in the structure) a bytes offset (depending on version)
% -- and a type.  This is an organized way of telling this program
% -- what to read:

[fieldnames,offsets,types] = getoffsetsX(ver);

%
for k=1:length(fieldnames)
	header = readandadd(fip,header,offsets(k),fieldnames{k},types{k});
end;


% -- Number of coils is actually not explicitly stored, so
% -- we just add it manually because it is useful.

ncoils = header.endrcvr - header.startrcvr + 1;
header = add2struct(header,'ncoils',ncoils);


% -- Read rhuser variables.
% -- Note that these are not contiguously stored (0-19, then 20-39 are.)
fseek(fip,216,-1);
rh019 = fread(fip,20,'float');
fseek(fip,1000,-1);
rh2048 = fread(fip,29,'float');
rhuser = [rh019(:); rh2048(:)];

% -- If requested, read full header in bytes.
%
if (nargout > 2)
  frewind(fip);
  [byteheader,nbytes] = fread(fip,header.rawhdrsize,'uint8');
  if (nbytes ~= header.rawhdrsize)
    tt=sprintf('Byte-header warning %d of %d bytes read.',nbytes,header.rawhdrsize); disp(tt);
  end;
end;

fclose(fip);

% -- Just dump out the full header information.
%
if (displayinfo == 1)
	header
end;


return;


%------------------------------------------------------------
function newstruct = readandadd(fid,oldstruct,offset,pname,ptype)
%
%	Internal function to read a header variable, then
%	add it to the structure.  This may be a bit slow, but
%	is designed to be simple to program.
%
%	INPUT:
%		fid = file pointer.
%		oldstruct = structure we are adding to.
%		offset = byte offset of parameter.
%		pname = name of parameter to use in structure (arbitrary).
%		ptype = type (int32, etc)

fseek(fid,offset,-1);
param = fread(fid,1,ptype);
newstruct = add2struct(oldstruct,pname,param);
return;


%------------------------------------------------------------
function newstruct = add2struct(oldstruct,newfield,newval)
%
%	Internal function to add a field to a structure.
%
%	INPUT:
%		oldstruct = structure we are adding to.
%		newfield = name of new field to add.
%		newval = value of field.
%
svals = struct2cell(oldstruct);
sfields=fieldnames(oldstruct);

svals = {svals{:} newval};
sfields = {sfields{:} newfield};
newstruct = cell2struct(svals,sfields,2);



