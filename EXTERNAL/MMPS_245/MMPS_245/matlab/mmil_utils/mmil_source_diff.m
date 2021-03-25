function varargout = mmil_source_diff(package, srcroot, destroot)
%function varargout = mmil_source_diff(package, [srcroot], [destroot])
%
% Purpose: 
%   To show the difference in source files for code management
%   This function recurses over all subdirectories of the source root.
%
% Parameters:
%   package - package name to analyze
%   srcroot - root of your new source
%             {default: ~/matlab/[package]}
%   destroot - root of where you want to copy your new source to
%             {default: ~mmildev/matlab/[package]}
%
% Output: 
%   - Status message is printed to the screen about differences
%   - An array of output objects is returned, each containing:
%      - a status message
%      - the result of a UNIX "diff" between source and dest files.
%
%
% Created:  08/13/07 by Ben Cipollini 
% Last Mod: 03/15/11 by Don Hagler
%

  if (~exist('package','var'))
    help(mfilename);
    return;
  end;

  % default src root to user's timesurfer location
  if (~exist('srcroot','var'))
    srcroot=fullfile('~','matlab',package);
  end;
  
  % default dest root to mmildev's timesurfer location
  if (~exist('destroot','var'))
    destroot=fullfile('~mmildev', 'matlab', package);
  end;
  
  % find source files & directories 
  %mmil_logstr('Checking %s', srcroot);
  
  srcall = dir(srcroot);
  srcdirs = srcall(find([srcall.isdir]==1));
  srcfiles = srcall(find([srcall.isdir]==0));
  diff = [];
  
  
  % diff the source & dest files
  for i=1:length(srcfiles)
  
    srcfile=srcfiles(i);
    srcfilepath = fullfile(srcroot, srcfile.name);
	
    destfilepath = fullfile(destroot, srcfile.name);
    destfile =  dir(destfilepath);

    %  skip files we don't care about
  	if (~isempty(findstr('~',srcfile.name))...
	    || (srcfile.name(1) == '.') ...
	  	|| ~isempty(findstr('.bak',srcfile.name)) ...
		  || ~isempty(findstr('.old.',srcfile.name)))
  	  continue;
	  
	  % in src but not dest
    elseif (isempty(destfile))
      diff(end+1).path = srcfilepath;
  	  diff(end).msg	= sprintf('not in dest: %s', srcfilepath);
      diff(end).diff = '[full file]';
      	  
    % different file timestamps
   	elseif ((isfield(srcfile, 'date') && datenum(srcfile.date) < datenum(destfile.date)) ...
           || (isfield(srcfile, 'datenum') && srcfile.datenum < destfile.datenum))
      diff(end+1).path = srcfilepath;
      diff(end).msg = sprintf('** DEST IS NEWER: %s', srcfilepath);
  	  [s,diff(end).diff] = mmil_unix(sprintf('diff %s %s', srcfilepath,destfilepath));
	  
	  
    % different file sizes
  	elseif (srcfile.bytes ~= destfile.bytes)
      diff(end+1).path = srcfilepath;
  	  diff(end).msg = sprintf('diff in src: %s', srcfilepath);
      [s,diff(end).diff] = mmil_unix(sprintf('diff %s %s', srcfilepath,destfilepath)); 
	  
	% No difference
	else
	  continue;
    end;
    
    % Log a message to the screen
    mmil_logstr(diff(end).msg);
  end ;

  % Make recursive calls for all subdirectories
  for i=1:length(srcdirs)
    srcdir = srcdirs(i);

    % Skip system directories
  	if (strcmp(srcdir.name, '.') ...
	     || strcmp(srcdir.name,'..'))
  	  continue;
	  
	% Skip directories with special prefixes
  	elseif (srcdir.name(1) == '_' ...
	      || ~isempty(findstr('old.',srcdir.name)))
      continue;
    end;
	
	%
  	diff2 = mmil_source_diff(package, ...
                             fullfile(srcroot, srcdir.name), ...
	                         fullfile(destroot, srcdir.name));

    diff  = [diff diff2];
  end;

  if (nargout==1)
    varargout{1} = diff;
  end;
