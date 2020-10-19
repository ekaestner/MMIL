function mmil_rm_all_empty_files(rootdir,dirstr,filestr)
%function mmil_rm_all_empty_files(rootdir,dirstr,filestr)
%
% Purpose: batch removal of empty files
%
% Required Input:
%   rootdir: root directory
%   dirstr: string to match subdirectories
%     e.g. 'MRIRAW*', 'MRIPROC*'
%
% Optional Input:
%   filestr: string to match files
%     {default = '*'}
%
% Example Usage: mmil_rm_all_empty_files(...
%   '/space/invaders/1/data/MMILDB/ADNI/Containers','MRIPROC*','*');
%
% Created:  08/02/07 by Anders Dale
% Last Mod: 09/09/12 by Don Hagler
%

%% todo: use varargin
%% todo: accept subdir
%% todo: accept test flag


if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('filestr','var') || isempty(filestr), filestr = '*'; end

dirlist = dir(sprintf('%s/%s',rootdir,dirstr));
for i = 1:length(dirlist)
  if dirlist(i).isdir
    dirname = sprintf('%s/%s',rootdir,dirlist(i).name);
    fprintf(1,'processing %s\n',dirname);
    filelist = dir(sprintf('%s/%s',dirname,filestr));
    for j = 1:length(filelist)
      if ~filelist(j).isdir
        if filelist(j).bytes == 0
          fname = sprintf('%s/%s',dirname,filelist(j).name);
          cmd = sprintf('rm %s',fname);
          disp(cmd)
          [s,w] = system(cmd);
          if s~=0
            fprintf(1,'error:\n');
            disp(w);
          end
        end
      end
    end
  end
end
