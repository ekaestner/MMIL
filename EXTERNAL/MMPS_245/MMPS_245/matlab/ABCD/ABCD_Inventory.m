function ABCD_Inventory(ProjID,varargin)
%function ABCD_Inventory(ProjID,[options])
%
% Purpose:
%   inventory of data and output from incoming through processing and QC
%
% Usage:
%  ABCD_Inventory(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%
% Optional Parameters:
%   'outdir': output directory
%     if empty, will place in /home/{user}/MetaData/{ProjID}
%     {default = []}
%   'outstem': output file stem
%     {default = 'inventory'}
%   'rootdir': root directory containing project directory
%     with incoming, unpack, orig, raw, proc, proc_dti, proc_bold dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'sites': cell array of site names
%     if empty, use all site dirs within incoming dir
%     {default = []}
%   'ContainerTypes': cell array of container types to check
%     {default = {'orig','pc','raw' 'proc' 'proc_dti' 'proc_bold'}}
%   'verbose': [0|1] display status messages and warnings
%     {default: 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  12/19/16 by Don Hagler
% Last Mod: 01/15/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(ProjID,varargin);

% check whether output file exists
if ~exist(parms.fname_out,'file') || parms.forceflag

  % compile struct with info about files in incoming
  inv_info = check_incoming(parms);

  % check files in unpack and orig
  inv_info = check_unpack(inv_info,parms);

  % check downstream containers
  inv_info = check_containers(inv_info,parms);

  % write summary file
  mmil_struct2csv(inv_info,parms.fname_out);

elseif parms.verbose
  fprintf('%s: output file %s exists... quitting\n',...
    mfilename,parms.fname_out);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'ProjID',ProjID,[],...
  ...
    'outdir',[],[],...
    'outstem','inventory',[],...
    'rootdir','/space/syn05/1/data/MMILDB',[],...
    'sites',[],[],...
    'ContainerTypes',{'orig','pc','raw' 'proc' 'proc_dti' 'proc_bold'},...
                     {'orig','pc','raw' 'proc' 'proc_dti' 'proc_bold' 'fsurf' 'fsico'},...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'info_tags',{'json_exists','json_readable','json_dicom',...
                 'SubjID','pGUID','EventName','SessionType',...
                 'StudyInstanceUID','SeriesInstanceUID',...
                 'StudyDate','SeriesDescription',...
                 'VisitID'},[],...
    'unpack_tags',{'unpacked','reconned','renamed',...
                   'unpack_error','recon_error','rename_error'},[],...
    'qc_types',{'proc' 'proc_dti' 'proc_bold'},[],...
    'qc_tags',{'proc_aQC','proc_dti_aQC','proc_bold_aQC'},[],...
  });

  parms.info_tags = cat(2,parms.info_tags,parms.unpack_tags);
  parms.info_tags = cat(2,parms.info_tags,parms.ContainerTypes);
  parms.info_tags = cat(2,parms.info_tags,parms.qc_tags);
  
  % set outdir and outstem and fname_out
  if isempty(parms.outdir)
    parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
  end;
  mmil_mkdir(parms.outdir);
  parms.outstem = sprintf('%s_%s',ProjID,parms.outstem);
  parms.fname_out = sprintf('%s/%s.csv',parms.outdir,parms.outstem);

  % set RootDirs
  parms.RootDirs = [];
  parms.RootDirs.incoming = sprintf('%s/%s/incoming',parms.rootdir,ProjID);
  parms.RootDirs.unpack = sprintf('%s/%s/unpack',parms.rootdir,ProjID);
  parms.RootDirs.orig = sprintf('%s/%s/orig',parms.rootdir,ProjID);
  parms.ntypes = length(parms.ContainerTypes);
  for i=1:parms.ntypes
    ContainerType = parms.ContainerTypes{i};
    parms.RootDirs.(ContainerType) = sprintf('%s/%s/%s',parms.rootdir,ProjID,ContainerType);
  end;
  parms.RootDirs = MMIL_Check_RootDirs(parms.RootDirs);

  % get list of sites
  if isempty(parms.sites)
    dlist = dir(sprintf('%s/*',parms.RootDirs.unpack));
    parms.sites = setdiff({dlist.name},{'.','..'});
  end;% check containers for each site

  parms.nsites = length(parms.sites);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inv_info = check_incoming(parms)
  inv_info = [];
  for s=1:parms.nsites
    site = parms.sites{s};
    if parms.verbose
      fprintf('%s: taking inventory of incoming for %s...\n',mfilename,site);
    end;
    tmp_info = check_incoming_site(site,parms);
    inv_info = cat(1,inv_info,tmp_info);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inv_info = check_incoming_site(site,parms)
  inv_info = [];
  sitepath = sprintf('%s/%s',parms.RootDirs.incoming,site);
  flist = dir(sprintf('%s/*.tgz',sitepath));
  for i=1:length(flist)
    fname = sprintf('%s/%s',sitepath,flist(i).name);
    [fpath,fstem] = fileparts(fname);
    iinfo = init_info(fstem,site,parms.info_tags);
    % check json file exists      
    fname_json = sprintf('%s/%s.json',sitepath,fstem);
    if exist(fname_json,'file')
      iinfo.json_exists = 1;
      % check json file is readable
      try
        jinfo = loadjson(fname_json);
        iinfo.json_readable = 1;
      catch ME
        iinfo.json_readable = 0;
      end
      if iinfo.json_readable
        % check if tgz contains dicoms with expected info
        if isfield(jinfo,'SeriesDescription') && isfield(jinfo,'SeriesInstanceUID')
          [iinfo.SubjID,iinfo.pGUID,iinfo.EventName,iinfo.SessionType] =...
            abcd_get_SubjID(fstem,jinfo);
          iinfo.StudyInstanceUID = jinfo.StudyInstanceUID;
          iinfo.SeriesInstanceUID = jinfo.SeriesInstanceUID;
          iinfo.StudyDate = jinfo.StudyDate;
          iinfo.SeriesDescription = jinfo.SeriesDescription;
          iinfo.json_dicom = 1;
        else
          iinfo.json_dicom = 0;
        end;
      end;
    else
      iinfo.json_exists = 0;
    end;
    inv_info = cat(1,inv_info,iinfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iinfo = init_info(fstem,site,tags)
  iinfo.site = site;
  iinfo.fstem = fstem;
  for i=1:length(tags)
    iinfo.(tags{i}) = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inv_info = check_unpack(inv_info,parms)
  if parms.verbose
    fprintf('%s: checking unpack dirs...\n',mfilename);
  end;
  for i=1:length(inv_info)
    % skip if for a non-dicom tgz (i.e., P-file)
    dicom_flag = inv_info(i).json_dicom;
    if isempty(dicom_flag) || ~dicom_flag, continue; end;

    sitepath = sprintf('%s/%s',parms.RootDirs.unpack,inv_info(i).site);
    fstem = inv_info(i).fstem;

    % check whether unpacked, reconned, renamed, errors
    for j=1:length(parms.unpack_tags)
      tag = parms.unpack_tags{j};
      fname_check = sprintf('%s/%s/%s.%s',sitepath,fstem,fstem,tag);
      if exist(fname_check,'file')
        inv_info(i).(tag) = 1;
      else      
        inv_info(i).(tag) = 0;
      end;
    end;

    % get VisitID from renamed file
    fname_rename = sprintf('%s/%s/%s.renamed',sitepath,fstem,fstem);
    if exist(fname_rename,'file')
      fid = fopen(fname_rename);
      origsubpath = fscanf(fid,'%s');
      fclose(fid);
      origpath = fileparts(origsubpath);
      [origrootpath,origdir] = fileparts(origpath);
      inv_info(i).VisitID = origdir;
    else
      continue;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inv_info = check_containers(inv_info,parms)
  if parms.verbose
    fprintf('%s: checking containers...\n',mfilename);
  end;
  for i=1:length(inv_info)
    VisitID = inv_info(i).VisitID;
    % skip if VisitID is missing
    if isempty(VisitID), continue; end;
    % check existence of each container
    for j=1:parms.ntypes
      ContainerType = parms.ContainerTypes{j};
      switch ContainerType
        case 'orig'
          ContainerPath = sprintf('%s/%s',parms.RootDirs.orig,VisitID);
          if exist(ContainerPath,'dir')
            inv_info(i).orig = 1;
          else
            inv_info(i).orig = 0;
          end;
        case 'pc'
          ContainerPath = sprintf('%s/%s',parms.RootDirs.pc,VisitID);          
          inv_info(i).pc = 0;
          % check for autoPC json files
          if exist(ContainerPath,'dir')
            dlist = dir(sprintf('%s/*',ContainerPath));
            dlist = dlist(find(~ismember({dlist.name},{'.','..'})));
            for k=1:length(dlist)
              flist = dir(sprintf('%s/%s/QC_PC.json',...
                          ContainerPath,dlist(k).name));
              if ~isempty(flist)
                inv_info(i).pc = 1;
                break;
              end;
            end;
          end;
        otherwise
          ContainerPath = ...
            MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
          if ~isempty(ContainerPath)
            inv_info(i).(ContainerType) = 1;
            % check for autoQC mat files
            if ismember(ContainerType,parms.qc_types)
              flist = dir(sprintf('%s/*auto_qcinfo.mat',ContainerPath));
              if ~isempty(flist)
                inv_info(i).([ContainerType '_aQC']) = 1;
              else
                inv_info(i).([ContainerType '_aQC']) = 0;
              end;
            end;
          else
            inv_info(i).(ContainerType) = 0;
          end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

