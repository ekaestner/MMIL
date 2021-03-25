function ts_write_mscript(parms,varargin)
% Last Mod: 09/15/12 by Don Hagler

inparms = mmil_args2parms(varargin, { ...
  'function','',[],...
  'cmdstr','',[],...
  'outfile','',[],...
  'outpath','',[],...
  'fileid','',[],...
  'filename','',[],...
  'hdr_flag',0,{0,1},...
});

if ~isstruct(parms), return; end
if isnumeric(inparms.fileid), inparms.fileid = num2str(inparms.fileid); end

if (~isfield(parms,'indata')  || isempty(parms.indata))  && (isfield(parms,'itype') && ~isempty(parms.itype))
  parms.indata = parms.itype;
end
if (~isfield(parms,'outdata') || isempty(parms.outdata)) && (isfield(parms,'otype') && ~isempty(parms.otype))
  parms.outdata = parms.otype;
end
  
try
  timesurfer_flag = isfield(parms,'datafile') && ~isempty(parms.datafile) && ...
                  any(parms.indata) && any(parms.outdata); 
catch
  timesurfer_flag = 0;
end

if ~isempty(inparms.outfile)
  [pathstr name ext] = fileparts(inparms.outfile);
  outfile = inparms.outfile;
  if ~exist(pathstr,'dir'),
      fprintf('making directory: %s\n',pathstr);
      unix(['mkdir -p ' pathstr]);
  end      
else
  if ~isempty(inparms.outpath)
    outfile = inparms.outpath;
    if ~exist(outfile,'dir'),
        fprintf('making directory: %s\n',outfile);
        unix(['mkdir -p ' outfile]);
    end            
  elseif isfield(parms,'rootoutdir')
    outfile = fullfile(parms.rootoutdir,'scripts');
    if ~exist(outfile,'dir'),
        fprintf('making scripts directory: %s\n',outfile);
        unix(['mkdir -p ' outfile]);
    end        
  else
    outfile = pwd;
  end
	if isfield(parms,'prefix')
		prefix = parms.prefix;
	else
		prefix = 'run';
	end
  if ~isempty(inparms.function)
    outfile = sprintf('%s/%s_%s%s.m',outfile,prefix,inparms.function,inparms.fileid);
  else
    outfile = sprintf('%s/%s_auto_options.m',outfile,prefix);
  end    
end
try parms = rmfield(parms,'previous'); end
fname = fieldnames(parms);

fprintf('writing matlab script: %s\n',outfile);
fid = fopen(outfile,'wt'); 
fprintf(fid,'tic\n');

if isfield(parms,'mmilclusterheader')
  fprintf(fid,'%s\n',parms.mmilclusterheader);
end

% if timesurfer_flag
if any(parms.indata) && ~isempty(parms.datafile)
  if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
  if inparms.hdr_flag
    for i = 1:length(parms.datafile)
      fprintf(fid,'%s.parms.filename{%d} = ''%s'';\n',parms.indata,i,parms.datafile{i});
    end
  else
    for i = 1:length(parms.datafile)
      fprintf(fid,'datafile{%d} = ''%s'';\n',i,parms.datafile{i});
    end
    fprintf(fid,'for i=1:length(datafile)\n');
    fprintf(fid,'\t%s(i) = getfield(load(datafile{i},''%s''),''%s'');\n',parms.indata,parms.indata,parms.indata);
    fprintf(fid,'end\n');
    fprintf(fid,'if i>1, %s = ts_combine_data(%s); end\n\n',parms.indata,parms.indata);
  end
end

for i = 1:length(fname)
    str = fname{i};
    val = parms.(str);
    valstr = '';
    if isnumeric(val)
        if length(val)==1
            valstr = num2str(val);
        else
            valstr = numarray2str(val);
        end
    elseif iscell(val)
        valstr = '';
        for j = 1:length(val)
            if isstr(val{j})
                if j==1
                    valstr=sprintf('''%s''',val{j}); 
                else
                    valstr = sprintf('%s,''%s''',valstr,val{j});
                end
            elseif isnumeric(val{j})
                if j==1, valstr=numarray2str(val{j}); else valstr=[valstr ',' numarray2str(val{j})]; end
            end
        end
        valstr = sprintf('{%s}',valstr);
    elseif isstruct(val)
        valstr = str;
    elseif isstr(val)
        valstr = sprintf('''%s''',val);
    end
    fprintf(fid,'parms.%s = %s;\n',str,valstr);
end

if timesurfer_flag
  fprintf(fid,'\nargs = mmil_parms2args(parms);\n');
  fprintf(fid,'%s = %s(%s,args{:});\n',parms.outdata,parms.function,parms.indata);
elseif ~isempty(inparms.function)
	if ~isempty(inparms.cmdstr)
		fprintf(fid,'%s\n',inparms.cmdstr);    
  else
    fprintf(fid,'\nargs = mmil_parms2args(parms);\n');    
  	fprintf(fid,'%s(args{:})\n',inparms.function);
	end
end
if any(parms.outdata)
  if any(parms.indata), fprintf(fid,'try parms.previous = %s.parms; end\n',parms.indata); end
  fprintf(fid,'if ~isfield(%s,''parms''), %s.parms = parms; end\n',parms.outdata,parms.outdata);  
  fprintf(fid,'save(''%s'',''%s'',''-v7.3'');\n',inparms.filename,parms.outdata);
end
fprintf(fid,'\ntoc\nexit\n');  % may need for cluster computing
fclose(fid);


function str = numarray2str(array)
str = '';
for i = 1:length(array)
    str = sprintf('%s %g',str,array(i));
end
str = sprintf('[%s]',str);

