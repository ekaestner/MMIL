function status = rc_calc_loose_svd(prefix,time0,time1);
%function status = rc_calc_loose_svd(prefix,time0,time1);
%
% Required Input:
%  prefix: RCSE prefix
%
% Optional Input:
%  time0 - start of time range (msec)
%    {default: 80}
%  time1 - end of time range (msec)
%    {default: 90}
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler 
%

status = 0;
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('time0','var'), time0 = 70; end;
if ~exist('time1','var'), time1 = 150; end;

% load matfiles
matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
if ~parms.loose_flag
  error('this function is only applicable when loose_flag=1');
end;
matfile = sprintf('matfiles/%s_ret_mapping.mat',prefix);
load(matfile);
matname=sprintf('matfiles/%s_results.mat',parms.prefix);
load(matname);
matfile = sprintf('matfiles/%s_dip_info.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

[num_tpoints num_sources] = size(S);
num_areas = retmap.num_areas;
num_locs = retmap.num_locs;
M = retmap.M;

sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((time0 - t_trigger)*sfreq/1000);
t1 = round((time1 - t_trigger)*sfreq/1000);
if t0 > t1, tmp=t0; t0=t1; t1=tmp; end;
if t0 < 1, t0 = 1; end;
if t1 > num_tpoints, t1 = num_tpoints; end;

svd_weights = zeros(3,num_areas*num_locs);
k=0; s=1;
for a=1:num_areas
  for theta=1:num_locs
    % get time course of normal and tangential components
    j = k+1; k = j+2;
    X = S(t0:t1,j:k);
    [U,SS,V] = svd(X,0);
    svd_weights(:,s) = diag(SS);
    s=s+1;
    dip = V(:,1); % is this correct?
    v_lh = M(a,theta).v_lh;
    v_rh = M(a,theta).v_rh;
    lh_dip_info(4:6,v_lh) = dip*ones(1,length(v_lh));
    rh_dip_info(4:6,v_rh) = dip*ones(1,length(v_rh));
  end;
end;

matfile = sprintf('matfiles/%s_dip_info_svd.mat',prefix);
save(matfile,'lh_dip_info','rh_dip_info','lh_dec_dips','rh_dec_dips','svd_weights');

status = 1;
return;
