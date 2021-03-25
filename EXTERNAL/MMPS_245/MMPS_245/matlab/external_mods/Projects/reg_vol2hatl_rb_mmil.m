function [M_atl_to_vol, min_cost, sf_opt]=reg_vol2hatl_rb(vol, volm, volstd, sampling, varargin)
%Multiscale Rigid body registration to atlas
%       
% [M_atl_to_vol, min_cost,sf_opt]=reg_vol2hatl_rb(vol, volm, volstd, sampling, 
%                                                 [bdispiter],[bsmooth], [scales],
%                                                 [tmin], [amin], [allowchanged],
%                                                 [M_reg], [thresh])
%
% Input:
%   vol: vol structure
%   volm: atlas mean volume
%   volstd: atlas std volume
%   sampling: sampling rate in each dim.
%   bdispiter: show info on each iteration
%     {default=false}
%   scales: multiple scales of translation or rotation
%     {default:[0 83 49 27 16 9 5 3 2 1]}
%     Be sure to put 0 in the first scale to initialize algorithm
%   bsmooth: smooth input volume
%     {default=true}
%   tmin: minimum translation (in mm)
%     {default=0.5}
%   amin: minimum in angle (in rad)
%     {default=0.5 degrees}
%   allowchanged: number of allowable parameter changes
%     {default=2}
%   M_reg: Initial Registration Matrix
%     {default=eye(4,4)}
%   thresh: values in vol below threshold will be ignored
%     important if volume has been resampled and includes zeros
%     threshold is applied after smoothing and interpolation
%     {default=20}
%
% Output:
%   M_atl_to_vol: registration matrix
%   min_cost: Minimization error
%   sf_opt: estimated Scaling factor
%
%
%  Early Mod: 04/21/08 by Don Hagler (added threshold)
%  Last Mod:  09/23/11 by Don Hagler (initialized sf_opt)
%

min_fract_points = 0.1; % need to have at least 10% of the sampled points
                        % to have an actual value in the input volume
                        % (i.e. after thresh)

dim =size(volm.imgs);
ind =find (volm.imgs>0);
inds = ind(1:sampling:length(ind));
[I J K] = ind2sub(dim, inds);
ni = length(inds);
vxl = ones(4, ni);
vxl(1,:)=I(:)';
vxl(2,:)=J(:)';
vxl(3,:)=K(:)';
lphpos = (volm.Mvxl2lph*vxl)';
mu=volm.imgs(inds);
sigma=volstd.imgs(inds);
clear volm volstd;

min_num_points = min_fract_points*size(lphpos,1);

bdispiter=false;
if nargin >= 5
  bdispiter = varargin{1};
end

bsmooth=true;
if nargin >=6
  bsmooth = varargin{2};
end

scales=[0 83 49 27 16 9 5 3 2 1];
if nargin >= 7
  scales = varargin{3};
end

tmin = 0.5;
if nargin >= 8
  tmin = varargin{4};
end

amin = (0.5)*(pi/180);
if nargin >= 9
  amin = varargin{5};
end

allowchanged=2;
if nargin >= 10
  allowchanged = varargin{6};
end

M_reg = eye(4,4);
if nargin >= 11
  M_reg = varargin{7};
end

thresh = 20;
if nargin >= 12
  thresh = varargin{8};
end

if (bsmooth)
  %disp    'Smoothing volume...'
  vols=vol_filter(vol, 1);
  clear vol;
else
  vols=vol;
  clear vol;
end

M_reg_opt = M_reg;
min_cost = 1e10;
sf=1;
sf_opt = sf;
si=0.;
lsi=length(scales);
for scale = scales
  si=si+1;
  if scale==0
    win = 0;
  else
    win = 1;
  end
  changed = 1;
  pass = 0;
  while changed
    pass = pass+1;
    changed = 0;
    M_reg_bak = M_reg_opt;
    
    for txi = -win:win
    for tyi = -win:win
    for tzi = -win:win
    for axi = -win:win
  	for ayi = -win:win
  	for azi = -win:win
		  if (sum([txi tyi tzi axi ayi azi]~=0)<=allowchanged)
        tx = txi*scale*tmin; ty = tyi*scale*tmin; tz = tzi*scale*tmin;
        ax = axi*scale*amin; ay = ayi*scale*amin; az = azi*scale*amin;
        M_reg = Mrotz(az)*Mroty(ay)*Mrotx(ax)*Mtrans(tx,ty,tz)*M_reg_bak;
        [vxlval inbound]= vol_getvxlsval(lphpos, vols, M_reg);
%        ind = find((inbound>0));
        ind = find(inbound>0 & vxlval>thresh);
        if length(ind)<min_num_points
          continue;
        end;
        sum_numer=sum((vxlval(ind).*mu(ind))./(sigma(ind).*sigma(ind)));
        sum_denom=sum((vxlval(ind).*vxlval(ind))./(sigma(ind).*sigma(ind)));
        tv=(sf*vxlval(ind)-mu(ind))./sigma(ind);
        sf=sum_numer/sum_denom;
        cost = tv'*tv/length(tv);
        str = 'scale=%d (%d) [%d %d %d %d %d %d ] sf=%f cost=%f min_cost=%f\n';
        if (bdispiter)
          fprintf(str,scale,pass,txi,tyi,tzi,axi,ayi,azi,sf, cost,min_cost);
        end  
        if (cost<min_cost)
          [mval maxldir] = max(abs(M_reg(:,1)));
          [mval maxpdir] = max(abs(M_reg(:,2)));
          [mval maxhdir] = max(abs(M_reg(:,3)));
          if (maxldir == 1) & (maxpdir == 2) & (maxhdir == 3)
            min_cost = cost;
            M_reg_opt = M_reg;
            sf_opt=sf;
            str = '*** scale=%d (%d) [%d %d %d %d %d %d]';
            str = [str 'sf=%f cost=%f min_cost=%f\n'];
            if (bdispiter)
              fprintf(str,scale,pass,txi,tyi,tzi,axi,ayi,azi,sf_opt, cost,min_cost);
            end
            changed = 1;
          end
        end
		  end
    end
  end
  end
  end
  end
  end
  end
end
M_atl_to_vol = M_reg_opt;



