function [M_v1_to_v2, min_cost] = rbreg_vol2vol_mi_mask_djh(vol1, vol2, volmask, varargin)
% Multiscale Rigid body registration by mutual information based on maskvolume
%       
% [M_v1_to_v2, min_cost] = rbreg_vol2vol_mi_mask_djh(vol1, vol2, volmask, [bdispiter], [mstep], ...
%                                             [[scales], [tmin], [amin], [allowchanged])
%
% Input:
%   vol1: template vol structure
%   vol2: regsitering vol structure
%   volmask: Mask volume
%   bdispiter: show ino ineach iteration
%   mstep:      sampling rate in volmask (default =1)   
%   scales: number of scales (default:[0 83 49 27 16 9 5 3 2 1]). Be sure
%   put 0 in the first scale to initialize algorithm
%   tmin: minimun in translation (in mm) default=0.05;
%   amin: minimum in angle (in rad) default= 0.05 degree
%   allowchanged: allowable parameter changes default=2;
%
% Output:
%   M_v1_to_v2: Registration matrix from Vol1 to Vol2
%   min_const:  Minimization error
%
% Last Mod: 07/02/08 by Don Hagler
%


bdispiter=false;
if nargin >= 4
  bdispiter = varargin{1};
end

mstep =1;
if nargin>=5
    mstep=varargin{2};
end

scales=[0 83 49 27 16 9 5 3 2 1];
if nargin >= 6
  scales = varargin{3};
end

tmin = 0.05;
if nargin >= 7
  tmin = varargin{4};
end

amin = (0.05)*(pi/180);
if nargin >= 8
  amin = varargin{5};
end

allowchanged=2;
if nargin >= 9
  allowchanged = varargin{6};
end


M_reg = eye(4,4);
M_reg_opt = M_reg;
min_cost = 1e10;
for scale = scales
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

      if 0
        if 0
%% todo: test thresholds
          minval1 = 1;
          maxval1 = 10^10;
          minval2 = 1;
          maxval2 = 10^10;
        else
          minval1 = 0;
          maxval1 = 10^10;
          minval2 = 0;
          maxval2 = 10^10;
        end;

		    err=vols_getMI_mask(vol1, vol2, volmask, eye(4,4), M_reg, mstep,...
          minval1,maxval1,...
          minval2,maxval2);
      else
        err=vols_getMI_mask(vol1, vol2, volmask, eye(4,4), M_reg, mstep);
      end;
		  cost = -err;
		  
		  str = 'scale=%d (%d) [%d %d %d %d %d %d] cost=%f min_cost=%f\n';
		  if (bdispiter)
		    fprintf(str,scale,pass,txi,tyi,tzi,axi,ayi,azi,cost,min_cost);
		  end  
          if (cost<min_cost)
              [mval maxldir] = max(abs(M_reg(:,1)));
              [mval maxpdir] = max(abs(M_reg(:,2)));
              [mval maxhdir] = max(abs(M_reg(:,3)));
              if (maxldir == 1) & (maxpdir == 2) & (maxhdir == 3)
                  min_cost = cost;
                  M_reg_opt = M_reg;
                  str = '*** scale=%d (%d) [%d %d %d %d %d %d]';
                  str = [str 'cost=%f min_cost=%f\n'];
                  if (bdispiter)
                      fprintf(str,scale,pass,txi,tyi,tzi,axi,ayi,azi,cost,min_cost);
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
M_v1_to_v2 = M_reg_opt;
