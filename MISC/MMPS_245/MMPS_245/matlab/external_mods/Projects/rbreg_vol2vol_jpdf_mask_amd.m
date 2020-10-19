function [M_v1_to_v2, min_cost, jpdf, bins1, bins2] = rbreg_vol2vol_jpdf_mask_amd(vol1, vol2,...
  volmask, Mreg0,bdispiter,mstep,scales,tmin,amin,allowchanged,jpdf,bins1,bins2,interpm)
%function [M_v1_to_v2, min_cost, jpdf, bins1, bins2] = rbreg_vol2vol_jpdf_mask(vol1, vol2,...
% volmask, Mreg0, [bdispiter],[mstep],[scales],[tmin],[amin],[allowchanged],[jpdf],...
% [bins1],[ bins2],[interpm])
%
% Multiscale Rigid body registration By Joint PDF based on maskvolume
%
% Input:
%   vol1: template vol structure
%   vol2: regsitering vol structure
%   volmask: Mask volume structure
%   Mreg0: initial registration matrix
%   bdispiter: show ino ineach iteration
%   mstep:      sampling rate in volmask (default =1)   
%   scales: number of scales (default:[0 83 49 27 16 9 5 3 2 1]).
%    (be sure to put 0 in the first scale to initialize algorithm)
%   tmin: minimun in translation (in mm) default=0.05;
%   amin: minimum in angle (in rad) default= 0.05 degree
%   allowchanged: allowable parameter changes default=2;
%   jpdf: precalculated joint probability distribution function
%   bins1: ?
%   bins2: ?
%   interpm :   0: Nearest  1:Linear(default) 2:cubic
%               3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%
% Output:
%   M_v1_to_v2: Registration matrix from Vol1 to Vol2
%   min_cost:  Minimization error

if ~exist('Mreg0','var') | isempty(Mreg0)
  Mreg0=eye(4);
end

if ~exist('bdispiter','var') | isempty(bdispiter)
  bdispiter=false;
end

if ~exist('mstep','var') | isempty(mstep)
  mstep =1;
end

if ~exist('interpm','var') | isempty(interpm)
  interpm=1;
end

if ~exist('scales','var') | isempty(scales)
  scales=[0 83 49 27 16 9 5 3 2 1];
end

if ~exist('tmin','var') | isempty(tmin)
  tmin = 0.05;
end

if ~exist('amin','var') | isempty(amin)
  amin = (0.05)*(pi/180);
end

if ~exist('allowchanged','var') | isempty(allowchanged)
  allowchanged=2;
end

if ~exist('jpdf','var') | isempty(jpdf)
  if bdispiter
    disp 'Calculating Joint PDF'
  end
  [sum_log10_val, jpdf, bins1, bins2, jentropy]=vols_jhist_mask_amd(vol1, vol2, volmask, 1, eye(4), Mreg0); % This should be pre-computed for MR/FDG-PET
  jpdf = max(0,smooth2(jpdf,11,11)); % Smooth jpdf function
else
  if bdispiter
    disp 'Using specified Joint PDF'
  end
  % Should determine scaling factors here?
end

min1=bins1(1);
dbr1=bins1(2)-bins1(1);
numbins1=length(bins1);
min2=bins2(1);
dbr2=bins2(2)-bins2(1);
numbins2=length(bins2);

M_reg = Mreg0;
%M_reg = eye(4);
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
		  err=jpdfregerr_mask(vol1, vol2, volmask, jpdf, numbins1, min1, dbr1, numbins2, min2, dbr2, M_reg, mstep,interpm);
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

