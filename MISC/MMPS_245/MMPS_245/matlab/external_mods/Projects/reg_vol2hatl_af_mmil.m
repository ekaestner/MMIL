function [Maffine, fmi, bconverged, optpnt]=reg_vol2hatl_af(vol, volm, volstd, sampling, M_atl_to_vol_rb, bsmooth, errfn, initp, wei, maxiter, ftol, mindp, thresh)
% Affine registration to head atlas By using  Simplex
%
%       [M_vol1_to_vol2_af, fmi, optpnt]=
%           reg_vol2hatl_af(vol, volm, volstd, sampling, M_volt_to_vol_rb, [bsmooth], [function_name], 
%                           , [initp], [wei], [maxiter], [ftol], [mindp], [thresh])
%   input:
%       vol: vol structure
%       volm: atlas mean volume
%       volstd: atlas std volume
%       Sampling: sampling rate in each dim.
%       M_atl_to_vol_rb: Rigid body registration matrix from atlas to vol
%       bsmooth  smooth volume (default =true);
%       lphpos : Sampling lph points
%       errfun: error function name (calmi_corr (default) or calmi_ss)
%       maxiter: max. iteration (default = 2000)
%       ftol: tolerance (default =0.001)
%       mindp: minimum distance between points (default = 0.01)
%       initp: Initial Point (default = [1 1 1 0 0.0 0.0])
%       wei: offset (default = [0.1 0.1 0.1 0.1 0.1 0.1])
%       thresh: values in vol below threshold will be ignored (default=20)
%
%   Output:
%       M_vol1_to_vol2_af: Affine transfrom from vol1 to vol2
%       fmi: minimum
%       optpnt: optimal point

dim =size(volm.imgs);
%sampling=[4 4 4];
[I J K] = ndgrid(1:sampling(1):dim(1), 1:sampling(2):dim(2), 1:sampling(3):dim(3));
vxl = [I(:) J(:) K(:) ones(size(I(:)))]';
lphpos = (volm.Mvxl2lph*vxl)';
ind = sub2ind(dim, I(:), J(:), K(:));
mu=volm.imgs(ind);
sigma=volstd.imgs(ind);
clear volm volstd;

if (~exist('bsmooth'))
    bsmooth=true;
end   

if (~exist('errfn'))
    errfn = 'calmi_hatl_mmil';
end

if ~exist('maxiter')
   maxiter = 2000;
end

if ~exist('ftol')
   ftol = 0.001;
end

if  ~exist('mindp')
   mindp = 0.01;
end

if  ~exist('thresh')
   thresh = 20;
end


dim = length(initp);
if (bsmooth)
  vol = vol_filter(vol, 1);
end

sf=1.;
coeff_reflect=1.;
coeff_expand=2.;
coeff_contract=0.5;
mi=zeros(dim+1,1);
pnts=zeros(dim+1,dim);
Mreg = M_atl_to_vol_rb;
bconverged =false;
for i=1:dim
    pnts(i,:)=initp;
    pnts(i,i)=initp(i)+wei(i);
end
pnts(dim+1,:)=initp;

for i=1:dim+1
  mi(i)=feval(errfn, vol, mu, sigma, lphpos, Mreg, pnts(i,:), sf, thresh);
end
    
i=0;
while(i<=maxiter)
    [Y, sind]=sort(mi);
    dp=zeros(dim+1,1);
    %Worst index
    wind=sind(dim+1);
    %Best index
    bind=sind(1);
    %2nd best index
    b2ind=sind(2);
    rtol=2.0*abs(mi(wind)-mi(bind))/(abs(mi(wind))+abs(mi(bind)));
    if (rtol<ftol)
        bconverged = true;
        break;
    end
    
    for j=1:dim+1
        dp(j)=norm(pnts(j,:)-pnts(bind,:));
    end
    
    if (max(dp)<mindp)
        bconverged = true;
        break;
    end
    
    %Calculate centroid without worst point
    cent=mean(pnts(sind(1:dim),:));
    
    %Reflection
    pntr=cent+coeff_reflect*(cent-pnts(wind,:));
    mir=feval(errfn, vol, mu, sigma, lphpos, Mreg, pntr, sf, thresh);
    
    %better than best
    if (mir<mi(bind))
        % expand
        pnte=cent+coeff_expand*(cent-pnts(wind,:));
        mie=feval(errfn, vol, mu, sigma, lphpos, Mreg, pnte, sf, thresh);
        % if better than the best, use this one
        if (mie<mi(bind))
            mi(wind)=mie;
            pnts(wind,:)=pnte;
        else
            mi(wind)=mir;
            pnts(wind,:)=pntr;
        end
    %better than worst
    elseif (mir <mi(wind))
         mi(wind)=mir;
         pnts(wind,:)=pntr;
    else
        pntcw=cent-coeff_contract*(cent-pnts(wind,:));
        micw=feval(errfn, vol, mu, sigma, lphpos, Mreg, pntcw, sf, thresh);
        if (micw<mi(wind))
            mi(wind)=micw;
            pnts(wind,:)=pntcw;
        else
            for j=1:dim+1
                if (j~=bind)
                    pnts(j,:)= coeff_contract*(pnts(j,:)+pnts(bind,:));
                    mi(j)=feval(errfn, vol, mu, sigma, lphpos, Mreg, pnts(j,:), sf, thresh);
                end        
            end
        end
    end
        
    i=i+1;
end

[Y, sind]=sort(mi);
bind=sind(1);
optpnt=pnts(bind,:);

Mt=eye(4,4);
Mt(1,4)=optpnt(1);
Mt(2,4)=optpnt(2);
Mt(3,4)=optpnt(3);

Ms=eye(4,4);
Ms(1,1)=optpnt(4);
Ms(2,2)=optpnt(5);
Ms(3,3)=optpnt(6);

Mrl=eye(4,4);
Mrl(2,2) = cos(optpnt(7));
Mrl(2,3) = sin(optpnt(7));
Mrl(3,2) = -sin(optpnt(7));
Mrl(3,3) = cos(optpnt(7));

Mrp=eye(4,4);
Mrp(1,1) = cos(optpnt(8));
Mrp(1,3) = sin(optpnt(8));
Mrp(3,1) = -sin(optpnt(8));
Mrp(3,3) = cos(optpnt(8));

Mrh=eye(4,4);
Mrh(1,1) = cos(optpnt(9));
Mrh(1,2) = sin(optpnt(9));
Mrh(2,1) = -sin(optpnt(9));
Mrh(2,2) = cos(optpnt(9));

Mss = eye(4,4);
Mss(1,2)=optpnt(10);
Mss(1,3)=optpnt(11);
Mss(2,3)=optpnt(12);

Maffine=Mt*Mrl*Mrp*Mrh*Ms*Mss*Mreg;

fmi=Y(1);
