function [Maffine, fmi, optpnt]=affine_reg_vol2vol(volt, vol, M_volt_to_vol_rb, bsmooth, lphpos, errfn, maxiter, ftol, mindp, initp, wei)
% Affine registration By using  Simplex
%
%       [M_vol1_to_vol2_af, fmi, optpnt]=
%           affine_reg_vol2vol(volt, vol, M_volt_to_vol_rb, [bsmooth], [lphpos], [function_name], 
%                           [maxiter], [ftol], [mindp], [initp], [wei]))
%   input:
%       volt: Target Volume 
%       vol: Volume to be registered
%       M_vol1_to_vol2_rb: Rigid body registration matrix from vol1 to vol2
%       bsmooth  smooth volume (default =true);
%       lphpos : Sampling lph points
%       errfun: error function name (calmi_corr (default) or calmi_ss)
%       maxiter: max. iteration (default = 2000)
%       ftol: tolerance (default =0.001)
%       mindp: minimum distance between points (default = 0.01)
%       initp: Initial Point (default = [1 1 1 0 0.0 0.0])
%       wei: offset (default = [0.1 0.1 0.1 0.1 0.1 0.1])
%
%   Output:
%       M_vol1_to_vol2_af: Affine transfrom from vol1 to vol2
%       fmi: minimum
%       optpnt: optimal point
%
%
% Created: ? by ?
% Last Mod: 04/01/10 by Don Hagler
%

if ~exist('lphpos')
    ith = vol_getBiModalThreshold(volt);
    range=[-70 -70 -70; 70 70 70; 2. 2. 2.];
    ln=ceil((range(2,1)-range(1,1))/range(3,1)+1);
    pn=ceil((range(2,2)-range(1,2))/range(3,2)+1);
    hn=ceil((range(2,3)-range(1,3))/range(3,3)+1);
    tsize=ln*pn*hn;
    numbins=ceil(tsize^(1/3));
    lphpos=zeros(tsize,4);
    count=1;
    lphcomv1=vol_getCOM(volt, true, ith);
    lphcomv1(3)=lphcomv1(3)+20;
    lphrange=[lphcomv1(1:3)';range];
     for i=1:ln
        for j=1:pn
            for k=1:hn
                lphpos(count,1:3)=lphcomv1(1:3)+(i-ln/2)*range(3,1)*volt.DirCol...
                    +(j-pn/2)*range(3,2)*volt.DirRow...
                    +(k-hn/2)*range(3,3)*volt.DirDep;
                lphpos(count,4)=1;
                count=count+1;
            end
        end
     end
end

if (~exist('bsmooth'))
    bsmooth=true;
end   


if (~exist('errfn'))
    errfn = 'calmi_corr';
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

initp=[1 1 1 0 0.0 0.0];
if ~exist('initp')
    initp=[1 1 1 0 0.0 0.0];
end

wei=[0.1 0.1 0.1 0.1 0.1 0.1];
if ~exist('wei')
   wei=[0.1 0.1 0.1 0.1 0.1 0.1];
end

if (bsmooth)
  volt = vol_filter(volt, 1);
  vol = vol_filter(vol, 1);
end

sf=1.;
dim=6;
coeff_reflect=1.;
coeff_expand=2.;
coeff_contract=0.5;
mi=zeros(dim+1,1);
pnts=zeros(dim+1,dim);
Mreg = M_volt_to_vol_rb;

for i=1:dim
    pnts(i,:)=initp;
    pnts(i,i)=initp(i)+wei(i);
end
pnts(dim+1,:)=initp;

for i=1:dim+1
    mi(i)=feval(errfn, volt,vol, lphpos, Mreg, pnts(i,:), sf);
end

%disp    'Registration...'
    
i=0;
ftol=0.001;
mindp=0.01;
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
        break;
    end
    
    for j=1:dim+1
        dp(j)=norm(pnts(j,:)-pnts(bind,:));
    end
    
    if (max(dp)<mindp)
        break;
    end
    
    %Calculate centroid without worst point
    cent=mean(pnts(sind(1:dim),:));
    
    %Reflection
    pntr=cent+coeff_reflect*(cent-pnts(wind,:));
    mir=feval(errfn, volt,vol, lphpos, Mreg, pntr, sf);
    
    %better than best
    if (mir<mi(bind))
        % expand
        pnte=cent+coeff_expand*(cent-pnts(wind,:));
        mie=feval(errfn, volt,vol, lphpos, Mreg, pnte, sf);
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
        micw=feval(errfn, volt,vol, lphpos, Mreg, pntcw, sf);
        if (micw<mi(wind))
            mi(wind)=micw;
            pnts(wind,:)=pntcw;
        else
            for j=1:dim+1
                if (j~=bind)
                    pnts(j,:)= coeff_contract*(pnts(j,:)+pnts(bind,:));
                    mi(j)=feval(errfn, volt,vol, lphpos, Mreg, pnts(j,:), sf);
                end        
            end
        end
    end
        
    
    %pnts;
    i=i+1;
end

[Y, sind]=sort(mi);
bind=sind(1);
optpnt=pnts(bind,:);
Msa=eye(4,4);
Msa(1,1)=optpnt(1);
Msa(2,2)=optpnt(2);
Msa(3,3)=optpnt(3);

Mss=eye(4,4);
Mss(1,2)=optpnt(4);
Mss(1,3)=optpnt(5);
Mss(2,3)=optpnt(6);

Maffine=Msa*Mss*Mreg;
fmi=Y(1);
