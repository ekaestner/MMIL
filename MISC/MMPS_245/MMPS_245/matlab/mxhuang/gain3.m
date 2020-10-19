function G = gain3(Rq,Re,R,s,layer,nmax);
%-------------------------------------------------------------------------------
% function G = gain3(Rq,Re,R,s,layer,nmax);
% computes potential in a three shell model
% Input:
% Rq   : dipole location(in meters)                     no_dipole x 3
% Re   : EEG sensors(in meters)                         no_sensor x 3
% R    : radii(in meters) of sphere from OUTERMOST to INNERMOST    no_layer x 1
% s: conductivity from OUTERMOST to INNERMOST           no_layer x 1
% layer: 1 -- scalp,                                    scalar
%        2 -- skull,
%        3 -- brain                  
% nmax : # of terms kept(OPTIONAL)                      scalar
%  
% Output:
% G    : EEG gain matrix                                no_sensor x (3*no_dipole)
%

% reference:
% (B.18)- B(.22) in ' Brain stimulation using electromagnetic sources:
% theoretic aspects', by Leon Heller and David B. van Hulsteyn
% Aside: to compute the EEG on the outermost layer for single or N shell
% model, call gainp_sph
% this function calls dlegpoly, rownorm, mx_dotprod
% the code is NOT optimized
% CCH, Dec 1995 at LANL
%-------------------------------------------------------------------------------


Ren = norm(Re(1,:));
% check input 
if R(1)~= max(R)
  error(' head radii must be specified from OUTERMOST to INNERMOST layer ! ')
end
if abs(Ren-R(layer)) > 1/1000*R(layer)
  error(' EEG sensors maybe not on the layer you specified !')
end
if size(Rq,2) ~= 3
  error('dipole location must have three columns !')
end

if nargin < 6,     % choose # terms
    nmax = fix(10/(1- (max(rownorm(Rq))/Ren)));
end

no_dipole = size(Rq,1);
G = zeros(size(Re,1),3*no_dipole);
s21 = s(2)-s(1); 
s32 = s(3)-s(2);
r21 = R(2)/R(1);
r31 = R(3)/R(1);
r32 = R(3)/R(2);
r21n =1;
r31n =1;
r32n =1;

D = zeros(nmax,1);
N = zeros(nmax,1);
for n=1:nmax,      % (B.19) - (B.22) in the reference
  r21n = r21n*r21;    % (R2/R1)^n
  r31n = r31n*r31;
  r32n = r32n*r32;
  D(n) = det([1           s21*r21n*r21        s32*r31n*r31;
            -(n+1)*r21n   n*s(2)+(n+1)*s(1)   n*s32*r32n*r32;
	    -(n+1)*r31n   -(n+1)*s21*r32n     n*s(3)+(n+1)*s(2)]);
	  
  if layer==1,
     N(n)= det([   1                n*s21*r21n*r21     n*s32*r31n*r31;
                (R(1)/R(2))^(n+1)   n*s(2)+(n+1)*s(1)  n*s32*r32n*r32;
                (R(1)/R(3))^(n+1)   -(n+1)*s21*r32n    n*s(3)+(n+1)*s(2)]);
  elseif layer==2,   	      
     N(n)= s(1)*det([ n           1                 n*s32*r31n*r31;
                   -(n+1)*r21n   (R(1)/R(2))^(n+1)  n*s32*r32n*r32;
                   -(n+1)*r31n   (R(1)/R(3))^(n+1)   n*s(3)+(n+1)*s(2)]);
  else	      
     N(n)= s(1)*det([ n          n*s21*r21n*r21       1 ; 
                   -(n+1)*r21n   n*s(2)+(n+1)*s(1)  (R(1)/R(2))^(n+1); 
                   -(n+1)*r31n    -(n+1)*s21*r32n    (R(1)/R(3))^(n+1)]);
  end     
end	    
      
n=[1:nmax]';       
nm1 = n-1;
common = (2*n+1).*N./(4*pi*s(1).*D.*R(1).^(n+1));

for i = 1:no_dipole,  % loop over all dipoles
  Rqn =  norm(Rq(i,:));
  commonnew = common.*Rqn.^nm1;
  cosx  = mx_dotprod(ones(size(Re,1),1)*Rq(i,:),Re)./(Rqn*Ren); 
  [P,dP] = dlegpoly(nmax,cosx);  % evaluate legendre poly and its derivative

  Gterm1 = P'*commonnew;
  Gterm2 = dP'*(commonnew./n ./Ren);
  % gain matrix,no_sensor x 3, related to the i-th dipole
  z= Re - mx_dotprod(ones(size(Re,1),1)*Rq(i,:),Re)./(sum(Rq(i,:).^2))*Rq(i,:);

  G(:,3*i-2:3*i) = Gterm1*Rq(i,:)/norm(Rq(i,:))+ [z(:,1).*Gterm2,z(:,2).*Gterm2,z(:,3).*Gterm2]; 
  
end


%keyboard
