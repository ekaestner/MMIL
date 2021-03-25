function sangle = solid_angle2(r,r1,r2,r3)
%SOLID_ANGLE2 - Solid angle of a viewed triangle
% function sangle = solid_angle2(r,r1,r2,r3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%         r : observation point (Nobs x 3) 
%        r1 : First Vertice of Triangles (Nt x 3) 
%        r2 : Second Vertice of Triangles (Nt x 3)
%        r3 : Third Vertice of Triangles (Nt x 3)
%
% Output:
%         sangle : Function Result (Nt x Nobs)
%
% Ref: Oosterom, Strackee, IEEE Trans. BME, pp125-126, 1983 
%
% Note: if r is on the triangle, sangle will be returned as 2*pi or -2*pi
%
% %%% John Ermer 04/04/00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
Nt = size(r1,1);
Nobs = size(r,1);
%
temp = reshape(repmat(r',Nt,1),3,Nt*Nobs)'; % (Nt*Nobs)x(3)
r1_r = repmat(r1,Nobs,1) - temp;            % (Nt*Nobs)x(3)
r2_r = repmat(r2,Nobs,1) - temp;            % (Nt*Nobs)x(3)
r3_r = repmat(r3,Nobs,1) - temp;            % (Nt*Nobs)x(3)
%
n1 = rownorm(r1_r);
n2 = rownorm(r2_r);
n3 = rownorm(r3_r);
%
%%%% compute sangle for all vector sets
%
sangle = 2*atan2(dotprod_fun(r1_r,crossprod_fun(r2_r,r3_r)),n1.*n2.*n3+...
                 n1.*dotprod_fun(r2_r,r3_r) + n2.*dotprod_fun(r1_r,r3_r) +...
                 n3.*dotprod_fun(r1_r,r2_r) );
sangle = reshape(sangle,Nt,Nobs);
%

