function G = sarvas5(L,P,Order);
[m,n] = size(L);
if(m~=3), % should be 3 x m
  % make it so, could be multiple dipoles
  L = reshape(L(:),3,length(L(:))/3);
end

G = 0;
for i = 1:length(P),
  G = G + sarvas(L,P(i),Order);
end

function G = sarvas(L,P,Order);
% SARVAS MEG Forward Model, spherical head
% function G = sarvas(L,P,order);
% L is 3 x nL, each column a source location
% P.sensor is 3 x nR,each column a sensor location
% P.orient is 3 x nR, the sensor orientation
% P.center is 3 x nR, the sphere center for each sensor
%  if P.center in nonexistant or null, then assumed to be
%  all zeros.
% Order is 
%  -1 current dipole
%   0 focal(magnetic) dipole
%   1 1st order multipole
% P.weight is optional scalar weight, e.g 1e-7
% If sensors are radially oriented, see RADIAL
   
if(~isfield(P,'center')), % user did not provide
  P.center = []; % initialize to null
end
if(isempty(P.center)), % user gave as null
  P.center = zeros(size(P.sensor));  % set to coordinate origin
end

P.sensor = P.sensor - P.center; % shift sensor coordinates

nR = size(P.sensor,2); % number of sensors
nL = size(L,2);  % number of source points

Rn2 = sum(P.sensor.^2,1); % distance to sensor squared
Rn = sqrt(Rn2); % distance

if(nR >= nL), % more sensors than dipoles
  if(Order == 1),
    G = zeros(nR,12*nL); % gain matrix
  else    
    G = zeros(nR,3*nL);  % gain matrix
  end
  
  for Li = 1:nL,
    Lmat = L(:,Li+zeros(1,nR)); % matrix of location repeated
    Lmat = Lmat - P.center; % each center shifted relative to its center
    D = P.sensor - Lmat;  % distance from souce to sensors
    Dn2 = sum(D.^2,1); % distance squared
    Dn = sqrt(Dn2);  % distance
    R_dot_D = sum(P.sensor .* D);  % dot product of sensor and distance
    R_dot_Dhat = R_dot_D ./ Dn;  % dot product of sensor and distance
    
    F = Dn2 .* Rn + Dn .* R_dot_D;  % Sarvas' function F
    
    GF_dot_o = Dn2 .* sum(P.sensor.*P.orient) ./ Rn + ...
      (2 * Rn + R_dot_Dhat) .* sum(D.*P.orient) + ...
      Dn .* sum((D+P.sensor).*P.orient);
    
    tempF = GF_dot_o ./ F.^2;
    
    if(Order == -1), % current dipole model
      temp = cross(Lmat,P.orient) ./ F([1 1 1],:) - ...
        cross(Lmat,P.sensor) .* tempF([1 1 1],:);
      G(:,Li*3+[-2 -1 0]) = temp';
    elseif(Order == 0) % magnetic dipole model
      temp = P.sensor .* tempF([1 1 1],:) - P.orient ./ F([1 1 1],:);
      G(:,Li*3+[-2 -1 0]) = temp';
    elseif(Order == 1),  % 1st order multipole
      % first the dipole
      temp_m = P.sensor .* tempF([1 1 1],:) - P.orient ./ F([1 1 1],:);
      % then the quadrupole
      temp1 = -(2*Rn + R_dot_Dhat + Dn);
      temp2 = -sum(D.*P.orient) ./ Dn;
      temp3 = -(2*sum(P.sensor.*P.orient)./Rn + sum((D+P.sensor).*P.orient)./Dn - ...
        sum(D.*P.orient).*R_dot_D./(Dn2.*Dn));
      
      GGpF_dot_o = temp1([1 1 1],:) .* P.orient + ...
        temp2([1 1 1],:).*P.sensor + temp3([1 1 1],:) .* D;
      temp1 = -(2*Rn + R_dot_Dhat);
      GpF = temp1([1 1 1],:) .* D - Dn([1 1 1],:) .* P.sensor;
      temp1 = 1 ./ F.^2;
      temp2 = 2*GF_dot_o./F;
      temp_q = temp1(ones(1,9),:) .* (kronmat(GGpF_dot_o,P.sensor) + ...
        kronmat(GpF,P.orient - temp2([1 1 1],:).*P.sensor));
      G(:,Li*12+[-11:0]) = [temp_m;temp_q]';
    end

  end
  
else  % more dipoles than sensors    
  if(Order == 1)
    G = zeros(12*nL,nR);  % 1st order multipole gain matrix transposed
  else
    G = zeros(3*nL,nR);  % gain matrix transposed
  end
  
  for Ri = 1:nR,
    Rmat = P.sensor(:,Ri+zeros(1,nL)); % matrix of sensor repeated
    Omat = P.orient(:,Ri+zeros(1,nL)); % orientations
    Lmat = L - P.center(:,Ri+zeros(1,nL)); % shift centers to this coordinate
    
    D = Rmat - Lmat;
    Dn2 = sum(D.^2,1); % distance squared
    Dn = sqrt(Dn2);  % distance
    R_dot_D = sum(Rmat .* D);  % dot product of sensor and distance
    R_dot_Dhat = R_dot_D ./ Dn;  % dot product of sensor and distance
    
    F = Dn2 * Rn(Ri) + Dn .* R_dot_D;  % Sarvas' function F
    
    GF_dot_o = Dn2 * sum(P.sensor(:,Ri).*P.orient(:,Ri)) / Rn(Ri) + ...
      (2 * Rn(Ri) + R_dot_D ./ Dn) .* sum(D.*Omat) + ...
      Dn .* sum((D+Rmat).*Omat);
    
    tempF = GF_dot_o ./ F.^2;
    
    if(Order == -1), % current dipole model
      temp = cross(Lmat,Omat) ./ F([1 1 1],:) - ...
        cross(Lmat,Rmat) .* tempF([1 1 1],:);
    elseif(Order == 0) % magnetic dipole model
      temp = Rmat .* tempF([1 1 1],:) - Omat ./ F([1 1 1],:);
    elseif(Order == 1),  % 1st order multipole
      % first the dipole
      temp_m = Rmat .* tempF([1 1 1],:) - Omat ./ F([1 1 1],:);
      % then the quadrupole
      temp1 = -(2*Rn(Ri) + R_dot_Dhat + Dn);
      temp2 = -sum(D.*Omat) ./ Dn;
      temp3 = -(2*sum(P.sensor(:,Ri).*P.orient(:,Ri))./Rn(Ri) + sum((D+Rmat).*Omat)./Dn - ...
        sum(D.*Omat).*R_dot_D./(Dn2.*Dn));
      
      GGpF_dot_o = temp1([1 1 1],:) .* Omat + ...
        temp2([1 1 1],:).*Rmat + temp3([1 1 1],:) .* D;
      temp1 = -(2*Rn(Ri) + R_dot_Dhat);
      GpF = temp1([1 1 1],:) .* D - Dn([1 1 1],:) .* Rmat;
      temp1 = 1 ./ F.^2;
      temp2 = 2*GF_dot_o./F;
      temp_q = temp1(ones(1,9),:) .* (kronmat(GGpF_dot_o,Rmat) + ...
        kronmat(GpF,Omat - temp2([1 1 1],:).*Rmat));
      temp = [temp_m;temp_q];
    end
    
    G(:,Ri) = temp(:);
  end
  
  G = G';
  
end

if(isfield(P,'weight')),
   G = P.weight*G;
end

return

function c = cross(a,b);
% fast and simple
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:); ...
    a(3,:).*b(1,:)-a(1,:).*b(3,:); ...
    a(1,:).*b(2,:)-a(2,:).*b(1,:)];
return

function k = kronmat(a,b);
% column by column, not element by matrix
k = [a([1 1 1],:) .* b; ...
    a([2 2 2],:) .* b; ...
    a([3 3 3],:) .* b];
return