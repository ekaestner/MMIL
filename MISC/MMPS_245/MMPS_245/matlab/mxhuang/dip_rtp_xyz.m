function [L_cart,T_cart]=dip_rtp_xyz(L_pol,T_pol);

% function [L_cart,T_cart]=dip_rtp_xyz(L_pol,T_pol);
% L_pol = [ro the phi];
% T_pol = 2 x Nt: [T_the(t) T_phi(t)]
% L_cart = [x y z];
% T_cart = 3 x Nt: [T_x(t) T_y(t) T_z(t)]

  L_cart = zeros(1,3);
  T_cart = zeros(3,size(T_pol,2));

  ro=L_pol(1);
  the=L_pol(2);
  phi=L_pol(3);
  Pthe=T_pol(1,:);
  Pphi=T_pol(2,:);

  x=ro*sin(the)*cos(phi);
  y=ro*sin(the)*sin(phi);
  z=ro*cos(the);       

  L_cart(1)=x;
  L_cart(2)=y;
  L_cart(3)=z;

  T_cart(1,:)=Pthe*sin(the+pi/2.0)*cos(phi)+Pphi*cos(phi+pi/2.0);
  T_cart(2,:)=Pthe*sin(the+pi/2.0)*sin(phi)+Pphi*sin(phi+pi/2.0);
  T_cart(3,:)=Pthe*cos(the+pi/2.0);
