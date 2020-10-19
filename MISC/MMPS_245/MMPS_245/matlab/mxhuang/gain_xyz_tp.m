function [G_the_phi,dip_rtp]=gain_xyz_tp(G_xyz,dip_xyz);
% function [G_the_phi,dip_rtp]=gain_xyz_tp(G_xyz,dip_xyz);
% G_xyz: M by 3p gain matrix in xyz system
% dip_xyz: p by 3 dipole locations in xyz
% G_the_phi: M by 2p gain in rtp system
% dip_rtp: dipole location in rtp (the and phi in degrees)

p=size(dip_xyz,1);
[m,n]=size(G_xyz);

if n~=3*p
   error('the num of columns in G_xyz has yo be 3*dipole_num');
end

dip_rtp=zeros(p,3);
G_the_phi=zeros(m,2*p);

for i=1:p,
   dip_rtp(i,:)=cart2rtp(dip_xyz(i,:));
   dip_rtp(i,2:3)=dip_rtp(i,2:3)/180*pi;
   unit_the=[1;0];
   unit_phi=[0;1];
   [dum,amp_the_xyz]=dip_rtp_xyz(dip_rtp(i,:),unit_the);
   [dum,amp_phi_xyz]=dip_rtp_xyz(dip_rtp(i,:),unit_phi);
   index_G_xyz=[-2:0]+3*i;
   G_xyz_tmp=G_xyz(:,index_G_xyz);
   index_G_the_phi=[-1:0]+2*i;
	G_the_phi(:,index_G_the_phi)=G_xyz_tmp*[amp_the_xyz amp_phi_xyz];	   
end
   
   
