function fmi=calmi_hatl_mmil(vol, mu, sigma, lphpos, Mreg, pnts, sf,thresh)
Mt=eye(4,4);
Mt(1,4)=pnts(1);
Mt(2,4)=pnts(2);
Mt(3,4)=pnts(3);

Ms=eye(4,4);
Ms(1,1)=pnts(4);
Ms(2,2)=pnts(5);
Ms(3,3)=pnts(6);

Mrl=eye(4,4);
Mrl(2,2) = cos(pnts(7));
Mrl(2,3) = sin(pnts(7));
Mrl(3,2) = -sin(pnts(7));
Mrl(3,3) = cos(pnts(7));

Mrp=eye(4,4);
Mrp(1,1) = cos(pnts(8));
Mrp(1,3) = sin(pnts(8));
Mrp(3,1) = -sin(pnts(8));
Mrp(3,3) = cos(pnts(8));

Mrh=eye(4,4);
Mrh(1,1) = cos(pnts(9));
Mrh(1,2) = sin(pnts(9));
Mrh(2,1) = -sin(pnts(9));
Mrh(2,2) = cos(pnts(9));

Mss = eye(4,4);
Mss(1,2)=pnts(10);
Mss(1,3)=pnts(11);
Mss(2,3)=pnts(12);

Mr=Mt*Mrl*Mrp*Mrh*Ms*Mss*Mreg;
[vxlval inbound]= vol_getvxlsval(lphpos, vol, Mr);
%ind = find (inbound>0);
ind = find(inbound>0 & vxlval>thresh);
sum_numer=sum((vxlval(ind).*mu(ind))./(sigma(ind).*sigma(ind)));
sum_denom=sum((vxlval(ind).*vxlval(ind))./(sigma(ind).*sigma(ind)));
tv=(sf*vxlval(ind)-mu(ind))./sigma(ind);
fmi = tv'*tv/length(tv);

min_fract_points = 0.1; % need to have at least this fraction of the sampled points
min_num_points = min_fract_points*size(lphpos,1);

if length(ind)<min_num_points
  fmi = 1e10;
end;


