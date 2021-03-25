function coordmat_out = ApplyMvxl2lphCoordlist(coordmat_in,M)

coordmat_out = cat(2,coordmat_in,ones(size(coordmat_in,1),1))*M(:,[2 1 3 4])';
coordmat_out = coordmat_out(:,1:3);
