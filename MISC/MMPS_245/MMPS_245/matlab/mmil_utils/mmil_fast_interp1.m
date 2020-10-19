function mati = fast_interp1(valmat,indmat,method,extrapval)

if ~exist('extrapval','var'), extrapval = NaN; end
dims = size(valmat);
wtmat = indmat-floor(indmat);
indmat0 = floor(indmat);
indmat1 = floor(indmat)+1;
outsideindvec = find(indmat0<1|indmat0>dims(1)|indmat1<1|indmat1>dims(1));
indmat0 = max(1,min(dims(1),floor(indmat)));
indmat1 = max(1,min(dims(1),floor(indmat)+1));
mati = (1-wtmat).*valmat(indmat0)+wtmat.*valmat(indmat1);
mati(outsideindvec) = extrapval;
% Should set values outside bounds to extrapval
