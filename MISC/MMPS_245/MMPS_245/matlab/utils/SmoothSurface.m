function surf_out = SmoothSurface(surf_in,niter,normcomp)

tmpvertices = surf_in.vertices;
for j = 1:3
  for smooth_iter = 1:niter
    tmpvertices(:,j) = single(smoothVertexValue(surf_in,double(tmpvertices(:,j))));
  end
end
dvert = tmpvertices-surf_in.vertices;
normproj = sum(dvert.*surf_in.normals,2);
dvert = dvert - (1-normcomp)*(normproj*ones(1,3)).*surf_in.normals;
surf_out = surf_in;
surf_out.vertices = surf_out.vertices + dvert;
surf_out.normals = computeNormals(surf_out);
