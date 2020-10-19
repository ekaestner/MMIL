// multiply matrix * ray, returns the resulting ray
template <class dataType>
inline void multiply( Ray3<dataType>& result, const Matrix4f& mat, const Ray3<dataType>& ray )
{
   // TODO: convert to seg, mul mat*seg (easy), convert back to ray.


   // Find the start and end point of the ray
	const Vec3<dataType>& rayStartPoint = ray->origin();
	Vec3<dataType> rayEndPoint = this->direction + this->origin;
	
	result.origin = mat * rayStartPoint;
	result.direction = mat * rayEndPoint;
	
	// make sure the direction is normalized.
	result.direction.normalize();
