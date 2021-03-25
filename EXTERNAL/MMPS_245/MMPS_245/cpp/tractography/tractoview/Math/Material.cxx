
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Material.h"

Material::Material() : _ambient(0.0f,0.0f,0.0f,1.0f),
						_diffuse(0.3f,0.3f,0.3f,1.0f),
						_specular(1.0f,1.0f,1.0f,1.0f),
						_shininess(0.3f),
						_emissive(0.0f,0.0f,0.0f,1.0f)
{
}

// give - a number between [0,1] where 0 == no shine.
void Material::setShininess( const float& shine )
{
	_shininess = shine;
}

// give - a number between [0,1] where 0 == no transparency.
void Material::setTransparency( const float& trans )
{
	//??? should all these be equal? or just the diffuse?
	_ambient[3] = trans;
	_diffuse[3] = trans;
	_specular[3] = trans;
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setAmbient( const float& red, const float& green, const float& blue )
{
	_ambient[0] = red;
	_ambient[1] = green;
	_ambient[2] = blue;
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setDiffuse( const float& red, const float& green, const float& blue )
{
	_diffuse[0] = red;
	_diffuse[1] = green;
	_diffuse[2] = blue;
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setSpecular( const float& red, const float& green, const float& blue )
{
	_specular[0] = red;
	_specular[1] = green;
	_specular[2] = blue;
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setEmissive( const float& red, const float& green, const float& blue )
{
	_emissive[0] = red;
	_emissive[1] = green;
	_emissive[2] = blue;
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setAmbient( const float* ambient )
{
	_ambient[0] = ambient[0];
	_ambient[1] = ambient[1];
	_ambient[2] = ambient[2];
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setDiffuse( const float* diffuse )
{
	_diffuse[0] = diffuse[0];
	_diffuse[1] = diffuse[1];
	_diffuse[2] = diffuse[2];
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setSpecular( const float* specular )
{
	_specular[0] = specular[0];
	_specular[1] = specular[1];
	_specular[2] = specular[2];
}

// give - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::setEmissive( const float* emissive )
{
	_emissive[0] = emissive[0];
	_emissive[1] = emissive[1];
	_emissive[2] = emissive[2];
}





// returns - a number between [0,1] where 0 == no shine.
void Material::getShininess( float& shine )
{
	shine = _shininess;
}

// returns - a number between [0,1] where 0 == no transparency.
void Material::getTransparency( float& trans )
{
	//??? should all these be equal? or just the diffuse?
	trans = _ambient[3];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getAmbient( float* ambient )
{
	ambient[0] = _ambient[0];
	ambient[1] = _ambient[1];
	ambient[2] = _ambient[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getDiffuse( float* diffuse )
{
	diffuse[0] = _diffuse[0];
	diffuse[1] = _diffuse[1];
	diffuse[2] = _diffuse[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getSpecular( float* specular )
{
	specular[0] = _specular[0];
	specular[1] = _specular[1];
	specular[2] = _specular[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getEmissive( float* emissive )
{
	emissive[0] = _emissive[0];
	emissive[1] = _emissive[1];
	emissive[2] = _emissive[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getAmbient( float& red, float& green, float& blue )
{
	red = _ambient[0];
	green = _ambient[1];
	blue = _ambient[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getDiffuse( float& red, float& green, float& blue )
{
	red = _diffuse[0];
	green = _diffuse[1];
	blue = _diffuse[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getSpecular( float& red, float& green, float& blue )
{
	red = _specular[0];
	green = _specular[1];
	blue = _specular[2];
}

// returns - three numbers [r,g,b] between [0,1] where 0 is no color
void Material::getEmissive( float& red, float& green, float& blue )
{
	red = _emissive[0];
	green = _emissive[1];
	blue = _emissive[2];
}
