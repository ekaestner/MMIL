
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#ifndef MATERIAL_INCLUDED
#define MATERIAL_INCLUDED

#include "ColorRGBA.h"
class Material
{
// Construction:
public:
	//: constructor
	Material();

// Methods:
public:
	// give - a number between [0,1] where 0 == no shine.
	void setShininess( const float& shine );

	// give - a number between [0,1] where 0 == no transparency.
	void setTransparency( const float& trans );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setAmbient( const float& red, const float& green, const float& blue );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setDiffuse( const float& red, const float& green, const float& blue );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setSpecular( const float& red, const float& green, const float& blue );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setEmissive( const float& red, const float& green, const float& blue );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setAmbient( const float* ambient );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setDiffuse( const float* diffuse );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setSpecular( const float* specular );

	// give - three numbers [r,g,b] between [0,1] where 0 is no color
	void setEmissive( const float* emissive );


	// returns - a number between [0,1] where 0 == no shine.
	void getShininess( float& shine );

	// returns - a number between [0,1] where 0 == no transparency.
	void getTransparency( float& trans );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getAmbient( float* ambient );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getDiffuse( float* diffuse );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getSpecular( float* specular );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getEmissive( float* emissive );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getAmbient( float& red, float& green, float& blue );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getDiffuse( float& red, float& green, float& blue );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getSpecular( float& red, float& green, float& blue );

	// returns - three numbers [r,g,b] between [0,1] where 0 is no color
	void getEmissive( float& red, float& green, float& blue );

//Aliases
public:
	//: Alias to the shininess component
	//  represents - a number between [0,1] where 0 == no shine.
	const float&		shininess() const;
	float&				shininess();

	//: Read only alias to the transparency component
	//  represents - a number between [0,1] where 0 == opaque.
	const float&		transparency() const;

	//: Alias to the ambient component
	// represents - three numbers [r,g,b] between [0,1] where 0 is no color
	const ColorRGBA&	ambient() const;
	ColorRGBA&			ambient();

	//: Alias to the diffuse component
	// represents - three numbers [r,g,b] between [0,1] where 0 is no color
	const ColorRGBA&	diffuse() const;
	ColorRGBA&			diffuse();

	//: Alias to the specular component
	// represents - three numbers [r,g,b] between [0,1] where 0 is no color
	const ColorRGBA&	specular() const;
	ColorRGBA&			specular();

	//: Alias to the emissive component
	// represents - three numbers [r,g,b] between [0,1] where 0 is no color
	const ColorRGBA&	emissive() const;
	ColorRGBA&			emissive();


protected:
    float		_shininess;
    ColorRGBA	_ambient; 
    ColorRGBA	_diffuse;
    ColorRGBA	_specular;
	ColorRGBA	_emissive;
    
     // TODO: add the flag for front, or back, or both sides...
};

////////////////////////////////////////////////////////////
//: Alias to the shininess component
//  represents - a number between [0,1] where 0 == no shine.
////////////////////////////////////////////////////////////
inline const float&	Material::shininess() const { return _shininess; }
inline float&		Material::shininess() { return _shininess; }

////////////////////////////////////////////////////////////
//: Read only alias to the transparency component
//  represents - a number between [0,1] where 0 == opaque.
////////////////////////////////////////////////////////////
inline const float&	Material::transparency() const { return _diffuse[3]; }

////////////////////////////////////////////////////////////
//: Alias to the ambient component
// represents - three numbers [r,g,b] between [0,1] where 0 is no color
////////////////////////////////////////////////////////////
inline const ColorRGBA&	Material::ambient() const { return _ambient; }
inline ColorRGBA&		Material::ambient() { return _ambient; }

////////////////////////////////////////////////////////////
//: Alias to the diffuse component
// represents - three numbers [r,g,b] between [0,1] where 0 is no color
////////////////////////////////////////////////////////////
inline const ColorRGBA&	Material::diffuse() const{ return _diffuse; }
inline ColorRGBA&		Material::diffuse(){ return _diffuse; }

////////////////////////////////////////////////////////////
//: Alias to the specular component
// represents - three numbers [r,g,b] between [0,1] where 0 is no color
////////////////////////////////////////////////////////////
inline const ColorRGBA&	Material::specular() const{ return _specular; }
inline ColorRGBA&		Material::specular(){ return _specular; }

////////////////////////////////////////////////////////////
//: Alias to the emissive component
// represents - three numbers [r,g,b] between [0,1] where 0 is no color
////////////////////////////////////////////////////////////
inline const ColorRGBA&	Material::emissive() const{ return _emissive; }
inline ColorRGBA&		Material::emissive(){ return _emissive; }

#endif
