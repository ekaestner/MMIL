
#ifndef NORMAL_DISK_COPY
#define NORMAL_DISK_COPY
#undef  NORMAL_DISK_COPY
#endif

/* make sure that max number of neighbor vertices is defined outside */

class mesh: public mri 
{

 public: 

        mesh(const string&, int, string, string);

	virtual ~mesh();

        void initializeSurface(float); 

        void createRepellingVolume(float, int); 

        void talairachForward(void); 

        void talairachBackward(void); 

        float smoothSurface(int); 

        void dilateSurface(float); 

        float fillContour(void); 

        float computeInsideIntensity(void); 

        void createBorderMap(float, int); 

        void exportFilling(string, string);

        void importSurface(string, string); 

        void exportSurface(int, int, string, string); 


	float computeSolidAngle(void); 

 protected: 

	string triInputDirectory; 

	int meshNumber;
 
	int nVertex, nFace;

        int LOCAL, GLOBAL; 

	int weightInside, weightBorder, weightOutside; 

	float coefSmoothLocal, coefSmoothGlobal; 

	float **talairach, **inverseTalairach; 

        int ***contourFilling;

        float ***borderMap;  

        int **icoFaceElement; 
	float **icoCoordinate; 
       
	struct structVertexType {

	  float r[3]; // Coordinates 

          float fNormal[3];  // z-axis in local coordinate system 
          float xTangent[3]; // x-axis in local coordinate system
          float yTangent[3]; // y-axis in local coordinate system
 
	  float force[3], repellingForce[3];

	  int vnum; // First-order neighboring vertices 
	  int v[6];

          int vnum2; // Second-order neighboring vertices 
          int v2[12]; 

	  int fnum; // Neighboring triangles 
	  int f[6];

          float curvature; // Local curvature based on 1st order neighbors  

          float curvature2; // Same based on 1st and 2nd order neighbors 

          float displacement; // surface position in the local system 
 
	} *vertex; 

	struct structFaceType { 
	  int nTri[3]; 
          float fNormal[3]; 
	} *face; 

	struct structTriFileInfoType { 
	  int nVertex; 
	  int nFace; 
	} *triFileInfo; 

	//float computeSolidAngle(void); 

	void computeRepellingForce(float);

	void addInertialForce(float, float);

	void localCoordinateSystem(void);

        void computeCurvature(void); 

	int getFilling(float, float, float); 

	float getBorderMap(float, float, float); 

	/* add optimization functions for the use in the derived classes */

	void simplexND(float**, float*, int); 

	float trial(float**, float*, float*, int, float, int); 

	virtual float costFunction(float*);  // abstract function 
};

