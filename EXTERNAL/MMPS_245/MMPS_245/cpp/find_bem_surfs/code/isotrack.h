

class isotrack: public mesh 
{

 public: 

        isotrack(string, string, int, string, string, int);

	virtual ~isotrack();

        void getIsotrackMatrix(void); 

        void rotateMesh(float, float, float, float, float, float); 

 protected: 

	string isotrackInputDirectory, isotrackInputFile; 
        
	int nIsotrack; // to be imported from the Isotrack input file 

        int nParameter; // 6 (i.e. 3+3) 

        int nMax; // number of iteration in the simplex method 

        float angularScale; // (rads)
  
        float translationalScale; // (mm) 
        
        float **isotrackCoordinate0; // a set to be computed 
    
        float **isotrackCoordinate; // a set to be read 

        float angle[3]; 

	float translation[3]; 

        float **transformationMatrix; // 4x4 matrix 

        float **transformationMatrix0; // 4x4 matrix 

        void computeMatrix(void); 

        virtual float costFunction(float*); 

};

