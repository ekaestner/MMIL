

class deterministic: public mesh
{

 public: 

        deterministic(const string&, int, string, string);

	virtual ~deterministic();

        void segmentSurface(string, string); 
                            
 protected: 

	struct structModelCoefficient {

	  int nMove; 
	  int nRest; 
	  int nRelax;
	  
	  float coefTangential; 
	  float coefNormal; 
	  float coefMRI; 
          float directionMRI; 
	  float coefRepelling; 
	  float intensityThreshould; 

	} *surface; 

	ifstream inputFile; 

	void addMriForce(float, float, float);

	float moveVertices(float, float, float, float, float, float);

	void readModelCoefficients(void); 

	float getNumber(void); 
 
};

