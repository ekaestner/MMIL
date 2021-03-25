

class deterministic: public mesh 
{

 public: 

        deterministic(string, string, int, string, string, int);

	virtual ~deterministic();

        void segmentSurface(string, string, int, int, int, 
			    float, float, float, float, float, float); 
                            
 protected: 

	void addMriForce(float, float, float);

	float moveVertices(float, float, float, float, float, float);
 
};

