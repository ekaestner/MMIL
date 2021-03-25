

class statistical: public mesh 
{

 public: 

        statistical(string, string, int, string, string, int);

	virtual ~statistical();

	/* make this function virtual to the same name as 
           in deterministic class, i.e. 
           segmentSurface or use different number of variables */ 

        void segmentSurfaceStat(float); 

        float adjustSurface(float, float, float, float, int, int); 

	/* add nTemplateMax or TemplateSize to this function */

        void exportTemplate(string, string, int); 

        void importMeanTemplate(string, string); 

        void importMeanSurface(string, string);   

        void importEof(string, string); 

        void moveEof(void); 

	virtual float costFunction(float*);  

        void optimizeCoefficient(void); 

	/* add function to free memory for the templates and eofs */ 

                         
 protected: 

        int nTemplateMax; 

        int nEof; 

        float stepTemplate; 

        float **meanSurface; 

	float **templateVal; 
	float **templateDev; 

	float *coefEof; 
	float ***eof; 

	/* make this function virtual as well ? */ 

	void computeShift(int, int);  

	float moveVerticesStat(float, float, float, float, float, float);

};

