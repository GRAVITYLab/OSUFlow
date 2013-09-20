#ifndef PATHLINE_LOADER_H
#define PATHLINE_LOADER_H

#include <list>
#include <map>
#include <string.h>
#include <OSUFlow.h>

typedef list<vtListTimeSeedTrace *> SeedTraceList4D;
typedef vtListTimeSeedTrace jcListTimeSeedTrace;

// save to file and delete immediately
class PathlineLoader
{
public:
    SeedTraceList4D traceList;
    VECTOR4 bmin, bmax;
public:
    PathlineLoader(const char *filename) {
    	printf("Loading trace file: %s\n", filename);
    	if (strcmp(filename+strlen(filename)-3, "raw")==0)
    		loadDataRaw(filename);
    	else //if (strcmp(filename+strlen(filename)-3,"out")==0)
    		loadDataFormal(filename);
    	printf("File loaded\n");
    }

    ~PathlineLoader()
    {

    	SeedTraceList4D::iterator it;
		for (it=traceList.begin(); it!=traceList.end(); ++it)
		{
			jcListTimeSeedTrace::iterator it1;
			for (it1=(*it)->begin(); it1!=(*it)->end(); ++it1)
				delete *it1;
			delete *it;
		}

    }

    class SeedMapComparer{ // simple comparison function
       public:
          bool operator()(const VECTOR4 x,const VECTOR4 y) const {
        	  return x[0]>y[0] ||
        			  x[0]==y[0] && x[1]>y[1] ||
        			  x[0]==y[0] && x[1]==y[1] && x[2]>y[2] ||
        			  x[0]==y[0] && x[1]==y[1] && x[2]==y[2] && x[3]>y[3];
          } // returns x>y
    };
    typedef std::map<VECTOR4, jcListTimeSeedTrace *, SeedMapComparer> SeedMap;


	void dump()
	{
		printf("Start dumping traces:\n");
    	SeedTraceList4D::iterator it;
    	for (it=traceList.begin(); it!=traceList.end(); ++it)
    	{
			jcListTimeSeedTrace::iterator it1;
			for (it1=(*it)->begin(); it1!=(*it)->end(); ++it1)
			{
				VECTOR4 &vec = **it1;
				printf("%f %f %f %f, ", vec[0], vec[1], vec[2], vec[3]);
			}
			printf(" <end trace>\n");
		}
		printf("End dump\n");
	}

	void connectTraces()
	{
		connectTraces(this->traceList);
	}

    static void connectTraces(SeedTraceList4D &traceList)
    {
    	printf("Original traces=%zu\n", traceList.size());
    	SeedMap seedMap;
    	SeedTraceList4D::iterator it;
    	for (it=traceList.begin(); it!=traceList.end(); ++it)
    	{
    		jcListTimeSeedTrace *shortTrace = *it;
    		if (shortTrace->size()<=1) continue;
    		jcListTimeSeedTrace::iterator it1;
    		it1 = shortTrace->end(); --it1;
    		VECTOR4 *p_start = *(shortTrace->begin());
    		VECTOR4 *p_end = *it1;

			//printf("Checking %f %f %f %f\n", (*p_start)[0], (*p_start)[1],(*p_start)[2],(*p_start)[3]);
    		SeedMap::iterator itSeedMap = seedMap.find(*p_start);
    		if (itSeedMap == seedMap.end()) {
    			seedMap.insert(std::pair<VECTOR4, jcListTimeSeedTrace *>(*p_end, shortTrace));
    		} else {
    			jcListTimeSeedTrace *longTrace = itSeedMap->second;
    			// remove first identical point
    			delete *(shortTrace->begin());
    			shortTrace->erase(shortTrace->begin());

    			longTrace->splice(longTrace->end(), *shortTrace);
    			seedMap.erase(itSeedMap);
    			seedMap.insert(std::pair<VECTOR4, jcListTimeSeedTrace *>(*p_end, longTrace));
    			//printf("merged\n");
    		}

    	}
    	assert(traceList.size());
    	for (it = traceList.begin(); it!= traceList.end(); ++it)
    		if ((*it)->size() == 0) {
    			SeedTraceList4D::iterator temp = it--;
    			delete *temp;
    			traceList.erase(temp);
    		}
    	printf("Merged traces=%zu\n", traceList.size());
    }

private:
	void loadDataRaw(const char *filename)
	{
	    FILE *fp;
        fp = fopen(filename, "rb");
        printf("Loading file: %s\n", filename);

    	fread(&bmin, 4, sizeof(float), fp);
    	fread(&bmax, 4, sizeof(float), fp);

        float vec[4];
        jcListTimeSeedTrace *trace = new jcListTimeSeedTrace;
        while(!feof(fp))
        {
        	fread(vec, sizeof(float), 4, fp);
        	if (isnan(vec[0]) ) {
        		traceList.push_back(trace);
        		trace = new jcListTimeSeedTrace;
#ifdef _DEBUG
    			printf("\n");
#endif
    			continue;
        	}
        	trace->push_back(new VECTOR4(vec[0], vec[1], vec[2], vec[3]));
#ifdef _DEBUG
			printf("%f %f %f %f, ", vec[0], vec[1], vec[2], vec[3]);
#endif
        }
       // printf("%d traces\n", (int)traceList.size());
        delete trace;
		fclose(fp);
	}

	void loadDataFormal(const char *filename)
	{
		FILE *fp = fopen(filename, "rb");

    	fread(&bmin, 4, sizeof(float), fp);
    	fread(&bmax, 4, sizeof(float), fp);

		std::vector<int> npts_ary;
		while(!feof(fp))
		{
			int n;
			fread(&n, sizeof(int),1, fp);
			if (n==-1) break;
			npts_ary.push_back(n);
			//printf("%d traces\n", n);
		}
		int i,j;
		for (i=0; i<npts_ary.size(); i++)
		{
			float vec[4];
			jcListTimeSeedTrace *ptrace = new jcListTimeSeedTrace;
			for (j=0; j<npts_ary[i]; j++)
			{
				float dum;
				fread(&vec, sizeof(float),4, fp);
				ptrace->push_back(new VECTOR4(vec[0], vec[1], vec[2], vec[3]));
#ifdef _DEBUG
				printf("%f %f %f %f, ", vec[0], vec[1], vec[2], vec[3]);
#endif
			}
#ifdef _DEBUG
			printf("\n");
#endif
			traceList.push_back(ptrace);

		}
		fclose(fp);
	}


};


#endif
