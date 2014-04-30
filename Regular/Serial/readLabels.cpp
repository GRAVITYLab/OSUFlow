#include <iostream>
#include <fstream> 
#include <string> 
#include <string.h> 
#include <stdlib.h> 

using namespace std;

int main(int argc, char** argv)
{
	if(argc<2)
	{
		cout<<"Usage: "<<argv[0]<<" expert-label.txt"<<endl;
	}
	ifstream infile;
	int startX,endX,startY,endY,startZ,endZ;

	infile.open(argv[1]);

	string line;

	getline(infile,line);

	cout<<"file content is "<<line<<endl;

	char * cstr = new char [line.length()+1];
	strcpy (cstr, line.c_str());
	char *pch;

	pch = strtok(cstr," ");
	endX = atoi(pch);
	pch = strtok(NULL," ");
	startX = atoi(pch);
	pch = strtok(NULL," ");
	endY = atoi(pch);
	pch = strtok(NULL," ");
	startY = atoi(pch);
	pch = strtok(NULL," ");
	endZ = atoi(pch);
	pch = strtok(NULL," ");
	startZ = atoi(pch);
	cout<<"From "<<startX<<" "<<startY<<" "<<startZ<<" To "<<endX<<" "<<endY<<" "<<endZ<<endl;

	for(int k=startZ; k<endZ; k++)
		for(int j=startY; j<endY; j++)
			for(int i=startX; i<endX; i++)
			{
				cout<<atoi(strtok(NULL," "))<<" ";
			}

	infile.close();
}