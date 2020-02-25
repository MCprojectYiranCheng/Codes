CC = g++

asa241.o: asa241.cpp asa241.hpp
	$(CC) -c asa241.cpp -std=c++17 

importanceSamplingAndVectorTools.o:  importanceSamplingAndVectorTools.cpp importanceSamplingAndVectorTools.hpp 
	$(CC) -c importanceSamplingAndVectorTools.cpp -std=c++17 

stratifiedMCTemplates.o:  stratifiedMCTemplates.cpp importanceSamplingAndVectorTools.cpp stratifiedMCTemplates.hpp importanceSamplingAndVectorTools.hpp 
	$(CC) -c stratifiedMCTemplates.cpp importanceSamplingAndVectorTools.cpp -std=c++17 

asianOptionNew.o: asianOptionNew.cpp importanceSamplingAndVectorTools.cpp asa241.cpp stratifiedMCTemplates.cpp stratifiedMCTemplates.hpp importanceSamplingAndVectorTools.hpp asa241.hpp 
	$(CC) -c asianOptionNew.cpp importanceSamplingAndVectorTools.cpp  asa241.cpp stratifiedMCTemplates.cpp -std=c++17 


simpleNormal.o: simpleNormal.cpp importanceSamplingAndVectorTools.cpp asa241.cpp stratifiedMCTemplates.cpp stratifiedMCTemplates.hpp importanceSamplingAndVectorTools.hpp asa241.hpp
	$(CC) -c simpleNormal.cpp importanceSamplingAndVectorTools.cpp asa241.cpp stratifiedMCTemplates.cpp -std=c++17



simpleNormal: simpleNormal.o stratifiedMCTemplates.o importanceSamplingAndVectorTools.o asa241.o
	$(CC) simpleNormal.o stratifiedMCTemplates.o importanceSamplingAndVectorTools.o asa241.o -o simpleNormal -std=c++17
	
asianOptionNew: asianOptionNew.o stratifiedMCTemplates.o importanceSamplingAndVectorTools.o asa241.o 
	$(CC) asianOptionNew.o stratifiedMCTemplates.o importanceSamplingAndVectorTools.o asa241.o  -o asianOptionNew -std=c++17

clean: 
	rm -f asianOptionNew simpleNormal
	rm -f *.o
	
