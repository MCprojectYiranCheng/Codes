CC = g++

asa241.o: asa241.cpp asa241.hpp
	$(CC) -c asa241.cpp -std=c++17 

gradientDescent.o: gradientDescent.cpp gradientDescent.hpp 
	$(CC) -c gradientDescent.cpp -std=c++17

stratifiedMCTemplates.o:  stratifiedMCTemplates.cpp gradientDescent.cpp stratifiedMCTemplates.hpp gradientDescent.hpp 
	$(CC) -c stratifiedMCTemplates.cpp gradientDescent.cpp -std=c++17 

asianOptionNew.o: asianOptionNew.cpp gradientDescent.cpp asa241.cpp stratifiedMCTemplates.cpp stratifiedMCTemplates.hpp gradientDescent.hpp asa241.hpp 
	$(CC) -c asianOptionNew.cpp gradientDescent.cpp asa241.cpp stratifiedMCTemplates.cpp -std=c++17 


simpleNormal.o: simpleNormal.cpp gradientDescent.cpp asa241.cpp stratifiedMCTemplates.cpp stratifiedMCTemplates.hpp gradientDescent.hpp asa241.hpp
	$(CC) -c simpleNormal.cpp gradientDescent.cpp asa241.cpp stratifiedMCTemplates.cpp -std=c++17



simpleNormal: simpleNormal.o stratifiedMCTemplates.o gradientDescent.o asa241.o
	$(CC) simpleNormal.o stratifiedMCTemplates.o gradientDescent.o asa241.o -o simpleNormal -std=c++17
	
asianOptionNew: asianOptionNew.o stratifiedMCTemplates.o gradientDescent.o asa241.o 
	$(CC) asianOptionNew.o stratifiedMCTemplates.o gradientDescent.o asa241.o  -o asianOptionNew -std=c++17

clean: 
	rm -f asianOptionNew simpleNormal
	rm -f *.o
	
