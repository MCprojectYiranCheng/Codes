CC = g++

asa241.o: asa241.cpp asa241.hpp
	$(CC) -c asa241.cpp -std=c++17 

gradientDescent.o: gradientDescent.cpp gradientDescent.hpp 
	$(CC) -c gradientDescent.cpp -std=c++17 

asianOption.o: asianOption.cpp gradientDescent.cpp asa241.cpp gradientDescent.hpp asa241.hpp 
	$(CC) -c asianOption.cpp gradientDescent.cpp asa241.cpp -std=c++17 

simpleNormal.o: simpleNormal.cpp gradientDescent.cpp asa241.cpp gradientDescent.hpp asa241.hpp
	$(CC) -c simpleNormal.cpp gradientDescent.cpp asa241.cpp -std=c++17



simpleNormal: simpleNormal.o gradientDescent.o asa241.o
	$(CC) simpleNormal.o  gradientDescent.o asa241.o  -o simpleNormal -std=c++17
	
asianOption: asianOption.o gradientDescent.o asa241.o
	$(CC) asianOption.o gradientDescent.o asa241.o -o asianOption -std=c++17
	
clean: 
	rm -f asianOption simpleNormal
	rm -f *.o
	
