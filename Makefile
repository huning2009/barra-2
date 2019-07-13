CXX = g++
CXXFLAGS = -std=c++11 -Wall -g

barra.exe: obj/main.o obj/marketdata.o obj/stylefactor.o obj/covariance.o
	${CXX} ${CXXFLAGS} -L/usr/local/lib -lgsl -lgslcblas -o barra.exe obj/main.o obj/marketdata.o obj/stylefactor.o obj/covariance.o
obj/main.o: main.cpp 
	${CXX} ${CXXFLAGS} -c main.cpp -o obj/main.o
obj/marketdata.o: source/marketdata.cpp
	${CXX} ${CXXFLAGS} -c source/marketdata.cpp -o obj/marketdata.o
obj/stylefactor.o: source/stylefactor.cpp
	${CXX} ${CXXFLAGS} -c source/stylefactor.cpp -o obj/stylefactor.o
obj/covariance.o: source/covariance.cpp
	${CXX} ${CXXFLAGS} -c source/covariance.cpp -o obj/covariance.o
clean:
	rm obj/*.o barra.exe


