CC:= clang++
OBJS:= main.o scene.o threadpool.o viewport.o xmlload.o tinyxml.o tinystr.o tinyxmlparser.o tinyxmlerror.o volumedata.o
INCLUDES:= -I ./tinyxml/ -I .
LDFLAGS:= -v -std=c++11  -stdlib=libstdc++  -framework Carbon -framework OpenGL -framework GLUT
CXXFLAGS:= -std=c++11  -Wall -c ${INCLDUES} -g

all: ${OBJS} 
	${CC} ${LDFLAGS} $^ -o raytracer

%.o: %.cpp
	${CC} ${CXXFLAGS} $*.cpp

tinyxml.o: tinyxml/tinyxml.cpp 
	${CC} ${CXXFLAGS} tinyxml/tinyxml.cpp
tinystr.o: tinyxml/tinystr.cpp 
	${CC} ${CXXFLAGS} tinyxml/tinystr.cpp
tinyxmlparser.o: tinyxml/tinyxmlparser.cpp
	${CC} ${CXXFLAGS} tinyxml/tinyxmlparser.cpp
tinyxmlerror.o: tinyxml/tinyxmlerror.cpp
	${CC} ${CXXFLAGS} tinyxml/tinyxmlerror.cpp


#%.o: $(wildcard tinyxml/*.cpp)
#	${CC} ${CXXFLAGS} $< -o $@

clean:
	rm -f *.o;
	rm raytracer
	