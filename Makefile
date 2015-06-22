CC:= g++
BIN:=./bin
SRC:= $(wildcard src/*.cpp) $(wildcard src/tinyxml/*.cpp)
OBJECT_FILES:=$(SRC:$(notdir %.cpp)=$(BIN)/$(notdir %.o))
EXECUTABLES:=$(BIN)/raytracer
INCLUDES:= -I ./tinyxml/ -I .
LDFLAGS:= -v -std=c++11  -stdlib=libstdc++  -pthread
CXXFLAGS:= -std=c++0x -Wall -c ${INCLDUES} -g

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJECT_FILES)
	@$(CC) $(LDFLAGS) -o $@ $^
	@echo "Build Successful"

$(OBJECT_FILES): $(BIN)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CXXFLAGS) -o $@ $<


clean:
	rm -f -r bin
