INCLUDES = ~/libs

LDLIBS = -litpp

% :: %.cpp
	g++ -I $(INCLUDES) $< -o $@ $(LDLIBS)
