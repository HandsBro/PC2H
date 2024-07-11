CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -O3
IFLAGS = -Iinclude

OBJ = main.o Graph.o Timer.o Index.o

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o 2hop $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ) 2hop
