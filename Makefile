CXX = g++
# CXXFLAGS = -Wall -Wextra -g -O2 -std=c++17 -pthread -fopenmp
CXXFLAGS = -Wall -Wextra -g -O3 -std=c++17
TARGET = aes

all:
	$(CXX) $(CXXFLAGS) aes.cpp main.cpp -o $(TARGET)