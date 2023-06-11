CXX = g++
CXXFLAGS = -Wall -Wextra -g -O3 -std=c++17 -pthread -fopenmp
TARGET = aes

all:
	$(CXX) $(CXXFLAGS) aes.cpp main.cpp -o $(TARGET)