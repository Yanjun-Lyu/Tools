# Default build all
.PHONY: all
all: den_profile coord_num_profile avg_density sp3OP OP test coord_num_profile_avg

# Compiler
CC = g++
# Compiler flags
CFLAGS = -std=c++17 -O3 -g

# Source files
COMMON_SRCS = trajectory.cpp helpers.cpp 
SPECIAL_SRCS = coord_num_profile.cpp den_profile.cpp avg_density.cpp sp3OP.cpp OP.cpp test.cpp coord_num_profile_avg.cpp

# Object files (specify file names from source files with pattern substitution)
COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)
SPECIAL_OBJS = $(SPECIAL_SRCS:.cpp=.o)

# Header files
COMMON_HEADERS = helpers.hpp trajectory.hpp



# Compile COMMON_SRCS files to COMMON_OBJS files without linking
trajectory.o: trajectory.cpp trajectory.hpp
	$(CC) $(CFLAGS) -c $< -o $@

helpers.o: helpers.cpp helpers.hpp
	$(CC) $(CFLAGS) -c $< -o $@


# Compile SPECIAL_SRCS files to SPECIAL_OBJS files without linking
den_profile.o: den_profile.cpp $(COMMON_HEADERS) 
	$(CC) $(CFLAGS) -c $< -o $@

coord_num_profile.o: coord_num_profile.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

avg_density.o: avg_density.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

sp3OP.o: sp3OP.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

OP.o: OP.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

test.o: test.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

coord_num_profile_avg.o: coord_num_profile_avg.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@


# Build executables (linking)
coord_num_profile: coord_num_profile.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@
	
den_profile: den_profile.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

avg_density: avg_density.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

sp3OP: sp3OP.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

OP: OP.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

test: test.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

coord_num_profile_avg: coord_num_profile_avg.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Clean rule
clean:
	rm -f $(COMMON_OBJS) $(SPECIAL_OBJS) den_profile coord_num_profile avg_density sp3OP OP test coord_num_profile_avg
