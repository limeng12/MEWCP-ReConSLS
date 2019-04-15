CC = g++
#CFLAGS = -std=c++11 -DNDEBUG -O3 -static
CFLAGS = -std=c++11 -DNDEBUG -O3 -pg
target = ReConSLS
all : $(target)

$(target) : main.cpp  
	$(CC) $(CFLAGS) main.cpp local_search.cpp -o $(target)

clean :
	-rm $(target)
