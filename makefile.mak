#We compile and execute the code :
output: main.cpp Film.cpp Film3D.cpp
	g++ main.cpp Film.cpp Film3D.cpp -o output
	./output

#This line will erase the output file :
clean:
	rm output
