// SimplexMethod.cpp : Defines the entry point for the console application.
//

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "data.h"
#include "simplex.h"
using namespace std;


int main(int argc, char** argv)
{
	ifstream in;
	in.open("input.txt");
	Data d;
	d.readUserData(in);
	in.close();

	ofstream out;
	out.open("output.txt");
	gomori::gomoriAlgorithm(&d, out);
	out.close();
	system("pause");
    
}

