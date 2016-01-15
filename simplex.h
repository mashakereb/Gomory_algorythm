#pragma once
#include <sstream>
#include "data.h"

class simplex {
	static void printSimplexTable(ostream&, Data*);
	static int findMinInTheVector(const vector<float>&);
	static void iterate(Data*, int, int);

	static void fillMFT(int, Data*);
	static bool checkIsOver(Data*);
	static void findLeadingElement(int&, int&, Data*);

	static void dualFillMFT(int, Data*);
	static bool dualCheckIsOver(Data*);
	static void dualFindLeadingElement(int&, int&, Data*);
public:
	static void simplexAlhorithm(Data*, ostream&);
	static void dualSimplexAlhorithm(Data*, ostream&);

};


class gomori {

	static int findMaxFractionalPart(vector<float>);
	static int findMaxFractionalPart_2(vector<float>, Data*);
	static inline float fractionalPart(float n);
	static void addRestriction(int, Data*);
	static void addRestriction_2(int, Data*);
	static inline float roundF(float);
	static bool checkIsOver(Data*);
	static bool checkIsOver_2(Data*);
public:
	static void gomoriAlgorithm(Data*, ostream&);
	
};