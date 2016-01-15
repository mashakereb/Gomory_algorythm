#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <set>
using namespace std;

class Data {
	friend class gomori;
	friend class simplex;
public:

	void readUserData(std::istream& in);
	void getCoefficients(string, vector<float>& );
	void convertation();

protected:
	bool isMaximization;
	bool isMixed = false;
	float max = 0, maxCoef = 0;
	vector <float> MFT;
	set <int> notIntegers;
	vector <int> base;

	//coefficients of main function
	vector<float> coefficients;

	//free members in restrictions
	vector<float> freeMembers;

	//matrix of coefficients of restrictions
	vector<vector<float>> system;

	//sign for each restriction ["=" -> 0, "<=" -> 1, ">=" -> -1]
	vector<int> signOfRestrictions;
	bool unresolved = false;

};