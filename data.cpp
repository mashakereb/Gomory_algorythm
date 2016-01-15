#include "data.h"

#define N 10
using std::string;

inline bool isNumber(char ch) {
	if (ch <= '9' && ch >= '0') return  true;
	return false;
}
//gets coefficients from string 
void Data::getCoefficients(string func, vector<float>& coeff)
{
	char ch = 0;
	float curr = 0;
	unsigned length = func.length();
	bool isNeg = false;
	for (unsigned i = 0; i < length; i++) {
		if (func[i] == 'x' && !isNumber(func[i + 1]))continue;
		if (func[i] == '-') isNeg = true;

		if (isNumber(func[i]) || func[i] == 'x') {

			string coef;
			if (func[i] == 'x') coef = "1";
			else {
				coef += func[i];
				while (func[++i] != 'x') {
					coef += func[i];
				}
			}

			if (isNeg) {
				coef = "-" + coef;
				isNeg = false;
			}

			string num;
			num += func[++i];
			while (isNumber(func[++i])) {
				num += func[i];
			}

			int curVar = stoi(num);
			while (coeff.size() + 1 < curVar) {

				coeff.push_back(0.0);
			}
			curr = stof(coef);
			maxCoef = (curr > maxCoef) ? curr : maxCoef;
			coeff.push_back(curr);
		}
	}
	while (coefficients.size() > coeff.size()) coeff.push_back(0.0);
}
//reads information from followed stream
// information must be structured in such form:
// function                f(x) = ax1 + bx2 + cx3 + ...  , where a, b, c... - function's coefficients
// mode                    (max or min)
// system of restrictions  a11x1 + a12x2 + ... + a1nxn <= b1  
//						   a21x1 + a22x2 + ... + a2nxn <= b1
//                         .................................
//                         am1x1 + am2x2 + ... + amnxn <= bm
//
//							END
// instead of "<=" can be used ">=" or "="
void Data::readUserData(std::istream& in){
	//
	string function;
	getline(in, function);
	getCoefficients(function, coefficients);

	string mode;
	getline(in, mode);
	if (mode == "min") isMaximization = false;
	else if (mode == "max") isMaximization = true;


	string st = "";
	getline(in, st);
	//reading system of restrictions
	while (!(in.eof() || st == "END")) {
		string num;
		int i = st.length();
		//the free member separation 
		while (isNumber(st[--i]) || st[i] == '.');
		num = st.substr(i + 1, st.length());
		float n = stof(num);
		while (st[i] != '=') { 
			if (st[i] == '-') n = -n;
			i--; 
		}
		freeMembers.push_back(n);
		//processing sign
		if (st[i - 1] == '<') signOfRestrictions.push_back(1);
		else if (st[i - 1] == '>') signOfRestrictions.push_back(-1);
		else signOfRestrictions.push_back(0);

		vector<float> vect;
		getCoefficients(st.substr(0, i), vect);
		system.push_back(vect);
		getline(in, st);
	}
	maxCoef *= 2;
	convertation();

}
//converts information after reading into standart form
void Data::convertation(){
	// standart mode - minimization
	if (isMaximization) {
		isMaximization = false;
		int s = coefficients.size();
		for (int i = 0; i < s; i++) {
			coefficients[i] = -coefficients[i];
		}
	}

	//standart sign - " <="
	//also here an artificial base is created
	int size = signOfRestrictions.size();
	for (int i = 0; i < size; i++) {
		if (signOfRestrictions[i] == 0) {
			coefficients.push_back(maxCoef);
			base.push_back(i);
			int s = system.size();
			for (int j = 0; j < s; j++) {
				if (j != i)system[j].push_back(0);
				else {
					system[j].push_back(1);
					base.push_back(coefficients.size());
				}
			}
		}
		else if (signOfRestrictions[i] > 0) {
			coefficients.push_back(0);
			coefficients.push_back(maxCoef);
			int s = system.size();
			for (int j = 0; j < s; j++) {
				if (j != i) {
					system[j].push_back(0);
					system[j].push_back(0);
				}
				else {
					system[j].push_back(1);
					system[j].push_back(1);
					base.push_back(coefficients.size());
				}
			}

		}
		else {
			coefficients.push_back(0);
			coefficients.push_back(maxCoef);
			int s = system.size();
			for (int j = 0; j < s; j++) {
				if (j != i) {
					system[j].push_back(0);
					system[j].push_back(0);
				}
				else {
					system[j].push_back(-1);
					system[j].push_back(1);
					base.push_back(coefficients.size());
				}
			}
		}
	}

	//we need the sum of elements in each column
		int s = coefficients.size();
		float free = 0; //the sum of free members
		for (int i = 0; i < freeMembers.size(); i++)
			free += freeMembers[i];

		float* sum = new float[s];
		for (int i = 0; i < s; i++)
			sum[i] = 0;

		for (int i = 0; i < system.size(); i++) {
			for (int j = 0; j < s; j++) {
				sum[j] += system[i][j];
			}
		}
		// after creating an artificial base we need to multiply coefficients
		for (int i = 0; i < s; i++)
			sum[i] *= maxCoef;
		free *= maxCoef;
		
		for (int i = 0; i < s; i++)
			coefficients[i] -= sum[i];

		max = -free; // the function value in this particular point
		MFT.assign(system.size(), 0.0); // vector of simplex coefficients
		
}