#include "simplex.h"
#define N 10
void simplex::simplexAlhorithm(Data* data, ostream& out) {

	printSimplexTable(out, data);

	int row, col;
	while (!checkIsOver(data)) {

		findLeadingElement(row, col, data);
		if (data->MFT[row] == INFINITY) data->unresolved = true;
		iterate(data, row, col);
		printSimplexTable(out, data);

	}
}
//returns the index of minimun in vector<>
int simplex::findMinInTheVector(const vector<float>& vec){
	int index = 0;
	float min = vec[0];
	int size = vec.size();
	for (int i = 1; i < size; i++)
		if (vec[i] < min) {
			min = vec[i];
			index = i;
		}
	return index;
}
//finds element with min simplex coefficient
void simplex::findLeadingElement(int& row, int& col, Data* data) {

	col = findMinInTheVector(data->coefficients);
	fillMFT(col, data);
	row = findMinInTheVector(data->MFT);
	data->base[row] = col + 1;
}
//fills simplex coefficients vector
void simplex::fillMFT(int n, Data* data) {
	data->MFT.assign(data->system.size(), 0.0);
	for (int i = 0; i < data->MFT.size(); i++) {
		if (data->system[i][n] > 0)
			data->MFT[i] = data->freeMembers[i] / data->system[i][n];
		else data->MFT[i] = INFINITY;
	}
}
//fills simplex coefficients vector in dual simplex method
void simplex::dualFillMFT(int n, Data* data) {
	data->MFT.assign(data->coefficients.size(), 0.0);

	for (int i = 0; i < data->MFT.size(); i++) {
		if (data->system[n][i] < 0 && abs(data->coefficients[i]) > 0)
			data->MFT[i] = abs(data->coefficients[i]) / abs(data->system[n][i]);
		else data->MFT[i] = INFINITY;
	}

}
void simplex::iterate(Data* data, int row, int col) {

	int s = data->system.size();
	int r = data->system[0].size();
	float lead = data->system[row][col];
	float curr;

	//filling simplex matrix
	for (int i = 0; i < s; i++) {
		curr = data->system[i][col];
		if (i != row) { 
			data->freeMembers[i] -= data->freeMembers[row] * curr / lead; 
			if (abs(data->freeMembers[i]) < 0.000001) data->freeMembers[i] = 0;
		}
		for (int j = 0; j < r; j++) {
			if (i == row) continue;
			else if (j == col) data->system[i][j] = 0;
			else data->system[i][j] -= (data->system[row][j] * curr / lead);
		}
	}

	//counting  value of function
	curr = data->coefficients[col];
	data->max -= (data->freeMembers[row] / lead * curr);
	data->freeMembers[row] /= lead;

	//filling coefficients' row
	for (int i = 0; i < r; i++) {
		data->coefficients[i] = (data->coefficients[i] - ((data->system[row][i] / lead) * curr));
	}

	//filling leading row
	for (int i = 0; i < r; i++) {
		data->system[row][i] /= lead;
	}
}
// simplex algorithm isn`t over if any of function coefficients a lees then zero
bool simplex::checkIsOver(Data* data) {
	if (data->unresolved) return true;
	for (int i = 0; i < data->coefficients.size(); i++) {

		if (abs(data->coefficients[i]) < 0.000001) data->coefficients[i] = 0; // prevents machine zero  errors

		if (data->coefficients[i] < 0) return false;
	}
	return true;
}

//prints the results into readable form
void simplex::printSimplexTable(ostream& out, Data* data){
	out.setf(ios::fixed);
	out.precision(3);
	
	int row = data->system.size();
	int col = data->system[0].size();

	out.width(N);
	out << "Base";
	out.width(N);
	out << "Plan";
	for (int j = 0; j < col; j++) {
		out.width(N - 1);
		out << "x";
		out << j + 1;
	}
	
	out << endl;
	for (int i = 0; i < row; i++) {
		out.width(N - 1);
		out << "x";
		out << data->base[i];
		out.width(N);
		out << data->freeMembers[i];
		for (int j = 0; j < col; j++) {
			out.width(N);
			out << data->system[i][j];
		}
		out << endl;

	}
	out.width(N);
	out << "F";
	out.width(N);
	out <<data-> max;
	int a = data->coefficients.size();
	for (int i = 0; i < a; i++) {
		out.width(N);
		out << data->coefficients[i];
	}
	out << endl << endl;

}

void simplex::dualSimplexAlhorithm(Data* data, ostream& out){

	printSimplexTable(out, data);
	int row, col;
	while (!dualCheckIsOver(data)) {
		dualFindLeadingElement(row, col, data);
		
		if (data->MFT[col] == INFINITY) data->unresolved = true;
		iterate(data, row, col);
		printSimplexTable(out, data);
	}
}

bool simplex::dualCheckIsOver(Data* data) {
	if (data->unresolved) return true;
	int s = data->freeMembers.size();
	for (int i = 0; i < s; i++) {
		if (data->freeMembers[i] < 0 ) return false;
	}
	return true;
}

void simplex::dualFindLeadingElement(int& row, int& col, Data* data){

	row = findMinInTheVector(data->freeMembers);

	dualFillMFT(row, data);

	col = findMinInTheVector(data->MFT);

	data->base[row] = col + 1;
}


int gomori::findMaxFractionalPart(vector<float> vec){
	unsigned s = vec.size();
	int ind = 0;
	float max = fractionalPart(vec[0]);
	for (unsigned i = 1; i < s; i++) {
		if (fractionalPart(vec[i]) > max) {
			max = fractionalPart(vec[i]); 
			ind = i;
		}
	}
	return ind;
}

float gomori::roundF(float n) {
	return n + 0.00005f; 
}

float gomori::fractionalPart(float n){
	if (abs(n) < 0.00001) return 0;
	
	if(n> 0) return (n - (int)roundF(n));
	if (fractionalPart(abs(n)) < 0.000001f) return 0;
	else return (n - (int)roundF(n - 1));
}

bool gomori::checkIsOver(Data* data) {
	unsigned s = data->freeMembers.size();
	for (unsigned i = 0; i < s; i++) {
		
		if (fractionalPart(data->freeMembers[i]) > 0.0001) 
			return false; 
	}
	return true;
}
// adds new restriction into simplex system 
// the restriction is based on followed line with index n 
void gomori::addRestriction(int n, Data* data){

	data->freeMembers.push_back(-fractionalPart(data->freeMembers[n]));

	vector<float> newRow;
	for (unsigned i = 0; i < data->system[n].size(); i++)
		newRow.push_back( -fractionalPart(data->system[n][i]));

	//adding new element into base
	newRow.push_back(1);  
	data->base.push_back(newRow.size()); // stores the number of base element

	//resizing of simplex system
	data->coefficients.push_back(0);
	for (unsigned i = 0; i < data->system.size(); i++)
		data->system[i].push_back(0);
	data->system.push_back(newRow);
}
void gomori::gomoriAlgorithm(Data* data, ostream& out){
	simplex::simplexAlhorithm(data, out);

	while (!checkIsOver(data)) {
		int n = findMaxFractionalPart(data->freeMembers);
		addRestriction(n, data);
		simplex::dualSimplexAlhorithm(data, out);
	}

}
