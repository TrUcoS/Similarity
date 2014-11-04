#include <iostream>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <stdlib.h>
#include <time.h>
#include <set>
#include <Eigen>

using namespace Eigen;

// Generates permutation matrix with size according to input (receives integer as input)
// Methodology: Generate unit matrix, swap entries randomly row-wise
PermutationMatrix<Dynamic, Dynamic> permutingIdentityMatrix(int rowCount){
	
	int random;
	srand(time(NULL));

	PermutationMatrix<Dynamic, Dynamic> perm(rowCount);
	perm.setIdentity();
	std::cout << "Insert number for random number generator\nd";
	
	int tmp;
	std::cin >> tmp;
	for (int i = rowCount - 1; i > 0; --i){
		random = (tmp+rand()) % rowCount;
		std::swap(perm.indices()[i], perm.indices()[rand() % (i + 1)]);
	}
	return perm;
}

// Permutes input matrix
MatrixXd permuteMatrixRows(MatrixXd original){

	// Use rows as dimension since rows are to be permuted
	MatrixXd perm = permutingIdentityMatrix(original.rows());
	// Multiplying from the left permutes the rows, from the right permutes the columns
	perm = perm * original;
	return perm;
}

// Returns index of first element in input vector
int rowIndexFirstEl(VectorXd input){ 
	int i = input.size();
	for (i = 0; i < input.size(); ++i){
		if (input[i] != 0){
			return i;
		}
	}
}

// returns minhash vector of a single input matrix
VectorXd minHash(MatrixXd input){
	int size = input.cols();
	VectorXd tmp(size);
	for (int i = 0; i < size; ++i){
		tmp(i) = rowIndexFirstEl(input.col(i));
	}
	VectorXd output = tmp;
	return output;
}

void randomProjection(std::vector<int> x, int dimension){
	
	int d = x.size();
	int k = dimension;

	float **matrix = new float*[k];
	float *matrixValues = new float[k*d];

	srand(time(0));
	double totalSum = 0;
	std::cout << "\n\nRandom projection:\n";
	for (int i = 0; i < k; ++i){
		
		matrix[i] = matrixValues + i*d;

		// fill matrix with random values
		for (int j = 0; j < d; ++j){
			*(matrixValues + j*sizeof(float)) = rand() % 1000;
			totalSum += *(matrixValues + j*sizeof(float));

		}

		// Normalize matrix rows
		for (int j = 0; j < d; ++j){
			*(matrixValues + j*sizeof(float)) /= totalSum;
			std::cout << *(matrixValues + j*sizeof(float)) << " " ;
		}
		// Set total sum back to zero for next row
		totalSum = 0;
		std::cout << "\n";
	}
	std::cout << "\n\n";
}



/*
float jaccardSimilarity(std::vector<int> x, std::vector<int> y){

	float jaccard_similarity = 0;
	std::set<int> a ( x.begin(), x.end() );
	std::set<int> b ( y.begin(), x.end() );

	std::set<int> all;//(a.size()+b.size());
	std::set<int>::iterator it;

	it=std::set_union(a.begin(), a.end(), b.begin(), b.end());

	return jaccard_similarity;
}

*/


// Input: two vectors of arbitrary type
// Output: Jaccard similarity coefficient
// Methodology: Check if element of vector x exists in y
// If yes, increase count
template<typename T>
float jaccardSimilarity(std::vector<T> x, std::vector<T> y){
 
	std::vector<T>::iterator it;
	int count = 0;

	// run for-loop over all elements of first
	for(int i=0; i<x.size();++i){
			
		// search for current element of first vector in second vector
		it = std::find(y.begin(), y.end(), x[i]);

		// if element of first vector does exist in second vector
		// increase count
		if(*it != 0  ){	
			count=count+1;
		}
	}

	int total_size = x.size()+y.size();

	// Calculate jaccard similarity
	float jaccard_similarity = ((float)count)/((float)total_size);
	
	std::cout<<"total_size: " << total_size << "\nintersecting elements: " << count<<"\n";
	
	return jaccard_similarity;
}



int main(){

	// create test vectors and fill with numbers
	std::vector<int> x;
	std::vector<int> y;

	std::vector<char> a;
	std::vector<char> b;

	a.push_back('a');
	a.push_back('b');
	a.push_back('c');
	a.push_back('d');

	b.push_back('a');
	b.push_back('k');

	// Generate seed for random integers
	srand(time(0));

	for (int i = 1; i < 1001; i++){
		x.push_back(i);
		double d = rand() % 10000;
		y.push_back(d);
	}

	std::cout << "Jaccard similarity of vectors is " << jaccardSimilarity(x, y) << "\n";
	std::cout << "Jaccard similarity of char vectors is " << jaccardSimilarity(a, b) << "\n";

	std::vector<int> v;
	v.push_back(1);
	v.push_back(1);
	v.push_back(1);
	v.push_back(1);

	randomProjection(v,5);

	int size = 5;
	MatrixXd m = MatrixXd::Random(size, size);
	Matrix4d m2;
	m2 << 0, 0, 0, 1, 
		  1, 0, 0, 0, 
		  0, 1, 0, 0, 
		  0, 0, 1, 0;

	std::cout << "Original: \n" << m2 << std::endl << std::endl;

	
	int tmp = 1;
	int permutations = 0;
	std::cout << "How many permutations?" << std::endl;
	std::cin >> permutations;

	MatrixXd perm;
	MatrixXd signatureMatrix = MatrixXd::Ones(m2.cols(), permutations);
	VectorXd hashedStuff;
	int i = 0;
	while (permutations > 0){
		i = permutations - 1;
		perm = permuteMatrixRows(m2);

		// Multiplying from the left permutes the rows, from the right permutes the columns
		std::cout << "Permuted:\n" << perm << std::endl << std::endl;
		hashedStuff = minHash(perm);

		std::cout << "Minhash-vector:\n" << hashedStuff << std::endl;
		signatureMatrix.col(i) = hashedStuff;
		--permutations;	
	}

	std::cout << "Minhash signature matrix:\n" << signatureMatrix;
	
	std::cin >> tmp;
		
	return 0;
}