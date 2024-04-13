#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <stdlib.h>
#include <ctime>
using namespace std;

double integrand(double c1, double c2, double c3, double c4, double c5, double xi, double x) {
	double change = x - xi;
	double value_inner =
		c1 * pow(change, 4)
		+ c2 * pow(change, 3)
		+ c3 * pow(change, 2)
		+ c4 * change
		+ c5;

	double proper_value = sqrt(value_inner); //Square root
	return proper_value;
}


int main() {
	double arc_length = 0, real_length;
	int n; // Number of Data points

	// Welcoming Prompt
	cout << "----------------------------------------------------\nPath Analyser\n----------------------------------------------------\n";
	cout << "Welcome to the Path Analyser tool, this program will take in the data from your text file and output:\n\n";
	cout << "1. Spline - A piecewise function that represents the path based on the sample points given\n";
	cout << "2. Path length - Calculated from the length of each component function of the spline\n";
	cout << "3. Curvature Equations - A piecewise function that can be used to calculate the curvature of any point in the spline\n";
	cout << "4. Curvature Calculator - You may input the point of interest into the terminal to calculate the curvature at that point\n";
	cout << "\n\nRequired input\n";
	cout << "Real road length (km): ";
	cin >> real_length;
	cout << "Number of data points: ";
	cin >> n;

	Eigen::VectorXd a(n), b(n), d(n), c(n), xi(n), m(n), h(n - 1), length(n - 1);
	Eigen::MatrixXd A(n, n), x(n, 2);//Coordinate 
	clock_t start, end;
	double runtime;
	start = clock();// Check time at the start of writing into file


	//Data initialisation
	ifstream inFile;
	inFile.open("data.txt");
	if (inFile.fail()) {
		cerr << "Output file could not be opened" << endl;
		exit(-1);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2; j++) {
			inFile >> x(i, j);
		}
	}
	inFile.close();
	for (int i = 0; i < n; i++) {
		xi(i) = x(i, 0);
		a(i) = x(i, 1);
	}

	// Solving for Coefficients of Natural Cubic Spline
	//Step Size
	for (int i = 0; i < n - 1; i++) {
		h[i] = xi[i + 1] - xi[i];
	}
	//Initialise matrix A
	for (int i = 0; i < n; i++) { //row
		for (int j = 0; j < n; j++) { //column
			if (i == 0 && j == 0) {
				A(i, j) = 1;
			}
			else if (i == n - 1 && j == n - 1) {
				A(i, j) = 1;
			}
			else {
				if (j == i - 1 && i != 0 && i != n - 1) {
					A(i, j) = h[i - 1];
				}
				else if (j == i && i != 0 && i != n - 1) {
					A(i, j) = 2 * (h[i - 1] + h[i]);
				}
				else if (j == i + 1 && i != 0 && i != n - 1) {
					A(i, j) = h[i];
				}
				else {
					A(i, j) = 0;
				}
			}
		}
	}
	//Initialise matrix m
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			m(0) = 0;
		}
		else if (i == n - 1) {
			m(i) = 0;
		}
		else {
			m(i) = ((3 / h[i]) * (a[i + 1] - a[i])) - ((3 / h[i - 1]) * (a[i] - a[i - 1]));
		}
	}
	// Solve for c 
	c = A.lu().solve(m);
	// Solve for b
	for (int i = 0; i < n - 1; i++) {
		b[i] = (1 / h[i]) * (a[i + 1] - a[i]) - (h[i] / 3) * (2 * c(i) + c(i + 1));
	}
	// Solve for d
	for (int i = 0; i < n - 1; i++) {
		d[i] = (1 / (3 * h[i])) * (c(i + 1) - c(i));
	}


	//Arc Length Simspons Rule
	for (int i = 0; i < n - 1; i++) {
		long double c1, c2, c3, c4, c5;
		c1 = 9 * pow(d[i], 2);
		c2 = 12 * d[i] * c[i];
		c3 = 6 * d[i] * b[i] + 4 * pow(c[i], 2);
		c4 = 4 * c[i] * b[i];
		c5 = pow(b[i], 2) + 1;

		double y0, y1, y2, x0 = xi[i], x2 = xi[i + 1], x1 = (x0 + x2) / 2;
		y0 = integrand(c1, c2, c3, c4, c5, x0, x0);
		y1 = integrand(c1, c2, c3, c4, c5, x0, x1);
		y2 = integrand(c1, c2, c3, c4, c5, x0, x2);


		length[i] = ((h[i] / 2) / 3) * (y0 + 4 * y1 + y2);
	}
	for (int i = 0; i < n - 1; i++) {
		arc_length += length[i];
	}

	//Curvature
	Eigen::MatrixXd coefficients(n - 1, 7);
	for (int i = 0; i < n - 1; i++) {
		coefficients(i, 0) = 6 * d[i];
		coefficients(i, 1) = 2 * c[i];
		coefficients(i, 2) = 9 * pow(d[i], 2);
		coefficients(i, 3) = 12 * d[i] * c[i];
		coefficients(i, 4) = 6 * d[i] * b[i] + 4 * pow(c[i], 2);
		coefficients(i, 5) = 4 * c[i] * b[i];
		coefficients(i, 6) = pow(b[i], 2) + 1;

	}

	// Write into text file
	ofstream outFile;
	outFile.open("file.txt");
	if (outFile.fail()) {
		cerr << "Output file could not be opened" << endl;
		exit(-1);
	}
	outFile << "Spline Equations: \n-----------------------------------------------------------------------\n";
	for (int i = 0; i < n - 1; i++) {
		outFile << "{" << xi[i] << " <= x <= " << xi[i + 1] << ": ";

		// Checking for 0 and writing into file.txt
		if (xi[i] == 0) {
			if (d[i] != 0)
				outFile << "( " << (d[i]) << " ) " << "(x) ^ 3 + ";
			if (c[i] != 0)
				outFile << "( " << (c[i]) << " ) " << "(x) ^ 2 + ";
			if (b[i] != 0)
				outFile << "( " << (b[i]) << " ) " << "(x) +";
			if (a[i] != 0)
				outFile << "( " << a[i] << " ) " << " }\n";
		}
		else {
			if (d[i] != 0)
				outFile << "( " << (d[i]) << " ) " << "(x-" << xi(i) << ") ^ 3 + ";
			if (c[i] != 0)
				outFile << "( " << (c[i]) << " ) " << "(x-" << xi[i] << ") ^ 2 + ";
			if (b[i] != 0)
				outFile << "( " << (b[i]) << " ) " << "(x-" << xi[i] << ") +";
			if (a[i] != 0)
				outFile << "( " << a[i] << " ) " << " }\n";
		}
		outFile << "Path Length: " << length[i] << " km" << "\n\n";
	}
	outFile << "-----------------------------------------------------------------------\nTotal Path Length: " << arc_length << " km" << "\n\n";
	outFile << "-----------------------------------------------------------------------\nCurvature: \n";
	for (int i = 0; i < n - 1; i++) {
		outFile << "\nCurvature Equation: ";
		outFile << "{" << xi[i] << " <= x <= " << xi[i + 1] << " = " << "\n|";

		// Checking for 0 and writing into file.txt
		if (xi[i] == 0) {
			if (coefficients(i, 0) != 0) {
				outFile << "( " << coefficients(i, 0) << " ) " << "(x)";
			}
			if (coefficients(i, 1) != 0) {
				outFile << " + " << "( " << coefficients(i, 1) << " ) " << "| / " << "[";
			}
			else
				outFile << ")| / " << "[";
			for (int n = 4, j = 2; j < 6; n--, j++) {
				if (coefficients(i, j) != 0) {
					outFile << "( " << coefficients(i, j) << " ) " << "(x)^" << n << " + ";
				}
			}
			outFile << coefficients(i, 6) << "] ^ 3/2\n";
			outFile << "\n";
		}
		else {
			if (coefficients(i, 0) != 0) {
				outFile << "( " << coefficients(i, 0) << " ) " << "(x-" << xi[i] << ")";
			}
			if (coefficients(i, 1) != 0) {
				outFile << " + " << "( " << coefficients(i, 1) << " ) " << "| / " << "[";
			}
			else
				outFile << ")| / " << "[";
			for (int n = 4, j = 2; j < 6; n--, j++) {
				if (coefficients(i, j) != 0) {
					outFile << "( " << coefficients(i, j) << " ) " << "(x-" << xi[i] << ")^" << n << " + ";
				}
			}
			outFile << coefficients(i, 6) << "] ^ 3/2\n";
			outFile << "\n";
		}

	}
	outFile.close();


	end = clock(); // Check time at the end of writing data into file
	runtime = ((end - start) / double(CLOCKS_PER_SEC)); //Calculate runtime
	cout << "\n----------------------------------------------------\nMetric:\n----------------------------------------------------\n";
	cout << "1.Path Length: " << arc_length << " km\n";
	cout << "2.Run Time:  " << runtime << " seconds \n";
	cout << "3.Absolute error: " << abs(real_length - arc_length) << " km\n";
	cout << "4.Relative error: " << (abs(real_length - arc_length) / real_length) * 100 << "%\n";


	double xval = 0, curvatureX = 0, numerator = 0, denominator = 0, check = 0;
	char repeat = 'y';
	cout << "\n----------------------------------------------------\nCurvature Calculator:\n----------------------------------------------------\n";
	do {
		check = 0;
		cout << "Please enter x value : ";
		cin >> xval;
		for (int i = 0; i < n - 1; i++) {
			// Check which polynomial to use, which interval x is in
			if (xval >= xi[i] && xval <= xi[i + 1]) {

				//Absolute Value 
				if ((coefficients(i, 0) * (xval - xi[i]) + coefficients(i, 1) < 0)) {
					numerator = -1 * (coefficients(i, 0) * (xval - xi[i]) + coefficients(i, 1));
				}
				else {
					numerator = (coefficients(i, 0) * (xval - xi[i]) + coefficients(i, 1)); // |f''(x)|
				}

				denominator = (coefficients(i, 2) * pow((xval - xi[i]), 4) // [1+(f'(x))^2]
					+ coefficients(i, 3) * pow((xval - xi[i]), 3)
					+ coefficients(i, 4) * pow((xval - xi[i]), 2)
					+ coefficients(i, 5) * (xval - xi[i])
					+ coefficients(i, 6));

				curvatureX = (numerator) / pow(denominator, 1.5);
				cout << "Curvature at x = " << xval << " is " << curvatureX << "\n\n";
				cout << "Would you like to calculate another point (Y/N)?: ";
				cin >> repeat;
				check = 1;
				break;
			}
		}
		if (check != 1) {
			cout << "Value is outside of interval, please try again\n\n";
		}

	} while (repeat == 'Y' || repeat == 'y');



}


