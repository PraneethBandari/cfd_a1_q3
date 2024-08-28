#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

int main()
{
	double L = 0.02;  // Plate thickness in meters
	int N = 50;       // Number of nodes
	double dx = L / (N - 1);  // Distance between nodes
	double T_left = 100;      // Temperature at the left face
	double T_right = 200;     // Temperature at the right face
	double k = 0.5;           // Thermal conductivity
	double Q = 1000000;       // Heat generation rate in W/mB3
	double accuracy = 1e-6;   // Convergence criterion

	// Generate node positions
	vector<double> x(N);
	for (int i = 0; i < N; i++) {
		x[i] = (L * i) / (N - 1);
	}

	// Analytical solution
	vector<double> T_analytical(N);
	for (int i = 0; i < N; i++) {
		double xi = x[i];
		T_analytical[i] = T_left + ((T_right - T_left) / L + (Q / (2 * k)) * L) * xi - (Q / (2 * k)) * xi * xi;
		cout << "Analytical temperature at node " << i + 1 << " is " << T_analytical[i] << "\n";
	}

	// Initialize FDM temperature
	vector<double> Temp(N, (T_right + T_left) / 2);
	vector<double> T_new(N);

	// Boundary conditions
	T_new[0] = T_left;
	T_new[N - 1] = T_right;

	double max_dt = 1;  // Difference between temperatures
	int iterations = 0;

	// Iterative solution using FDM
	while (max_dt > accuracy) {
		max_dt = 0;
		for (int i = 1; i < (N - 1); i++) {
			T_new[i] = (Temp[i + 1] + Temp[i - 1] + (Q * dx * dx / k)) / 2;
		}

		// Apply boundary conditions
		T_new[0] = T_left;
		T_new[N - 1] = T_right;

		// Check convergence
		for (int i = 0; i < N; i++) {
			double dt = fabs(Temp[i] - T_new[i]);
			if (dt > max_dt) {
				max_dt = dt;
			}
			Temp[i] = T_new[i];
		}

		iterations++;
	}

	cout << "Number of iterations: " << iterations << "\n";

	// Output results
	for (int i = 0; i < N; i++) {
		cout << "Analytical temp at node " << i + 1 << " is " << T_analytical[i] << ", FDM temp is " << Temp[i] << "\n";
	}

	ofstream myfile("fdmq3sol.txt");

	if (myfile.is_open()) {

		for (int i = 0; i < N; i++) {
			myfile<< i + 1 << ":" << Temp[i] << endl;
		}
		myfile.close();
		cout << "Temperature values have been written to fdmq3sol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}
	ofstream outfile("fdmq3analyticalsol.txt");

	if (outfile.is_open()) {

		for (int i = 0; i < N; i++) {
			outfile << i + 1 << ":" << T_analytical[i] << endl;
		}
		outfile.close();
		cout << "Temperature values have been written to fdmq3analyticalsol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}


	return 0;
}
