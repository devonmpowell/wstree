/***************************************************

	wstree.cpp

	A 3D watershed transform implementation with
	emphasis on the heirarchy of voids

	by Devon Powell

***************************************************/

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <functional>
#include "HDF_IO.hh"

using namespace std;

// other
void init();
void release();
void read_hdf5(string filename, string fieldname);
void write_hdf5();

vector<ulong> argsort(const vector<double> &field);

typedef unsigned long ulong;


// data arrays
ulong nx, ny, nz, ntot;
vector<double> field;
vector<ulong> inds_sorted;


int main(int argc, char **argv) {

	// we want two arguments
	/*if(argc != 3) {
		printf("Requires one input file and one output file!\n");
		exit(0);
	}

	printf("\n"); 

	*/

	// set up the simulation and allocate memory
	
	//readInput(argv[1]);
	init(); 
	read_hdf5("data/dset128.hdf5", "RHO");


	printf("Testing argsort code...\n");

	vector<double> unsorted {1.0, 0.9, 0.5, 3.0, 4.0, 1.6, 1.5};
/*	vector<int> indices(unsorted.size());

	int i_tmp = 0;
	generate(indices.begin(), indices.end(), [&] { return i_tmp++; });


	for(int i = 0; i < unsorted.size(); ++i) {
		printf("%d\t%f\n", indices[i], unsorted[i]);
	}
		
	printf("Sorting...\n");
*/
	vector<ulong> inds_sorted = argsort(unsorted);

	
	for(int i = 0; i < unsorted.size(); ++i) {
		printf("%d\t%f\n", inds_sorted[i], unsorted[inds_sorted[i]]);
	}



	write_hdf5();

	// free up pointers
	release();
	
	return 0;    
}


void read_hdf5(string filename, string fieldname) {

	vector<int> dims;

	HDFGetDatasetExtent(filename, fieldname, dims);

	for(int i = 0; i < dims.size(); ++i) {
		printf("%d ", dims[i]);
	}
	printf("\n");


	vector<double> data;

	HDFReadDataset(filename, fieldname, data);
	printf("Read %d data values.\n", data.size());
	printf("Sorting...\n");
	sort(data.begin(), data.end());
	printf("...done.\n");

	for(int i = 0; i < 50; ++i) {
		printf("%f ", data[i]);
	}
	printf("\n");

}

void write_hdf5() {
/*
	char* outfilename = "output/test.hdf5";

	printf("Writing to file: %s...\n", outfilename);

	cudaMemcpy(pos_h, pos_d, npart*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_h, vel_d, npart*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(phi_h, phi_d, ngrid*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(rho_h, rho_d, ngrid*sizeof(float), cudaMemcpyDeviceToHost);
	cudaErrCheck("Copy device arrays to host");

	vector<float> pos_v(pos_h, pos_h + npart);
	vector<float> vel_v(vel_h, vel_h + npart);
	vector<float> phi_v(phi_h, phi_h + ngrid);
	vector<float> rho_v(rho_h, rho_h + ngrid);

	HDFCreateFile(outfilename);
	HDFWriteDataset(outfilename, "POS", pos_v);
	HDFWriteDataset(outfilename, "VEL", vel_v);
	HDFWriteDataset(outfilename, "PHI", phi_v);
	HDFWriteDataset(outfilename, "RHO", rho_v);

	printf("...done.\n");
*/
}

void init() {

	printf("Initializing...\n");


}

void release() {

	//free(pos_h); free(vel_h); free(phi_h); free(rho_h);
	//cudaFree(pos_d); cudaFree(vel_d); cudaFree(phi_d); cudaFree(rho_d);

	return;
}

// helper function for argsort
bool indcmp(int a, int b, vector<double> &field) {
	return field[a] < field[b];
}
vector<ulong> argsort(const vector<double> &field) {
	vector<ulong> indices(field.size());
	int i_tmp = 0;
	generate(indices.begin(), indices.end(), [&] { return i_tmp++; });
	sort(indices.begin(), indices.end(), bind(indcmp, placeholders::_1, placeholders::_2, field));
	return indices;
}

