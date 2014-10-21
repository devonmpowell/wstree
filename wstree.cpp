/***************************************************

	wstree.cpp

	A 3D watershed transform implementation with
	emphasis on the heirarchy of voids

	by Devon Powell

***************************************************/

#include <cstdio>
#include <cmath>
#include <vector>
#include "HDF_IO.hh"

using namespace std;

// other
void init();
void release();
void write_hdf5();

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



//	write_hdf5();

	// free up pointers
	release();
	
	return 0;    
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

	printf("Initializing...");


}

void release() {

	//free(pos_h); free(vel_h); free(phi_h); free(rho_h);
	//cudaFree(pos_d); cudaFree(vel_d); cudaFree(phi_d); cudaFree(rho_d);

	return;
}


