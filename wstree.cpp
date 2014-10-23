/***************************************************

	wstree.cpp

	A 3D watershed transform implementation with
	emphasis on the heirarchy of voids

	by Devon Powell

***************************************************/

#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdint>
#include <cfloat>
#include "HDF_IO.hh"

// custom macros and typedefs
using namespace std;
typedef uint32_t uint;
#define uint_max UINT32_MAX

// a POD struct for watershed info
typedef struct {
//	uint zone_id;
	uint subzone_0, subzone_1;
	double vol;
	double mass;
	double fmin;
//	double favg;

	double barrier;
//	uint barrier_neighbor;

} Watershed;

// forward declarations
void read_hdf5(string filename, string fieldname);
void write_hdf5(string filename, string fieldname);
void argsort();
void watershed();
void tree();
void write_tree(string filename);

// global data arrays
uint nx, ny, nz, ntot;
vector<double> field;
vector<uint> inds_sorted;
uint nzones;
vector<uint> zones; 

vector<Watershed> watersheds;

int main(int argc, char **argv) {

	// we want four arguments
	if(argc != 5) {
		printf("-------------------------------------\n");
		printf(" Usage:\n");
		printf("   ./wstree input.hdf5 input_field output.hdf5 output_field\n");
		printf(" Example:\n");
		printf("   ./wstree data/dset128.hdf5 RHO output/ws128.hdf5 WS\n");
		printf("-------------------------------------\n");
		return 0;
	}

	printf("-------------------------------------\n");
	printf(" Reading field %s from file %s...", argv[2], argv[1]);
	read_hdf5(argv[1], argv[2]);
	printf(" done.\n");
	printf("   Read %d data values.\n", (int)field.size());
	printf("   nx = %d, ny = %d, nz = %d for %d total.\n", nx, ny, nz, ntot);
	printf("-------------------------------------\n");
	printf(" Sorting array indices...");
	argsort();
	printf(" done.\n");
	printf("-------------------------------------\n");
	printf(" Running watershed transform...");
	watershed();
	printf(" done.\n");
	printf("   Found %u distinct zones.\n", nzones);
	printf("-------------------------------------\n");
	printf(" Writing field %s to file %s...", argv[4], argv[3]);
	write_hdf5(argv[3], argv[4]);
	printf(" done.\n");
	printf("-------------------------------------\n");

	printf(" Building heirarchy...");
	tree();

	//for(uint z = 0; z < 20; ++z) {
		//printf(" Watershed %u:\n", z);
		//printf("   vol = %f\tmass = %f\tfavg = %f\tfmin = %f\n", watersheds[z].vol, watersheds[z].mass, watersheds[z].favg, watersheds[z].fmin);
	//}
	printf("nzones = %u\n", nzones);

	write_tree("output/tree128.dat");

	printf(" done.\n");
	printf("-------------------------------------\n");

	return 0;    
}

void write_tree(string filename) {

	FILE* file = fopen(filename.c_str(), "w");
	printf("\tWriting to file: %s...", filename.c_str());

	if(file != NULL) {

		fprintf(file, "###########################################################################\n");
		fprintf(file, "#    TREE											                     #\n");
		fprintf(file, "###########################################################################\n");

/*typedef struct {
//	uint zone_id;
	uint subzone_0, subzone_1;
	double vol;
	double mass;
	double fmin;
//	double favg;

	double barrier;
//	uint barrier_neighbor;

} Watershed;
*/


		fprintf(file, "# sub0\tsub1\tvol\tmass\tfmin\tbar\n");

		for(uint z = 0; z < nzones; ++z) {
			Watershed ws = watersheds[z];
			//fprintf(file, "%011u\t%011u\t%.11f\t%.11f\t%.11f\t%.11f\n",
			fprintf(file, "%u\t%u\t%lf\t%lf\t%lf\t%lf\n",
					ws.subzone_0, ws.subzone_1, ws.vol, ws.mass, ws.fmin, ws.barrier);
		}

		fclose(file);
	}
	else {
		printf("\n\tFailed to open file.\n");
	}

}

void tree() {

	// fill in leaf information
	watersheds.resize(nzones);
	for(uint z = 0; z < nzones; ++z) {
		watersheds[z].subzone_0 = uint_max;
		watersheds[z].subzone_1 = uint_max;
		watersheds[z].vol = 0.0;
		watersheds[z].mass = 0.0;
		watersheds[z].fmin = DBL_MAX;
		watersheds[z].barrier = DBL_MAX;
	}
	for(uint i = 0; i < ntot; ++i) {
		uint z = zones[i];
		double f = field[i];
		// TODO: use physical units for volume and mass
		watersheds[z].vol += 1.0;
		watersheds[z].mass += f;
		if(f < watersheds[z].fmin)
			watersheds[z].fmin = f;
	}
//	for(uint z = 0; z < nzones; ++z)
//		watersheds[z].favg = watersheds[z].mass/watersheds[z].vol;
	
	// find barrier saddle points and merge zones
	for(uint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		uint ind_flat = inds_sorted[ind_uns];
		uint z = zones[ind_flat];
		double f0 = field[ind_flat];

		// get 3D indices from ind_flat
		uint ix = ind_flat/(ny*nz);
		uint iy = (ind_flat - ix*ny*nz)/nz;
		uint iz = ind_flat - ix*ny*nz - iy*nz;


		// TODO: test this barrier-finding method more rigorously...

		// inspect the 26 neighboring cells
		double fnmin = DBL_MAX;
		uint znmin = uint_max;
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					uint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					uint z_neighbor = zones[tmp_flat];
					double f_neighbor = field[tmp_flat];

					// check the neighboring watershed
					// TODO: check this...
					if(z != z_neighbor && f_neighbor < fnmin) {	
						fnmin = f_neighbor;
						znmin = z_neighbor;
					}
				}
			}
		}
		// merge the two zones
		if(znmin != uint_max) {

			// create a new watershed from the union of the two
			Watershed ws_new;
			uint z0, z1;
			if(watersheds[z].fmin < watersheds[znmin].fmin) {
				z0 = z; z1 = znmin;
			}
			else {
				z0 = znmin; z1 = z;
			}
			watersheds[z0].barrier = fnmin;
			watersheds[z1].barrier = fnmin;

			ws_new.subzone_0 = z0;
			ws_new.subzone_1 = z1;
			ws_new.vol = watersheds[z0].vol + watersheds[z1].vol;
			ws_new.mass = watersheds[z0].mass + watersheds[z1].mass;
			ws_new.fmin = watersheds[z0].fmin;
			ws_new.barrier = DBL_MAX;//fnmin;

			// flood the two subzones with the new zone
			// TODO: make the loop bounds smarter
			for(uint i = 0; i < ntot; ++i) {
				if(zones[i] == z0 || zones[i] == z1)
					zones[i] = nzones;
			}

			// add the newest parent zone
			watersheds.push_back(ws_new);
			++nzones;
		}
	}

	return;
}


void watershed() {

	nzones = 0;
	zones.assign(ntot, uint_max); // unassigned zones use uint_max

	for(uint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		uint ind_flat = inds_sorted[ind_uns];


		// get 3D indices from ind_flat
		uint ix = ind_flat/(ny*nz);
		uint iy = (ind_flat - ix*ny*nz)/nz;
		uint iz = ind_flat - ix*ny*nz - iy*nz;

		// iterate over the 26 neighboring cells
		double f0 = field[ind_flat];
		//double dmin = f0; 
		double grad_max = 0.0; 
		uint zmin = zones[ind_flat];// = uint_max;
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					uint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f0 - field[tmp_flat])/sqrt(ox*ox + oy*oy + oz*oz);

					//if(zones[tmp_flat] < zmin) { // This disambiguation has an inherent bias towards deeper voids!
					//if(field[tmp_flat] < dmin) { // finds the neighboring cell with the lowest density 
					if(grad > grad_max) { // finds the largest gradient to a neighboring cell 
					
						//dmin = field[tmp_flat];
						grad_max = grad;
						
						zmin = zones[tmp_flat];
						zones[ind_flat] = zmin;
					}
				}
			}
		}
		if(zmin == uint_max) {
			zones[ind_flat] = nzones++;
		}
	}
	return;
}

// read in a multidimensional hdf5 file
void read_hdf5(string filename, string fieldname) {
	vector<int> dims;
	HDFGetDatasetExtent(filename, fieldname, dims);
	nx = dims[0]; ny = dims[1]; nz = dims[2];
	ntot = nx*ny*nz;
	HDFReadDataset(filename, fieldname, field);
	return;
}

void write_hdf5(string filename, string fieldname) {
	HDFCreateFile(filename);
	uint dims[3] = {nx, ny, nz};
	HDFWriteDataset3D(filename, fieldname, dims, zones);
	return;
}

bool indcmp(uint a, uint b) { // helper function for argsort
	return field[a] < field[b];
}
void argsort() {
	inds_sorted.resize(ntot);
	uint itmp = 0;
	generate(inds_sorted.begin(), inds_sorted.end(), [&] { return itmp++; });
	sort(inds_sorted.begin(), inds_sorted.end(), indcmp);
	return;
}

