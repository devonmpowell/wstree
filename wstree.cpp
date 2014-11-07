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
#include <cmath>
#include "HDF_IO.hh"

// custom macros and typedefs
using namespace std;
typedef uint32_t luint;
#define LUINT_MAX UINT32_MAX

// a POD struct for watershed info
typedef struct {
	//luint subzone0, subzone1;
	//vector<luint> subws;
	luint parent;
	double vol;
	double mass;
	double fmin;
	double barrier;
	//luint core; // flat index of zone minimum
	luint depth;
} Watershed;

// forward declarations
void read_hdf5(string filename, string fieldname);
void write_hdf5(string filename, string fieldname);
void argsort();
void watershed();
void tree();
void write_tree(string filename);

// global data arrays
luint nx, ny, nz, ntot;
luint nzones;
vector<double> field;
vector<luint> inds_sorted;
vector<luint> zones; 
vector<Watershed> watersheds;

int main(int argc, char **argv) {

	setbuf(stdout, NULL);

	// we want four arguments
	if(argc != 5) {
		printf("-------------------------------------\n");
		printf(" Usage:\n");
		printf("   ./wstree input.hdf5 input_field output.hdf5 output.txt\n");
		printf(" Example:\n");
		printf("   ./wstree data/dset128.hdf5 RHO output/ws128.hdf5 output/tree128.txt\n");
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
	write_hdf5(argv[3], "WS");
	printf(" done.\n");
	printf("-------------------------------------\n");

	printf(" Building heirarchy...");
	tree();

	//for(luint z = 0; z < 20; ++z) {
		//printf(" Watershed %u:\n", z);
		//printf("   vol = %f\tmass = %f\tfavg = %f\tfmin = %f\n", watersheds[z].vol, watersheds[z].mass, watersheds[z].favg, watersheds[z].fmin);
	//}

	write_tree(argv[4]);

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
		fprintf(file, "# parent\ttdepth\tvol\tmass\tfmin\tbar\n");

		for(luint z = 0; z < nzones; ++z) {
			Watershed ws = watersheds[z];
			fprintf(file, "%u\t%u\t%lf\t%lf\t%lf\t%lf\n",
					ws.parent, ws.depth, ws.vol, ws.mass, ws.fmin, ws.barrier);
		}

		fclose(file);
	}
	else {
		printf("\n\tFailed to open file.\n");
	}

}
/*
// recursively floods a region of the box
void rfill(luint ix, luint iy, luint iz, luint zsearch, luint zfill, luint depth) {

	for(int ox = -1; ox <= 1; ++ox) {
		for(int oy = -1; oy <= 1; ++oy) {
			for(int oz = -1; oz <= 1; ++oz) {
	
				luint ixn = (ix + ox + nx)%nx;
				luint iyn = (iy + oy + ny)%ny;
				luint izn = (iz + oz + nz)%nz;
				luint tmp_flat = ixn*ny*nz +iyn*nz + izn;

				//printf("%u ", depth);
	
				if(ox == 0 && oy == 0 && oz == 0) {
					zones[tmp_flat] = zfill;
				}
				else if(zones[tmp_flat] == zsearch)	 {
					rfill(ixn, iyn, izn, zsearch, zfill, depth + 1);
				}
			}
		
		}
	}
	return;
}*/


void tree() {

	// fill in leaf information
	watersheds.resize(nzones);
	for(luint z = 0; z < nzones; ++z) {
		//watersheds[z].subws.push_back(z);
		watersheds[z].parent = LUINT_MAX;
		watersheds[z].vol = 0.0;
		watersheds[z].mass = 0.0;
		watersheds[z].fmin = DBL_MAX;
		watersheds[z].barrier = DBL_MAX;
	}
	for(luint i = 0; i < ntot; ++i) {
		luint z = zones[i];
		double f = field[i];
		// TODO: use physical units for volume and mass
		watersheds[z].vol += 1.0;
		watersheds[z].mass += f;
		if(f < watersheds[z].fmin)
			watersheds[z].fmin = f;
	}
	
	
	// find barrier saddle points and merge zones
	for(luint ind_uns = 0; ind_uns < ntot; ++ind_uns) {


		// get the flattened array index
		luint ind_flat = inds_sorted[ind_uns];
		luint z = zones[ind_flat];
		double f0 = field[ind_flat];

		// get 3D indices from ind_flat
		luint ix = ind_flat/(ny*nz);
		luint iy = (ind_flat - ix*ny*nz)/nz;
		luint iz = ind_flat - ix*ny*nz - iy*nz;

		// inspect the 26 neighboring cells
		double fnmin = DBL_MAX;
		luint znmin = LUINT_MAX;
		double grad_min = DBL_MAX; 
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					luint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					luint z_neighbor = zones[tmp_flat];
					double f_neighbor = field[tmp_flat];

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f_neighbor - f0)/sqrt(ox*ox + oy*oy + oz*oz);

					// check for a neighboring watershed with the shallowest gradient
					if(z != z_neighbor && grad < grad_min) {	
						grad_min = grad;
						fnmin = f_neighbor;
						znmin = z_neighbor;
					}
				}
			}
		}
		if(znmin != LUINT_MAX) {
			// we have found a saddle point
			
			// find the global parent void of z0
			luint z0 = z;
			luint gp0 = z0; 
			while(watersheds[gp0].parent < LUINT_MAX) {
				gp0 = watersheds[gp0].parent;
			}

			// find the global parent void of z1
			luint z1 = znmin;
			luint gp1 = z1; 
			while(watersheds[gp1].parent < LUINT_MAX) {
				gp1 = watersheds[gp1].parent;
			}

			if(gp0 == gp1) continue; // this saddle point is between voids that have already merged

			if(watersheds[gp0].fmin < watersheds[gp1].fmin) {
				watersheds[gp1].barrier = fnmin;
				watersheds[gp1].parent = gp0;
				watersheds[gp0].vol += watersheds[gp1].vol;
				watersheds[gp0].mass += watersheds[gp1].mass;
			} else {
				watersheds[gp0].barrier = fnmin;
				watersheds[gp0].parent = gp1;
				watersheds[gp1].vol += watersheds[gp0].vol;
				watersheds[gp1].mass += watersheds[gp0].mass;

			}
		}
	}

	// post-processing to get tree depth 
	for(luint z = 0; z < nzones; ++z) {
		luint depth = 0;
		luint parent = watersheds[z].parent;

		while(parent < LUINT_MAX) {
			//watersheds[parent].vol += watersheds[z].vol;
			//watersheds[parent].mass += watersheds[z].mass;
			parent = watersheds[parent].parent;
			++depth;
		}
		watersheds[z].depth = depth;
	}

	return;
}


void watershed() {

	nzones = 0;
	zones.assign(ntot, LUINT_MAX); // unassigned zones use LUINT_MAX

	for(luint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		luint ind_flat = inds_sorted[ind_uns];

		// get 3D indices from ind_flat
		luint ix = ind_flat/(ny*nz);
		luint iy = (ind_flat - ix*ny*nz)/nz;
		luint iz = ind_flat - ix*ny*nz - iy*nz;

		// iterate over the 26 neighboring cells
		double f0 = field[ind_flat];
		double grad_max = 0.0; 
		luint zmin = zones[ind_flat];// = LUINT_MAX;
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					luint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f0 - field[tmp_flat])/sqrt(ox*ox + oy*oy + oz*oz);

					if(grad > grad_max) { // finds the largest gradient to a neighboring cell 
						grad_max = grad;
						zmin = zones[tmp_flat];
						zones[ind_flat] = zmin;
					}
				}
			}
		}
		if(zmin == LUINT_MAX) {
			//Watershed wsn;
			//wsn.parent = LUINT_MAX;
			//wsn.core = ind_flat;
			//wsn.vol = 1.0;
			//wsn.mass = f0;
			//wsn.fmin = f0;
			//wsn.barrier = DBL_MAX;
			//watersheds.push_back(wsn);
			zones[ind_flat] = nzones++;
		}
		//else {
			//watersheds[zmin].vol += 1.0;
			//watersheds[zmin].mass += f0;
		//}
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
	luint dims[3] = {nx, ny, nz};
	HDFWriteDataset3D(filename, fieldname, dims, zones);
	return;
}

void argsort() {
	inds_sorted.resize(ntot);
	luint itmp = 0;
	generate(inds_sorted.begin(), inds_sorted.end(), [&] { return itmp++; });
	sort(inds_sorted.begin(), inds_sorted.end(), [&](luint a, luint b) { return field[a] < field[b]; });
	return;
}

