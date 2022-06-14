/****************************************************************************/
/* This is the beginning of this program that calculates properties that 
   are related to the size of a nanocrystal (e.g., volume of the NC) */

#include "nc.h"

/****************************************************************************/
// Initiates the main functions 

int main(int argc, char *argv[]) {
	FILE *pf;
	int i, nNewAtoms, nFinalAtoms;

	/****************************************************************************/
	// Defines the variables that are used throughout this function 
	param params;
	atom *atoms, *cutAtoms, *noDanglingAtoms, *passivatedAtoms;

	/****************************************************************************/ 
	// Reads input parameters that determine the size of the system
	writeCurrentTime(stdout);
	writeSeparation(stdout);
	readInput(&params);
	if (params.buildNewNC && params.buildNewBulk) params.nAtoms = 6000000;

	/****************************************************************************/
	// Allocates the input dependent memory 
	// TODO: clean this up - only need 2 total
	atoms = (atom *) calloc(2*params.nAtoms, sizeof(atom));
	cutAtoms = (atom *) calloc(2*params.nAtoms, sizeof(atom));
	noDanglingAtoms = (atom *) calloc(2*params.nAtoms, sizeof(atom));
	passivatedAtoms = (atom *) calloc(2*params.nAtoms, sizeof(atom));  	
	for (i = 0; i < 2*params.nAtoms; i++) {
		atoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		atoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		cutAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		cutAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		noDanglingAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		noDanglingAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		passivatedAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		passivatedAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));				
	}

	/****************************************************************************/
  	// Read in the starting nanocrystal positions to atoms or build a new nanocrystal 
  	if (params.buildNewNC) {
  		params.nAtoms = buildNewNanocrystal(atoms, params);
  		params.nAtomTypes = assignAtomTypes(atoms, params.nAtoms);
  	}
  	else if (params.nHalfLayersToAdd) {
  		readConf(atoms, &params);
  		params.nAtoms = addAllNewLayers(atoms, params);
  	}
  	else {
  		readConf(atoms, &params);
  	}

  	/****************************************************************************/
  	// Calculate all the original nanocrystals size and bond related statistics 
  	calcSizeStatistics(atoms, params);
	calcBondStatistics(atoms, params);

	/****************************************************************************/
	// Check if any singly bonded semiconductor atoms exist in the 
	// original nanocrystal and remove them and write new conf
	if (params.remDanglingAtoms) {
		nNewAtoms = cutDanglingAtoms(noDanglingAtoms, atoms, params.nMaxBonds, params.nMinBonds, params.nAtoms);
		pf = fopen("removedDanglingAtoms_conf.xyz", "w");
		writeXYZFile(noDanglingAtoms, nNewAtoms, 1.0, pf);
		fclose(pf);
	}
	
	/****************************************************************************/
	// Cut the desired number of layers and write final conf
	for (i = 0; i < params.nLayersToCut; i++) {
		if (i) deepCopyAllAtoms(noDanglingAtoms, cutAtoms, nNewAtoms, params.nMaxBonds);
		nNewAtoms = cutOneLayer(cutAtoms, noDanglingAtoms, nNewAtoms, params.nMaxBonds);
		if (i == (params.nLayersToCut-1)) {
			pf = fopen("cut_conf.xyz", "w");
			writeXYZFile(cutAtoms, nNewAtoms, 1.0, pf);
			fclose(pf);
		} 
	}

	/****************************************************************************/
	// Passivate the nanocyrstal by replacing the outer layer of semiconductor atoms with passivation atoms	
	if (params.passivate) {
		if (params.nLayersToCut) {
			nFinalAtoms = passivateNanocrystal(passivatedAtoms, cutAtoms, nNewAtoms, params.nMaxBonds);
		}
		else if (params.remDanglingAtoms) { 
			nFinalAtoms = passivateNanocrystal(passivatedAtoms, noDanglingAtoms, nNewAtoms, params.nMaxBonds);
		}
		else { 
			nFinalAtoms = passivateNanocrystal(passivatedAtoms, atoms, params.nAtoms, params.nMaxBonds);
		}
		if (params.makeChiralNC) {
			addChiralCenters(passivatedAtoms, nFinalAtoms, params.nMaxBonds, params.chiralNanocrystalParams);
		}
	}

	/****************************************************************************/
	// Update params and calculate size and bond related statistics for the final nanocyrstal
	if (params.passivate) {
		params.nAtoms = nFinalAtoms;
		calcSizeStatistics(passivatedAtoms, params);
		calcBondStatistics(passivatedAtoms, params);
	} 
	else if (params.nLayersToCut) { 
		params.nAtoms = nNewAtoms;
		calcSizeStatistics(cutAtoms, params);
		calcBondStatistics(cutAtoms, params);
	} 
	else if (params.remDanglingAtoms) { 
		params.nAtoms = nNewAtoms;
		calcSizeStatistics(noDanglingAtoms, params);
		calcBondStatistics(noDanglingAtoms, params);
	} 
	else { // Just calculate the size and bond statistics of the original nanocrystal
	  	calcSizeStatistics(atoms, params);
		calcBondStatistics(atoms, params);
	}

  	/****************************************************************************/
  	// Free dynamically allocated memory
	for (i = 0; i < params.nAtoms; i++) {
		free(atoms[i].neighborList); free(atoms[i].neighborPos);
		free(cutAtoms[i].neighborList); free(cutAtoms[i].neighborPos);
		free(noDanglingAtoms[i].neighborList); free(noDanglingAtoms[i].neighborPos);
		free(passivatedAtoms[i].neighborList); free(passivatedAtoms[i].neighborPos);
	}
	free(cutAtoms); free(noDanglingAtoms);
	free(passivatedAtoms); free(atoms); 

	// Print the time and exit the program
	writeSeparation(stdout);
	writeCurrentTime(stdout);

  	return 0;
}

/****************************************************************************/
