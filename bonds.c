/****************************************************************************/
/* This file deals with bond related properties of nanocrystals */

#include "nc.h"

/****************************************************************************/
// Calculates the bond related statistics of nanocrystals 

void calcBondStatistics(atom *atoms, param params) {
	FILE *pf;
	int i, j, k, pairIndex = 0;
	int nBonds, *nBondedAtoms, *underCoordinatedFlag;
	double *aveBondLength, *bondLengthSTD;
	char atomSymbol1[10], atomSymbol2[10];

    // Make sure all atom types are assigned properly
	params.nAtomTypes = assignAtomTypes(atoms, params.nAtoms);

	// Allocate memory
	aveBondLength = (double *) calloc(params.nAtomTypes*params.nAtomTypes, sizeof(double));
	bondLengthSTD = (double *) calloc(params.nAtomTypes*params.nAtomTypes, sizeof(double));
	nBondedAtoms = (int *) calloc(params.nMaxBonds+1, sizeof(int));
	underCoordinatedFlag = (int *) calloc(params.nAtomTypes, sizeof(int));

	// Open the file that will contain the output/ resulting statistics
	pf = fopen("nc_bonds.dat", "w");
	writeCurrentTime(pf);  
	writeSeparation(pf);
	writeSystemInfo(atoms, params, pf);
	writeSeparation(pf);
	fprintf(pf, "Nanocrystal Bond Statistics [Angstroms]\n\n");

	// Calculate the nearest neighbor lists for each atom and store the
	// positions in atoms[i].neighborPos[j] and index in atoms[i].neighborList[j]
	calcNearestNeighborLists(atoms, params.nMaxBonds, params.nAtoms);
	writeNearestNeighbors(atoms, params, pf);

	// Write out all bonds
	writeAllBondsWithPositions(atoms, params);

	// Determine and print the number of atoms that have a certain number of bonds
	for (i = 0; i < params.nAtoms; i++) {
		nBonds = 0;
		// Count the number of bonds atom i has
		for (j = 0; j < params.nMaxBonds; j++) {
			if (atoms[i].neighborPos[j].mag > EPS) {
				nBonds++;
			}
		}
		// Increment the index of nBondedAtoms for which atom i has that many bonds 
		for (j = 0; j < params.nMaxBonds+1; j++) {
			if (nBonds == j) {
				nBondedAtoms[nBonds]++;
				if (isAPassivationSymbol(atoms[i].symbol)) {
					if (nBonds != 1) {
						underCoordinatedFlag[atoms[i].type]++;
					}
				}
				else if (nBonds != params.nMaxBonds) {
					underCoordinatedFlag[atoms[i].type]++;
				}
				break;
			}
		}
	}
	fprintf(pf, "Number of atoms with specific number of bonds:\n\n");
	for (i = 0; i < params.nMaxBonds+1; i++) {
		fprintf(pf, "Number of atoms with %2d bonds = %d\n", i, nBondedAtoms[i]);
	}
	for (i = 0; i < params.nAtomTypes; i++) {
		if (underCoordinatedFlag[i]) {
			for (j = 0; j < params.nAtoms; j++) { 
				if (atoms[j].type == i) {
					fprintf(pf, "\nThere are %d undercoordinated  %s atoms!\n", underCoordinatedFlag[i], atoms[j].symbol);
					break;
				}
			}
		}
	}
	writeSeparation(pf);

	// Loop over all atom type pairs and calculate the average bond length between
	// each atom type pairs - will be zero if two atom types are never nearest neighbors
	fprintf(pf, "Average bond lengths for all unique atom type pairs:\n\n");
	for (i = 0; i < params.nAtomTypes; i++) {
		for (k = 0; k < params.nAtoms; k++) if (atoms[k].type == i) {
			strcpy(atomSymbol1, atoms[k].symbol);
			break;
		}
		for (j = i; j < params.nAtomTypes; j++) {
			for (k = 0; k < params.nAtoms; k++) if (atoms[k].type == j) {
				strcpy(atomSymbol2, atoms[k].symbol);
				break;
			}
			aveBondLength[pairIndex] = retAverageBondLength(atoms, atomSymbol1, atomSymbol2, params);
			bondLengthSTD[pairIndex] = retSTDBondLength(aveBondLength[pairIndex], atoms, atomSymbol1, atomSymbol2, params);
			if (aveBondLength[pairIndex] > EPS) { // only print bond pairs for which at least 1 bond exists
				fprintf(pf, "%s-%s average bond length = %.4f with std = %.6f\n\n", atomSymbol1, atomSymbol2, 
					aveBondLength[pairIndex], bondLengthSTD[pairIndex]);
			}
			pairIndex++;
		}
	}
	writeSeparation(pf);

	// Close file that contains bond related statistics
	fclose(pf);

	// Free dynamically allocated memory
	free(bondLengthSTD); free(aveBondLength); 
	free(nBondedAtoms);

	return;
}


/****************************************************************************/
// Builds the nearest neighbor list for each atom: 
// It stores the atoms of up to nMaxBonds nearest neighbors in the 
// neighborPos and neighborList fields of the atom struct

void calcNearestNeighborLists(atom *atoms, int nMaxBonds, int nAtoms) {
	int i, j, k, numBonds;
	double idealBondLength = 0.0;
	vector diff;

	for (i = 0; i < nAtoms; i++) {
		if (atoms[i].pos.mag > 900000.0) numBonds = 0; // used when removing atoms
		else {
			numBonds = 0;
			for (j = 0; j < nAtoms; j++) {
				if (i == j) continue;
				diff = retSubtractedVectors(atoms[i].pos, atoms[j].pos);
				idealBondLength = retIdealBondLength(atoms[i].symbol, atoms[j].symbol);
				if (diff.mag < idealBondLength*1.1) {
					atoms[i].neighborPos[numBonds] = diff;
					atoms[i].neighborList[numBonds] = j;
					numBonds++;
				}
				if (numBonds > nMaxBonds) { // error check
					writeSeparation(stdout);
					printf("PROGRAM EXITING: An atom has more than %d nearest neighbors!!\n", nMaxBonds);
					writeSeparation(stdout);
					exit(EXIT_FAILURE);
				}
			}
		}
		if (numBonds < nMaxBonds) for (k = numBonds; k < nMaxBonds; k++) { // zero rest of neighbors
			atoms[i].neighborPos[k] = retZeroVector();
			atoms[i].neighborList[k] = -1; // not an atom index			
		}
	}

	return;
}

/****************************************************************************/
// Returns 1 (True) if atoms are nearest neighbors or 0 (False) if they are not

int areAtomsNearestNeighbors(atom atom1, atom atom2) {
	vector diff;
	double idealBondLength;
	int nearestNeighbors;

	diff = retSubtractedVectors(atom1.pos, atom2.pos);
	idealBondLength = retIdealBondLength(atom1.symbol, atom2.symbol);
	if (diff.mag < idealBondLength*1.1) nearestNeighbors = 1;
	else nearestNeighbors = 0;

	return nearestNeighbors;
}

/****************************************************************************/
// Calculates the average bond length for a given atom type pair

double retAverageBondLength(atom *atoms, char *atom1, char *atom2, param params) {
	int i, k, nBonds = 0;
	double totalBondLength = 0.0;

	for (i = 0; i < params.nAtoms; i++) {
		if (! strcmp(atoms[i].symbol, atom1)) {
			for (k = 0; k < params.nMaxBonds; k++) {
				if (atoms[i].neighborPos[k].mag > EPS) {
					if (! strcmp(atoms[atoms[i].neighborList[k]].symbol, atom2)) {
						totalBondLength += atoms[i].neighborPos[k].mag;
						nBonds++;
					}
				}
			}
		}
	}
	
	// Return average if a nearest neighbor bond existed, else return 0.0
	if (nBonds) return totalBondLength / (double)(nBonds);
	else return 0.0;
}

/****************************************************************************/
// Calculates the standard deviation of the bond lengths for a given atom type pair

double retSTDBondLength(double aveBondLength, atom *atoms, char *atom1, char *atom2, param params) {
	int i, k, nBonds = 0;
	double variance = 0.0;

	for (i = 0; i < params.nAtoms; i++) {
		if (! strcmp(atoms[i].symbol, atom1)) {
			for (k = 0; k < params.nMaxBonds; k++) {
				if (atoms[i].neighborPos[k].mag > EPS) {
					if (! strcmp(atoms[atoms[i].neighborList[k]].symbol, atom2)) {
						variance += sqr(atoms[i].neighborPos[k].mag - aveBondLength);
						nBonds++;
					}
				}
			}
		}
	}
	
	// Return sqrt(variance/(N-1)) if a nearest neighbor bond existed, else return 0.0
	if (nBonds > 1) return sqrt(variance / (double)(nBonds-1));		
	else if (nBonds == 1) return 9999.99; // only 1 bond - return non-sense 
	else return 0.0; 
}

/****************************************************************************/
// Returns 

int allAtomsHaveTwoNeighbors(atom *atoms, int nAtoms, int nMaxBonds) {
	int i, j, numBonds;

	for (i = 0; i < nAtoms; i++) {
		numBonds = nMaxBonds;
		for (j = 0; j < nMaxBonds; j++) if (atoms[i].neighborPos[j].mag < EPS) numBonds -= 1; 
		if (numBonds < 2) return 0; // an atom has less than 2 bonds -> return false/0
	}

	// all neighbors have nMaxBonds nearest neighbors so return true/ 1
	return 1;
}

/****************************************************************************/
// Returns the ideal bond length  of the given atom pair in Angstroms 
// currently Wurtzite crystal structure is assumed!

double retIdealBondLength(char *atom1, char *atom2) {

    // TODO: dipti -- add bonds for InP, and P-P2??
	if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "Se"))) return 2.6326; //3.0; //2.8; //2.6326;
	else if (! strcmp(atom1, "Se") && (! strcmp(atom2, "Cd"))) return 2.6326; //3.0; //2.8; //2.6326;
	else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "S")))  return 2.6326; //2.5292; // 2.7;
	else if (! strcmp(atom1, "S") && (! strcmp(atom2, "Cd")))  return 2.6326; //2.5292; // 2.7;
	else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "Te"))) return 2.81;
	else if (! strcmp(atom1, "Te") && (! strcmp(atom2, "Cd"))) return 2.81;
	else if (! strcmp(atom1, "Zn") && (! strcmp(atom2, "Se"))) return 2.6326; //2.46; // 2.6; 
	else if (! strcmp(atom1, "Se") && (! strcmp(atom2, "Zn"))) return 2.6326; //2.46; // 2.6; 
	else if (! strcmp(atom1, "Zn") && (! strcmp(atom2, "S")))  return 2.6326; //2.34; // 2.6; 
	else if (! strcmp(atom1, "S") && (! strcmp(atom2, "Zn")))  return 2.6326; //2.34; // 2.6; 
	else if (! strcmp(atom1, "Zn") && (! strcmp(atom2, "Te"))) return 2.64;
	else if (! strcmp(atom1, "Te") && (! strcmp(atom2, "Zn"))) return 2.64;
	else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "P1"))) return 2.6326*0.6;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "Cd"))) return 2.6326*0.6;
    else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "P3"))) return 2.6326*0.6;
	else if (! strcmp(atom1, "P3") && (! strcmp(atom2, "Cd"))) return 2.6326*0.6;
	else if (! strcmp(atom1, "Zn") && (! strcmp(atom2, "P1"))) return 2.46*0.6;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "Zn"))) return 2.46*0.6;
    else if (! strcmp(atom1, "Zn") && (! strcmp(atom2, "P3"))) return 2.46*0.6;
	else if (! strcmp(atom1, "P3") && (! strcmp(atom2, "Zn"))) return 2.46*0.6;
	else if (! strcmp(atom1, "Se") && (! strcmp(atom2, "P2"))) return 2.6326*0.5;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "Se"))) return 2.6326*0.5;
	else if (! strcmp(atom1, "S") && (! strcmp(atom2, "P2")))  return 2.5292*0.5;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "S")))  return 2.5292*0.5;
	else if (! strcmp(atom1, "Te") && (! strcmp(atom2, "P2"))) return 2.81*0.5;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "Te"))) return 2.81*0.5;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "P2"))) return EPS;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "P1"))) return EPS;
	else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "Se1"))) return 2.6326;
    else if (! strcmp(atom1, "Se1") && (! strcmp(atom2, "Cd"))) return 2.6326;
	else if (! strcmp(atom1, "Se1") && (! strcmp(atom2, "P2"))) return 2.6326*0.5;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "Se1"))) return 2.6326*0.5;
	else if (! strcmp(atom1, "Cdz") && (! strcmp(atom2, "Sez"))) return 2.62059;
	else if (! strcmp(atom1, "Sez") && (! strcmp(atom2, "Cdz"))) return 2.62059;
	else if (! strcmp(atom1, "Cdz") && (! strcmp(atom2, "P1"))) return 2.62059*0.6;
	else if (! strcmp(atom1, "Sez") && (! strcmp(atom2, "Cd"))) return 2.6326;
	else if (! strcmp(atom1, "Cd") && (! strcmp(atom2, "Sez"))) return 2.6326;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "Cdz"))) return 2.62059*0.6;
	else if (! strcmp(atom1, "Sez") && (! strcmp(atom2, "P2"))) return 2.62059*0.5;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "Sez"))) return 2.62059*0.5;
	else if (! strcmp(atom1, "Ar") && (! strcmp(atom2, "Ar"))) return 3.71655;
	else if (! strcmp(atom1, "In") && (! strcmp(atom2, "As"))) return 2.6228;
	else if (! strcmp(atom1, "As") && (! strcmp(atom2, "In"))) return 2.6228;
        else if (! strcmp(atom1, "In") && (! strcmp(atom2, "P"))) return 2.5228;
        else if (! strcmp(atom1, "P") && (! strcmp(atom2, "In"))) return 2.5228;
	else if (! strcmp(atom1, "Ga") && (! strcmp(atom2, "As"))) return 2.4480;
	else if (! strcmp(atom1, "As") && (! strcmp(atom2, "Ga"))) return 2.4480;
	else if (! strcmp(atom1, "In") && (! strcmp(atom2, "P1"))) return 2.6228*0.6;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "In"))) return 2.6228*0.6;
	else if (! strcmp(atom1, "Ga") && (! strcmp(atom2, "P1"))) return 2.4480*0.5;
	else if (! strcmp(atom1, "P1") && (! strcmp(atom2, "Ga"))) return 2.4480*0.5;
	else if (! strcmp(atom1, "As") && (! strcmp(atom2, "P2"))) return 2.6228*0.6;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "As"))) return 2.6228*0.6;	
	else if (! strcmp(atom1, "P") && (! strcmp(atom2, "P2"))) return 2.5228*0.6;
	else if (! strcmp(atom1, "P2") && (! strcmp(atom2, "P"))) return 2.5228*0.6;	
	else if (! strcmp(atom1, atom2)) return EPS; // currently one atom type cannot bind to itself
	else {
		return EPS;
		writeSeparation(stdout);
		printf("WARNING: An atom type pair was found for which no ideal bond length exists\n");
		printf("%s and %s bond length is not defined!\n", atom1, atom2);
		printf("An ideal bond length can be added to the program in retIdealBondLength in bonds.c\n");
		writeSeparation(stdout);
	}

	return 0.0;
}

/****************************************************************************/
