/****************************************************************************/
/* This file builds new nanocrystals */

#include "nc.h"

/****************************************************************************/

int buildNewNanocrystal(atom *atoms, param params) {
	FILE *pf;
	int i, nAtoms = 0;
	atom *bulkAtoms;

	// Dynamically allocate memory
	bulkAtoms = (atom *) calloc(params.nAtoms, sizeof(atom));
	for (i = 0; i < params.nAtoms; i++) {
	 	bulkAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
	 	bulkAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
	}

	// Build or read in the bulk structure that will be cut
	if (params.buildNewBulk) {
		printf("Building a new bulk crystal:\n\n");
		nAtoms = buildBulkCrystal(bulkAtoms, params);
	}
	else {
		printf("Cutting new nanocrystal from the supplied bulk crystal configuration file:\n\n");
		readConf(bulkAtoms, &params);
		nAtoms = params.nAtoms; 
		// Make sure the read in bulk crystal is centered properly
		translateAllAtoms(bulkAtoms, nAtoms, params.mainNC.ncCenter);		
	}
	printf("The total number of bulk crystal atoms = %d\n", nAtoms); fflush(stdout);

	// Cut a rough nanocrystal (QD or NR or NPL or attachedNCs) from a bulk crystal
	nAtoms = cutRoughNanocrystalFromBulk(bulkAtoms, nAtoms, params.mainNC, params);
	printf("The total number of nanocrystal atoms (with dangling atoms) = %d\n", nAtoms); fflush(stdout);

	// Remove the dangling atoms from the rough nanocrystal self-consistently
	nAtoms = cutDanglingAtoms(atoms, bulkAtoms, params.nMaxBonds, params.nMinBonds, nAtoms);
	printf("\nThe total number of nanocrystal atoms (without dangling atoms) = %d\n\n", nAtoms); fflush(stdout);

	// Make a core/shell nanocrystal 
	if (params.nCores) {
		changeAllCoreAtomTypes(atoms, nAtoms, params);
	}

	// Write the final configuration
	pf = fopen("nanocrystal_unpassivated.xyz", "w");
	writeXYZFile(atoms, nAtoms, 1.0, pf);
	fclose(pf);

	// Free dynamically allocated memory
	for (i = 0; i < params.nAtoms; i++) {
	 	free(bulkAtoms[i].neighborList);
	 	free(bulkAtoms[i].neighborPos);
	}
	free(bulkAtoms);

	return nAtoms;
}

/****************************************************************************/

int cutRoughNanocrystalFromBulk(atom *atoms, int nAtoms, nanostructure nc, param params) {
	double skinDiam, skinLen;
	
	if (params.passivate) { 
		skinDiam = 2.0*retIdealBondLength(params.mainNC.ncAtomSymbol1, params.mainNC.ncAtomSymbol2);
		skinLen = 1.0*skinDiam;
	}
	else {
		skinDiam = 0.5*retIdealBondLength(params.mainNC.ncAtomSymbol1, params.mainNC.ncAtomSymbol2);
		skinLen = 2.0*skinDiam;
	}

	if (! strcmp(nc.ncType, "QD")) {
		nAtoms = cutQDFromBulk(atoms, nAtoms, params.nMaxBonds, nc.ncSize.mag+skinDiam);
	}
	else if (! strcmp(nc.ncType, "NR")) { 
		nAtoms = cutNRFromBulk(atoms, nAtoms, params.nMaxBonds, nc.ncSize.mag+skinDiam, nc.ncSize.z+skinLen);
	}
	else if (! strcmp(nc.ncType, "NPL")) {
		nAtoms = cutNPLFromBulk(atoms, nAtoms, params.nMaxBonds, nc.ncSize.x+skinLen, 
				nc.ncSize.y+skinLen, nc.ncSize.z+skinLen);
	}
	else if (! strcmp(nc.ncType, "attachedNCs")) {
		nAtoms = cutAttachedNCsFromBulk(atoms, nAtoms, params);
	}
	else {
		writeSeparation(stdout);
		printf("Invalid ncType! ncType = QD or NR or NPL or attachedNCs\n");
		printf("The program is exiting due to this error!\n\n");
		writeSeparation(stdout);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}

	return nAtoms;
}

/****************************************************************************/
// 1) duplicate bulk atoms (deep copy)
// 2) center one copy to ncOriginQD1
// 3) cutRoughNanocrystal() (might need to have new ncType for hexagonal cylinder)
// 4) cutDanglingAtoms()
// 5) print out configuration of the nanocrystal with it centered to the origin (before translating it back)
// 6) repeat 2-4 for second copy/ second nanocrystal
// 7) translate back both nanocrystals to make their center ncOriginQD[i]
// 8) join their atom lists and remove duplicants 
// 9) print configurations (possibly add in cores to each nanocrystal)

int cutAttachedNCsFromBulk(atom *atoms, int nOrigAtoms, param params) {
	FILE *pf;
	int i, iNC, nFinalAtoms, nCoreAtoms[params.nAttachedNCs]; 
	char fileName[50];
	atom *tmpAtoms, *finalAtoms;

	// Dynamically allocate memory
	tmpAtoms = (atom *) calloc(nOrigAtoms, sizeof(atom));
	finalAtoms = (atom *) calloc(2*nOrigAtoms, sizeof(atom));
	for (i = 0; i < nOrigAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		finalAtoms[i].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		finalAtoms[i].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		finalAtoms[i+nOrigAtoms].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		finalAtoms[i+nOrigAtoms].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
	}

	// Create each nanocrystal one at a time
	nFinalAtoms = 0; 
	for (iNC = 0; iNC < params.nAttachedNCs; iNC++) {
		// Deep copy all the original atoms into tmpAtoms
		deepCopyAllAtoms(tmpAtoms, atoms, nOrigAtoms, params.nMaxBonds);
		// Translate tmpAtoms to the origin of this nanocrystal
		translateAllAtoms(tmpAtoms, nOrigAtoms, params.attachedNC[iNC].ncCenter);
		// Use cutRoughNanocrystalFromBulk function to create specific NC shape and size specified in params.attachedNC[iNC]
		nCoreAtoms[iNC] = cutRoughNanocrystalFromBulk(tmpAtoms, nOrigAtoms, params.attachedNC[iNC], params);
		// Remove the dangling atoms from the rough nanocrystal self-consistently and store results in finalAtoms
		nCoreAtoms[iNC] = cutDanglingAtoms(&(finalAtoms[nFinalAtoms]), tmpAtoms, params.nMaxBonds, params.nMinBonds, nCoreAtoms[iNC]);
		// Write the final configuration
		printf("The number of atoms in attached NC %d = %d\n", iNC, nCoreAtoms[iNC]); fflush(stdout);
		sprintf(fileName, "core_%d_unpass_conf.xyz", iNC);
		pf = fopen(fileName, "w");
		writeXYZFile(&(finalAtoms[nFinalAtoms]), nCoreAtoms[iNC], 1.0, pf);
		fclose(pf);
		// Translate the new finalAtoms back to their original position
		translateAllAtoms(&(finalAtoms[nFinalAtoms]), nCoreAtoms[iNC], retScaledVector(params.attachedNC[iNC].ncCenter, -1.0));
		// Remove repeated atoms from finalAtomList
		nFinalAtoms += nCoreAtoms[iNC];
		nFinalAtoms = removeDuplicatedAtoms(finalAtoms, nFinalAtoms, params.nMaxBonds);
	}

	// Deep copy finalAtoms into the first nFinalAtoms spots of atoms
	calcNearestNeighborLists(finalAtoms, params.nMaxBonds, nFinalAtoms);
	deepCopyAllAtoms(atoms, finalAtoms, nFinalAtoms, params.nMaxBonds);

	// Free dynamically allocated memory
	for (i = 0; i < nOrigAtoms; i++) {
		free(finalAtoms[i].neighborList); free(finalAtoms[i+nOrigAtoms].neighborList);
		free(finalAtoms[i].neighborPos); free(finalAtoms[i+nOrigAtoms].neighborPos);
		free(tmpAtoms[i].neighborPos);
		free(tmpAtoms[i].neighborList);
	}
	free(tmpAtoms); free(finalAtoms);

	return nFinalAtoms;
}

/****************************************************************************/

int cutQDFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double qdDiameter) {
	int i, nNewAtoms = 0;
	atom *tmpAtoms;

	// Dynamically allocate memory
	tmpAtoms = (atom *) calloc(nAtoms, sizeof(atom));
	for (i = 0; i < nAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}

	// Copy atom information to tmpAtoms if within QD circle
	for (i = 0; i < nAtoms; i++) {
		if (retVectorMagnitude(atoms[i].pos) < 0.5*qdDiameter) {
			deepCopySingleAtom(tmpAtoms, nNewAtoms, atoms, i, nMaxBonds);
			nNewAtoms++; 
		} 
		else {
			continue;
		}
	}

	// Deep copy tmpAtoms into the first nNewAtoms spots of atoms
	deepCopyAllAtoms(atoms, tmpAtoms, nNewAtoms, nMaxBonds);

	// Free dynamically allocated memory
	for (i = 0; i < nAtoms; i++) {
		free(tmpAtoms[i].neighborList);
		free(tmpAtoms[i].neighborPos);
	}
	free(tmpAtoms);

	return nNewAtoms;
}

/****************************************************************************/

int cutNRFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double nrDiameter, double nrLength) {
	int i, nNewAtoms = 0;
	double xyProjectionMag;
	atom *tmpAtoms;

	// Dynamically allocate memory
	tmpAtoms = (atom *) calloc(nAtoms, sizeof(atom));
	for (i = 0; i < nAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}

	// Copy atom information to tmpAtoms if within the NR cylinder
	for (i = 0; i < nAtoms; i++) {
		xyProjectionMag = 2.0*sqrt(sqr(atoms[i].pos.x) + sqr(atoms[i].pos.y));
		if (xyProjectionMag < (nrDiameter) && fabs(atoms[i].pos.z) < (0.5*nrLength)) {
			deepCopySingleAtom(tmpAtoms, nNewAtoms, atoms, i, nMaxBonds);
			nNewAtoms++; 
		} 
		else {
			continue;
		}
	}

	// Deep copy tmpAtoms into the first nNewAtoms spots of atoms
	deepCopyAllAtoms(atoms, tmpAtoms, nNewAtoms, nMaxBonds);

	// Free dynamically allocated memory
	for (i = 0; i < nAtoms; i++) {
		free(tmpAtoms[i].neighborList);
		free(tmpAtoms[i].neighborPos);
	}
	free(tmpAtoms);

	return nNewAtoms;
}

/****************************************************************************/

int cutNPLFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double nplLength, double nplWidth, double nplThickness) {
	int i, nNewAtoms = 0;
	atom *tmpAtoms;

	// Dynamically allocate memory
	tmpAtoms = (atom *) calloc(nAtoms, sizeof(atom));
	for (i = 0; i < nAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}

	// Copy atom information to tmpAtoms if within the NPL box
	for (i = 0; i < nAtoms; i++) {
		if (fabs(atoms[i].pos.x) < (0.5*nplLength) && fabs(atoms[i].pos.y) < (0.5*nplWidth)
												&& fabs(atoms[i].pos.z) < (0.5*nplThickness)) {
			deepCopySingleAtom(tmpAtoms, nNewAtoms, atoms, i, nMaxBonds);
			nNewAtoms++; 
		} 
		else {
			continue;
		}
	}

	// Deep copy tmpAtoms into the first nNewAtoms spots of atoms
	deepCopyAllAtoms(atoms, tmpAtoms, nNewAtoms, nMaxBonds);

	// Free dynamically allocated memory
	for (i = 0; i < nAtoms; i++) {
		free(tmpAtoms[i].neighborList);
		free(tmpAtoms[i].neighborPos);
	}
	free(tmpAtoms);

	return nNewAtoms;
}

/****************************************************************************/

int buildBulkCrystal(atom *atoms, param params) {
	int i, ix, iy, iz, nAtomsInUnitCell, nCells[NDIM], nAtoms = 0;
	vector origin, requiredSize, unitCell, *primitiveVectors; // the primitive lattice vectors
        FILE *pf;

	// Dynamically allocate memory
	primitiveVectors = (vector *) calloc(NDIM, sizeof(vector));

	// Get information on the unit cell 
	nAtomsInUnitCell = fillPrimitiveLatticeVectors(primitiveVectors, params.mainNC.ncCrystalStructure, params.mainNC.ncAtomSymbol1, params.mainNC.ncAtomSymbol2);
	printf("\nThe number of atoms in a unit cell = %d\n", nAtomsInUnitCell); fflush(stdout);

	// Get the required size of the bulk crystal 
	requiredSize = retRequiredSizeVector(primitiveVectors, params.mainNC, params); 
	unitCell = retAddedVectors(primitiveVectors[0], retAddedVectors(primitiveVectors[1], primitiveVectors[2]));
	printf("\nRequired size of bulk crystal to be cut:\n");
	writeVector(requiredSize, stdout);
	printf("Unit cell vector is:\n");
	writeVector(unitCell, stdout);
	nCells[0] = (int)(requiredSize.x/unitCell.x*2.0) + 40;
	nCells[1] = (int)(requiredSize.y/unitCell.y*2.0) + 40;
	nCells[2] = (int)(requiredSize.z/unitCell.z*2.0) + 40;
	for (i = 0; i < NDIM; i++) if (nCells[i] % 2) nCells[i]++; // must be even -> [0]=x [1]=y [2]=z
	printf("Bulk crystal size in x = %.3f\n", unitCell.x * (double)(nCells[0]));
	printf("Bulk crystal size in y = %.3f\n", unitCell.y * (double)(nCells[1]));
	printf("Bulk crystal size in z = %.3f\n", unitCell.z * (double)(nCells[2]));
	// Build a bulk crystal until system is at least the requiredSize
	origin = retZeroVector();
	for (iz = 0; iz < nCells[2]; iz++) {
		for (iy = 0; iy < nCells[1]; iy++) {
			for (ix = 0; ix < nCells[0]; ix++) {
				nAtoms += addUnitCellAtoms(atoms, nAtoms, primitiveVectors, origin, params);
				origin = retAddedVectors(primitiveVectors[0], origin);
			}
			origin = retAddedVectors(retScaledVector(primitiveVectors[0], (double)(-nCells[0])), origin);
			origin = retAddedVectors(primitiveVectors[1], origin);
		}
		origin = retAddedVectors(retScaledVector(primitiveVectors[1], (double)(-nCells[1])), origin);
		origin = retAddedVectors(primitiveVectors[2], origin);
	}

	// Center the bulk slab
	vector avePosition = retZeroVector();
	for (i = 0; i < nAtoms; i++) avePosition = retAddedVectors(atoms[i].pos, avePosition);
	avePosition = retScaledVector(avePosition, (1.0/((double)(nAtoms))));
	for (i = 0; i < nAtoms; i++) atoms[i].pos = retSubtractedVectors(atoms[i].pos, avePosition);
	// TODO: make centerNanocrystal more general (e.g., have it be unit cell or dislocation centered as well as a specific atomSymbol)
	centerNanocrystal(atoms, nAtoms, params.centerAtomSymbol);

	// Make sure atoms types are correctly assigned
	assignAtomTypes(atoms, nAtoms);

	// Add in stacking faults
	if (params.nStackingFaults) {
		addStackingFaults(atoms, nAtoms, primitiveVectors, params);
	}

        // Dipti: print bulk structure
        // pf = fopen("bulk.xyz", "w");
        // fprintf(pf, "%d\n", nAtoms);
        // for (i = 0; i < nAtoms; i++) {
        //     fprintf(pf, "%s ", atoms[i].symbol);
        //     fprintf(pf, "% .8f % .8f % .8f\n", atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
        // }
        // fclose(pf);

	// Free dynamically allocated memory
	free(primitiveVectors);

	return nAtoms;
}

/*****************************************************************************/
// Returns a vector with the minimum size + 4 nm that the bulk crystal 
// must be built to given the desired final nanocrystal size

vector retRequiredSizeVector(vector *primitiveVectors, nanostructure nc, param params) {
	int iNC;
	double extraSize = 40.0;
	vector requiredSize, tmpVector, paddingVector;
	paddingVector.x = paddingVector.y = paddingVector.z = extraSize;

	// Get size NC type specific size requirements for the bulk crystal
	if (! strcmp(nc.ncType, "QD")) {
		requiredSize.x = nc.ncSize.mag;
		requiredSize.y = nc.ncSize.mag;
		requiredSize.z = nc.ncSize.mag;
	} 
	else if (! strcmp(nc.ncType, "NR")) {
		requiredSize.x = nc.ncSize.mag;
		requiredSize.y = nc.ncSize.mag;
		requiredSize.z = nc.ncSize.z;
	} 
	else if (! strcmp(nc.ncType, "NPL")) {
		requiredSize.x = nc.ncSize.x;
		requiredSize.y = nc.ncSize.y;
		requiredSize.z = nc.ncSize.z;
	}  
	else if (! strcmp(nc.ncType, "attachedNCs")) {
		requiredSize = retZeroVector();
		for (iNC = 0; iNC < params.nAttachedNCs; iNC++) {
			tmpVector = retRequiredSizeVector(primitiveVectors, params.attachedNC[iNC], params);
			tmpVector = retSubtractedVectors(tmpVector, paddingVector);
			tmpVector = retAddedVectors(tmpVector, params.attachedNC[iNC].ncCenter);
			requiredSize = retAddedVectors(requiredSize, tmpVector);
		}
	} 
	else {
		writeSeparation(stdout);
		printf("Invalid ncType! ncType = QD or NR or NPL or attachedNCs\n");
		printf("The program is exiting due to this error!\n\n");
		exit(EXIT_FAILURE);
	}

	// Add some extra padding to be safe
	requiredSize = retAddedVectors(requiredSize, paddingVector);
	requiredSize.mag = retVectorMagnitude(requiredSize);

	return (requiredSize);
}

/*****************************************************************************/
// Change the core atoms from mainNC.ncAtomSymbol 1,2 to 
// coreNC.ncAtomSymbol 1,2 for all nCores specified in input.par

void changeAllCoreAtomTypes(atom *atoms, int nAtoms, param params) {
	int i, iNC, jNC;
	double xyProjectionMag;
	vector diff;

	for (iNC = 0; iNC < params.nCores; iNC++) {
		if (! strcmp(params.coreNC[iNC].ncType, "QD")) {
			for (i = 0; i < nAtoms; i++) {
				diff = retSubtractedVectors(params.coreNC[iNC].ncCenter, atoms[i].pos);
				if (diff.mag < (0.5*params.coreNC[iNC].ncSize.mag)) {
					if ( ! changeSingleCoreAtomType(&(atoms[i]), params.mainNC, params.coreNC[iNC]) ) {
						continue;
					}
					else {
						for (jNC = 0; jNC < iNC; jNC++) {
							changeSingleCoreAtomType(&(atoms[i]), params.coreNC[jNC], params.coreNC[iNC]);
						}
					}
				}
			}
		} 
		else if (! strcmp(params.coreNC[iNC].ncType, "NR")) {
			for (i = 0; i < nAtoms; i++) {
				diff = retSubtractedVectors(params.coreNC[iNC].ncCenter, atoms[i].pos);
				xyProjectionMag = 2.0*sqrt(sqr(diff.x) + sqr(diff.y));
				if (xyProjectionMag < (params.coreNC[iNC].ncSize.mag) && fabs(diff.z) < (0.5*params.coreNC[iNC].ncSize.z)) {
					if ( ! changeSingleCoreAtomType(&(atoms[i]), params.mainNC, params.coreNC[iNC]) ) {
						continue;
					}
					else {
						for (jNC = 0; jNC < iNC; jNC++) {
							changeSingleCoreAtomType(&(atoms[i]), params.coreNC[jNC], params.coreNC[iNC]);
						}
					}				
				}
			}
		} 
		else if (! strcmp(params.coreNC[iNC].ncType, "NPL")) {
			for (i = 0; i < nAtoms; i++) {
				diff = retSubtractedVectors(params.coreNC[iNC].ncCenter, atoms[i].pos);
				if (fabs(diff.x) < (0.5*params.coreNC[iNC].ncSize.x) && fabs(diff.y) < (0.5*params.coreNC[iNC].ncSize.y)
						&& fabs(diff.z) < (0.5*params.coreNC[iNC].ncSize.z)) {
					if ( ! changeSingleCoreAtomType(&(atoms[i]), params.mainNC, params.coreNC[iNC]) ) {
						continue;
					}
					else {
						for (jNC = 0; jNC < iNC; jNC++) {
							changeSingleCoreAtomType(&(atoms[i]), params.coreNC[jNC], params.coreNC[iNC]);
						}						
					}
				}
			}
		}
	}

	// Make sure atoms types are correctly assigned
	assignAtomTypes(atoms, nAtoms);

	return;
}

/*****************************************************************************/
// Change the core atoms from mainNC.ncAtomSymbol 1,2 to coreNC.ncAtomSymbol 1,2
// Returns 0 if atom changed and returns 1 if not

int changeSingleCoreAtomType(atom *atoms, nanostructure mainNC, nanostructure coreNC) {

	if (! strcmp(atoms[0].symbol, mainNC.ncAtomSymbol1)) {
		strcpy(atoms[0].symbol, coreNC.ncAtomSymbol1);
	}
	else if (! strcmp(atoms[0].symbol, mainNC.ncAtomSymbol2)) {
		strcpy(atoms[0].symbol, coreNC.ncAtomSymbol2);
	}
	else {
		return 1;
	}

	return 0;
}

/*****************************************************************************/
// Translate all atoms to make the centerAtomSymbol that 
// was closest to the center be located at the origin

void centerNanocrystal(atom *atoms, int nAtoms, char *centerAtomSymbol) {
	int i;
	double minMag = 100000.0;
	vector shiftVector = retZeroVector();

	printf("\nPutting %s at the origin of the nanocrystal\n\n", centerAtomSymbol);

	// Find the atom closest to the center
	for (i = 0; i < nAtoms; i++) {
		if (! strcmp(atoms[i].symbol, centerAtomSymbol)) {
			if (retVectorMagnitude(atoms[i].pos) < minMag) {
				minMag = retVectorMagnitude(atoms[i].pos);
				shiftVector = atoms[i].pos;
			}
		}
	} 

	// Translate all atoms equally so that the new atoms[iMin].pos.mag = 0.0
	for (i = 0; i < nAtoms; i++) atoms[i].pos = retSubtractedVectors(atoms[i].pos, shiftVector);

	return;
}

/****************************************************************************/
// Adds in stacking faults the atoms (by shifting specific layers) 
// assumes that the atoms have the Wurtzite crystal structure

void addStackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, param params) {
	int i, j, sfOrder[params.nStackingFaults];
	int counter = 0, minPosition = -100;
	int flag;

	// Initialize the sfOrder list
	for (i = 0; i < params.nStackingFaults; i++) sfOrder[i] = -1;

	// TODO: add stacking faults in correct order so user input can be in any order
	// Determine correct order to add in the type 1 and type 2 stacking faults 
	// all type 1s first followed by lowest to highest type 2s
	//for (i = 0; i < params.nStackingFaults; i++) {
	// 	for (j = i; j < params.nStackingFaults; j++) {
	// 		if (j == sfOrder[j]) continue;
	// 		else if (params.stackFault[j].type == 1 && params.stackFault[j].position > minPosition) {
	// 			sfOrder[0] = j;
	// 			counter++;
	// 		}
	// 	}
	// }

	// Add in the stacking faults one by one 
	// flag used to determine if normal shifting (flag = 0) or if opposite shifting (flag = 1)
	// is needed to add the tacking faults correctly when multiple stacking faults are present
	flag = 0;
	for (i = 0; i < params.nStackingFaults; i++) {
	 	flag = flag % 2;
	 	if (params.stackFault[i].type == 1) {
	 		addType1StackingFaults(atoms, nAtoms, primitiveVectors, params.stackFault[i].position, flag);
	 		flag++;
	 		minPosition = params.stackFault[i].position;
	 	}
	 	else if (params.stackFault[i].type == 2) {
	 		if (params.stackFault[i].position < minPosition) flag = 0; // if below type 1s then flag = 0
	 		addType2StackingFaults(atoms, nAtoms, primitiveVectors, params.stackFault[i].position, flag);
	 	}
	 	else {
	 		writeSeparation(stdout);
	 		printf("WARNING: No stacking fault added due to incorrect stackingFaultType!\n");
	 		printf("stackingFaultType = 1 or 2 only!\n");
	 		writeSeparation(stdout);
	 	}
	}

	return;
}

/****************************************************************************/
// Creates a type 2 stacking fault at location specified by position
// For example, position= 0 shifts all the atoms above the origin
// position must be an integer -> -n,-n+1,...,0,...,n 
// flag should = 1 if adding a type 2 sf above one type 1 sf -> caller's responsibility

void addType2StackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, int position, int flag) {
	int i;
	double shiftY = 2.48; // original -> should be positive when adding just 1 SF, 2.48 A for CdSe
	double shiftAboveZ = 0.5*(double)(position)*primitiveVectors[2].z;

	// Shift atoms to the negative y direction if position is even
	if (! (abs(position) % 2)) shiftY *= -1.0;
	
	// Reverse shifting direction if flag != 0
	if (flag) shiftY *= -1.0;

	printf("Adding in a type 2 stacking fault at stacking fault position = %d and flag = %d\n\n", position, flag);

	// Shift select wurtzite atoms in the y-direction to generate type 2 stacking fault
	for (i = 0; i < nAtoms; i++) {
		if (atoms[i].pos.z > shiftAboveZ+EPS) {
			atoms[i].pos.y += shiftY; 
		}
	}

	return;
}

/****************************************************************************/
// Creates a type 1 stacking fault at location specified by position
// For example, position = 0 shifts the atoms directly above the origin
// then alternates between not shifting and shifting the next vertical Cd-Se pairs
// position must be an integer -> -n,-n+1,...,0,...,n 
// flag should be 1 when putting in second type 1 sf above another type 1 -> 0 otherwise

void addType1StackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, int position, int flag) {
	int i, k;
	double shiftY = 2.48;
	double shift = 0.5*(double)(position)*primitiveVectors[2].z;
	double normalizedZ;

	printf("Adding in a type 1 stacking fault at stacking fault position = %d and flag = %d\n\n", position, flag);

	// Shift every other verticle pair of wurtzite atoms in the y-direction to generate type 1 stacking fault
	for (i = 0; i < nAtoms; i++) {
		normalizedZ = (atoms[i].pos.z-EPS-shift) / primitiveVectors[2].z;
		while (normalizedZ > 1.0) normalizedZ -= 1.0;
		if (! (abs(position + flag) % 2)) {
			if (normalizedZ < 0.5 && normalizedZ > 0.0) atoms[i].pos.y -= shiftY; // shift left for even positions
		} else {
			if (normalizedZ < 0.5 && normalizedZ > 0.0) atoms[i].pos.y += shiftY; // shift right for even positions
		} 
	}

	return;
}

/****************************************************************************/

int addUnitCellAtoms(atom *atoms, int nOrigAtoms, vector *primitiveVectors, vector origin, param params) {
	int i, nAtomsInUnitCell = 4;

	// wurtzite: Cd (1/6, 1/6, 0.4375) Cd (-1/6, -1/6, -0.0625) Se (1/6, 1/6, 0.0625) Se (-1/6, -1/6, -0.4375)
	// zincblende: Cd (-1/8, -1/8, -1/8) Se (1/8, 1/8, 1/8)

	if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) {
		nAtomsInUnitCell = 4;
		strcpy(atoms[nOrigAtoms].symbol, params.mainNC.ncAtomSymbol1);
		atoms[nOrigAtoms].pos.x = 1.0/6.0; atoms[nOrigAtoms].pos.y = 1.0/6.0; atoms[nOrigAtoms].pos.z = 0.4375;
		strcpy(atoms[nOrigAtoms+1].symbol, params.mainNC.ncAtomSymbol1);
		atoms[nOrigAtoms+1].pos.x = -1.0/6.0; atoms[nOrigAtoms+1].pos.y = -1.0/6.0; atoms[nOrigAtoms+1].pos.z = -0.0625; 
		strcpy(atoms[nOrigAtoms+2].symbol, params.mainNC.ncAtomSymbol2);
		atoms[nOrigAtoms+2].pos.x = 1.0/6.0; atoms[nOrigAtoms+2].pos.y = 1.0/6.0; atoms[nOrigAtoms+2].pos.z = 0.0625;
		strcpy(atoms[nOrigAtoms+3].symbol, params.mainNC.ncAtomSymbol2);
		atoms[nOrigAtoms+3].pos.x = -1.0/6.0; atoms[nOrigAtoms+3].pos.y = -1.0/6.0; atoms[nOrigAtoms+3].pos.z = -0.4375;		
	} else if (! strcmp(params.mainNC.ncCrystalStructure, "zincblende")) {
		nAtomsInUnitCell = 2;
		strcpy(atoms[nOrigAtoms].symbol, params.mainNC.ncAtomSymbol1);
		atoms[nOrigAtoms].pos.x = -0.125; atoms[nOrigAtoms].pos.y = -0.125; atoms[nOrigAtoms].pos.z = -0.125;
		strcpy(atoms[nOrigAtoms+1].symbol, params.mainNC.ncAtomSymbol2);
		atoms[nOrigAtoms+1].pos.x = 0.125; atoms[nOrigAtoms+1].pos.y = 0.125; atoms[nOrigAtoms+1].pos.z = 0.125; 
	}

	// Scale and shift the vectors
	for (i = 0; i < nAtomsInUnitCell; i++) {
		atoms[nOrigAtoms+i].pos.x *= (primitiveVectors[0].x + primitiveVectors[1].x + primitiveVectors[2].x);
		atoms[nOrigAtoms+i].pos.y *= (primitiveVectors[0].y + primitiveVectors[1].y + primitiveVectors[2].y);
		atoms[nOrigAtoms+i].pos.z *= (primitiveVectors[0].z + primitiveVectors[1].z + primitiveVectors[2].z);
		atoms[nOrigAtoms+i].pos.mag = retVectorMagnitude(atoms[i].pos);
  	}
  	for (i = 0; i < nAtomsInUnitCell; i++) atoms[nOrigAtoms+i].pos = retAddedVectors(atoms[nOrigAtoms+i].pos, origin);

	return nAtomsInUnitCell;
}

/****************************************************************************/
// 

int addAllNewLayers(atom *atoms, param params) {
	FILE *pf;
	int iAtom, iAtomType, iNeighbor, iCenterAtom, iTetrohedron, iHalfLayer, wurtziteLayerFlag;
	int iDuplicate;
	int nAtoms, nFinalAtoms, nOrigAtomTypes, nUniqueTetrohedrons, nOrigAtoms = params.nAtoms;
	int nAtomsInNewLayer, nAddedAtoms;
	int nMaxAtoms = 500000;
	double latticeConstant, c, diameter, centerAtomZPos[params.nAtomTypes];
	double piOverThree = PIE/3.0;
	double diameterIncrease;
	char newAtomSymbol[100], fileName[100];	
	vector *atomPositions, *tetrohedronVectors, maxDimensions, scaledMaxDimensions;
	vector newTestPosition, plusXShift, minusXShift;
	atom *tmpAtoms;
	atom *finalAtoms;

	// Write beginning of function
	fprintf(stdout, "Adding %d half layers to the original configuration\n", params.nHalfLayersToAdd);

	// Make sure all atom types are assigned properly
	params.nAtomTypes = nOrigAtomTypes = assignAtomTypes(atoms, nOrigAtoms);
	if (nOrigAtomTypes != 2) inputError("The seed confFile must have only 1 metal and 1 non-metal atom types!");

	// Calculate the nearest neighbor list for the original structure
	calcNearestNeighborLists(atoms, params.nMaxBonds, nOrigAtoms);

	// Allocate memory for the positions of all atoms 
	atomPositions = (vector *) calloc(nMaxAtoms, sizeof(vector));

	// Copy all atom positions into vector only atomPositions and 
	// calculate the size statistics for the entire nanocrystal
	for (iAtom = 0; iAtom < params.nAtoms; iAtom++) atomPositions[iAtom] = atoms[iAtom].pos;
	calcNanocrystalDimensions(&maxDimensions, atomPositions, params.nAtoms);
	diameter = retMaxDistance(atomPositions, params.nAtoms);

	// Determine the total number of unique tetrohedron vectors
	if (! strcmp(params.mainNC.ncCrystalStructure, "zincblende")) {
		nUniqueTetrohedrons = nOrigAtomTypes*params.nMaxBonds;
		latticeConstant = c = retLatticeConstant(params.mainNC.ncCrystalStructure, atoms[0].symbol, atoms[atoms[0].neighborList[0]].symbol);
	}
	else if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) {
		nUniqueTetrohedrons = 2*nOrigAtomTypes*params.nMaxBonds;
		latticeConstant = retLatticeConstant(params.mainNC.ncCrystalStructure, atoms[0].symbol, atoms[atoms[0].neighborList[0]].symbol);
		c = sqrt(8.0/3.0)*latticeConstant;
	}
	else {
		inputError("Only wurtzite and zincblende crystal structures are allowed");
	}
	fprintf(stdout, "The number of original atoms = %d\n", nOrigAtoms);
	fprintf(stdout, "The number of original atom types = %d\n", nOrigAtomTypes);
	fprintf(stdout, "The number of unique tetrohedrons = %d\n", nUniqueTetrohedrons); fflush(stdout);

	// Allocate memory for centerAtoms and tetrohedronVectors
	tetrohedronVectors = (vector *) calloc(nUniqueTetrohedrons, sizeof(vector));
	tmpAtoms = (atom *) calloc(nMaxAtoms, sizeof(atom));
	finalAtoms = (atom *) calloc(nMaxAtoms, sizeof(atom));
	for (iAtom = 0; iAtom < nMaxAtoms; iAtom++) {
		finalAtoms[iAtom].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		finalAtoms[iAtom].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
		tmpAtoms[iAtom].neighborList = (int *) calloc(params.nMaxBonds, sizeof(int));
		tmpAtoms[iAtom].neighborPos = (vector *) calloc(params.nMaxBonds, sizeof(vector));
	}
	deepCopyAllAtoms(finalAtoms, atoms, nOrigAtoms, params.nMaxBonds);

	// Store the unique sets of 4 vectors that make up the tetrohedron for all unique atom types and bonding environments  	
	double atomMagNearestOrigin;
	iTetrohedron = 0;
	for (iAtomType = 0; iAtomType < nOrigAtomTypes; iAtomType++) {
		// Find the index of the atom nearest the origin, TODO: have this in a separate function
		atomMagNearestOrigin = 100000.0;
		for (iAtom = 0; iAtom < nOrigAtoms; iAtom++) {
			if ((atoms[iAtom].isMetalAtomFlag == iAtomType) && (atoms[iAtom].pos.mag < (atomMagNearestOrigin + EPS))
					&& ! isAtomASurfaceAtom(atoms, iAtom, params.nMaxBonds)) {
				iCenterAtom = iAtom;
				atomMagNearestOrigin = atoms[iAtom].pos.mag;
				centerAtomZPos[iAtomType] = atoms[iAtom].pos.z; // only used for wurtzite structures
			}
		}
		fprintf(stdout, "\nThe index of the center atom in the original structure = %d\n", iCenterAtom);
		fprintf(stdout, "The z coordinate of the center atom = % .8f\n", centerAtomZPos[iAtomType]); 
		fprintf(stdout, "Center atom: %s %d % .8f % .8f % .8f % .8f\n", atoms[iCenterAtom].symbol, atoms[iCenterAtom].isMetalAtomFlag,
							atoms[iCenterAtom].pos.x, atoms[iCenterAtom].pos.y, atoms[iCenterAtom].pos.z, atoms[iCenterAtom].pos.mag);
		fflush(stdout);
		// Store the neighbor list for the atom closest to the origin
		for (iNeighbor = 0; iNeighbor < params.nMaxBonds; iNeighbor++) {
			tetrohedronVectors[iTetrohedron] = atoms[iCenterAtom].neighborPos[iNeighbor];
			tetrohedronVectors[iTetrohedron] = retScaledVector(tetrohedronVectors[iTetrohedron], -1.0);
			if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) { // wurtzite has two bonding environments per atom type
				tetrohedronVectors[iTetrohedron + params.nMaxBonds] = retXYPlaneRotatedVector(tetrohedronVectors[iTetrohedron], piOverThree);
			}
			fprintf(stdout, "Center atom neighbor %d: %d % .8f % .8f % .8f % .8f\n",  
				iNeighbor, atoms[iCenterAtom].neighborList[iNeighbor],
				tetrohedronVectors[iTetrohedron].x, tetrohedronVectors[iTetrohedron].y, 
				tetrohedronVectors[iTetrohedron].z, tetrohedronVectors[iTetrohedron].mag); fflush(stdout);
			iTetrohedron++;
		}
		if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) { // account for wurtzite having two bonding environments per atom type
			iTetrohedron += params.nMaxBonds;
		}		
	}

	// Print tetrohedron vectors
	pf = fopen("tetrahedronVectors.dat", "w");
	fprintf(pf, "%d\n", nUniqueTetrohedrons);
	for (iTetrohedron = 0; iTetrohedron < nUniqueTetrohedrons; iTetrohedron++) {
		writeVector(tetrohedronVectors[iTetrohedron], pf);
	}
	fclose(pf);

	// Add in a tetrohedron of atoms for all surface atoms and remove duplicants for each layer
	nFinalAtoms = nOrigAtoms;
	int nTmpAtoms;
	vector tmpVector;
	int increaseNCSeparation = 1;
	for (iHalfLayer = 0; iHalfLayer < params.nHalfLayersToAdd; iHalfLayer++) {
		nAtomsInNewLayer = 0; // TODO: have this next for loop in a separate function
		for (iAtom = 0; iAtom < nFinalAtoms; iAtom++) {
			if (! strcmp(params.mainNC.ncType, "attachedNCs") && fabs(finalAtoms[iAtom].pos.x) > latticeConstant/sqrt(3.0)+0.2) {
				continue;			
			}
			// if (! strcmp(params.mainNC.ncType, "attachedNCs") && increaseNCSeparation) {
			// 	if (fabs(finalAtoms[iAtom].pos.x) < latticeConstant/sqrt(3.0)-0.2) {

			// 	}
			// }
			else if (isAtomASurfaceAtom(finalAtoms, iAtom, params.nMaxBonds) && newIsAtomASurfaceAtom(finalAtoms, iAtom, nFinalAtoms)) {
				// Determine which set of tetrohedron vectors will be used to add the new atoms
				if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) {
					wurtziteLayerFlag = getZLayer(finalAtoms[iAtom].pos.z, centerAtomZPos[finalAtoms[iAtom].isMetalAtomFlag], c);
					iTetrohedron = (2*finalAtoms[iAtom].isMetalAtomFlag + wurtziteLayerFlag)*params.nMaxBonds;
				}
				else { 
					iTetrohedron = finalAtoms[iAtom].isMetalAtomFlag*params.nMaxBonds;
				}
				// Determine if adding new metal atoms or new chalcogen
				if ( isAtomAMetal(finalAtoms[iAtom].symbol) ) {
				 	strcpy(newAtomSymbol, params.newLayer[iHalfLayer].nonMetalAtomSymbol);
				}
				else {
					strcpy(newAtomSymbol, params.newLayer[iHalfLayer].metalAtomSymbol);
				}
				// Add in the new atoms
				nAddedAtoms = addTetrohedronAtoms(&(tmpAtoms[nAtomsInNewLayer]), newAtomSymbol, 
													finalAtoms[iAtom].pos, &(tetrohedronVectors[iTetrohedron]), 
													params.nMaxBonds);
				nAtomsInNewLayer += nAddedAtoms;
			}
		}

		// Make sure the corner atoms were added when adding NPL layers:
		if (! strcmp(params.mainNC.ncType, "NPL") && ! strcmp(params.mainNC.ncCrystalStructure, "zincblende")) {
			nTmpAtoms = nAtomsInNewLayer;
			assignAtomTypes(tmpAtoms, nTmpAtoms);
			for (iAtom = 0; iAtom < nTmpAtoms; iAtom++) atomPositions[iAtom] = tmpAtoms[iAtom].pos;
			calcNanocrystalDimensions(&maxDimensions, atomPositions, nTmpAtoms);
			fprintf(stdout, "\nAdding the atoms within the box:\n");
			for (iAtom = 0; iAtom < nTmpAtoms; iAtom++) {
				if ( isAtomAMetal(tmpAtoms[iAtom].symbol) ) {
				 	strcpy(newAtomSymbol, params.newLayer[iHalfLayer].nonMetalAtomSymbol);
				}
				else {
					strcpy(newAtomSymbol, params.newLayer[iHalfLayer].metalAtomSymbol);
				}
				iTetrohedron = tmpAtoms[iAtom].isMetalAtomFlag*params.nMaxBonds;
				nAddedAtoms = addTetrohedronAtoms(&(tmpAtoms[nAtomsInNewLayer]), newAtomSymbol, 
														tmpAtoms[iAtom].pos, &(tetrohedronVectors[iTetrohedron]), 
														params.nMaxBonds);
				nAtomsInNewLayer += nAddedAtoms;
			}
			fprintf(stdout, "nAtomsInNewLayer before NPL Bulk cutting = %d\n", nAtomsInNewLayer);	
			nAtomsInNewLayer = cutNPLFromBulk(tmpAtoms, nAtomsInNewLayer, params.nMaxBonds, maxDimensions.x+0.1, maxDimensions.y+0.1, maxDimensions.z+0.1);
			fprintf(stdout, "nAtomsInNewLayer after  NPL Bulk cutting = %d\n", nAtomsInNewLayer);
		}
		else if (! strcmp(params.mainNC.ncType, "attachedNCs")) {
			fprintf(stdout, "nAtomsInNewLayer before NPL Bulk cutting = %d\n", nAtomsInNewLayer);
			nAtomsInNewLayer = cutNPLFromBulk(tmpAtoms, nAtomsInNewLayer, params.nMaxBonds, 
										2.0*latticeConstant/sqrt(3.0), 10000., 10000.);
			fprintf(stdout, "nAtomsInNewLayer after  NPL Bulk cutting = %d\n", nAtomsInNewLayer);

			// Add in atoms that are closer to origin than 
			nTmpAtoms = nAtomsInNewLayer;
			//for (iAtom = 0; iAtom < nTmpAtoms; iAtom++) atomPositions[iAtom] = tmpAtoms[iAtom].pos;
			//calcNanocrystalDimensions(&maxDimensions, atomPositions, nTmpAtoms);
			assignAtomTypes(tmpAtoms, nTmpAtoms);
			for (iAtom = 0; iAtom < nTmpAtoms; iAtom++) {
			 	// Determine which set of tetrohedron vectors will be used to add the new atoms
			 	if (! strcmp(params.mainNC.ncCrystalStructure, "wurtzite")) {
			 		wurtziteLayerFlag = getZLayer(tmpAtoms[iAtom].pos.z, centerAtomZPos[tmpAtoms[iAtom].isMetalAtomFlag], c);
			 		iTetrohedron = (2*tmpAtoms[iAtom].isMetalAtomFlag + wurtziteLayerFlag)*params.nMaxBonds;
			 	}
			 	else { 
			 		iTetrohedron = tmpAtoms[iAtom].isMetalAtomFlag*params.nMaxBonds;
			 	}
			 	if ( isAtomAMetal(tmpAtoms[iAtom].symbol) ) {
			 	 	strcpy(newAtomSymbol, params.newLayer[iHalfLayer].nonMetalAtomSymbol);
			 	}
			 	else {
			 		strcpy(newAtomSymbol, params.newLayer[iHalfLayer].metalAtomSymbol);
			 	}
			 	for (iNeighbor = 0; iNeighbor < params.nMaxBonds; iNeighbor++) {
			 		tmpVector = retAddedVectors(tmpAtoms[iAtom].pos, tetrohedronVectors[iTetrohedron + iNeighbor]);
			 		if (fabs(tmpVector.x) < latticeConstant/sqrt(3.0)+0.2) {
			 			strcpy(tmpAtoms[nAtomsInNewLayer].symbol, newAtomSymbol);
			 			tmpAtoms[nAtomsInNewLayer].pos = tmpVector;
			 			nAtomsInNewLayer++;
			 		}
			 	}
			}
			plusXShift = retZeroVector();
			minusXShift = retZeroVector();
			nAddedAtoms = 0;
			assignAtomTypes(tmpAtoms, nAtomsInNewLayer);
			sprintf(fileName, "tmpAtoms_%d_%d.xyz", iHalfLayer, 0);
			pf = fopen(fileName, "w");
			writeXYZFile(tmpAtoms, nAtomsInNewLayer, 1.0, pf);
			fclose(pf);
			for (iDuplicate = 0; iDuplicate < iHalfLayer+1; iDuplicate++) {
				plusXShift.x = sqrt(3.0)*latticeConstant*(double)(1+iDuplicate);
				minusXShift.x = -plusXShift.x;
				deepCopyAllAtoms(&(tmpAtoms[nAtomsInNewLayer  + nAddedAtoms]), tmpAtoms, nAtomsInNewLayer, params.nMaxBonds);
				translateAllAtoms(&(tmpAtoms[nAtomsInNewLayer + nAddedAtoms]), nAtomsInNewLayer, plusXShift);
				nAddedAtoms += nAtomsInNewLayer;
				deepCopyAllAtoms(&(tmpAtoms[nAtomsInNewLayer  + nAddedAtoms]), tmpAtoms, nAtomsInNewLayer, params.nMaxBonds);
				translateAllAtoms(&(tmpAtoms[nAtomsInNewLayer + nAddedAtoms]), nAtomsInNewLayer, minusXShift);
				nAddedAtoms += nAtomsInNewLayer;
				sprintf(fileName, "tmpAtoms_%d_%d.xyz", iHalfLayer, iDuplicate+1);
				pf = fopen(fileName, "w");
				writeXYZFile(tmpAtoms, nAtomsInNewLayer+nAddedAtoms, 1.0, pf);
				fclose(pf);
			}
			nAtomsInNewLayer += nAddedAtoms;
			fprintf(stdout, "nAtomsInNewLayer after  duplicating atoms %d times = %d\n", iDuplicate, nAtomsInNewLayer);
			nAtomsInNewLayer = cutNPLFromBulk(tmpAtoms, nAtomsInNewLayer, params.nMaxBonds, 
											0.5*maxDimensions.x+0.1, maxDimensions.y+0.1, maxDimensions.z+0.1);
			fprintf(stdout, "nAtomsInNewLayer after  NPL x y and z cutting = %d\n", nAtomsInNewLayer);
		}

		// Add the new atoms stored in tmpAtoms to finalAtoms and remove all duplicants
		fprintf(stdout, "\nAdding new layer of: %s %s\n", params.newLayer[iHalfLayer].metalAtomSymbol, 
														params.newLayer[iHalfLayer].nonMetalAtomSymbol);
		fprintf(stdout, "nAtomsInNewLayer before removing duplicants = %d\n", nAtomsInNewLayer);
		nAtomsInNewLayer = removeDuplicatedAtoms(tmpAtoms, nAtomsInNewLayer, params.nMaxBonds);
		fprintf(stdout, "nAtomsInNewLayer after  removing duplicants = %d\n", nAtomsInNewLayer);
		fprintf(stdout, "nFinalAtoms before copying and removing duplicants = %d\n", nFinalAtoms);
		deepCopyAllAtoms(&(finalAtoms[nFinalAtoms]), tmpAtoms, nAtomsInNewLayer, params.nMaxBonds);
		nFinalAtoms = removeDuplicatedAtoms(finalAtoms, nFinalAtoms+nAtomsInNewLayer, params.nMaxBonds);
		params.nAtomTypes = assignAtomTypes(finalAtoms, nFinalAtoms);
		fprintf(stdout, "nFinalAtoms after  copying and removing duplicants = %d\n", nFinalAtoms); fflush(stdout);
	}

	// Remove dangling atoms, TODO: check if necessary or not
	if (params.remDanglingAtoms) {
		fprintf(stdout, "\nRemoving dangling atoms from the finalAtoms configuration\n"); fflush(stdout);
		deepCopyAllAtoms(tmpAtoms, finalAtoms, nFinalAtoms, params.nMaxBonds);
		pf = fopen("dangling_atoms_conf.xyz", "w");
		writeXYZFile(tmpAtoms, nFinalAtoms, 1.0, pf);
		fclose(pf);
		nFinalAtoms = cutDanglingAtoms(finalAtoms, tmpAtoms, params.nMaxBonds, params.nMinBonds, nFinalAtoms);
	}

	// Print final atom configuration
	pf = fopen("final_unpassivated_conf.xyz", "w");
	writeXYZFile(finalAtoms, nFinalAtoms, 1.0, pf);
	fclose(pf);

	// Calculate size and bond statistics of the final structure
	params.nAtomTypes = assignAtomTypes(finalAtoms, nFinalAtoms);
	params.nAtoms = nFinalAtoms;
	calcSizeStatistics(finalAtoms, params);
	calcBondStatistics(finalAtoms, params);

	// Free dynamically allocated memory
	free(tetrohedronVectors);
	for (iAtom = 0; iAtom < nMaxAtoms; iAtom++) {
		free(finalAtoms[iAtom].neighborList); free(finalAtoms[iAtom].neighborPos);
		free(tmpAtoms[iAtom].neighborList); free(tmpAtoms[iAtom].neighborPos); 
	}
	free(finalAtoms);
	free(tmpAtoms);
	free(atomPositions);

	// Print the time and exit the program
	writeSeparation(stdout);
	writeCurrentTime(stdout);
	exit(EXIT_SUCCESS);
	
	return nFinalAtoms;
}

/****************************************************************************/
// returns 0 if 
// returns 1 if ...

int getZLayer(double testAtomZPos, double centerAtomZPos, double cLatticeConstant) {
	int sameNeighborAsCenterAtom = 0;
	double layerNum;

	// layerNum should be close to 0 or 1 
	layerNum = fmod( ( 2.0*(testAtomZPos - centerAtomZPos) / cLatticeConstant ), 1.999);

	if (fabs(layerNum) > 0.5) {
		return sameNeighborAsCenterAtom = 1;
	}

	return sameNeighborAsCenterAtom;
}

/****************************************************************************/
// Function to add the 4 atoms that make up the tetrohedron surrounding a given atom

int addTetrohedronAtoms(atom *atoms, char *newAtomSymbol, vector origAtomPos, vector *tetrohedronVectors, int nNewAtoms) {
	int iTetrohedron;

	for (iTetrohedron = 0; iTetrohedron < nNewAtoms; iTetrohedron++) {
		strcpy(atoms[iTetrohedron].symbol, newAtomSymbol);  
		atoms[iTetrohedron].pos = retAddedVectors(origAtomPos, tetrohedronVectors[iTetrohedron]);
	}

 	return nNewAtoms;
} 

/****************************************************************************/
// Function

int addCornerAtoms(atom *atoms, char *newAtomSymbol, vector maxDimensions) {
	int iAtom, nNewAtoms = 0;

	// Add in the atom symbol for all atoms and (+-xmax, +-ymax, +-zmax)
	for (iAtom = 0; iAtom < 8; iAtom++) { 
		strcpy(atoms[iAtom].symbol, newAtomSymbol);
		atoms[iAtom].pos.x = maxDimensions.x;
		atoms[iAtom].pos.y = maxDimensions.y;
		atoms[iAtom].pos.z = maxDimensions.z;
		if (iAtom == 1) atoms[iAtom].pos.x *= -1.0;
		else if (iAtom == 2) atoms[iAtom].pos.y *= -1.0;
		else if (iAtom == 3) atoms[iAtom].pos.z *= -1.0;
		else if (iAtom == 4) {
			atoms[iAtom].pos.x *= -1.0; atoms[iAtom].pos.y *= -1.0;
		}
		else if (iAtom == 5) {
			atoms[iAtom].pos.x *= -1.0; atoms[iAtom].pos.z *= -1.0;
		}
		else if (iAtom == 6) {
			atoms[iAtom].pos.y *= -1.0; atoms[iAtom].pos.z *= -1.0;
		}
		else if (iAtom == 7) atoms[iAtom].pos = retScaledVector(atoms[iAtom].pos, -1.0);
		atoms[iAtom].pos.mag = retVectorMagnitude(atoms[iAtom].pos);
		nNewAtoms++;
	}

	return nNewAtoms;
}

/****************************************************************************/
// 

double retLatticeConstant(char *ncCrystalStructure, char *atomSymbol1, char *atomSymbol2) {
	//double latticeConstant = EPS;
	double latticeConstant = -1.0;

	// Lattice constant scaling, the lattice constants are defined in nc.h
	if (! strcmp(atomSymbol1, "Cd")) {
		if (! strcmp(atomSymbol2, "Se") && ! strcmp(ncCrystalStructure, "wurtzite")) latticeConstant = ACdSeWZ;
		else if (! strcmp(atomSymbol2, "Se") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = ACdSeZB;
		else if (! strcmp(atomSymbol2, "S") && ! strcmp(ncCrystalStructure, "wurtzite"))    latticeConstant = ACdSWZ;
		else if (! strcmp(atomSymbol2, "S") && ! strcmp(ncCrystalStructure, "zincblende"))  latticeConstant = ACdSZB;
		else if (! strcmp(atomSymbol2, "Te") && ! strcmp(ncCrystalStructure, "wurtzite"))   latticeConstant = ACdTeWZ;
		else if (! strcmp(atomSymbol2, "Te") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = ACdTeZB;	
	}
	else if (! strcmp(atomSymbol1, "Zn")) {
		if (! strcmp(atomSymbol2, "Se") && ! strcmp(ncCrystalStructure, "wurtzite")) latticeConstant = AZnSeWZ;
		else if (! strcmp(atomSymbol2, "Se") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = AZnSeZB;
		else if (! strcmp(atomSymbol2, "S") && ! strcmp(ncCrystalStructure, "wurtzite"))    latticeConstant = AZnSWZ;
		else if (! strcmp(atomSymbol2, "S") && ! strcmp(ncCrystalStructure, "zincblende"))  latticeConstant = AZnSZB;
		else if (! strcmp(atomSymbol2, "Te") && ! strcmp(ncCrystalStructure, "wurtzite"))   latticeConstant = AZnTeWZ;
		else if (! strcmp(atomSymbol2, "Te") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = AZnTeZB;
	}
	else if (! strcmp(atomSymbol1, "In")) {
		if (! strcmp(atomSymbol2, "As") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = AInAsZB;
                else if (! strcmp(atomSymbol2, "P") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = AInPZB;
	}
	else if (! strcmp(atomSymbol1, "Ga")) {
		if (! strcmp(atomSymbol2, "As") && ! strcmp(ncCrystalStructure, "zincblende")) latticeConstant = AGaAsZB;
	}
	else if ( isAtomAMetal(atomSymbol2) ) {
		latticeConstant = retLatticeConstant(ncCrystalStructure, atomSymbol2, atomSymbol1);
	}
	else if (latticeConstant < 0.0) {
		printf("Allowed crystal structure are, ncCrystalStructure = wurtzite or zincblende\n");
		printf("Allowed atomSymbol1/2 are only Cd, Zn, In, Ga\n");
		printf("Allowed atomSymbol2/1 are only Se, S, Te or As\n");
		printf("The program is exiting!!!\n");
		exit(EXIT_FAILURE);
	}

	return latticeConstant;
}

/****************************************************************************/

int fillPrimitiveLatticeVectors(vector *primitiveVectors, char *ncCrystalStructure, char *atomSymbol1, char *atomSymbol2) {
	int i, numAtomsInUnitCell;
	double latticeConstant = -1.0;

	// Fill in the primitive vectors
	if (! strcmp(ncCrystalStructure, "wurtzite")) {
		primitiveVectors[0].x = 1.0; primitiveVectors[0].y = 0.0; primitiveVectors[0].z = 0.0;
		primitiveVectors[1].x = 0.5; primitiveVectors[1].y = 0.5*sqrt(3.0); primitiveVectors[1].z = 0.0;
		primitiveVectors[2].x = 0.0; primitiveVectors[2].y = 0.0; primitiveVectors[2].z = sqrt(8.0/3.0);
		numAtomsInUnitCell = 4; // 4 atoms in the wurtzite unit cell (e.g., 2 Cd and 2 Se)
	} 
	else if (! strcmp(ncCrystalStructure, "zincblende")) {
		primitiveVectors[0].x = 0.0; primitiveVectors[0].y = 0.5; primitiveVectors[0].z = 0.5;
		primitiveVectors[1].x = 0.5; primitiveVectors[1].y = 0.0; primitiveVectors[1].z = 0.5;
		primitiveVectors[2].x = 0.5; primitiveVectors[2].y = 0.5; primitiveVectors[2].z = 0.0;
		numAtomsInUnitCell = 2; // 2 atoms in the zincblende unit cell (e.g., 1 Cd and 1 Se)
	} 
	else {
		printf("Only wurtzite and zincblende nanocrystals are currently implemented!\n");
		printf("The program is exiting!!!\n");
		exit(EXIT_FAILURE);
	}

	// Lattice constant scaling, the lattice constants are defined in nc.h
	latticeConstant = retLatticeConstant(ncCrystalStructure, atomSymbol1, atomSymbol2);

	// Scale and print the resulting vectors
	for (i = 0; i < NDIM; i++) primitiveVectors[i] = retScaledVector(primitiveVectors[i], latticeConstant);
	printf("The primitive lattice vectors are (%s %s%s):\n\n", ncCrystalStructure, atomSymbol1, atomSymbol2);
	for (i = 0; i < NDIM; i++) writeVector(primitiveVectors[i], stdout);
	fflush(stdout);

	return numAtomsInUnitCell; 
}

/****************************************************************************/
