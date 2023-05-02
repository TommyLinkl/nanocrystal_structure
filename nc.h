/****************************************************************************/
//
//
//
/****************************************************************************/

#ifndef NC_H
#define NC_H

/****************************************************************************/
/* These are the library functions that are used */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <time.h>
#include <assert.h>
#include "vector.h"
#include "atom.h"

/****************************************************************************/
/* Macro definitions: common unit conversion and multiplication schemes
  that help make the program more readable */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define AUTONM    0.05291772108
#define AUTOANG   0.5291772108
#define ANGTOAU   1.889725988579
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define HBAR      1.0
#define NDIM      3
#define MASS      1.0
#define EPS       1e-10

// All lattice constants are in Angstroms
#define ACdSeWZ 4.2989773 // Cd-Se bond length = 2.62 A, Eg =  eV
#define ACdSeZB 6.0583000 // Cd-Se bond length = 2.62 A, Eg =  
#define ACdSWZ  4.1302281 // Cd-S  bond length = 2.52 A, Eg = 2.50  eV
#define ACdSZB  5.8179695 // Cd-S  bond length = 2.52 A, Eg = 2.50  eV
#define ACdTeWZ 4.60      // TODO: find better value! Cd-Te bond length =  A, Eg =  eV 
#define ACdTeZB 6.48      // Cd-Te bond length = 2.81 A, Eg = 1.474 eV 
// www.semiconductors.co.uk/propiivi5410.htm (all Zn taken from here, 300 K) and above comments are
#define AZnSeWZ 3.98      // c = 6.53, c/a = 1.641,   Zn-Se bond length = 2.45 A, Eg = 
#define AZnSeZB 5.67      //                          Zn-Se bond length = 2.46 A, Eg = 2.8215 eV
#define AZnSWZ  3.811     // c = 6.234, c/a = 1.636,  Zn-S  bond length = 2.34 A, Eg = 3.911  eV
//#define AZnSWZ  4.1302281
#define AZnSZB  5.421     //                          Zn-S  bond length = 2.34 A, Eg = 3.68   eV
#define AZnTeWZ 4.27      // c = 6.99, c/a = 1.637,   Zn-Te bond length = 2.62 A, Eg = 
#define AZnTeZB 6.10      //                          Zn-Te bond length = 2.64 A, Eg = 2.394  eV 
#define AZnOWZ  3.2495    // c = 5.2069, c/a = 1.602, Zn-O  bond length = 1.95 A, Eg = 3.4    eV
#define AInAsZB 6.0570    // In-As bond length = 2.62 A, Eg  =   eV
#define AGaAsZB 5.6535    // Ga-As bond length = 2.45 A, Eg  =   eV
// added by djasras from Fu and Zunger PRB 66, 1642 (1997) 
#define AInPZB  5.826     // In-P bond length is 2.5228 A

/****************************************************************************/
/* Structure declarations */

typedef struct  chiralNC_ {
  int index;
  int type;
  char chiralAtomSymbol[4];
} chiralNC;

typedef struct stackingFault_ {
  int type; // only type 1 or 2 currently allowed
  int position; // indexed layer(s) to be shifted
} stackingFault;

typedef struct ncLayer_ {
  int index, nHalfLayers;
  char metalAtomSymbol[4];
  char nonMetalAtomSymbol[4];
} ncLayer;

typedef struct nanostructure_ {
  int index;
  int nAtoms, nAtomTypes, nMaxBonds;
  int nStackingFaults, nCores;
  char ncType[100], ncCrystalStructure[100];
  char ncAtomSymbol1[4], ncAtomSymbol2[4], centerAtomSymbol[4];
  vector ncSize, ncCenter;
//stackingFault stackFault[5];
//atom *atoms; // pointer to the first atom
} nanostructure;

typedef struct param_ {
  int nAtoms, nAtomTypes, buildNewNC, buildNewBulk;
  int nMaxBonds, nMinBonds;
  int nLayersToCut, passivate, remDanglingAtoms, makeChiralNC;
  int nHalfLayersToAdd;
  int nStackingFaults, nCores, nAttachedNCs, nChiralCenters;
  int chiralFlag, doNotCutCdSurfaceFlag, allCdSurfaceFlag;
  char confFile[100], confFileUnits[100];
  char ncType[100];
  char centerAtomSymbol[10];
  char ncAtomSymbol1[4], ncAtomSymbol2[4];
  nanostructure mainNC, coreNC[10], attachedNC[10]; // maximum of 10 cores
  ncLayer newLayer[100]; // maximum of 100 new half layers
  stackingFault stackFault[5]; // maximum of 5 stacking faults
  chiralNC chiralNanocrystalParams; // only 1 chiral nc type allowed currently
} param;

/****************************************************************************/
/* Function declarations - public interface */

/* Functions that initialize the system - in read.c */
void readInput(param *params);
void readConf(atom *atoms, param *params);
int fillNanocrystalStructure(nanostructure *nc, FILE *pf);
int assignAtomTypes(atom *atoms, int nAtoms);
int isNewAtomType(atom *atoms, int currIndex);
int isAtomAMetal(char *atomSymbol);

/* Functions that calculate nanocrystal size related properties - size.c */
void calcSizeStatistics(atom *atoms, param params);
void calcNanocrystalDimensions(vector *maxDimensions, vector *atomPositions, int numVectors);
double retNanocrystalVolume(vector *atomPositions, int numVectors, char *ncType);
double retNanocrystalArea(vector *atomPositions, int numVectors, char *ncType);
double retMaxPlaneProjectedLength(vector *vectors, int numVectors, char *plane);
double retMaxDistance(vector *vectors, int numVectors);
double retMaxLength(vector *vectors, int numVectors, char axis);

/* Functions that calculate nanocrystal size related properties - bonds.c */
void calcBondStatistics(atom *atoms, param params);
void calcNearestNeighborLists(atom *atoms, int nMaxBonds, int nAtoms);
int areAtomsNearestNeighbors(atom atom1, atom atom2);
double retAverageBondLength(atom *atoms, char *atom1, char *atom2, param params);
double retSTDBondLength(double aveBondLength, atom *atoms, char *atom1, char *atom2, param params);
int allAtomsHaveTwoNeighbors(atom *atoms, int numAtoms, int nMaxBonds);
double retIdealBondLength(char *atom1, char *atom2);
double retAtomMass(char *atomSymbol);

/* Functions that build new nanocrystals (QDs, NRs or NPLs) - build.c */
int buildNewNanocrystal(atom *atoms, param params);
int buildBulkCrystal(atom *atoms, param params);
vector retRequiredSizeVector(vector *primitiveVectors, nanostructure nc, param params);
void centerNanocrystal(atom *atoms, int nAtoms, char *centerAtomSymbol);
int addUnitCellAtoms(atom *atoms, int nOrigAtoms, vector *primitiveVectors, vector origin, param params);
void addStackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, param params);
void addType2StackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, int position, int flag);
void addType1StackingFaults(atom *atoms, int nAtoms, vector *primitiveVectors, int position, int flag);
int fillPrimitiveLatticeVectors(vector *primitiveVectors, char *ncCrystalStructure, char *atomSymbol1, char *atomSymbol2);
int cutRoughNanocrystalFromBulk(atom *atoms, int nAtoms, nanostructure nc, param params);
int cutQDFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double qdDiameter);
int cutNRFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double nrDiameter, double nrLength);
int cutNPLFromBulk(atom *atoms, int nAtoms, int nMaxBonds, double nplLength, double nplWidth, double nplThickness);
int cutAttachedNCsFromBulk(atom *atoms, int nOrigAtoms, param params);
void changeAllCoreAtomTypes(atom *atoms, int nAtoms, param params);
int changeSingleCoreAtomType(atom *atoms, nanostructure mainNC, nanostructure coreNC);
int addAllNewLayers(atom *atoms, param params);
int addTetrohedronAtoms(atom *atoms, char *newAtomSymbol, vector origAtomPos, vector *tetrohedronVectors, int nNewAtoms);
int addCornerAtoms(atom *atoms, char *newAtomSymbol, vector maxDimensions);
int getZLayer(double testAtomZPos, double centerAtomZPos, double cLatticeConstant);
double retLatticeConstant(char *ncCrystalStructure, char *atomSymbol1, char *atomSymbol2);

/* Functions that remove the surface atoms of a nanocrystal - cut.c */
int cutOneLayer(atom *newAtoms, atom *origAtoms, int nOrigAtoms, int nMaxBonds);
int cutDanglingAtoms(atom *newAtoms, atom *atoms, int nMaxBonds, int nMinBonds, int nOrigAtoms);
int hasMinNumBonds(atom *atoms, int atomIndex, int nMaxBonds, int nMinBonds);
int removeAtomsWithMissingBonds(atom *atoms, int nMaxBonds, int nMinBonds, int nOrigAtoms);
int isAtomASurfaceAtom(atom *atoms, int atomIndex, int nMaxBonds);
int newIsAtomASurfaceAtom(atom *atoms, int atomIndex, int nAtoms);
int isAPassivationSymbol(char *atomSymbol);

/* Functions that passivate the surface atoms of a nanocrystal - passivate.c */
int passivateNanocrystal(atom *passivatedAtoms, atom *origAtoms, int nOrigAtoms, int nMaxBonds);
int isCrystalPassivated(atom *atoms, int nAtoms);

/* Functions that make a chiral nanocrystal - chiral.c */
int addChiralCenters(atom *atoms, int nAtoms, int nMaxBonds, chiralNC chiralParams);

/* Functions that write the output - write.c */
void writeInput(atom *atoms, param params, FILE *pf);
void writeConfFile(atom *atoms, int nAtoms, double unitChange, FILE *pf);
void writeXYZFile(atom *atoms, int nAtoms, double unitChange, FILE *pf);
void writeLammpsConfFile(atom *atoms, vector ncDimensions, int nAtoms, int nAtomTypes, FILE *pf);
void writeNearestNeighbors(atom *atoms, param params, FILE *pf);
void writeAllBondsWithPositions(atom *atoms, param params);
void writeSizeResults(vector *maxDimensions, atom *atoms, param params, FILE *pf);
void writeIntegerArray(int *intArray, int arrayLength, char *fileName);
void writeVectorShort(vector vect, FILE *pf);
void writeVector(vector vect, FILE *pf);
void writeCurrentTime(FILE *pf);
void writeSeparation(FILE *pf);
void writeSystemInfo(atom *atoms, param params, FILE *pf);

/* Functions that handle errors/ perform error checks - errorHandling.c */
void memoryError(const char *errorMessage);
void inputError(const char *errorMessage);

/****************************************************************************/

#endif

/****************************************************************************/
