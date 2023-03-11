# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# ifndef SPKMEANS_H_
# define SPKMEANS_H_

extern int NUMS_PER_POINT;
extern int POINTS_NUM;
extern int CLUSTERS_NUM;

/*General Functions*/
int * createDupArray(int, int);
void copyDoublePFromTo(double *, double *, int);

/*Matrix Handling Functions*/
double ** createMatrix(int, int);
void copyMatFromTo(double **, double **, int, int);
double ** doubleMatLength(double **, int, int);
void zeroMat (double **, int, int);
void divideMat(double **, int, int, int *);
double ** transposeMat(double **, int, int);
double ** idMat(int);
double ** mulSqMats(double **, double **, int);
void printMat(double **, int, int);
void freeMat(double **);

/*File Processing Functions*/
double ** processFile(char *);
double * processFirstLine(FILE *);
char * getLine(FILE *);
int checkGotLine(char *, int *);
void enterAfter(char *, FILE *);
double processNum(char *, int *);
double ** processPoints(double **, FILE *);
double * processLine(FILE *);

/*Kmeans Functions*/
int itterPoints(double **, double **);
int closestCluster(double **, double *);
double distance(double *, double *, int);
void addPoint(double *, double *);
int didPointsConverge(double **, double **, int, int);

/*Jacobi Functions*/
double ** itterRots(double **, int);
int * maxAbsVal(double **, int);
double ** createRotMat(double **, int, int, int);
int didRotConverge(double **, double **, int);
double calcOff(double **, int);
double ** buildJacobiRet(double **, double **, int);

# endif