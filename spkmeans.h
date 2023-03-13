# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

int MEM_ERR = 0;
static int MAX_R_ITTER = 100;
static double EPS_R = 0.00001;

/*General Functions*/
int * createDupArray(int, int);
void copyDoublePFromTo(double *, double *, int);

/*Matrix Handling Functions*/
double ** createMat(int, int);
void copyMatFromTo(double **, double **, int, int);
double ** doubleMatLength(double **, int, int);
void zeroMat (double **, int, int);
void divideMat(double **, int, int, int *);
double ** transposeMat(double **, int, int);
double ** cutMat(double **, int, int);
double ** sortMat(double **, int, int);
double ** idMat(int);
void subtractMatrix(double **, double **, int n);
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

/*Jacobi Functions*/
double ** itterRots(double **, int);
int * maxAbsVal(double **, int);
double ** createRotMat(double **, int, int, int);
int didRotConverge(double **, double **, int);
double calcOff(double **, int);
double ** buildJacobiRet(double **, double **, int);
int eigenHur(double *, int);

/*wam, ddg, gl Functions*/
double squaredEuclideanDistance(double *, double *, int);
double ** calcWeightedAdjencyMatrix(double **, int ,int);
double findTheDegreeOfaVertex(double *, int);
double ** calcDiagonalDegreeMatrix(double **, int);
void calcLaplasianMatrix(double **, double **, int);

/*kmeans Functions*/
int itterPoints(double **, double **);
int closestCluster(double **, double *);
double distance(double *, double *, int);
void addPoint(double *, double *);
int didPointsConverge(double **, double **, int, int);

/*******************
 * General Functions
 *******************/

/*
Func creates an integer pointer (array) of a given length and value for all cells.
param length: the amount of cells rquired.
param num: the numbr which is to be repaeated in every cell.
return: an int * pointer of size length in which each cell contains num.
*/
int * createDupArray(int length, int num) {
    int i = 0, *dupArray = (int *) malloc((size_t) (length * sizeof(int)));
    if (dupArray == NULL) {return NULL;}
    for (; i < length; i++) {dupArray[i] = num;}
    return dupArray;
}

/*
Func copies a double * pointer into another double *  until a given index.
param from: the array storing the copied information.
param to: the array which the information will be copied to.
param length: the size of the arrays.
*/
void copyDoublePFromTo(double *from, double *to, int length) {
    int i = 0;
    for (; i < length; i++) {to[i] = from[i];}
}

/***************************
 * Matrix Handling Functions
 ***************************/

/*
Func creates double ** pointer to a contigous block (matrix) of memory.
param rows: the ammount of required rows.
param cols: the amount of required cols.
return: the double ** pointer to the matrix.
*/
double ** createMat(int rows, int cols) {
    double *mat, **mat_loc;
    int i = 0;
    mat = (double *) malloc((size_t) (rows * cols * sizeof(double)));
    if (mat == NULL) {return NULL;}
    mat_loc = (double **) malloc((size_t) (rows * sizeof(double *)));
    if (mat_loc == NULL) {
        free(mat);
        return NULL;
    }
    for (; i < rows; i++) {mat_loc[i] = mat + i * cols;}
    return mat_loc;
}

/*
Func copies a double ** matrix (contigous block) into another double ** matrix.
param from: the matrix storing the copied information.
param to: the matrix which the information will be copied to.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
*/
void copyMatFromTo(double **from, double **to, int rows, int cols) {
    int i = 0;
    for (; i < rows; i++) {copyDoublePFromTo(from[i], to[i], cols);}
}

/*
Func doubles the length of a given matrix while preserving the infomratin.
param points: the given matrix.
param length: the amount of rows in the matrix.
return: A pointer to the doubled length matrix.
*/
double ** doubleMatLength(double **points, int rows, int cols) {
    double **tmp = createMat(rows * 2, cols);
    if (tmp == NULL) {return NULL;}
    copyMatFromTo(points, tmp, rows, cols);
    freeMat(points);
    return tmp;
}

/*
Func sets all the values of a matrix to zero.
param matrix: the matrix we would like to zero.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
*/
void zeroMat(double **mat, int rows, int cols) {
    int i = 0, j;
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {mat[i][j] = 0;}}
}

/*
Func divides each element in a matrix according to the value matching it's row.
param matrix: the matrix we would like to divide.
param divider: an array of ints signifying which row should be divided by which number.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
*/
void divideMat(double **matrix, int rows, int cols, int *divider) {
    int i = 0, j;
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {matrix[i][j] /= divider[i];}}
}

/*
Func transposes a matrix.
param mat: the matrix we would like to zero.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
return: the transposed matrix.
*/
double ** transposeMat(double **mat, int rows, int cols) {
    int i = 0, j;
    double **tposed = createMat(cols, rows);
    if (tposed == NULL) {return NULL;}
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {tposed[j][i] = mat[i][j];}}
    return tposed;
}

/*
Func creates a new matrix from the first param rows rows and param cols columns of mat.
param mat: the matrix we would like to cut.
param rows: the amount of rows in the new matrix.
param cols: the amount of columns in the new matrix.
return: the new cut matrix.
*/
double ** cutMat(double **mat, int rows, int cols) {
    double **cutMat = createMat(rows, cols);
    if (cutMat == NULL) {return NULL;}
    copyMatFromTo(mat, cutMat, rows, cols);
    return cutMat;
}

/*
Func sorts a matrix using it's first row.
param mat: the matrix we would like to sort.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
return: the sorted matrix.
*/
double ** sortMat(double **mat, int rows, int cols) {
    int i, j;
    double *tmp, **sorted, **toSort = transposeMat(mat, rows, cols);
    if (toSort == NULL) {return NULL;}
    for (i = 1; i < cols; i++) {
        tmp = toSort[i];
        j = i;
        while (j > 0 && toSort[j - 1][0] > tmp[0]) {
            toSort[j] = toSort[j - 1];
            j -= 1;
        }
        toSort[j] = tmp;
    }
    sorted = transposeMat(toSort, cols, rows);
    freeMat(toSort);
    return sorted;
}

/*Func creates and returns an identity matrix of deimension dim*/
double ** idMat(int dim) {
    double **mat = createMat(dim, dim);
    if (mat == NULL) {return NULL;}
    zeroMat(mat, dim, dim);
    for (; dim > 0; dim--) {mat[dim - 1][dim - 1] = 1;}
    return mat;
}

/*subtracts mat2 from mat1 (the result is in mat1, both of the matrices are of the same dimension [n]).*/
void subtractMatrix(double ** mat1, double ** mat2, int n) {
    int i = 0, j;
    for (; i < n; i++) {for (j = 0;  j < n; j++) {mat1[i][j] = mat1[i][j] - mat2[i][j];}}
}

/*
Func implements matrix multiplication for teo square matrices - A x B = P.
param mat1: the left matrix in the multiplication process (A).
param mat2: the right matrix in the multiplication process (B).
return: the product of the multiplication process (P).
*/
double ** mulSqMats(double **mat1, double **mat2, int dim) {
    int i = 0, j, cur;
    double **prod = createMat(dim, dim);
    if (prod == NULL) {return NULL;}
    zeroMat(prod, dim, dim);
    for (; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            for (cur = 0; cur < dim; cur++) {prod[i][j] += mat1[i][cur] * mat2[cur][j];}
        }
    }
    return prod;
}

/*
Func prints a given matrix.
param mat: the matrix we would like to print.
param rows: the amount of rows in the matrix.
param cols: the amount of columns in the matrix.
*/
void printMat(double **mat, int rows, int cols) {
    int i = 0, j;
    if (mat == NULL) {printf("An error has accured\n");}
    else {
        for (; i < rows; i++) {
            for (j = 0; j < cols - 1; j++) {printf("%.4f,", mat[i][j]);}
            printf("%.4f", mat[i][cols - 1]);
            printf("\n");
        }
    }
}

/*Func frees a matrix created using createMat*/
void freeMat(double **mat) {
    free(mat[0]);    
    free(mat);
}

/******************
 * Jacobi Functions 
 ******************/

/*
Func rotates a given symetrical matrix until it's diagonal/converged/max iteration has passed.
param mat: a symetrical matrix to be diagonalized.
param dim: the dimension of the matrix.
return: The first line is the eigen values, from there on it's the product of the rotation matrices.
*/
double ** itterRots(double **mat, int dim) {
    int *maxDim, exCode;
    double **rotated, **curRot, **prev, **product = idMat(dim);
    if (product == NULL) {return NULL;}
    prev = idMat(dim);
    if (prev == NULL) {
        freeMat(product);
        return NULL;
    }
    rotated = mat;
    exCode = 0;
    do {
        freeMat(prev);
        prev = rotated;
        maxDim = maxAbsVal(mat, dim);
        if (maxDim == NULL) {
            exCode = 1;
            break;
        }
        if (mat[*maxDim][*(maxDim + 1)] == 0) {break;}
        curRot = createRotMat(mat, dim, *maxDim, *(maxDim + 1));
        if (curRot == NULL) {
            exCode = 2;
            break;
        }
        rotated = mulSqMats(mat, curRot, dim);
        if (rotated == NULL) {
            exCode = 3;
            break;
        }
        product = mulSqMats(product, curRot, dim);
        if (product == NULL) {
            exCode = 4;
            break;
        }
        freeMat(curRot);
    } while (--MAX_R_ITTER > 0 && !didRotConverge(prev, rotated, dim));
    if (exCode != 1) {free(maxDim);}
    if (exCode != 2 && exCode != 0) {freeMat(curRot);}
    if (exCode == 3) {freeMat(product);}
    if (exCode == 4) {freeMat(rotated);}
    if (exCode == 2 || exCode == 1) {free(prev);}
    if (exCode != 0) {return NULL;}
    return buildJacobiRet(rotated, product, dim);
}

/*
Func funds the location of the off axis element with the largest absolute value.
param mat: the matrix we want to find the element from.
param dim: the dimension of the matrix.
return: the row and column of the element.
*/
int * maxAbsVal(double **mat, int dim) {
    int i = 0, j;
    int *maxDim = calloc((size_t) 2, sizeof(int));
    if (maxDim == NULL) {return NULL;}
    for (; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if ((i != j) && (fabs(mat[i][j]) >= fabs(mat[*maxDim][*(maxDim + 1)]))) {
                *maxDim = i;
                *(maxDim + 1) = j;
            }
        }
    }
    return maxDim;
}

/*
Func creates a rotation matrix for the jacobi algorithm.
param mat: the matrix we want to rotate.
param dim: the dimension of the matrix.
param maxRow: the row of the off axis element with the largest absolute value.
param maxCol: the column of the off axis element with the largest absolute value.
return: the rotation matrix.
*/
double ** createRotMat(double **mat, int dim, int maxRow, int maxCol) {
    double **rot, t, c, s, theta = (mat[maxCol][maxCol] - mat[maxRow][maxRow]) / (2 * mat[maxRow][maxCol]);
    t = (theta >= 0 ? 1:-1) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / (sqrt(1 + pow(t, 2)));
    s = t * c;
    rot = idMat(dim);
    if (rot == NULL) {return NULL;}
    rot[maxRow][maxRow] = rot[maxCol][maxCol] = c;
    rot[maxRow][maxCol] = rot[maxCol][maxRow] = s;
    if (maxRow < maxCol) {rot[maxCol][maxRow] *= -1;}
    else {rot[maxRow][maxCol] *= -1;}
    return rot;
}

/*
Func checks for the convergence of the rotation process.
param org: the matrix before the rotation.
param new: the matrix after the rotation.
param dim: the dimension of the matrices.
return: 1 if converged, 0 otherwise.
*/
int didRotConverge(double **org, double **new, int dim) {return calcOff(org, dim) - calcOff(new, dim) <= EPS_R;}

/*
Func claculates the off parameter for convergence checking.
param mat: the matrix we want to calculate the off of.
param dim: the dimension of the matrix.
return: the off value.
*/
double calcOff(double **mat, int dim) {
    int i = 0, j;
    double off = 0;
    for (; i < dim; i++) {for (j = 0; j < dim; j++) {if (i != j) {off += pow(mat[i][j], 2);}}}
    return off;
}

/*
Func builds the matrix the jacobi algorithm returns from the given values.
param vals: the matrix after the diagonlization process.
param vectors: the product of the rotation matrices.
param dim: the dimension of the matrices.
return: the jacobi algorithm's wanted return value.
*/
double ** buildJacobiRet(double **vals, double **vectors, int dim) {
    int i = 0;
    double **jacobi = createMat(dim + 1, dim);
    if (jacobi == NULL) {return NULL;}
    for (; i < dim; i++) {jacobi[0][i] = vals[i][i];}
    copyMatFromTo(vectors, jacobi + 1, dim ,dim);
    return jacobi;
}

/*
Func finds the largest gap between 2 eigen values (numbers in the first line) and returns it's index + 1.
param row: the first row of the matrix (the row in which we want to find the lragest gap).
param length: the length of the row.
return: the index of the gap + 1.
*/
int eigenHur(double *row, int length) {
    double curGap, maxGap = 0;
    int maxInd = 0, i = 0;
    for (; i < (length - 1) / 2; i++) {
        curGap = fabs(row[i] - row[i + 1]);
        if (maxGap < curGap) {
            maxInd = i;
            maxGap = curGap;
        }
    }
    return maxInd + 1;
}

/************************
 * wam, ddg, gl Functions
 ************************/

/* gets two points (as lists of doubles - p1, p2), same dimension (d). returns their squared euclidean distance.*/
double squaredEuclideanDistance(double *p1, double *p2, int d) {
    int i = 0;
    double result = 0.0;
    for (; i < d; i++) {result += pow((p1[i] - p2[i]), 2);}
    return result;
}

/*gets a list of data points (dataPoints), the list length (n), and their dimension (d). clacs their weighted and returns it.*/
double ** calcWeightedAdjencyMatrix(double ** dataPoints, int n ,int d) {
    int i = 0, j;
    double ** weightedMatrix = createMat(n, n);
    if (weightedMatrix == NULL) {return NULL;}
    for (; i < n; i++) {
        for (j = 0; j < n; j++){
            if (i == j) {weightedMatrix[i][j] = 0.0;}
            else {
                weightedMatrix[i][j] = exp(-0.5 * (squaredEuclideanDistance(dataPoints[i], dataPoints[j], d)));
            }
        }
    }
    return weightedMatrix;
}

/*gets the vertex's row from the weighted adjency matrix (vertexRow), and the number of data points (n). returns the vertex's degree as it define. */
double findTheDegreeOfaVertex(double * vertexRow, int n) {
    int i = 0;
    double degree = 0.0;
    for (; i < n; i++) {degree = degree + vertexRow[i];}
    return degree;
}

/*gets the weighted Adjency Matrix and the matrix dimension (n), returns te diagonal degree matrix as it describe.*/
double ** calcDiagonalDegreeMatrix(double ** weightedAdjencyMatrix, int n) {
    int i = 0;
    double ** diagonalMatrix = createMat(n, n);
    if (diagonalMatrix == NULL) {return NULL;}
    zeroMat(diagonalMatrix, n, n);
    for (; i < n; i++) {diagonalMatrix[i][i] = findTheDegreeOfaVertex(weightedAdjencyMatrix[i], n);}
    return diagonalMatrix;
}

/*puts the Laplasian matrix in the diagonal matrix.*/
void calcLaplasianMatrix(double ** weightedMatrix, double ** diagonalMatrix, int n) {
    subtractMatrix(diagonalMatrix, weightedMatrix, n);
}
