# include "spkmeans.h"

static int MAX_P_ITTER = 300;
static int MAX_R_ITTER = 100;
static double MAX_DIS = 2147483647.0;
static double EPS_P = 0;
static double EPS_R = 0.00001;

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
double ** createMatrix(int rows, int cols) {
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
    double **tmp = createMatrix(rows * 2, cols);
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
    double **tposed = createMatrix(cols, rows);
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
    double **cutMat = createMatrix(rows, cols);
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
    double **mat = createMatrix(dim, dim);
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
    double **prod = createMatrix(dim, dim);
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
 * Kmeans Functions
 ******************/

/*
Func itterates over the points given and sorts them into the different clusters.
returns: 1 on success, zero on failure.
*/
int itterPoints(double **clusters, double **points) {
    int i, closest, *nums_entered;
    double **means = createMatrix(CLUSTERS_NUM, NUMS_PER_POINT);
    if (means == NULL) {return 0;}
    copyMatFromTo(clusters, means, NUMS_PER_POINT, CLUSTERS_NUM);
    do {
        i = 0;
        nums_entered = createDupArray(CLUSTERS_NUM, 0);
        if (nums_entered == NULL) {
            freeMat(means);
            return 0;
        }
        copyMatFromTo(means, clusters, NUMS_PER_POINT, CLUSTERS_NUM);
        zeroMat(means, CLUSTERS_NUM, NUMS_PER_POINT);
        for (; i < POINTS_NUM; i++) {
            closest = closestCluster(clusters, points[i]);
            addPoint(means[closest], points[i]);
            nums_entered[closest]++; 
        }
        divideMat(means, CLUSTERS_NUM, NUMS_PER_POINT, nums_entered);
        free(nums_entered);
    } while (--MAX_P_ITTER > 0 && !didPointsConverge(clusters, means, NUMS_PER_POINT, CLUSTERS_NUM));
    copyMatFromTo(means, clusters, NUMS_PER_POINT, CLUSTERS_NUM);
    freeMat(means);
    return 1;
}

/*Func returns the index of the closest cluster to a point.*/
int closestCluster(double **clusters, double *point) {
    int min_index = 0;
    double min_dis = MAX_DIS , cur_dis;
    int i = 0;
    for (; i < CLUSTERS_NUM; i++) {
        cur_dis = distance(clusters[i], point, NUMS_PER_POINT);
        if (cur_dis < min_dis) {
            min_index = i;
            min_dis = cur_dis;
        }
    }
    return min_index;
}

/*Func returns the distance between two given data points.*/
double distance(double *point1, double *point2, int length) {
    double sum = 0;
    int i = 0;
    for (; i < length; i++) {sum += pow(point1[i] - point2[i], 2);}
    return sqrt(sum);
}

/*Func adds a point to a given cluster-mean.*/
void addPoint(double *cluster, double *point)  {
    int i = 0;
    for (; i < NUMS_PER_POINT; i++) {cluster[i] += point[i];}
}

/*
Func checks for convergence of the means vector to the clusters vector.
return: 1 if converged, 0 if didn't converge.
*/
int didPointsConverge(double **points1, double **points2, int row_length, int num_rows) {
    int i = 0;
    for (;i < num_rows; i++) {if (distance(points1[i], points2[i], row_length) >= EPS_P) {return 0;}}
    return 1;
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
    double **rot, theta, t, c, s;
    int i = 0;
    theta = (mat[maxCol][maxCol] - mat[maxRow][maxRow]) / (2 * mat[maxRow][maxCol]);
    t = (theta >= 0 ? 1:-1) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / (sqrt(1 + pow(t, 2)));
    s = t * c;
    rot = createMatrix(dim, dim);
    if (rot == NULL) {return NULL;}
    zeroMat(rot, dim, dim);
    for (; i < dim; i++) {rot[i][i] = 1;}
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
    double **jacobi = createMatrix(dim + 1, dim);
    if (jacobi == NULL) {return NULL;}
    for (; i < dim; i++) {jacobi[0][i] = vals[i][i];}
    copyMatFromTo(vectors, jacobi + 1, dim ,dim);
    return jacobi;
}

int eigenHur(double **mat, int dim) {
    double curGap, maxGap = 0;
    int maxInd = 0, i = 0;
    for (; i < (dim - 1) / 2; i++) {
        curGap = fabs(mat[0][i] - mat[0][i + 1]);
        if (maxGap < curGap) {
            maxInd = i;
            maxGap = curGap;
        }
    }
    return maxInd + 1;
}

/************************
 * wam, ddg, gl functions
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
    double ** weightedMatrix = createMatrix(n, n);
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
    double ** diagonalMatrix = createMatrix(n, n);
    if (diagonalMatrix == NULL) {return NULL;}
    zeroMat(diagonalMatrix, n, n);
    for (; i < n; i++) {diagonalMatrix[i][i] = findTheDegreeOfaVertex(weightedAdjencyMatrix[i], n);}
    return diagonalMatrix;
}

/*puts the Laplasian matrix in the diagonal matrix.*/
void calcLaplasianMatrix(double ** weightedMatrix, double ** diagonalMatrix, int n) {
    subtractMatrix(diagonalMatrix, weightedMatrix, n);
}

/********************
 * Python-C functions
 ********************/

/*converts PyObject matrix to matrix in C. [r],[c]: rows and columns number.
double ** convertMatrixToC(PyObject * matrix, int r, int c) {
    if matrix == NULL{
        return NULL;
    }
    double ** cMat = createMatrix(r, c);
    double value;  
    for (Py_ssize_t i = 0; i < r; ++i){
        row = PyList_GetItem(matrix, i);
        for (Py_ssize_t j = 0; j < pointLen; ++j) {
            PyObject* node = PyList_GET_ITEM(row, j);
            value = PyFloat_AsDouble(node);
            cMat[i][j] = value;
        }
    }
    return cMat;
}*/

/*converts a list in C to a list as PyObject
PyObject* convertListToPy(double cArray[], size_t size) {
    PyObject *l = PyList_New(size);
    size_t i;
    for (Py_ssize_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(cArray[i]));
    }
    return l;
}*/

/*converts a matrix in C to a matrix as PyObject. frees the c matrix at the end.
PyObject* convertMatrixToPy(double ** cMatrix, size_t r, size_t c) {
    pyObject* pyMatrix = PyList_New(r);
    for (Py_ssize_t i = 0; i < r; i++){
        PyList_SetItem(pyMatrix, i,  convertListToPy(cMatrix[i], c));
    }

    freeMat(cMatrix);
    return pyMatrix;
}*/

/*returns the wanted matrix due to mode. 2- wam. 3 - ddg. 4 - gl.
PyObject* getTheWantedMatrix(PyObject *self, PyObject *args, int mode) {
    int d;
    int n;
    PyObject* dataPointsInput;
    PyObject* point;
    if(!PyArg_ParseTuple(args, "O" ,&dataPointsInput)) {
        return NULL;
    }
    if (!PyList_CheckExact(dataPointsInput)) {
        PyErr_SetString(PyExc_RuntimeError, "Received non-list type object.");
        return NULL;
    }

    n = PyList_GET_SIZE(dataPointsInput);
    point = PyList_GetItem(dataPointsInput, 0);
    d = PyList_GET_SIZE(point);

    double ** dataPoints = convertMatrixToC(PyObject * dataPointsInput, int n, int d);
    double ** weightMatrix = calcWeightedAdjencyMatrix(dataPoints, n , d);
    freeMat(dataPoints);
    if (mode == 2) {return convertMatrixToPy(weightMatrix, n, n);}
    double ** diagonalMatrix = calcDiagonalDegreeMatrix(weightMatrix, n);
    freeMat(weightMatrix);
    if (mode == 3) {return convertMatrixToPy(diagonalMatrix, n, n);}
    calcLaplasianMatrix(eightedMatrix, diagonalMatrix, n);
    if (mode == 4) {return convertMatrixToPy(diagonalMatrix, n, n);}   
}*/

/*
static PyObject* wam(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 2);
}

static PyObject* ddg(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 3);
}

static PyObject* gl(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 4);
}
*/
