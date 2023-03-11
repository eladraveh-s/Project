# include "spkmeans.h"

static int MAX_P_ITTER = 300;
static int MAX_R_ITTER = 100;
static double MAX_DIS = 2147483647.0;
static double EPS_P = 0;
static double EPS_R = 0.00001;

/*******************
 * General FUnctions
 *******************/

/*
Func creates an integer pointer (array) of a given length and value for all cells.
param length: the amount of cells rquired.
param num: the numbr which is to be repaeated in every cell.
return: an int * pointer of size length in which each cell contains num.
*/
int * createDupArray(int length, int num) {
    int *dupArray = (int *) malloc((size_t) (length * sizeof(int)));
    int i = 0;
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
    double *mat = (double *) malloc((size_t) (rows * cols * sizeof(double)));
    double **mat_loc = (double **) malloc((size_t) (rows * sizeof(double *)));
    int i = 0;
    for (; i < rows; i++) {mat_loc[i] = mat + i * cols;}
    return mat_loc;
}

/*
Func copies a double ** matrix (contigous block) into another double ** matrix.
param from: the matrix storing the copied information.
param to: the matrix which the information will be copied to.
param rows: the amount of rows in the matrix.
param cols: tha amount of columns in th matrix.
*/
void copyMatFromTo(double **from, double **to, int rows, int cols) {
    int i = 0;
    for (; i < rows; i++) {copyDoublePFromTo(from[i], to[i], cols);}
}

/*
Func doubles the length of a given matrix while preserving the infomratin.
param points: the given matrix.
param length: the amount of rows in the matrix.
*/
double ** doubleMatLength(double **points, int rows, int cols) {
    double **tmp = createMatrix(rows * 2, cols);
    if (tmp == NULL) {return NULL;}
    copyMatFromTo(points, tmp, rows, cols);
    freeMat(points);
    return tmp;
}

/*Func sets all the values of a matrix to zero*/
void zeroMat(double **mat, int rows, int cols) {
    int i = 0, j;
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {mat[i][j] = 0;}}
}

/*
Func divides each element in a matrix according to the value matching it's row.
param matrix: the matrix we would like to divide.
param divider: an array of ints signifying which row should be divided by which number.
*/
void divideMat(double **matrix, int rows, int cols, int *divider) {
    int i = 0, j;
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {matrix[i][j] /= divider[i];}}
}

double ** transposeMat(double **mat, int rows, int cols) {
    int i = 0, j;
    double **tposed = createMatrix(rows, cols);
    if (tposed == NULL) {return NULL;}
    for (; i < rows; i++) {for (j = 0; j < cols; j++) {tposed[i][j] = mat[j][i];}}
    return tposed;
}

double ** idMat(int dim) {
    double **mat = createMatrix(dim, dim);
    if (mat == NULL) {return NULL;}
    zeroMat(mat, dim, dim);
    for (; dim > 0; dim--) {mat[dim - 1][dim - 1] = 1;}
    return mat;
}

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

/*Func prints a given matrix.*/
void printMat(double **mat, int rows, int cols) {
    int i = 0, j;
    for (; i < rows; i++) {
        for (j = 0; j < cols - 1; j++) {printf("%.4f,", mat[i][j]);}
        printf("%.4f", mat[i][cols - 1]);
        printf("\n");
    }
}

/*Func frees a matrix created using createMat*/
void freeMat(double **mat) {
    free(mat[0]);    
    free(mat);
}

/*subtracts mat2 from mat1 (the result is in mat1, both of matrix are the same dimension [n]).*/
void subtractMatrix(double ** mat1, double ** mat2, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0;  j < n; j++){
            mat1[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
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
        if (nums_entered == NULL) {return 0;}
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

double ** itterRots(double **mat, int dim) {
    double **rotated, **product ,**curRot;
    int *maxDim;
    product = idMat(dim);
    rotated = mat;
    do {
        mat = rotated;
        maxDim = maxAbsVal(mat, dim);
        if (mat[*maxDim][*(maxDim + 1)] == 0) {break;}
        curRot = createRotMat(mat, dim, *maxDim, *(maxDim + 1));
        rotated = mulSqMats(mat, curRot, dim);
        product = mulSqMats(product, curRot, dim);
    } while (--MAX_R_ITTER > 0 && !didRotConverge(mat, rotated, dim));
    freeMat(curRot);
    free(maxDim);
    return buildJacobiRet(rotated, product, dim);
}

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

double ** createRotMat(double **mat, int dim, int maxRow, int maxCol) {
    double **rot, theta, t, c, s;
    int i = 0;
    theta = (mat[maxCol][maxCol] - mat[maxRow][maxRow]) / (2 * mat[maxRow][maxCol]);
    t = (theta >= 0 ? 1:-1) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / (sqrt(1 + pow(t, 2)));
    s = t * c;
    rot = createMatrix(dim, dim);
    zeroMat(rot, dim, dim);
    for (; i < dim; i++) {rot[i][i] = 1;}
    rot[maxRow][maxRow] = rot[maxCol][maxCol] = c;
    rot[maxRow][maxCol] = rot[maxCol][maxRow] = s;
    if (maxRow < maxCol) {rot[maxCol][maxRow] *= -1;}
    else {rot[maxRow][maxCol] *= -1;}
    return rot;
}

int didRotConverge(double **org, double **new, int dim) {
    return calcOff(org, dim) - calcOff(new, dim) <= EPS_R;
}

double calcOff(double **mat, int dim) {
    int i = 0, j;
    double off = 0;
    for (; i < dim; i++) {for (j = 0; j < dim; j++) {if (i != j) {off += pow(mat[i][j], 2);}}}
    return off;
}

double ** buildJacobiRet(double **vals, double **vectors, int dim) {
    int i = 0;
    double **jacobi = createMatrix(dim + 1, dim);
    for (; i < dim; i++) {jacobi[0][i] = vals[i][i];}
    copyMatFromTo(vectors, jacobi + 1, dim ,dim);
    return jacobi;
}

/******************
 * wam, ddg, gl functions
 ******************/

/* gets two points (as lists of doubles - p1, p2), same dimension (d). returns their squared euclidean distance.*/
double SquaredEuclideanDistance(double *p1, double *p2, int d){
    double result = 0.0;
    for (int i = 0 ; i < d; i ++){
        result = result + (p1[i]-p2[i])**2;
    }
    return result;
}

/*gets a list of data points (dataPoints), the list length (n), and their dimension (d). clacs their weighted and returns it.*/
double ** calcWeightedAdjencyMatrix(double ** dataPoints, int n ,int d){
    double ** weightedMatrix = createMatrix(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                weightedMatrix[i][j] = 0.0;
            }
            else{
                weightedMatrix[i][j] = exp(-0.5*(SquaredEuclideanDistance(dataPoints[i], dataPoints[j], d)));
            }
        }
    }
    return weightedMatrix;
}

/*gets the vertex's row from the weighted adjency matrix (vertexRow), and the number of data points (n). returns the vertex's degree as it define. */
double findTheDegreeOfaVertex(double * vertexRow, int n){
    double degree = 0.0;
    for (int i; i < n; i++){
        degree = degree + vertexRow[i];
    }
    return degree;
}

/*gets the weighted Adjency Matrix and the matrix dimension (n), returns te diagonal degree matrix as it describe.*/
double ** calcDiagonalDegreeMatrix(double ** weightedAdjencyMatrix, int n){
    double ** diagonalMatrix = createMatrix(n, n);
    zeroMat(diagonalMatrix, n, n);
    for (int i = 0; i < n; i++){
        diagonalMatrix[i][i] = findTheDegreeOfaVertex(weightedAdjencyMatrix[i]);
    }
    return diagonalMatrix;
}

/*puts the Laplasian matrix in the diagonal matrix.*/
void calcLaplasianMatrix(double ** weightedMatrix, double ** diagonalMatrix, int n){
    subtractMatrix(diagonalMatrix, weightedMatrix, n);
}

/******************
 * Python-C functions
 ******************/

/*converts PyObject matrix to matrix in C. [r],[c]: rows and columns number.*/
double ** convertMatrixToC(PyObject * matrix, int r, int c){
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
}

/*
converts list in C to list as PyObject.
*/
PyObject* convertListToPy(double cArray[], size_t size) {
    PyObject *l = PyList_New(size);
    size_t i;
    for (Py_ssize_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(cArray[i]));
    }
    return l;
}

/*converts matrix in C to matrix as PyObject. frees the c matrix at the end.*/
PyObject* convertMatrixToPy(double ** cMatrix, size_t r, size_t c){
    pyObject* pyMatrix = PyList_New(r);
    for (Py_ssize_t i = 0; i < r; i++){
        PyList_SetItem(pyMatrix, i,  convertListToPy(cMatrix[i], c));
    }

    freeMat(cMatrix);
    return pyMatrix;
}

/*returns the wanted matrix due to mode. 2- wam. 3 - ddg. 4 - gl.*/
PyObject* getTheWantedMatrix(PyObject *self, PyObject *args, int mode){
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
    if mode == 2 {
        return convertMatrixToPy(weightMatrix, n, n);
    }
    double ** diagonalMatrix = calcDiagonalDegreeMatrix(weightMatrix, n);
    freeMat(weightMatrix);
    if mode == 3 {
        return convertMatrixToPy(diagonalMatrix, n, n);
    }
    calcLaplasianMatrix(double ** weightedMatrix, double ** diagonalMatrix, int n);
    if mode == 4 {
        return convertMatrixToPy(diagonalMatrix, n, n);
    }   
}

static PyObject* wam(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 2);
}

static PyObject* ddg(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 3);
}

static PyObject* gl(PyObject *self, PyObject *args) {
    return getTheWantedMatrix(PyObject *self, PyObject *args, 4);
}
