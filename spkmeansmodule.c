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
