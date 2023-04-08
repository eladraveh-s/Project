# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "spkmeans.h"

int CLUSTERS_NUM = 0;
static int MAX_P_ITTER = 300;
static double MAX_DIS = 2147483647.0;
static double EPS_P = 0;
int COLS = 0;
int ROWS = 0;

/******************
 * Kmeans Functions
 ******************/

/*
Func itterates over the points given and sorts them into the different clusters.
returns: 1 on success, zero on failure.
*/
int itterPoints(double **clusters, double **points) {
    int i, closest, *nums_entered;
    double **means = createMat(CLUSTERS_NUM, COLS);
    if (means == NULL) {return 0;}
    copyMatFromTo(clusters, means, COLS, CLUSTERS_NUM);
    do {
        i = 0;
        nums_entered = createDupArray(CLUSTERS_NUM, 0);
        if (nums_entered == NULL) {
            freeMat(means);
            return 0;
        }
        copyMatFromTo(means, clusters, COLS, CLUSTERS_NUM);
        zeroMat(means, CLUSTERS_NUM, COLS);
        for (; i < ROWS; i++) {
            closest = closestCluster(clusters, points[i]);
            addPoint(means[closest], points[i]);
            nums_entered[closest]++; 
        }
        divideMat(means, CLUSTERS_NUM, COLS, nums_entered);
        free(nums_entered);
    } while (--MAX_P_ITTER > 0 && !didPointsConverge(clusters, means, COLS, CLUSTERS_NUM));
    copyMatFromTo(means, clusters, COLS, CLUSTERS_NUM);
    freeMat(means);
    return 1;
}

/*Func returns the index of the closest cluster to a point.*/
int closestCluster(double **clusters, double *point) {
    double cur_dis, min_dis = MAX_DIS;
    int i = 0, min_index = 0;
    for (; i < CLUSTERS_NUM; i++) {
        cur_dis = distance(clusters[i], point, COLS);
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
    for (; i < COLS; i++) {cluster[i] += point[i];}
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

/********************
 * C-Python Functions 
 ********************/

/*converts PyObject matrix to matrix in C. [r],[c]: rows and columns number.*/
double ** convertMatrixToC(PyObject * matrix, int r, int c) {
    PyObject *row, *node;
    Py_ssize_t i, j;
    double value, **cMat = createMat(r, c);
    if (matrix == NULL) {return NULL;}
    if (cMat == NULL) {return NULL;}
    for (i = 0; i < r; ++i) {
        row = PyList_GetItem(matrix, i);
        for (j = 0; j < c; ++j) {
            node = PyList_GET_ITEM(row, j);
            value = PyFloat_AsDouble(node);
            cMat[i][j] = value;
        }
    }
    return cMat;
}

/*converts a list in C to a list as PyObject*/
PyObject* convertListToPy(double cArray[], size_t size) {
    Py_ssize_t i = 0;
    PyObject *l = PyList_New(size);
    for (; i != size; ++i) {PyList_SET_ITEM(l, i, PyFloat_FromDouble(cArray[i]));}
    return l;
}

/*converts a matrix in C to a matrix as PyObject and frees the c matrix at the end.*/
PyObject* convertMatrixToPy(double ** cMatrix, size_t r, size_t c) {
    Py_ssize_t i = 0;
    PyObject *pyMatrix = PyList_New(r);
    for (; i < r; i++) {PyList_SetItem(pyMatrix, i,  convertListToPy(cMatrix[i], c));}
    freeMat(cMatrix);
    return pyMatrix;
}

/*returns the wanted matrix according to the mode: 1- jacobi, 2- wam, 3 - ddg, 4 - gl, 5 - mat4spk.*/
PyObject * getTheWantedMatrix(PyObject *self, PyObject *args, int mode) {
    int d, n, k;
    double **dataPoints, **weightedMatrix, **diagonalMatrix;
    PyObject *dataPointsInput, *point;
    if (mode == 5) {if (!PyArg_ParseTuple(args, "Oi" ,&dataPointsInput, &k)) {return NULL;}}
    else {if (!PyArg_ParseTuple(args, "O" ,&dataPointsInput)) {return NULL;}}
    if (!PyList_CheckExact(dataPointsInput)) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred (Received non-list type object)");
        return NULL;
    }

    n = PyList_GET_SIZE(dataPointsInput);
    point = PyList_GetItem(dataPointsInput, 0);
    d = PyList_GET_SIZE(point);

    dataPoints = convertMatrixToC(dataPointsInput, n, d);
    if (dataPoints == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    if (mode == 1) {
        weightedMatrix = itterRots(dataPoints, n);
        if (weightedMatrix == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        return convertMatrixToPy(weightedMatrix, n + 1, n);
    }
    weightedMatrix = calcWeightedAdjencyMatrix(dataPoints, n , d);
    if (weightedMatrix == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    freeMat(dataPoints);
    if (mode == 2) {return convertMatrixToPy(weightedMatrix, n, n);}
    diagonalMatrix = calcDiagonalDegreeMatrix(weightedMatrix, n);
    if (diagonalMatrix == NULL) {
        freeMat(weightedMatrix);
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    if (mode > 3) {calcLaplasianMatrix(weightedMatrix, diagonalMatrix, n);}
    if (diagonalMatrix == NULL) {
        freeMat(weightedMatrix);
        PyErr_SetString(PyExc_RuntimeError, "An error has accured.");
        return NULL;
    }
    if (mode == 5) {
        freeMat(weightedMatrix);
        weightedMatrix = itterRots(diagonalMatrix, n);
        freeMat(diagonalMatrix);
        if (weightedMatrix == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        diagonalMatrix = sortMat(weightedMatrix, n + 1, n);
        if (diagonalMatrix == NULL) {
            freeMat(weightedMatrix);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        if (!k) {k = eigenHur(diagonalMatrix[0], n);}
    }
    freeMat(weightedMatrix);
    return convertMatrixToPy(diagonalMatrix, n + (mode == 5 ? 1:0), mode == 5 ? k:n);
}

static PyObject * spk(PyObject *self, PyObject *args) {
    double **points, **clusters;
    PyObject *pointsInput, *clustersInput, *point;

    if(!PyArg_ParseTuple(args, "OO" ,&pointsInput, &clustersInput)) {return NULL;}
    if (!(PyList_CheckExact(pointsInput) && PyList_CheckExact(clustersInput))) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred (Received non-list type object)");
        return NULL;
    }

    ROWS = PyList_GET_SIZE(pointsInput);
    point = PyList_GetItem(pointsInput, 0);
    COLS = PyList_GET_SIZE(point);
    CLUSTERS_NUM = PyList_GET_SIZE(clustersInput);

    points = convertMatrixToC(pointsInput, ROWS, COLS);
    if (points == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    clusters = convertMatrixToC(clustersInput, CLUSTERS_NUM, COLS);
    if (clusters == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        free(points);
        return NULL;
    }
    if (itterPoints(clusters, points)) {return convertMatrixToPy(clusters, CLUSTERS_NUM, COLS);}
    else {return NULL;}
}

static PyObject * jacobi(PyObject *self, PyObject *args) {return getTheWantedMatrix(self, args, 1);}

static PyObject * wam(PyObject *self, PyObject *args) {return getTheWantedMatrix(self, args, 2);}

static PyObject * ddg(PyObject *self, PyObject *args) {return getTheWantedMatrix(self, args, 3);}

static PyObject * gl(PyObject *self, PyObject *args) {return getTheWantedMatrix(self, args, 4);}

static PyObject * mat4spk(PyObject *self, PyObject *args) {return getTheWantedMatrix(self, args, 5);}

static PyMethodDef spkOptions[] = {
    {"spk", (PyCFunction) spk, METH_VARARGS,           
      PyDoc_STR("spk(matrix points, matrix clusters). Clustering given data points to given number of clusters, based on given centroids. Returns the new centroids after the clustering.")},
    {"jacobi", (PyCFunction) jacobi, METH_VARARGS,           
      PyDoc_STR("jacobi(matrix points). Returns a matrix in which the firs line contains the eigen values of the given symmetrical matrix, and the rest are it's eigen vectors.")},
    {"wam", (PyCFunction) wam, METH_VARARGS,           
      PyDoc_STR("wam(matrix points). Returns a diagonal degree matrix based upon the given points.")},
    {"ddg", (PyCFunction) ddg, METH_VARARGS,           
      PyDoc_STR("ddg(matrix points). Returns a weighted adjacency matrix based upon the given points.")},
    {"gl", (PyCFunction) gl, METH_VARARGS,           
      PyDoc_STR("gl(matrix points). Returns a graph laplacian matrix based upon the given points.")},
    {"mat4spk", (PyCFunction) mat4spk, METH_VARARGS,           
      PyDoc_STR("mat4spk(matrix points). Returns the U matrix.")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    spkOptions /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m = PyModule_Create(&spkmeansmodule);
    if (!m) {return NULL;}
    return m;
}
