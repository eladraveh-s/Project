# include "spkmeans.h"

int LINE_SIZE = 10;
int NUMS_PER_POINT = 0;
int POINTS_NUM = 0;
int CLUSTERS_NUM = 0;

/***************************
 * File Processing Functions
 ***************************/

/*
Func processes the points file.
param path: the path of the points file.
return: the processed points matrix.
*/
double ** processFile(char *path) {
    double *point, **points;
    FILE *file = fopen(path, "r");
    if (file == NULL) {return NULL;}
    point = processFirstLine(file); 
    if (point == NULL) {
        fclose(file);
        return NULL;
    }
    points = createMatrix(10, NUMS_PER_POINT);
    if (points == NULL) {
        free(point);
        fclose(file);
        return NULL;
    }
    copyDoublePFromTo(point, points[0], NUMS_PER_POINT);
    points = processPoints(points, file);
    free(point);
    fclose(file);
    return points;
}

/*
Func processes the first line of the file to NUMS_PER_POINT and suspected LINE_SIZE.
param file: the pointer to the file object.
return: double *point - the first point in the file.
*/
double * processFirstLine(FILE *file) {
    char *line;
    int *indp, size = 2;
    double *point = (double *) malloc((size_t) (size * sizeof(double)));
    if (point == NULL) {return NULL;}
    line = getLine(file);
    if (line == NULL) {
        free(point);
        return NULL;
    }
    indp = (int *) malloc(sizeof(int));
    if (indp == NULL) {
        free(point);
        free(line);
        return NULL;
    }
    *indp = 0;
    while (line[(*indp)] != '\0') {
        point[NUMS_PER_POINT++] = processNum(line, indp);
        if (size == NUMS_PER_POINT) {
            point = (double *) realloc(point, (size_t) ((size *= 2) * sizeof(double)));
            if (point == NULL) {
                free(indp);
                free(line);
                return NULL;
            }
        }
        (*indp)++;
    }
    free(indp);
    free(line);
    return point;
}

/*
Func recieves a single line from the given file.
param file: the pointer to the file object.
return: given char * of line.
*/
char * getLine(FILE *file) {
    char *pointer, *line;
    int *indp = (int *) malloc(sizeof(int));
    if (indp == NULL) {return NULL;}
    *indp = 0;
    line = (char *) malloc((size_t) (LINE_SIZE * sizeof(char)));
    if (line == NULL) {
        free(indp);
        return NULL;
    }
    pointer = fgets(line, LINE_SIZE, file);
    if (pointer == NULL) {
        free(line);
        free(indp);
        return NULL;
    }
    while (!checkGotLine(line, indp)) {
        line = (char *) realloc(line, (size_t) ((LINE_SIZE * 2) * sizeof(char)));
        if (line == NULL) {
            free(pointer);
            free(indp);
            return NULL;
        }
        enterAfter(line + *indp, file);
        LINE_SIZE *= 2;
    }
    free(indp);
    return line;
}

/*
Func checks for the success of the memory allocatin needed to recieve a line.
param line: points to the recieved line - saved as a char array.
param indp: a pointer to the last checked index of the array.
return: 1 on success, 0 on failure.
*/
int checkGotLine(char *line, int *indp) {
    for (; *indp < LINE_SIZE; (*indp)++) {if (line[*indp] == '\0') {return line[*indp - 1] == '\n';}}
    return 0;
}

/*
Func adds the rest of the line to the information recieved if the line size was too small.
param line: points to the end of the recieved line.
param file: the pointer to the file object.
*/
void enterAfter(char *line, FILE *file) {line = fgets(line, LINE_SIZE, file);}

/*
Func computes the numerical value of a double given as a string.
param line: the line from which we try to read the number.
param indp: the index of the first char which is a part of the number.
return: the numerical value of a double given as a string.
*/
double processNum(char *line, int *indp) {
    int size = 8, i = 0, start = *indp;
    char *num_str = (char *) malloc((size_t) (size * sizeof(char)));
    double num;
    /*if (num_str == NULL) {return NULL;}*/
    while (line[*indp] != ',' && line[*indp] != '\n') {
        num_str[i] = line[(*indp)++];
        if (size == ++i) {num_str = (char *) realloc(num_str, (size_t) ((size *= 2) * sizeof(char)));}
    }
    num_str = (char *) realloc(num_str, (size_t) ((*indp - start + 1) * sizeof(char)));
    num_str[*indp - start] = '\0';
    num = atof(num_str);
    free(num_str);
    return num;
}

/*
Func processes the points given in the file to the given pointers.
param clusters: a pointer to the (empty) clusters matrix.
param nums_entered: an integer array signifyign how many numbers went into each centroid.
param file: the pointer to the file object.
return: double ** matrix containing all of the points in the file.
*/
double ** processPoints(double **points, FILE *file) {
    int length = 10;
    double *curPoint;
    while ((curPoint = processLine(file)) != NULL) {
        if (POINTS_NUM++ == length - 1) {
            points = doubleMatLength(points, length, NUMS_PER_POINT);
            if (points == NULL) {
                free(curPoint);
                return NULL;
            }
            length *= 2;
        }
        copyDoublePFromTo(curPoint, points[POINTS_NUM], NUMS_PER_POINT);
        free(curPoint);
    }
    POINTS_NUM++;
    return points;
}

/*
Func processes a line of input from the given file.
param file: the pointer to the file object.
return: double * point - a numerical represntation of the line.
*/
double * processLine(FILE *file) {
    double *point;
    int *indp, i = 0;
    char *line = getLine(file);
    if (line == NULL) {return NULL;}
    point = (double *) malloc((size_t) (NUMS_PER_POINT * sizeof(double)));
    if (point == NULL) {
        free(line);
        return NULL;
    }
    indp = (int *) malloc(sizeof(int));
    if (indp == NULL) {
        free(line);
        free(point);
        return NULL;
    }
    *indp = 0;
    for (; i < NUMS_PER_POINT; i++) {
        point[i] = processNum(line, indp);
        (*indp)++;
    }
    free(line);
    free(indp);
    return point;
}

int main(int argc, char *argv[]) {
    double **points, **mat, **helpMat;
    if (argc != 3) {
        printf("The arguments given don't match the requierments!\n");
        return 0;
    }
    if (strcmp(argv[1], "jacobi") && strcmp(argv[1], "wam") && strcmp(argv[1], "ddg") && strcmp(argv[1], "gl")) {
        printf("Invalid goal!\n");
        return 0;
    }
    points = processFile(argv[2]);
    if (points == NULL) {
        printf("An error has accured\n");
        return 0;
    }
    if (!strcmp(argv[1], "jacobi")) {mat = itterRots(points, POINTS_NUM);}
    else {
        mat = calcWeightedAdjencyMatrix(points, POINTS_NUM, NUMS_PER_POINT);
        if (strcmp(argv[1], "wam")) {
            helpMat = calcDiagonalDegreeMatrix(mat, POINTS_NUM);
            if (!strcmp(argv[1], "gl")) {calcLaplasianMatrix(mat, helpMat, POINTS_NUM);}
            if (mat != NULL) {freeMat(mat);}
            mat = helpMat;
        }
    }
    printMat(mat, POINTS_NUM + (strcmp(argv[1], "jacobi") ? 0:1), POINTS_NUM);
    freeMat(points);
    if (mat != NULL) {freeMat(mat);}
    else {return 0;}
    return 1;
}
