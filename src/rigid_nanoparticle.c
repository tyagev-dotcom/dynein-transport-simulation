#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

#ifdef _WIN32
#define OS_WINDOWS
#include <windows.h>
#else
#define OS_LINUX
#include <unistd.h> // For sysconf
#include <stdint.h>
#endif

// ----------------------------
// Global Parameters
// ----------------------------



#define numMotors 6              // Number of grafted motors
#define numTrials 10            // Number of nanoparticle simulation runs
#define radiusMT  (12.5e-9)         // Microtubule radius
#define radiusNC  20e-9             // Nanoparticle radius
#define effTemp 1900                // Effective temperature for energy scaling
#define dynein_height 60e-9         // Dynein height
#define dyenin_diameter 10e-9       // Dynein diameter
#define dyenin_diameter_path 12e-9  // Dynein path width
#define dyenin_pers_length 2e-6    // Persistence length of dynein (semi-flexible polymer)
#define meshLimit_walk 102          // Size of walk mesh array
#define meshLimit_att 49            // Size of attachment mesh array
#define dataSize 80000              // Max number of iterations (steps/attach/detach)
#define dataIntervalSize 500        // Max size for time interval analysis

// Gradient Descent Parameters
#define Angle_LEARNING_RATE 1e19
#define Step_Learning_Rate 1e3
#define MAX_GRADIENT_ITER 100
#define GRADIENT_THRESHOLD 1e-27


//Dynein Gradient decsent variables
#define DYENIN_N 10
#define min_delta_E  1e-27 // Tunable threshold for energy change
#define  DYNEIN_MAX_ITER  30
#define  DYNEIN_LEARNING_RATE 1e16


// ----------------------------
// Global Variables
// ----------------------------

double PolymerLength = 4.256e-8; //4.256e-8
double aPolymer = 0.75e-9;
double Na;                    // Will be calculated at runtime
double R0_Globe;              // Will be calculated at runtime
double PI = 3.14159265358979323846;
double attEnergy = 3.2899e-20; // Dynein to MT attachment energy gain
double Temper = 298;          // Physical temperature (K)
double kb = 1.38e-23;         // Boltzmann constant
double fs = 3.6e-12;          // Stalling force

// Simulation parameters

double springRelax;           // Typically R0/sqrt(3)
double tauStep = 0.01088;     // Average step duration [s]
double MTBSxShift = 8e-9;     // MT Binding Site spacing in X

// Binding site angular spacing (13 MTBS around MT)
double MTBSphiShift = (2 * 3.14159265358979323846 / 13); // There are 13 MTBS around the MT
double dF;
double q0;


int numThreads;   // number of available threads for the simulation 
struct ThreadData* threadData;
pthread_mutex_t mutex, task_mutex;
pthread_t* threads;

#ifdef _WIN32
LONG trial = 0;
#define ATOMIC_FETCH_AND_ADD(ptr, val) InterlockedExchangeAdd(ptr, val)
#else
int trial = 0;
#define ATOMIC_FETCH_AND_ADD(ptr, val) __sync_fetch_and_add(ptr, val)
#endif

// ----------------------------
// Mesh Arrays (Attachment / Walk / Probabilities)
// ----------------------------

// Attachment mesh points
double Xmesh_att[] = { -2.213867e-08 ,-2.245691e-08 ,-2.312864e-08 ,-2.400000e-08 ,-2.487136e-08 ,-2.554309e-08 ,-2.586133e-08 ,-1.413867e-08 ,-1.445691e-08 ,-1.512864e-08 ,-1.600000e-08 ,-1.687136e-08 ,-1.754309e-08 ,-1.786133e-08 ,-6.138671e-09 ,-6.456905e-09 ,-7.128644e-09 ,-8.000000e-09 ,-8.871356e-09 ,-9.543095e-09 ,-9.861329e-09 ,1.861329e-09 ,1.543095e-09 ,8.713559e-10 ,0 ,-8.713559e-10 ,-1.543095e-09 ,-1.861329e-09 ,9.861329e-09 ,9.543095e-09 ,8.871356e-09 ,8.000000e-09 ,7.128644e-09 ,6.456905e-09 ,6.138671e-09 ,1.786133e-08 ,1.754309e-08 ,1.687136e-08 ,1.600000e-08 ,1.512864e-08 ,1.445691e-08 ,1.413867e-08 ,2.586133e-08 ,2.554309e-08 ,2.487136e-08 ,2.400000e-08 ,2.312864e-08 ,2.245691e-08 ,2.213867e-08 };
double Ymesh_att[] = { -1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 };

// Stepping mesh points and probabilities
double Xmesh_walk[] = { -6.077684e-08,  -5.277684e-08,  -4.893779e-08,  -4.000000e-08,  -3.693779e-08,  -3.122316e-08,  -2.843905e-08,  -2.800000e-08,  -2.706221e-08,  -2.477684e-08,  -2.306221e-08,  -2.093779e-08,  -2.077684e-08,  -2.043905e-08,  -2.000000e-08,  -1.956095e-08,  -1.922316e-08,  -1.906221e-08,  -1.693779e-08,  -1.677684e-08,  -1.600000e-08,  -1.556095e-08,  -1.506221e-08,  -1.293779e-08,  -1.277684e-08,  -1.243905e-08,  -1.200000e-08,  -1.156095e-08,  -1.122316e-08,  -1.106221e-08,  -8.937794e-09,  -8.776835e-09,  -8.439048e-09,  -8.000000e-09,  -7.560952e-09,  -7.223165e-09,  -7.062206e-09,  -4.937794e-09,  -4.776835e-09,  -4.439048e-09,  -4.000000e-09,  -3.223165e-09,  -3.062206e-09,  -9.377940e-10,  9.377940e-10,  3.062206e-09,  3.223165e-09,  3.560952e-09,  4.000000e-09,  4.439048e-09,  4.776835e-09,  4.937794e-09,  7.062206e-09,  7.223165e-09,  7.560952e-09,  8.000000e-09,  8.439048e-09,  8.776835e-09,  8.937794e-09,  1.106221e-08,  1.122316e-08,  1.156095e-08,  1.200000e-08,  1.243905e-08,  1.277684e-08,  1.293779e-08,  1.506221e-08,  1.600000e-08,  1.677684e-08,  1.693779e-08,  1.906221e-08,  2.000000e-08,  2.043905e-08,  2.306221e-08,  2.356095e-08,  2.477684e-08,  2.706221e-08,  2.800000e-08,  3.293779e-08,  3.506221e-08,  4.000000e-08,  4.306221e-08,  4.400000e-08,  5.906221e-08,  6.077684e-08,  5.277684e-08,  4.893779e-08,  3.693779e-08,  3.122316e-08,  2.843905e-08,  2.093779e-08,  2.077684e-08,  1.956095e-08,  1.922316e-08,  1.556095e-08,  -3.560952e-09,  -2.356095e-08,  -3.293779e-08,  -3.506221e-08,  -4.306221e-08,  -4.400000e-08,  -5.906221e-08 };
double Ymesh_walk[] = { -5.140000e-09,  -5.140000e-09,  -6.205000e-09,  0,  -6.205000e-09,  5.140000e-09,  -2.905000e-09,  0,  6.205000e-09,  -5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  0,  2.905000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  0,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  0,  2.905000e-09,  -6.205000e-09,  -2.905000e-09,  5.140000e-09,  -6.205000e-09,  0,  6.205000e-09,  -6.205000e-09,  0,  -6.205000e-09,  0,  -6.205000e-09,  5.140000e-09,  5.140000e-09,  6.205000e-09,  6.205000e-09,  -5.140000e-09,  2.905000e-09,  6.205000e-09,  5.140000e-09,  -2.905000e-09,  -5.140000e-09,  -2.905000e-09,  2.905000e-09,  2.905000e-09,  -6.205000e-09,  6.205000e-09,  6.205000e-09,  0,  6.205000e-09 };
double xy_WalkProb[] = { 2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  7.832898e-03,  7.832898e-03,  5.221932e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  5.221932e-03,  5.221932e-03,  2.610966e-03,  5.221932e-03,  7.832898e-03,  2.088773e-02,  7.832898e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  2.610966e-02,  3.655352e-02,  1.305483e-02,  1.566580e-02,  1.044386e-02,  1.566580e-02,  7.832898e-03,  4.177546e-02,  7.571802e-02,  5.221932e-03,  1.044386e-02,  1.827676e-02,  7.832898e-03,  1.566580e-02,  2.610966e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  1.305483e-02,  7.832898e-03,  1.305483e-02,  7.832898e-03,  8.093995e-02,  6.527415e-02,  1.044386e-02,  2.088773e-02,  1.305483e-02,  2.610966e-02,  7.832898e-03,  4.177546e-02,  2.872063e-02,  7.832898e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  5.221932e-03,  1.044386e-02,  2.349869e-02,  2.610966e-03,  2.610966e-03,  2.610966e-03,  7.832898e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  7.832898e-03,  7.832898e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  1.305483e-02,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03 };


// ----------------------------
// Structures
// ----------------------------

struct MotorData {
    double motors[numMotors][6]; // [attached?, phi, theta, x_MT, y_MT, unused]
    double distance2MT;
    double IterationDeltaX, IterationDeltaY;
    double TotalDeltaX, TotalDeltaY;
};

struct BalanceMovment {
    struct MotorData NewLocations;
    double Energy;
};

struct DeltaTime {
    double Probability[3]; // [stepRate, SumDetach, SumAttach]
    double attachRate[numMotors][meshLimit_att];
    double detachRate[numMotors];
    double AttachProb[numMotors];
};

struct SimulationResults {
    double** IterationMotorsCount;
    double* TrailXtotal;
    double* TrailYtotal;
    double* TrailTimetotal;
    double** Yvelocity, ** Xvelocity, ** XvelocityPlus, ** XvelocityMinus;
    double** XIntervalMSD, ** YvelocityPlus, ** YvelocityMinus;
    double** HelicalPitch;
};

struct SingleTrialResults {
    double IterationMotorsCount[dataSize];
    double TrailXtotal, TrailYtotal, TrailTimetotal;
    double Yvelocity[dataIntervalSize];
    double Xvelocity[dataIntervalSize];
    double XvelocityPlus[dataIntervalSize];
    double XvelocityMinus[dataIntervalSize];
    double XIntervalMSD[dataIntervalSize];
    double YvelocityPlus[dataIntervalSize];
    double YvelocityMinus[dataIntervalSize];
    double HelicalPitch[dataIntervalSize];
};

struct VariableStatistics {
    double mean, STD, SEM;
};

struct Statistics {
    struct VariableStatistics X, Y, Time, Motors;
    struct VariableStatistics VelocityX, VelocityXplus, VelocityXminus;
    struct VariableStatistics VelocityY, VelocityYplus, VelocityYminus;
    struct VariableStatistics HelicalPitch;
};

struct ThreadData {
    int* PositiveTotal;
    int* NegativeTotal;
    double TimeIntervalsLength;
    struct SimulationResults* Results;
    int* trial;
};

// ----------------------------
// Function Declarations
// ----------------------------

void allocateSimulationResults(struct SimulationResults* Results);
void freeSimulationResults(struct SimulationResults* Results);
void spreadMotors(struct MotorData* Motors, unsigned int* seed);
void rotateMotor(struct MotorData* Motors, double xRot, double yRot, double zRot);
int attMotors(struct MotorData Motors);
double energy(struct MotorData Motors);
_Bool excludedVolume(int ChoosenMotor, double yNextLocation, double xNextLocation, double motors[numMotors][6]);
_Bool isPathClear(int steppingMotor, double xCurr, double yCurr, double xNext, double yNext, double motors[numMotors][6]);
struct DeltaTime deltaTime(struct MotorData Motors);
int monteCarloSelect(double* values, int size, unsigned int* seed);
int selectAttachedMotor(struct MotorData);
void copyMotorDataStruct(struct MotorData* dest, struct MotorData* src);
double CalculateEnergyGradient(struct MotorData Motors, int axis);
struct BalanceMovment GradientDescent(struct MotorData Motors);
void isLimiteisExceeded(double* ymtRot, double* DistShift, struct MotorData Motors, double* xRot, double* yRot);
void SingleMotorAlign(struct MotorData* Motors, int whoattach);
void initializeThreads();
void cleanupThreads();
void assignThreadRanges(int numItems, int numThreads, int* Thread_active, int* ActionsPerThread);
double newRand(unsigned int* seed);
double mean_non_zero_doubles(const double* list, int size);
double mean_non_zero_doubles_2d(double** array, int rows, int cols);
double std_non_zero_doubles(const double* list, int size, double mean);
double std_non_zero_doubles_2d(double** array, int rows, int cols, double mean);
double sem_non_zero_doubles_2d(double** array, int rows, int cols, double std_dev);
double compute_msd_and_b_value(double** Results_XIntervalMSD, double timeIntervalLength, double timeMean, double* outMSDlist);
void* simulate_trial(void* arg);

// Main entry point
int main(void);

// ----------------------------
// Utility Functions
// ----------------------------

// Generate a pseudo-random number in the range [0, 1]
double newRand(unsigned int* seed) {
    *seed = (*seed * 1103515245 + 12345) & 0x7fffffff;  // 31-bit RNG
    return (double)(*seed) / 0x7fffffff;
}

// Initialize thread-related resources
void initializeThreads() {
    threadData = (struct ThreadData*)malloc(numThreads * sizeof(struct ThreadData));
    threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));

    if (!threadData || !threads) {
        fprintf(stderr, "Error: Unable to allocate memory.\n");
        exit(EXIT_FAILURE);
    }

    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&task_mutex, NULL);
}

// Clean up thread-related memory and synchronization objects
void cleanupThreads() {
    free(threadData);
    free(threads);
    pthread_mutex_destroy(&mutex);
}

// Allocate memory for all fields of SimulationResults
void allocateSimulationResults(struct SimulationResults* Results) {
    Results->TrailYtotal = calloc(numTrials, sizeof(double));
    Results->TrailXtotal = calloc(numTrials, sizeof(double));
    Results->TrailTimetotal = calloc(numTrials, sizeof(double));

    Results->YvelocityPlus = calloc(numTrials, sizeof(double*));
    Results->YvelocityMinus = calloc(numTrials, sizeof(double*));
    Results->Yvelocity = calloc(numTrials, sizeof(double*));
    Results->Xvelocity = calloc(numTrials, sizeof(double*));
    Results->XvelocityPlus = calloc(numTrials, sizeof(double*));
    Results->XvelocityMinus = calloc(numTrials, sizeof(double*));
    Results->IterationMotorsCount = calloc(numTrials, sizeof(double*));
    Results->HelicalPitch = calloc(numTrials, sizeof(double*));
    Results->XIntervalMSD = calloc(numTrials, sizeof(double*));

    for (int i = 0; i < numTrials; i++) {
        Results->YvelocityPlus[i] = calloc(dataIntervalSize, sizeof(double));
        Results->YvelocityMinus[i] = calloc(dataIntervalSize, sizeof(double));
        Results->Yvelocity[i] = calloc(dataIntervalSize, sizeof(double));
        Results->Xvelocity[i] = calloc(dataIntervalSize, sizeof(double));
        Results->XvelocityPlus[i] = calloc(dataIntervalSize, sizeof(double));
        Results->XvelocityMinus[i] = calloc(dataIntervalSize, sizeof(double));
        Results->HelicalPitch[i] = calloc(dataIntervalSize, sizeof(double));
        Results->IterationMotorsCount[i] = calloc(dataSize, sizeof(double));
        Results->XIntervalMSD[i] = calloc(dataIntervalSize, sizeof(double));
    }
}

// Free all dynamically allocated memory in SimulationResults
void freeSimulationResults(struct SimulationResults* Results) {
    for (int i = 0; i < numTrials; i++) {
        free(Results->YvelocityPlus[i]);
        free(Results->YvelocityMinus[i]);
        free(Results->Yvelocity[i]);
        free(Results->Xvelocity[i]);
        free(Results->XvelocityPlus[i]);
        free(Results->XvelocityMinus[i]);
        free(Results->IterationMotorsCount[i]);
        free(Results->HelicalPitch[i]);
        free(Results->XIntervalMSD[i]);
    }

    free(Results->Yvelocity);
    free(Results->Xvelocity);
    free(Results->XvelocityPlus);
    free(Results->XvelocityMinus);
    free(Results->YvelocityPlus);
    free(Results->YvelocityMinus);
    free(Results->HelicalPitch);
    free(Results->XIntervalMSD);
    free(Results->IterationMotorsCount);
    free(Results->TrailYtotal);
    free(Results->TrailXtotal);
    free(Results->TrailTimetotal);
}

// Spread motors uniformly over the nanoparticle surface and attach the first one
void spreadMotors(struct MotorData* Motors, unsigned int* seed) {
    double phi, teta;

    for (int i = 0; i < numMotors; i++) {
        phi = newRand(seed) * 2 * PI;
        teta = newRand(seed) * PI;

        Motors->motors[i][0] = 0;  // Not attached
        Motors->motors[i][1] = radiusNC * sin(teta) * cos(phi); // x
        Motors->motors[i][2] = radiusNC * sin(teta) * sin(phi); // y
        Motors->motors[i][3] = radiusNC * cos(teta);            // z
        Motors->motors[i][4] = 0;
        Motors->motors[i][5] = 0;
    }

    // Attach the first motor
    Motors->motors[0][0] = 1;

    // Initialize NP distance from MT and align first motor to (0, 0, -R_NP)
    Motors->distance2MT = springRelax + radiusNC;
    SingleMotorAlign(Motors, 0);

    Motors->IterationDeltaX = 0;
    Motors->IterationDeltaY = 0;
    Motors->TotalDeltaX = 0;
    Motors->TotalDeltaY = 0;
}

// Copy contents from one MotorData structure to another
void copyMotorDataStruct(struct MotorData* dest, struct MotorData* src) {
    for (int i = 0; i < numMotors; i++) {
        for (int j = 0; j < 6; j++) {
            dest->motors[i][j] = src->motors[i][j];
        }
    }
}

// Rotate each motor in the MotorData structure around X, Y, and Z axes
void rotateMotor(struct MotorData* Motors, double xRot, double yRot, double zRot) {
    double x1, y1, z1, x2, y2, z2;
    for (int i = 0; i < numMotors; i++) {
        // X-axis rotation
        x1 = Motors->motors[i][1];
        y1 = Motors->motors[i][2] * cos(xRot) - Motors->motors[i][3] * sin(xRot);
        z1 = Motors->motors[i][2] * sin(xRot) + Motors->motors[i][3] * cos(xRot);

        // Y-axis rotation
        x2 = x1 * cos(yRot) + z1 * sin(yRot);
        y2 = y1;
        z2 = -x1 * sin(yRot) + z1 * cos(yRot);

        // Z-axis rotation
        Motors->motors[i][1] = x2 * cos(zRot) - y2 * sin(zRot);
        Motors->motors[i][2] = x2 * sin(zRot) + y2 * cos(zRot);
        Motors->motors[i][3] = z2;
    }
}

// ----------------------------
// Energy and Attachment Calculation
// ----------------------------

// Count how many motors are currently attached
int attMotors(struct MotorData Motors) {
    int numAttachedMotors = 0;
    for (int i = 0; i < numMotors; i++) {
        if (Motors.motors[i][0]) {
            numAttachedMotors++;
        }
    }
    return numAttachedMotors;
}




// Utility: normalize 3D vector
void normalize(double v[3]) {
    double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len > 1e-12) {
        v[0] /= len; v[1] /= len; v[2] /= len;
    }
}

// Plane basis from dynein and polymer vectors
void define_plane_basis(double V_dyenin[3], double V_polymer[3], double e1[3], double e2[3]) {
    e1[0] = V_dyenin[0]; e1[1] = V_dyenin[1]; e1[2] = V_dyenin[2];
    normalize(e1);

    double dot = e1[0] * V_polymer[0] + e1[1] * V_polymer[1] + e1[2] * V_polymer[2];
    double proj[3] = { dot * e1[0], dot * e1[1], dot * e1[2] };
    double perp[3] = { V_polymer[0] - proj[0], V_polymer[1] - proj[1], V_polymer[2] - proj[2] };
    normalize(perp);
    e2[0] = perp[0]; e2[1] = perp[1]; e2[2] = perp[2];
}

// Project 3D point onto 2D plane coordinates
void project_point(double point[3], double e1[3], double e2[3], double result[2]) {
    result[0] = point[0] * e1[0] + point[1] * e1[1] + point[2] * e1[2];
    result[1] = point[0] * e2[0] + point[1] * e2[1] + point[2] * e2[2];
}

// Compute zPrime in 2D
double compute_zPrime2D(double P_attach2D[2], double V_polymer2D[2]) {
    return (1.0 / radiusNC) * (P_attach2D[0] * V_polymer2D[0] + P_attach2D[1] * V_polymer2D[1]);
}

double compute_total_energy(double angles[DYENIN_N - 1], double P_base2D[2], double P_attach2D[2], double phi0) {
    double kappa = kb * Temper * dyenin_pers_length;
    double phi = phi0;
    double pos[2] = { P_base2D[0], P_base2D[1] };  // Start at actual dynein base
    double SEGMENT_LENGTH = dynein_height / DYENIN_N;

    for (int i = 0; i < DYENIN_N; i++) {
        if (i > 0) phi += angles[i - 1];
        pos[0] += SEGMENT_LENGTH * cos(phi);
        pos[1] += SEGMENT_LENGTH * sin(phi);
    }

    // Polymer linker vector from dynein tip to nanoparticle attach
    double V_polymer[2] = {
        P_attach2D[0] - pos[0],
        P_attach2D[1] - pos[1]
    };
    double linker_length = sqrt(V_polymer[0] * V_polymer[0] + V_polymer[1] * V_polymer[1]);
    double M = linker_length * linker_length;

    double Ebend = 0.0;
    for (int i = 0; i < DYENIN_N - 1; i++) {
        Ebend += (kappa / SEGMENT_LENGTH) * (1.0 - cos(angles[i]));
    }

    double zPrime = fabs(compute_zPrime2D(P_attach2D, V_polymer));

    double Estretch = -attEnergy
        - kb * Temper * (
            log((9.0 * zPrime * pow(aPolymer, 3)) / (2 * PI * pow(R0_Globe, 4)))
            - 1.5 * (M / (R0_Globe * R0_Globe))
            );

    return Ebend + Estretch;
}

double minimize_total_energy(double P_attach[3], double P_base[3], double P_tip[3], double angles[DYENIN_N - 1]) {

    double V_dyenin[3] = { P_tip[0] - P_base[0],
                           P_tip[1] - P_base[1],
                           P_tip[2] - P_base[2] };

    double V_polymer[3] = { P_attach[0] - P_tip[0],
                            P_attach[1] - P_tip[1],
                            P_attach[2] - P_tip[2] };

    double e1[3], e2[3];
    define_plane_basis(V_dyenin, V_polymer, e1, e2);

    double P_attach2D[2], P_tip2D[2], P_base2D[2];
    project_point(P_attach, e1, e2, P_attach2D);
    project_point(P_tip, e1, e2, P_tip2D);
    project_point(P_base, e1, e2, P_base2D);

    double V_base_tip_2D[2] = { P_tip2D[0] - P_base2D[0], P_tip2D[1] - P_base2D[1] };
    double phi0 = atan2(V_base_tip_2D[1], V_base_tip_2D[0]);

    for (int i = 0; i < DYENIN_N - 1; i++) angles[i] = 0.0;

    double initial_E = compute_total_energy(angles, P_base2D, P_attach2D, phi0);
    double E_prev = initial_E;

    // === Energy history ===
    double energy_history[DYNEIN_MAX_ITER + 1];
    int energy_count = 0;
    energy_history[energy_count++] = E_prev;


    for (int iter = 0; iter < DYNEIN_MAX_ITER; iter++) {
        double grad[DYENIN_N - 1];
        double delta = 1e-6;

        for (int j = 0; j < DYENIN_N - 1; j++) {
            double saved = angles[j];
            angles[j] += delta;
            double E2 = compute_total_energy(angles, P_base2D, P_attach2D, phi0);
            grad[j] = (E2 - E_prev) / delta;
            angles[j] = saved;
        }

        for (int j = 0; j < DYENIN_N - 1; j++) {
            angles[j] -= DYNEIN_LEARNING_RATE * grad[j];
        }

        double E_new = compute_total_energy(angles, P_base2D, P_attach2D, phi0);

        energy_history[energy_count++] = E_new;


        if (fabs(E_new - E_prev) < min_delta_E)
            break;
        E_prev = E_new;
    }

    double final_E = compute_total_energy(angles, P_base2D, P_attach2D, phi0);
    if (final_E > initial_E) {
        final_E = initial_E;
    }


    return final_E;
}
// Calculate total system energy from all attached motors
double energy_discrete(struct MotorData Motors) {
    double result = 0, result0 = 0;
    for (int i = 0; i < numMotors; i++) {
        if (Motors.motors[i][0]) {

            double x = Motors.motors[i][1], y = Motors.motors[i][2], z = Motors.motors[i][3];
            double Xmt = Motors.motors[i][4], Ymt = Motors.motors[i][5];

            double zMT = -Motors.distance2MT - ((radiusMT + dynein_height) - sqrt(pow((radiusMT + dynein_height), 2) - pow((Ymt * (radiusMT + dynein_height) / radiusMT), 2)));

            double angles[DYENIN_N - 1];


            // Example configuration
            double P_attach[3] = { x, y, z };
            double P_tip[3] = { Xmt, Ymt * (radiusMT + dynein_height) / radiusMT, zMT };

            double P_base[3] = { Xmt, Ymt, zMT - (sqrt(dynein_height * dynein_height - pow((P_tip[1] - Ymt),2))) };


            result += minimize_total_energy(P_attach, P_base, P_tip, angles);


            if (isnan(result)) result = 0;


        }
    }

    return result;
}
double energy(struct MotorData Motors)
{
    double result = 0, zMT, M, zPrime, R0;
    double x, y, z, Xmt, Ymt;

    for (int i = 0; i < numMotors; i++)
    {

        if (Motors.motors[i][0])
        {   // Convert the angle phi and teta angles into cartesian coordinate
            x = Motors.motors[i][1];
            y = Motors.motors[i][2];
            z = Motors.motors[i][3];
            Xmt = Motors.motors[i][4];
            Ymt = Motors.motors[i][5];
            R0 = R0_Globe; //R0 = R0_list[i] ? R0_list[i] : R0_Globe;


            zMT = -Motors.distance2MT - ((radiusMT + dynein_height) - sqrt(pow((radiusMT + dynein_height), 2) - pow((Ymt * (radiusMT + dynein_height) / radiusMT), 2))); //Z axis distance from NP center to end of polymer
            M = pow((Xmt - x), 2) + pow(((Ymt * (radiusMT + dynein_height) / radiusMT) - y), 2) + pow((zMT - z), 2); //Total length of polymer
            zPrime = (1 / radiusNC) * (x * (Xmt - x) + y * ((Ymt * (radiusMT + dynein_height) / radiusMT) - y) + z * (zMT - z));

            if (zPrime < 0)
                zPrime *= -1;

            result += -attEnergy - kb * Temper * (log((9 * zPrime * pow(aPolymer, 3)) / (2 * PI * pow(R0, 4))) - 1.5 * (M / pow(R0, 2)));

            if (isnan(result))
            {
                result = 0;
            }
        }
    }
    return result;
}
// Check and reset invalid limits after rotation and shifting
void isLimiteisExceeded(double* ymtRot, double* DistShift, struct MotorData Motors, double* xRot, double* yRot) {
    for (int i = 0; i < numMotors; i++) {
        if (Motors.motors[i][0]) {
            if (fabs(fmod((*ymtRot + asin(Motors.motors[i][5] / radiusMT)), 2 * PI)) > PI / 2) {
                *ymtRot = 0;
            }
        }
    }
    if ((*DistShift + Motors.distance2MT) > (springRelax + radiusNC) || (*DistShift + Motors.distance2MT) < radiusNC) {
        *DistShift = 0;
    }
}

// Align a single attached motor to the bottom of the nanoparticle
void SingleMotorAlign(struct MotorData* Motors, int whoattach) {
    double SinglemotorYrotate = atan(-Motors->motors[whoattach][1] / Motors->motors[whoattach][3]);
    double SinglemotorXrotate = atan(Motors->motors[whoattach][2] / Motors->motors[whoattach][3]);

    rotateMotor(Motors, SinglemotorXrotate, SinglemotorYrotate, 0);

    if (Motors->motors[whoattach][3] > 0) {
        rotateMotor(Motors, PI, 0, 0); // Flip 180 degrees if necessary
    }
}

// Perform gradient descent to minimize energy and compute bead balance
struct BalanceMovment GradientDescent(struct MotorData Motors) {
    double xRot = 0, yRot = 0, zRot = 0, xmtShift = 0, ymtRot = 0, disShift = 0;
    double gradx = 0, grady = 0, gradz = 0, gradxmt = 0, gradymt = 0, graddis = 0;
    double prevEnergy, currentEnergy, momentum_x = 0, momentum_y = 0, momentum_z = 0;
    double momentum_xmt = 0, momentum_ymt = 0, momentum_dis = 0;
    double beta = 0.9; // Momentum decay
    double DeltaX = 0, DeltaY = 0;
    struct BalanceMovment BeadBalance;

    int nAttached = attMotors(Motors);

    if (nAttached == 0) {
        DeltaX = 0;
        DeltaY = 0;
        currentEnergy = dF;
    }
    else if (nAttached == 1) {
        for (int i = 0; i < numMotors; i++) {
            if (Motors.motors[i][0]) {
                DeltaX = Motors.motors[i][4];
                DeltaY = asin(Motors.motors[i][5] / radiusMT);
                Motors.motors[i][4] = 0;
                Motors.motors[i][5] = 0;
                Motors.distance2MT = springRelax + radiusNC;
                SingleMotorAlign(&Motors, i);
                currentEnergy = dF;
            }
        }
    }
    else {
        currentEnergy = energy(Motors);

        for (int iter = 0; iter < MAX_GRADIENT_ITER; iter++) {
            gradx = CalculateEnergyGradient(Motors, 1);
            grady = CalculateEnergyGradient(Motors, 2);
            gradz = CalculateEnergyGradient(Motors, 3);
            gradxmt = CalculateEnergyGradient(Motors, 4);
            gradymt = CalculateEnergyGradient(Motors, 5);
            graddis = CalculateEnergyGradient(Motors, 6);

            // Apply momentum
            momentum_x = beta * momentum_x + (1 - beta) * gradx;
            momentum_y = beta * momentum_y + (1 - beta) * grady;
            momentum_z = beta * momentum_z + (1 - beta) * gradz;
            momentum_xmt = beta * momentum_xmt + (1 - beta) * gradxmt;
            momentum_ymt = beta * momentum_ymt + (1 - beta) * gradymt;
            momentum_dis = beta * momentum_dis + (1 - beta) * graddis;

            // Update parameters
            xRot = -Angle_LEARNING_RATE * momentum_x;
            yRot = -Angle_LEARNING_RATE * momentum_y;
            zRot = -Angle_LEARNING_RATE * momentum_z;
            xmtShift = -Step_Learning_Rate * momentum_xmt;
            ymtRot = -Angle_LEARNING_RATE * 1e-2 * momentum_ymt;
            disShift = -Step_Learning_Rate * 1e-2 * momentum_dis;

            isLimiteisExceeded(&ymtRot, &disShift, Motors, &xRot, &yRot);
            rotateMotor(&Motors, xRot, yRot, zRot);

            for (int i = 0; i < numMotors; i++) {
                Motors.motors[i][4] += xmtShift;
                Motors.motors[i][5] = radiusMT * sin(ymtRot + asin(Motors.motors[i][5] / radiusMT));
            }
            Motors.distance2MT += disShift;

            prevEnergy = currentEnergy;
            currentEnergy = energy(Motors);

            DeltaX -= xmtShift;
            DeltaY -= ymtRot;

            if (fabs(prevEnergy - currentEnergy) < GRADIENT_THRESHOLD) {

                break;
            }
        }
    }

    Motors.IterationDeltaX = DeltaX;
    Motors.IterationDeltaY = DeltaY;
    Motors.TotalDeltaX += DeltaX;
    Motors.TotalDeltaY += DeltaY;

    BeadBalance.NewLocations = Motors;
    BeadBalance.Energy = energy_discrete(Motors);
    return BeadBalance;
}

// Compute numerical gradient of energy with respect to a parameter (rotation or translation)
double CalculateEnergyGradient(struct MotorData Motors, int axis) {
    double energyPlus, energyMinus = energy(Motors), h = 1e-6;
    struct MotorData Motorscopy;
    copyMotorDataStruct(&Motorscopy, &Motors);
    Motorscopy.distance2MT = Motors.distance2MT;

    if (axis == 1) rotateMotor(&Motorscopy, h, 0, 0);
    else if (axis == 2) rotateMotor(&Motorscopy, 0, h, 0);
    else if (axis == 3) rotateMotor(&Motorscopy, 0, 0, h);
    else if (axis == 4) {
        h = 1e-12;
        for (int i = 0; i < numMotors; i++) Motorscopy.motors[i][4] += h;
    }
    else if (axis == 5) {
        for (int i = 0; i < numMotors; i++) Motorscopy.motors[i][5] = radiusMT * sin(h + asin(Motorscopy.motors[i][5] / radiusMT));
    }
    else if (axis == 6) {
        h = 1e-12;
        Motorscopy.distance2MT += h;
    }

    energyPlus = energy(Motorscopy);
    return (energyPlus - energyMinus) / h;
}

// Check if proposed position violates excluded volume constraint
_Bool excludedVolume(int ChoosenMotor, double yNextLocation, double xNextLocation, double motors[numMotors][6]) {
    for (int i = 0; i < numMotors; i++) {
        if (i != ChoosenMotor && motors[i][0]) {
            double dx = xNextLocation - motors[i][4];
            double dy = yNextLocation - motors[i][5];
            double dz = sqrt(radiusMT * radiusMT - yNextLocation * yNextLocation) - sqrt(radiusMT * radiusMT - motors[i][5] * motors[i][5]);
            double distance = sqrt(dx * dx + dy * dy + dz * dz);
            if (distance <= dyenin_diameter) return false;
        }
    }
    return true;
}

// Check if the stepping path is clear from obstructions
_Bool isPathClear(int steppingMotor, double xCurr, double yCurr, double xNext, double yNext, double motors[numMotors][6]) {
    int numSteps = 100;
    for (int i = 0; i < numMotors; i++) {
        if (i != steppingMotor && motors[i][0]) {
            double xM = motors[i][4];
            double yM = motors[i][5];
            double zM = sqrt(radiusMT * radiusMT - yM * yM);
            for (int step = 0; step <= numSteps; step++) {
                double t = (double)step / numSteps;
                double xStep = xCurr + t * (xNext - xCurr);
                double yStep = yCurr + t * (yNext - yCurr);
                double zStep = sqrt(radiusMT * radiusMT - yStep * yStep);
                double distance = sqrt(pow(xM - xStep, 2) + pow(yM - yStep, 2) + pow(zM - zStep, 2));
                if (distance <= dyenin_diameter_path) return false;
            }
        }
    }
    return true;
}

// Assign a work range to each thread
void assignThreadRanges(int numItems, int numThreads, int* Thread_active, int* ActionsPerThread) {
    *Thread_active = (numItems < numThreads) ? numItems : numThreads;
    *ActionsPerThread = (numItems < numThreads) ? 1 : numItems / numThreads;
}

// Select a value based on Monte Carlo probability distribution
int monteCarloSelect(double* values, int size, unsigned int* seed) {
    double total = 0.0;
    for (int i = 0; i < size; i++) total += values[i];
    double randNum = newRand(seed) * total;
    double cumulativeSum = 0.0;
    for (int i = 0; i < size; i++) {
        cumulativeSum += values[i];
        if (randNum < cumulativeSum) return i;
    }
    return size - 1; // Fallback
}

// Randomly select an attached motor
int selectAttachedMotor(struct MotorData Motors) {
    int numAttached = attMotors(Motors);
    if (numAttached == 0) return -1;
    int* attachedMotors = (int*)malloc(numAttached * sizeof(int));
    if (!attachedMotors) return -1;
    int count = 0;
    for (int i = 0; i < numMotors; i++) {
        if (Motors.motors[i][0]) attachedMotors[count++] = i;
    }
    int selectedMotor = attachedMotors[rand() % count];
    free(attachedMotors);
    return selectedMotor;
}

// Compute stepping, detachment, and attachment probabilities
struct DeltaTime deltaTime(struct MotorData Motors) {
    struct DeltaTime dTime = { 0 };
    struct BalanceMovment BeadBalance;
    int NumMotorsAttach = attMotors(Motors);
    double i_Energy = (NumMotorsAttach == 1) ? dF : energy_discrete(Motors);

    dTime.Probability[0] = (double)NumMotorsAttach / tauStep;

    // Detachment rates
    for (int i = 0; i < numMotors; i++) {
        if (Motors.motors[i][0]) {
            if (NumMotorsAttach == 1) {
                dTime.detachRate[i] = 1;
            }
            else {
                Motors.motors[i][0] = 0;
                BeadBalance = GradientDescent(Motors);
                dTime.detachRate[i] = 1 / (q0 * (1 + exp((BeadBalance.Energy - i_Energy) / (kb * Temper))));
                Motors.motors[i][0] = 1;
            }
            dTime.Probability[1] += dTime.detachRate[i];
        }
    }

    // Attachment rates
    for (int i = 0; i < numMotors; i++) {
        if (!Motors.motors[i][0] && Motors.motors[i][3] < 0) {
            for (int k = 0; k < meshLimit_att; k++) {
                double PhiMeshShift = asin(Ymesh_att[k] / radiusMT) - fmod(Motors.TotalDeltaY, MTBSphiShift);
                Motors.motors[i][0] = 1;
                Motors.motors[i][4] = Xmesh_att[k] - fmod(Motors.TotalDeltaX, MTBSxShift);
                Motors.motors[i][5] = radiusMT * sin(PhiMeshShift);

                if (excludedVolume(i, Motors.motors[i][5], Motors.motors[i][4], Motors.motors) && fabs(PhiMeshShift) <= PI / 2) {
                    BeadBalance = GradientDescent(Motors);
                    double attachRate = 1 / (q0 * (1 + exp((BeadBalance.Energy - i_Energy) / (kb * Temper))));
                    dTime.attachRate[i][k] = attachRate;
                    dTime.AttachProb[i] += attachRate;
                    dTime.Probability[2] += attachRate;
                }

                Motors.motors[i][0] = 0;
                Motors.motors[i][4] = 0;
                Motors.motors[i][5] = 0;
            }
        }
    }

    return dTime;
}


// Evaluate the probability of each possible walk step
void walkAlgoritm(int whoWalks, struct MotorData Motors, double walkProb[meshLimit_walk]) {
    int NumMotorsAttach = attMotors(Motors);
    double i_Energy = (NumMotorsAttach == 1) ? dF : energy_discrete(Motors);
    double TempX = Motors.motors[whoWalks][4];
    double TempY = Motors.motors[whoWalks][5];

    for (int k = 0; k < meshLimit_walk; k++) {
        Motors.motors[whoWalks][4] += Xmesh_walk[k];
        Motors.motors[whoWalks][5] -= Ymesh_walk[k] * cos(asin(Motors.motors[whoWalks][5] / radiusMT));

        if (Motors.motors[whoWalks][5] < -radiusMT || Motors.motors[whoWalks][5] > radiusMT) {
            walkProb[k] = 0;
        }
        else {
            struct BalanceMovment BeadBalance = GradientDescent(Motors);
            double prob = 1 / (1 + exp((-fs * Xmesh_walk[k] + BeadBalance.Energy - i_Energy) / (kb * effTemp))) * xy_WalkProb[k];
            walkProb[k] = prob;
        }

        Motors.motors[whoWalks][4] = TempX;
        Motors.motors[whoWalks][5] = TempY;
    }
}


// ----------------------------
// Statitstics Functions 
// ----------------------------

// Compute the slope of the MSD curve in log-log scale (b-value)
double compute_b_value(int N, double TimeIntervalsLength, double msd_values[]) {
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (int i = 0; i < N; i++) {
        double x = log(TimeIntervalsLength * (i + 1));
        double y = log(msd_values[i]);
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }
    return (N * sum_xy - sum_x * sum_y) / (N * sum_x2 - sum_x * sum_x);
}

// Compute mean of non-zero entries (1D)
double mean_non_zero_doubles(const double* list, int size) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < size; i++) {
        if (list[i] != 0.0) {
            sum += list[i];
            count++;
        }
    }
    return (count > 0) ? sum / count : 0.0;
}

// Compute mean of non-zero entries (2D)
double mean_non_zero_doubles_2d(double** array, int rows, int cols) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (array[i][j] != 0.0) {
                sum += array[i][j];
                count++;
            }
        }
    }
    return (count > 0) ? sum / count : 0.0;
}

// Compute standard deviation of non-zero values (1D)
double std_non_zero_doubles(const double* list, int size, double mean) {
    double variance_sum = 0.0;
    int count = 0;
    for (int i = 0; i < size; i++) {
        if (list[i] != 0.0) {
            double diff = list[i] - mean;
            variance_sum += diff * diff;
            count++;
        }
    }
    return (count > 1) ? sqrt(variance_sum / (count - 1)) : 0.0;
}

// Compute standard deviation of non-zero values (2D)
double std_non_zero_doubles_2d(double** array, int rows, int cols, double mean) {
    double sum_squared_diff = 0.0;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (array[i][j] != 0.0) {
                double diff = array[i][j] - mean;
                sum_squared_diff += diff * diff;
                count++;
            }
        }
    }
    return (count > 1) ? sqrt(sum_squared_diff / (count - 1)) : 0.0;
}

// Compute standard error of the mean (SEM) from std and non-zero count
double sem_non_zero_doubles_2d(double** array, int rows, int cols, double std_dev) {
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (array[i][j] != 0.0) count++;
        }
    }
    return (count > 0) ? std_dev / sqrt(count) : 0.0;
}

// Calculate mean, STD, and SEM for simulation results
void calculate_statistics(struct SimulationResults Results, struct Statistics* Stat) {
    Stat->X.mean = mean_non_zero_doubles(Results.TrailXtotal, numTrials);
    Stat->Y.mean = mean_non_zero_doubles(Results.TrailYtotal, numTrials);
    Stat->Time.mean = mean_non_zero_doubles(Results.TrailTimetotal, numTrials);

    Stat->VelocityX.mean = mean_non_zero_doubles_2d(Results.Xvelocity, numTrials, dataIntervalSize);
    Stat->HelicalPitch.mean = mean_non_zero_doubles_2d(Results.HelicalPitch, numTrials, dataIntervalSize);
    Stat->VelocityY.mean = mean_non_zero_doubles_2d(Results.Yvelocity, numTrials, dataIntervalSize);
    Stat->VelocityXplus.mean = mean_non_zero_doubles_2d(Results.XvelocityPlus, numTrials, dataIntervalSize);
    Stat->VelocityXminus.mean = mean_non_zero_doubles_2d(Results.XvelocityMinus, numTrials, dataIntervalSize);
    Stat->VelocityYplus.mean = mean_non_zero_doubles_2d(Results.YvelocityPlus, numTrials, dataIntervalSize);
    Stat->VelocityYminus.mean = mean_non_zero_doubles_2d(Results.YvelocityMinus, numTrials, dataIntervalSize);
    Stat->Motors.mean = mean_non_zero_doubles_2d(Results.IterationMotorsCount, numTrials, dataSize);

    Stat->X.STD = std_non_zero_doubles(Results.TrailXtotal, numTrials, Stat->X.mean);
    Stat->Y.STD = std_non_zero_doubles(Results.TrailYtotal, numTrials, Stat->Y.mean);
    Stat->Time.STD = std_non_zero_doubles(Results.TrailTimetotal, numTrials, Stat->Time.mean);

    Stat->Motors.STD = std_non_zero_doubles_2d(Results.IterationMotorsCount, numTrials, dataSize, Stat->Motors.mean);
    Stat->HelicalPitch.STD = std_non_zero_doubles_2d(Results.HelicalPitch, numTrials, dataIntervalSize, Stat->HelicalPitch.mean);
    Stat->VelocityX.STD = std_non_zero_doubles_2d(Results.Xvelocity, numTrials, dataIntervalSize, Stat->VelocityX.mean);
    Stat->VelocityY.STD = std_non_zero_doubles_2d(Results.Yvelocity, numTrials, dataIntervalSize, Stat->VelocityY.mean);
    Stat->VelocityXplus.STD = std_non_zero_doubles_2d(Results.XvelocityPlus, numTrials, dataIntervalSize, Stat->VelocityXplus.mean);
    Stat->VelocityXminus.STD = std_non_zero_doubles_2d(Results.XvelocityMinus, numTrials, dataIntervalSize, Stat->VelocityXminus.mean);
    Stat->VelocityYplus.STD = std_non_zero_doubles_2d(Results.YvelocityPlus, numTrials, dataIntervalSize, Stat->VelocityYplus.mean);
    Stat->VelocityYminus.STD = std_non_zero_doubles_2d(Results.YvelocityMinus, numTrials, dataIntervalSize, Stat->VelocityYminus.mean);

    Stat->Motors.SEM = sem_non_zero_doubles_2d(Results.IterationMotorsCount, numTrials, dataSize, Stat->Motors.STD);
    Stat->HelicalPitch.SEM = sem_non_zero_doubles_2d(Results.HelicalPitch, numTrials, dataIntervalSize, Stat->HelicalPitch.STD);
    Stat->VelocityX.SEM = sem_non_zero_doubles_2d(Results.Xvelocity, numTrials, dataIntervalSize, Stat->VelocityX.STD);
    Stat->VelocityY.SEM = sem_non_zero_doubles_2d(Results.Yvelocity, numTrials, dataIntervalSize, Stat->VelocityY.STD);
    Stat->VelocityXplus.SEM = sem_non_zero_doubles_2d(Results.XvelocityPlus, numTrials, dataIntervalSize, Stat->VelocityXplus.STD);
    Stat->VelocityXminus.SEM = sem_non_zero_doubles_2d(Results.XvelocityMinus, numTrials, dataIntervalSize, Stat->VelocityXminus.STD);
    Stat->VelocityYplus.SEM = sem_non_zero_doubles_2d(Results.YvelocityPlus, numTrials, dataIntervalSize, Stat->VelocityYplus.STD);
    Stat->VelocityYminus.SEM = sem_non_zero_doubles_2d(Results.YvelocityMinus, numTrials, dataIntervalSize, Stat->VelocityYminus.STD);

    double sqrtTrials = sqrt((double)numTrials);
    Stat->X.SEM = Stat->X.STD / sqrtTrials;
    Stat->Y.SEM = Stat->Y.STD / sqrtTrials;
    Stat->Time.SEM = Stat->Time.STD / sqrtTrials;
}

// Compute MSD values and extract b-value from log-log linear regression
double compute_msd_and_b_value(double** Results_XIntervalMSD, double timeIntervalLength, double timeMean, double* outMSDlist) {
    double XMSDMeanlist[dataIntervalSize] = { 0 };
    int MSDcount[dataIntervalSize] = { 0 };
    int MSDmaxTime = (int)(timeMean / timeIntervalLength);
    if (MSDmaxTime >= dataIntervalSize) MSDmaxTime = dataIntervalSize - 1;

    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j <= MSDmaxTime; j++) {
            if (Results_XIntervalMSD[i][j] != 0.0) {
                XMSDMeanlist[j] += Results_XIntervalMSD[i][j];
                MSDcount[j]++;
            }
        }
    }

    for (int j = 0; j <= MSDmaxTime; j++) {
        if (MSDcount[j] > 0) {
            XMSDMeanlist[j] /= MSDcount[j];
        }
    }

    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j <= MSDmaxTime; j++) {
            if (Results_XIntervalMSD[i][j] != 0.0 && MSDcount[j] > 0) {
                double diff = Results_XIntervalMSD[i][j] - XMSDMeanlist[j];
                outMSDlist[j] += diff * diff;
            }
        }
    }

    for (int j = 0; j <= MSDmaxTime; j++) {
        if (MSDcount[j] > 0) {
            outMSDlist[j] /= MSDcount[j];
        }
    }

    return compute_b_value(MSDmaxTime, timeIntervalLength, outMSDlist);
}


// ----------------------------
//       Simulation body
// ----------------------------
void* simulate_trial(void* arg) {

    struct ThreadData* data = (struct ThreadData*)arg; // Cast input to ThreadData
    int* PositiveTotal = data->PositiveTotal;
    int* NegativeTotal = data->NegativeTotal;
    struct SimulationResults* Results = data->Results;
    double TimeIntervalsLength = data->TimeIntervalsLength;



    double walkProb[meshLimit_walk], OneMotorWalkProb[meshLimit_walk];
    int OneMotorFlag = 0, PositiveTotal_copy = 0, NegativeTotal_copy = 0;
    _Bool legalmove;
    struct DeltaTime dTime;
    struct BalanceMovment beadMovement;
    struct MotorData Motors;
    int current_trial;


    while (1) {
        int current_trial = ATOMIC_FETCH_AND_ADD(&trial, 1);
        if (current_trial >= numTrials) {
            break;
        }


        int	 whoWalks, whoAttach, whoDetach, AttachLocation, Step, IndCount = 0, NumMotorsAttach, TimeintervalCount = 0, PitchCountTiral = 0;
        double StepTime, timeinterval = 0, XstepInterval = 0, YstepInterval = 0, Pitchinterval = 0, PitchLastX = 0, PhiMeshShift;
        struct SingleTrialResults* Results_Copy = calloc(1, sizeof(struct SingleTrialResults));

#ifdef _WIN32
#include <windows.h>
        DWORD thread_id = GetCurrentThreadId();
#else
#include <unistd.h>
#include <pthread.h>
        uintptr_t thread_id = (uintptr_t)pthread_self();
#endif

        unsigned int seed = (unsigned int)(time(NULL) + (current_trial * 7919) + thread_id);

        //Spread motors 
        spreadMotors(&Motors, &seed);
        NumMotorsAttach = 1;

        printf("Start NP %d\n", current_trial);

        while (NumMotorsAttach && IndCount < dataSize)
        {


            dTime = deltaTime(Motors);// the function is auto defined as int. Thats why it shows error at the moment
            //printf("the time is %e, %e, %e\n", OneMotorWalkProb[0], dTime.Probability[1], dTime.Probability[2]);

            StepTime = 1 / (dTime.Probability[0] + dTime.Probability[1] + dTime.Probability[2]);
            // Perform one of the following actions:  walk,detach, attach or do nothing

            legalmove = false;

            while (!legalmove) {

                legalmove = true;
                Results_Copy->TrailTimetotal += StepTime; //Add to the total movment time 
                timeinterval += StepTime; //Add to the interval counting time

                // Perform one of the following actions:  walk,detach, attach or do nothing
                switch (monteCarloSelect(dTime.Probability, 3, &seed))
                {
                case 0: //walk
                    //pick random motor.
                {

                    whoWalks = selectAttachedMotor(Motors);

                    if (NumMotorsAttach == 1) {
                        if (OneMotorFlag == 0) {
                            walkAlgoritm(whoWalks, Motors, OneMotorWalkProb);
                            OneMotorFlag = 1;
                        }
                        Step = monteCarloSelect(OneMotorWalkProb, meshLimit_walk, &seed);

                    }
                    else {

                        //Get a probability list of all the possible steps for the detected motor
                        walkAlgoritm(whoWalks, Motors, walkProb);
                        // Choose Step with MonteCarlo
                        Step = monteCarloSelect(walkProb, meshLimit_walk, &seed);

                        //legalmove = excludedVolume(whoWalks, Motors.motors[whoWalks][5] - Ymesh_walk[Step] * cos(asin(Motors.motors[whoWalks][5] / radiusMT)), Motors.motors[whoWalks][4] + Xmesh_walk[Step], Motors.motors);
                        if (!(isPathClear(whoWalks, Motors.motors[whoWalks][4], Motors.motors[whoWalks][5], Motors.motors[whoWalks][4] + Xmesh_walk[Step]
                            , Motors.motors[whoWalks][5] - Ymesh_walk[Step] * cos(asin(Motors.motors[whoWalks][5] / radiusMT)), Motors.motors) && excludedVolume(whoWalks, Motors.motors[whoWalks][5] - Ymesh_walk[Step] * cos(asin(Motors.motors[whoWalks][5] / radiusMT)),
                                Motors.motors[whoWalks][4] + Xmesh_walk[Step], Motors.motors))) {

                            legalmove = false;
                        }


                    }

                    if (legalmove) {
                        // Take the Chosen step 
                        Motors.motors[whoWalks][4] += Xmesh_walk[Step];
                        Motors.motors[whoWalks][5] -= Ymesh_walk[Step] * cos(asin(Motors.motors[whoWalks][5] / radiusMT)); // Change the motor postion on the y-axis with resepect to the NP reference plane
                    }
                    break;

                }
                case 1: //detach
                    //pick random motor
                {
                    whoDetach = monteCarloSelect(dTime.detachRate, numMotors, &seed);
                    //detach
                    Motors.motors[whoDetach][0] = 0;
                    Motors.motors[whoDetach][4] = 0;
                    Motors.motors[whoDetach][5] = 0;

                    break;
                }

                case 2: //attach

                {
                    // Pick the attaching motor

                    whoAttach = monteCarloSelect(dTime.AttachProb, numMotors, &seed);
                    //Pick the attaching location 
                    AttachLocation = monteCarloSelect(dTime.attachRate[whoAttach], meshLimit_att, &seed);
                    //attach
                    Motors.motors[whoAttach][0] = 1;

                    PhiMeshShift = asin(Ymesh_att[AttachLocation] / radiusMT) - fmod(Motors.TotalDeltaY, MTBSphiShift);//The angle that the NP has shifted from the original mesh
                    Motors.motors[whoAttach][0] = 1; // Temporaly attach the selected motor
                    Motors.motors[whoAttach][4] = Xmesh_att[AttachLocation] - fmod(Motors.TotalDeltaX, MTBSxShift); // Mesh shift on X axis
                    Motors.motors[whoAttach][5] = radiusMT * sin(PhiMeshShift);// Mesh shift on Y axis

                    break;
                }
                default:
                    printf("That should never happened");
                    break;
                }

                if (legalmove) {
                    //balance the np
                    beadMovement = GradientDescent(Motors); // caclulate the movment of the NP in order to reach equilibrium 

                    Motors = beadMovement.NewLocations;

                    Pitchinterval += beadMovement.NewLocations.IterationDeltaY;
                    XstepInterval += beadMovement.NewLocations.IterationDeltaX;
                    YstepInterval += beadMovement.NewLocations.IterationDeltaY;
                }
                //Calculation for Y total,X total,X positive and X negative velocities
                if (timeinterval >= TimeIntervalsLength) {

                    if ((TimeintervalCount) < dataIntervalSize) {

                        Results_Copy->Xvelocity[TimeintervalCount] = XstepInterval / timeinterval;
                        Results_Copy->Yvelocity[TimeintervalCount] = YstepInterval / timeinterval;
                        Results_Copy->XIntervalMSD[TimeintervalCount] = beadMovement.NewLocations.TotalDeltaX;

                        if (XstepInterval > 0) {
                            Results_Copy->XvelocityPlus[TimeintervalCount] = XstepInterval / timeinterval;
                        }
                        else {
                            Results_Copy->XvelocityMinus[TimeintervalCount] = XstepInterval / timeinterval;
                        }
                        //Calculation for Y positive and negative velocity in MC intervals 
                        if (YstepInterval > 0) {
                            Results_Copy->YvelocityPlus[TimeintervalCount] += YstepInterval / timeinterval;
                        }
                        else {
                            Results_Copy->YvelocityMinus[TimeintervalCount] += YstepInterval / timeinterval;
                        }

                        timeinterval = 0;
                        XstepInterval = 0;
                        YstepInterval = 0;
                        TimeintervalCount++;
                    }
                }

            }

            NumMotorsAttach = attMotors(Motors);

            //Calculation for mean motors attached
            if (NumMotorsAttach) {
                Results_Copy->IterationMotorsCount[IndCount] = (double)NumMotorsAttach; //Create a list of the amount of motors attached in every round
            }

            //Calculation for forward or backward steps
            if (beadMovement.NewLocations.IterationDeltaX > 0) {
                (PositiveTotal_copy)++;
            }
            else {
                (NegativeTotal_copy)++;
            }
            //Calculation for Helical pitch 
            if (Pitchinterval > (2 * PI)) {
                Results_Copy->HelicalPitch[PitchCountTiral] = fabs(beadMovement.NewLocations.TotalDeltaX - PitchLastX);
                PitchLastX = beadMovement.NewLocations.TotalDeltaX;
                PitchCountTiral++;

                Pitchinterval = 0;
            }

            IndCount++;
        }
        Results_Copy->TrailXtotal = beadMovement.NewLocations.TotalDeltaX; //Create a list of the X length
        Results_Copy->TrailYtotal = beadMovement.NewLocations.TotalDeltaY; //Create a list of the Y length

        // Update the shared walkProb array with a single lock
        pthread_mutex_lock(&mutex);
        *PositiveTotal += PositiveTotal_copy;
        *NegativeTotal += NegativeTotal_copy;
        Results->TrailTimetotal[current_trial] = Results_Copy->TrailTimetotal;
        Results->TrailXtotal[current_trial] = Results_Copy->TrailXtotal;
        Results->TrailYtotal[current_trial] = Results_Copy->TrailYtotal;
        for (int i = 0; i < dataSize; i++) {
            Results->IterationMotorsCount[current_trial][i] = Results_Copy->IterationMotorsCount[i];
        }
        for (int i = 0; i < dataIntervalSize; i++) {
            Results->HelicalPitch[current_trial][i] = Results_Copy->HelicalPitch[i];
            Results->Xvelocity[current_trial][i] = Results_Copy->Xvelocity[i];
            Results->Yvelocity[current_trial][i] = Results_Copy->Yvelocity[i];
            Results->XIntervalMSD[current_trial][i] = Results_Copy->XIntervalMSD[i];
            Results->XvelocityPlus[current_trial][i] = Results_Copy->XvelocityPlus[i];
            Results->XvelocityMinus[current_trial][i] = Results_Copy->XvelocityMinus[i];
            Results->YvelocityPlus[current_trial][i] = Results_Copy->YvelocityPlus[i];
            Results->YvelocityMinus[current_trial][i] = Results_Copy->YvelocityMinus[i];
        }

        pthread_mutex_unlock(&mutex);
        free(Results_Copy);

    }

    pthread_exit(NULL);
    return NULL;
}


// ----------------------------
//       Main function
// ----------------------------

int main(void)
{
    setvbuf(stdout, NULL, _IONBF, 0);
    int PositiveTotal = 0, NegativeTotal = 0, OneMotorFlag = 0;
    int Thread_active, ActionsPerThread, start, end;
    struct SimulationResults Results;
    struct Statistics Stat = { 0 };
    double cpu_time_used, b_value, TimeIntervalsLength = 0.27, elapsed_seconds;
    double  XMSDlist[dataIntervalSize] = { 0 }, OneMotorWalkProb[meshLimit_walk];
    time_t  start_time, end_time;


    // Calculate q0.
    Na = (PolymerLength) / (aPolymer);  // Compute Na at runtime
    R0_Globe = sqrt(Na) * aPolymer;  // Compute R0_Globe at runtime
    springRelax = R0_Globe / 1.7321; // 1.7321 ~= sqrt (3)
    dF = -attEnergy - kb * Temper * (log(9 * springRelax * pow(aPolymer, 3) / (2 * PI * pow(R0_Globe, 4))) - 3 * (pow(springRelax, 2) / (2 * pow(R0_Globe, 2))));
    q0 = 1 / (1 + exp(-dF / (kb * Temper)));

#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    numThreads = sysinfo.dwNumberOfProcessors;

#elif defined(__linux__)
    long numCores = sysconf(_SC_NPROCESSORS_ONLN);
    numThreads = (int)numCores;
#endif

    start_time = time(NULL);  // Record the start time
    initializeThreads();
    allocateSimulationResults(&Results);
    printf("The simulation of %d motors as started \n", numMotors);
    // Parallel NP Running calculation using threads

    //assignThreadRanges(numTrials, numThreads, &Thread_active, &ActionsPerThread);

    for (int t = 0; t < numThreads; t++) {

        threadData[t] = (struct ThreadData){
            // Pass by value
       .PositiveTotal = &PositiveTotal,
       .NegativeTotal = &NegativeTotal,
       .Results = &Results,
       .TimeIntervalsLength = TimeIntervalsLength
        };

        pthread_create(&threads[t], NULL, simulate_trial, &threadData[t]);

    }

    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }


    // Claculate Statistics like mean,STD and SEM of all the system variables 
    calculate_statistics(Results, &Stat);

    // Claculate the MSD
    b_value = compute_msd_and_b_value(Results.XIntervalMSD, TimeIntervalsLength, Stat.Time.mean, XMSDlist);

    //save data
    FILE* file = fopen("output-NP.txt", "w");  // Open the file in write mode
    if (file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    // End time
    end_time = time(NULL);  // Record the end time

    // Calculate the elapsed time
    elapsed_seconds = difftime(end_time, start_time);
    // Print the runtime
    printf("Time taken: %f seconds\n", elapsed_seconds);



    fprintf(file, "number of particles = %d\n", numTrials);
    fprintf(file, "This is a Simulation for NP with %d motors\n\n", numMotors);
    fprintf(file, "%-25s %-20s %-20s %-20s\n", "Parameter", "Mean", "+/- SEM", "STD");  // Header
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Motors connected", Stat.Motors.mean, Stat.Motors.SEM, Stat.Motors.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Time [s]", Stat.Time.mean, Stat.Time.SEM, Stat.Time.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "X [nm]", Stat.X.mean * 1e9, Stat.X.SEM * 1e9, Stat.X.STD * 1e9);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Y [Rad]", Stat.Y.mean, Stat.Y.SEM, Stat.Y.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "X Velocity [nm/s]", Stat.VelocityX.mean * 1e9, Stat.VelocityX.SEM * 1e9, Stat.VelocityX.STD * 1e9);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "X Velocity > 0 [nm/s]", Stat.VelocityXplus.mean * 1e9, Stat.VelocityXplus.SEM * 1e9, Stat.VelocityXplus.STD * 1e9);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "X Velocity < 0 [nm/s]", Stat.VelocityXminus.mean * 1e9, Stat.VelocityXminus.SEM * 1e9, Stat.VelocityXminus.STD * 1e9);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Y Velocity [Rad/s]", Stat.VelocityY.mean, Stat.VelocityY.SEM, Stat.VelocityY.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Y Velocity > 0 [Rad/s]", Stat.VelocityYplus.mean, Stat.VelocityYplus.SEM, Stat.VelocityYplus.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Y Velocity < 0 [Rad/s]", Stat.VelocityYminus.mean, Stat.VelocityYminus.SEM, Stat.VelocityYminus.STD);
    fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "<H> [nm]", Stat.HelicalPitch.mean, Stat.HelicalPitch.SEM, Stat.HelicalPitch.STD);
    fprintf(file, "%-25s %-20.2e \n", "Mean Helical pitch [nm]", 2 * PI * Stat.VelocityX.mean / Stat.VelocityY.mean);

    fprintf(file, "%-25s %15.2e | %15.2e\n", "Fractions of forward & backward steps",
        (double)PositiveTotal / (double)(PositiveTotal + NegativeTotal),
        (double)NegativeTotal / (double)(PositiveTotal + NegativeTotal));

    fprintf(file, "This is the MSD LIST \n");
    for (int i = 0; i < dataIntervalSize; i++) {
        if (XMSDlist[i]) {
            fprintf(file, "%e,", XMSDlist[i]);
        }
    }
    fprintf(file, " \n This is the b value %e \n", b_value);


    fclose(file);  // Close the file

    freeSimulationResults(&Results);

    cleanupThreads();
    return 0;
}
