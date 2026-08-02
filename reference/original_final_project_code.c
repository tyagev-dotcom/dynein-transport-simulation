# define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "mtwister.h"
#ifdef _WIN32
#define OS_WINDOWS
#include <windows.h>
#else
#define OS_LINUX
#include <unistd.h> // For sysconf
#endif


// define global parameters 
#define WINDING false
#define numMotors 13 // # of grafted motor
#define numTrials 1000 // # of NP runs
#define radiusMT  (12.5e-9) //microtubule radius 
#define radiusNC  20e-9 
#define effTemp 1900
#define meshLimit_att 49 // size of attach array (e.g., Xmesh_att)
#define dynein_height 86e-9 
#define dyenin_diameter 10e-9
#define dyenin_diameter_path 12e-9
#define meshLimit_walk 102 // size of walk array (e.g., Xmesh_walk)
#define dataSize 80000 // The maximum number of iteration (steps/attach/detach)
#define dataIntervalSize 500 // The maximum number of interval counting


// Constants for Gradient Descent
#define Angle_LEARNING_RATE 1e19 // Reduced learning rate for Gradient Descent
#define Step_Learning_Rate 1e3
#define MAX_GRADIENT_ITER 100  // Maximum iterations for Gradient Descent
#define GRADIENT_THRESHOLD 1e-27 // Stopping criterion for Gradient Descent

// Define arrays

double PolymerLength = 4.256e-8; //4.256e-8
double aPolymer = 0.75e-9;
double Na;  // Compute Na at runtime
double R0_Globe;  // Compute R0_Globe at runtime
double PI = 3.14159265358979323846;
double attEnergy = 3.2899e-20;  //energy gain from motor attachment to microtubule
double Temper = 298; // real/physical temperature 
double kb = 1.38e-23; //boltzmann constant
double fs = 3.6e-12; //stalling force
double springRelax = 5.65e-9 / 1.7321; // 1.7321 ~= sqrt (3)
double tauStep = 0.01088; // The average time for a step: <L_x> / <V_x> = ~8.8 / 800 [sec]  
double MTBSxShift = 8e-9; // Distand between MTBS on the x axis 
double MTBSphiShift = (2 * 3.14159265358979323846 / 13); // There are 13 MTBS around the MT
double dF;
double q0;
// Optinal: change q0 to deal with the different rate of detachment
//q0 = q0 * some factor (for yeast dynein, the factor = 19 since the characterstic time for detachment is longer...)
int numThreads;                  // Number of threads based on available cores
struct ThreadData* threadData;   // Dynamically allocated global ThreadData
pthread_mutex_t mutex;           // Mutex for thread synchronization
pthread_t* threads;            // Pointer to array of thread handles


double stepRate;

double Xmesh_att[] = { -2.213867e-08 ,-2.245691e-08 ,-2.312864e-08 ,-2.400000e-08 ,-2.487136e-08 ,-2.554309e-08 ,-2.586133e-08 ,-1.413867e-08 ,-1.445691e-08 ,-1.512864e-08 ,-1.600000e-08 ,-1.687136e-08 ,-1.754309e-08 ,-1.786133e-08 ,-6.138671e-09 ,-6.456905e-09 ,-7.128644e-09 ,-8.000000e-09 ,-8.871356e-09 ,-9.543095e-09 ,-9.861329e-09 ,1.861329e-09 ,1.543095e-09 ,8.713559e-10 ,0 ,-8.713559e-10 ,-1.543095e-09 ,-1.861329e-09 ,9.861329e-09 ,9.543095e-09 ,8.871356e-09 ,8.000000e-09 ,7.128644e-09 ,6.456905e-09 ,6.138671e-09 ,1.786133e-08 ,1.754309e-08 ,1.687136e-08 ,1.600000e-08 ,1.512864e-08 ,1.445691e-08 ,1.413867e-08 ,2.586133e-08 ,2.554309e-08 ,2.487136e-08 ,2.400000e-08 ,2.312864e-08 ,2.245691e-08 ,2.213867e-08 };
double Ymesh_att[] = { -1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 ,-1.240886e-08 ,-1.028730e-08 ,-5.809040e-09 ,0 ,5.809040e-09 ,1.028730e-08 ,1.240886e-08 };

// The stepping is specific for cytoplasmic mammalian dynein!!!
double Xmesh_walk[] = { -6.077684e-08,  -5.277684e-08,  -4.893779e-08,  -4.000000e-08,  -3.693779e-08,  -3.122316e-08,  -2.843905e-08,  -2.800000e-08,  -2.706221e-08,  -2.477684e-08,  -2.306221e-08,  -2.093779e-08,  -2.077684e-08,  -2.043905e-08,  -2.000000e-08,  -1.956095e-08,  -1.922316e-08,  -1.906221e-08,  -1.693779e-08,  -1.677684e-08,  -1.600000e-08,  -1.556095e-08,  -1.506221e-08,  -1.293779e-08,  -1.277684e-08,  -1.243905e-08,  -1.200000e-08,  -1.156095e-08,  -1.122316e-08,  -1.106221e-08,  -8.937794e-09,  -8.776835e-09,  -8.439048e-09,  -8.000000e-09,  -7.560952e-09,  -7.223165e-09,  -7.062206e-09,  -4.937794e-09,  -4.776835e-09,  -4.439048e-09,  -4.000000e-09,  -3.223165e-09,  -3.062206e-09,  -9.377940e-10,  9.377940e-10,  3.062206e-09,  3.223165e-09,  3.560952e-09,  4.000000e-09,  4.439048e-09,  4.776835e-09,  4.937794e-09,  7.062206e-09,  7.223165e-09,  7.560952e-09,  8.000000e-09,  8.439048e-09,  8.776835e-09,  8.937794e-09,  1.106221e-08,  1.122316e-08,  1.156095e-08,  1.200000e-08,  1.243905e-08,  1.277684e-08,  1.293779e-08,  1.506221e-08,  1.600000e-08,  1.677684e-08,  1.693779e-08,  1.906221e-08,  2.000000e-08,  2.043905e-08,  2.306221e-08,  2.356095e-08,  2.477684e-08,  2.706221e-08,  2.800000e-08,  3.293779e-08,  3.506221e-08,  4.000000e-08,  4.306221e-08,  4.400000e-08,  5.906221e-08,  6.077684e-08,  5.277684e-08,  4.893779e-08,  3.693779e-08,  3.122316e-08,  2.843905e-08,  2.093779e-08,  2.077684e-08,  1.956095e-08,  1.922316e-08,  1.556095e-08,  -3.560952e-09,  -2.356095e-08,  -3.293779e-08,  -3.506221e-08,  -4.306221e-08,  -4.400000e-08,  -5.906221e-08 };
double Ymesh_walk[] = { -5.140000e-09,  -5.140000e-09,  -6.205000e-09,  0,  -6.205000e-09,  5.140000e-09,  -2.905000e-09,  0,  6.205000e-09,  -5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  0,  2.905000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  -5.140000e-09,  -2.905000e-09,  0,  2.905000e-09,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  0,  5.140000e-09,  6.205000e-09,  -6.205000e-09,  0,  2.905000e-09,  -6.205000e-09,  -2.905000e-09,  5.140000e-09,  -6.205000e-09,  0,  6.205000e-09,  -6.205000e-09,  0,  -6.205000e-09,  0,  -6.205000e-09,  5.140000e-09,  5.140000e-09,  6.205000e-09,  6.205000e-09,  -5.140000e-09,  2.905000e-09,  6.205000e-09,  5.140000e-09,  -2.905000e-09,  -5.140000e-09,  -2.905000e-09,  2.905000e-09,  2.905000e-09,  -6.205000e-09,  6.205000e-09,  6.205000e-09,  0,  6.205000e-09 };
double xy_WalkProb[] = { 2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  7.832898e-03,  7.832898e-03,  5.221932e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  5.221932e-03,  5.221932e-03,  2.610966e-03,  5.221932e-03,  7.832898e-03,  2.088773e-02,  7.832898e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  2.610966e-02,  3.655352e-02,  1.305483e-02,  1.566580e-02,  1.044386e-02,  1.566580e-02,  7.832898e-03,  4.177546e-02,  7.571802e-02,  5.221932e-03,  1.044386e-02,  1.827676e-02,  7.832898e-03,  1.566580e-02,  2.610966e-03,  2.610966e-03,  1.044386e-02,  2.610966e-03,  1.305483e-02,  7.832898e-03,  1.305483e-02,  7.832898e-03,  8.093995e-02,  6.527415e-02,  1.044386e-02,  2.088773e-02,  1.305483e-02,  2.610966e-02,  7.832898e-03,  4.177546e-02,  2.872063e-02,  7.832898e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  5.221932e-03,  1.044386e-02,  2.349869e-02,  2.610966e-03,  2.610966e-03,  2.610966e-03,  7.832898e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  5.221932e-03,  2.610966e-03,  7.832898e-03,  7.832898e-03,  1.044386e-02,  2.610966e-03,  5.221932e-03,  1.305483e-02,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03,  2.610966e-03 };


struct MotorData {
	double motors[numMotors][6]; // Every motors have 5 values [attach/detach (1/0),phi angle,teta angle,x on MT,y on MT]
	double distance2MT; // distance of the center of the NP from the surface of MT
	double IterationDeltaX;
	double IterationDeltaY;
	double TotalDeltaX;
	double TotalDeltaY;

};

struct BalanceMovment {
	struct MotorData NewLocations;
	double Energy;

};

struct DeltaTime {
	double Probability[3]; // List of probabilities : [stepRate,SumDetach,SumAttach]
	double attachRate[numMotors][meshLimit_att];
	double detachRate[numMotors];
	double AttachProb[numMotors];
};

struct SimulationResults {
	double** IterationMotorsCount;
	double* TrailXtotal;
	double* TrailYtotal;
	double* TrailTimetotal;
	double** Yvelocity;
	double** Xvelocity;
	double** XvelocityAbs;
	double** XvelocityPlus;
	double** XvelocityMinus;
	double** XIntervalMSD;
	double* YvelocityPlus;
	double* YvelocityMinus;
	double** HelicalPitch;
	int* StepsNum;

};

struct VariableStatistics {
	double mean;
	double STD;
	double SEM;
};

struct Statistics {
	struct VariableStatistics X;
	struct VariableStatistics Y;
	struct VariableStatistics Time;
	struct VariableStatistics Motors;
	struct VariableStatistics VelocityX;
	struct VariableStatistics VelocityXplus;
	struct VariableStatistics VelocityXminus;
	struct VariableStatistics VelocityXabs;
	struct VariableStatistics VelocityY;
	struct VariableStatistics VelocityYplus;
	struct VariableStatistics VelocityYminus;
	struct VariableStatistics Steps;
	struct VariableStatistics angle;
	struct VariableStatistics HelicalPitch;


};

// Define the structure for a linked list node
typedef struct MotorPair {
	int motor1; // Index of the first motor
	int motor2; // Index of the second motor
	double angle_sum; // Sum of angles
	double last_angle; // Last angle calculated
	struct MotorPair* next; // Pointer to the next node
} MotorPair;

// Shared data structures
struct ThreadData {
	struct MotorData Motors;
	struct DeltaTime* dTime;
	double i_Energy;
	double q0;
	int start;
	int end;
	double* walkProb;
	int whoWalks;
};



double compute_b_value(int N, double TimeIntervalsLength, double msd_values[]);
void allocateSimulationResults(struct SimulationResults* Results);
void freeSimulationResults(struct SimulationResults* Results);
void spreadMotors(struct MotorData* Motors, MTRand* rand);
void rotateMotor(struct MotorData* Motors, double xRot, double yRot, double zRot);
int attMotors(struct MotorData Motors);
double energy(struct MotorData Motors);
_Bool excludedVolume(int ChoosenMotor, double yNextLocation, double xNextLocation, double motors[numMotors][6]);
_Bool isPathClear(int steppingMotor, double xCurr, double yCurr, double xNext, double yNext, double motors[numMotors][6]);
struct DeltaTime deltaTime(struct MotorData Motors);
int monteCarloSelect(double* values, int size, MTRand* rand);
int selectAttachedMotor(struct MotorData);
void copyMotorDataStruct(struct MotorData* dest, struct MotorData* src);
double CalculateEnergyGradient(struct MotorData Motors, int axis);
struct BalanceMovment GradientDescent(struct MotorData Motors);
void isLimiteisExceeded(double* ymtRot, double* DistShift, struct MotorData Motors, double* xRot, double* yRot);
void SingleMotorAlign(struct MotorData* Motors, int whoattach);
double normalizeAngle(double angle);
double calculateAngle(double Motor1[6], double Motor2[6]);
MotorPair* createNode(int motor1, int motor2);
void addMotorPair(MotorPair** head, int motor1, int motor2);
void removeMotorPair(MotorPair** head, int motor1, int motor2, double* angle_mean, int* coil_count, int* angle_count);
void updateMotorAngles(MotorPair* head, double Motors[][6]);
void initializeThreads();
void cleanupThreads();
void assignThreadRanges(int numItems, int numThreads, int* Thread_active, int* ActionsPerThread);
void main(void);

MotorPair* motorPairs = NULL;


void initializeThreads() {
	

	threadData = (struct ThreadData*)malloc(numThreads * sizeof(struct ThreadData));
	threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));

	if (!threadData || !threads) {
		fprintf(stderr, "Error: Unable to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	pthread_mutex_init(&mutex, NULL);
}


void cleanupThreads() {
	free(threadData);
	free(threads);
	pthread_mutex_destroy(&mutex);
}

// Function to allocate memory for the results
void allocateSimulationResults(struct SimulationResults* Results) {
	Results->TrailYtotal = calloc(numTrials, sizeof(double));
	Results->TrailXtotal = calloc(numTrials, sizeof(double));
	Results->TrailTimetotal = calloc(numTrials, sizeof(double));
	Results->StepsNum = calloc(numTrials, sizeof(int));
	Results->YvelocityPlus = calloc(numTrials, sizeof(double));
	Results->YvelocityMinus = calloc(numTrials, sizeof(double));

	// Allocate memory for 2D arrays
	Results->Yvelocity = calloc(numTrials, sizeof(double*));
	Results->Xvelocity = calloc(numTrials, sizeof(double*));
	Results->XvelocityAbs = calloc(numTrials, sizeof(double*));
	Results->XvelocityPlus = calloc(numTrials, sizeof(double*));
	Results->XvelocityMinus = calloc(numTrials, sizeof(double*));
	Results->IterationMotorsCount = calloc(numTrials, sizeof(double*));
	Results->HelicalPitch = calloc(numTrials, sizeof(double*));
	Results->XIntervalMSD = calloc(numTrials, sizeof(double*));


	for (int i = 0; i < numTrials; i++) {
		Results->Yvelocity[i] = calloc(dataIntervalSize, sizeof(double));
		Results->Xvelocity[i] = calloc(dataIntervalSize, sizeof(double));
		Results->XvelocityAbs[i] = calloc(dataIntervalSize, sizeof(double));
		Results->XvelocityPlus[i] = calloc(dataIntervalSize, sizeof(double));
		Results->XvelocityMinus[i] = calloc(dataIntervalSize, sizeof(double));
		Results->HelicalPitch[i] = calloc(dataIntervalSize, sizeof(double));
		Results->IterationMotorsCount[i] = calloc(dataSize, sizeof(double));
		Results->XIntervalMSD[i] = calloc(dataIntervalSize, sizeof(double));

	}
}

void freeSimulationResults(struct SimulationResults* Results) {
	for (int i = 0; i < numTrials; i++) {
		free(Results->Yvelocity[i]);
		free(Results->Xvelocity[i]);
		free(Results->XvelocityAbs[i]);
		free(Results->XvelocityPlus[i]);
		free(Results->XvelocityMinus[i]);
		free(Results->IterationMotorsCount[i]);
		free(Results->HelicalPitch[i]);
		free(Results->XIntervalMSD[i]);

	}

	free(Results->Yvelocity);
	free(Results->Xvelocity);
	free(Results->XvelocityAbs);
	free(Results->XvelocityPlus);
	free(Results->XvelocityMinus);
	free(Results->YvelocityPlus);
	free(Results->YvelocityMinus);
	free(Results->StepsNum);
	free(Results->HelicalPitch);
	free(Results->XIntervalMSD);

	free(Results->IterationMotorsCount);
	free(Results->TrailYtotal);
	free(Results->TrailXtotal);
	free(Results->TrailTimetotal);
}


/* This function returns an array of motors that spread randomly on the surface of the NP */
void spreadMotors(struct MotorData* Motors, MTRand* rand)
{
	double phi, teta;

	for (int i = 0; i < numMotors; i++)
	{

		phi = genRand(rand) * 2 * PI;
		teta = genRand(rand) * PI;
		Motors->motors[i][0] = 0;
		Motors->motors[i][1] = radiusNC * sin(teta) * cos(phi);
		Motors->motors[i][2] = radiusNC * sin(teta) * sin(phi);
		Motors->motors[i][3] = radiusNC * cos(teta);
		Motors->motors[i][4] = 0;
		Motors->motors[i][5] = 0;

	}
	//attach first motor
	Motors->motors[0][0] = 1;
	// Rotate the NP so the first motor will be attached at coordinate (0,0,-R(NP))
	//reset Distance 
	Motors->distance2MT = springRelax + radiusNC;
	SingleMotorAlign(Motors, 0);
	Motors->IterationDeltaX = 0;
	Motors->IterationDeltaY = 0;
	Motors->TotalDeltaX = 0;
	Motors->TotalDeltaY = 0;

	return;
}



// Count the number of attached motors in the real array "motors"
int attMotors(struct MotorData Motors)
{
	int numAttachedMotors = 0;
	for (int i = 0; i < numMotors; i++)
	{
		if (Motors.motors[i][0]) {

			numAttachedMotors++;
		}
	}
	return numAttachedMotors;
}

// Function to calculate the shortest distance and midpoint
double calculateDistanceAndMidpoint(double motor1[6], double motor2[6]) {
	double distance, midpoint_x, midpoint_y, midpoint_z;
	// Dot product of the two points
	double dotProduct = motor1[1] * motor2[1] + motor1[2] * motor2[2] + motor1[3] * motor2[3];

	// Magnitude of the radius squared
	double rSquared = radiusNC * radiusNC;

	// Central angle (in radians) using arccos of the dot product
	double theta = acos(dotProduct / rSquared);

	// Shortest distance along the surface of the sphere
	distance = radiusNC * theta;

	// Calculate the midpoint vector (vector sum)
	midpoint_x = (motor1[1] + motor2[1]) / 2;
	midpoint_y = (motor1[2] + motor2[2]) / 2;
	midpoint_z = (motor1[3] + motor2[3]) / 2;

	// Normalize the midpoint to lie on the sphere's surface
	double magnitude = sqrt(midpoint_x * midpoint_x +
		midpoint_y * midpoint_y +
		midpoint_z * midpoint_z);
	motor1[1] = (midpoint_x / magnitude) * radiusNC;
	motor2[1] = (midpoint_x / magnitude) * radiusNC;
	motor1[2] = (midpoint_y / magnitude) * radiusNC;
	motor2[2] = (midpoint_y / magnitude) * radiusNC;
	motor1[3] = (midpoint_z / magnitude) * radiusNC;
	motor2[3] = (midpoint_z / magnitude) * radiusNC;

	return distance;
}



double energy(struct MotorData Motors)
{
	double result = 0, zMT, M, zPrime, R0;
	double x, y, z, Xmt, Ymt;
	double winding, angle, distance, New_Na, R0_list[numMotors] = { 0 };
	MotorPair* temp = motorPairs;
	if (WINDING) {
		while (temp != NULL) {
			if (Motors.motors[temp->motor1] && Motors.motors[temp->motor2]) {
				angle = normalizeAngle(calculateAngle(Motors.motors[temp->motor1], Motors.motors[temp->motor2]));
				winding = temp->angle_sum + normalizeAngle(angle - temp->last_angle);
				if (fabs(winding) >= 2 * PI) {

					distance = calculateDistanceAndMidpoint(Motors.motors[temp->motor1], Motors.motors[temp->motor2]);
					New_Na = Na - 1.5 - distance / (2 * aPolymer);
					if (New_Na > 0) {
						R0_list[temp->motor1] = R0_Globe * sqrt(New_Na / Na);
						R0_list[temp->motor2] = R0_Globe * sqrt(New_Na / Na);
					}
					else {
						R0_list[temp->motor1] = 1e-12;
						R0_list[temp->motor2] = 1e-12;
					}

				}
			}

			temp = temp->next;
		}
	}
	for (int i = 0; i < numMotors; i++)
	{

		if (Motors.motors[i][0])
		{   // Convert the angle phi and teta angles into cartesian coordinate
			x = Motors.motors[i][1];
			y = Motors.motors[i][2];
			z = Motors.motors[i][3];
			Xmt = Motors.motors[i][4];
			Ymt = Motors.motors[i][5];
			R0 = R0 = R0_list[i] ? R0_list[i] : R0_Globe; //R0 = R0_list[i] ? R0_list[i] : R0_Globe;


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

void copyMotorDataStruct(struct MotorData* dest, struct MotorData* src) {
	for (int i = 0; i < numMotors; i++) {
		for (int j = 0; j < 6; j++) {
			dest->motors[i][j] = src->motors[i][j];
		}
	}
}

// Rotate a single motor by given angles along X, Y, and Z axes 
void rotateMotor(struct MotorData* Motors, double xRot, double yRot, double zRot) {
	double x1, y1, z1, x2, y2, z2;
	for (int i = 0; i < numMotors; i++) {
		// Rotate around X-axis
		x1 = Motors->motors[i][1];
		y1 = Motors->motors[i][2] * cos(xRot) - Motors->motors[i][3] * sin(xRot);
		z1 = Motors->motors[i][2] * sin(xRot) + Motors->motors[i][3] * cos(xRot);

		// Rotate around Y-axis
		x2 = x1 * cos(yRot) + z1 * sin(yRot);
		y2 = y1;
		z2 = -x1 * sin(yRot) + z1 * cos(yRot);

		// Rotate around Z-axis
		Motors->motors[i][1] = x2 * cos(zRot) - y2 * sin(zRot);
		Motors->motors[i][2] = x2 * sin(zRot) + y2 * cos(zRot);
		Motors->motors[i][3] = z2;
	}
}

void isLimiteisExceeded(double* ymtRot, double* DistShift, struct MotorData Motors, double* xRot, double* yRot) {
	rotateMotor(&Motors, *xRot, *yRot, 0);
	for (int i = 0; i < numMotors; i++) {
		if (Motors.motors[i][0]) {
			if (abs(fmod((*ymtRot + asin(Motors.motors[i][5] / radiusMT)), 2 * PI)) > PI / 2) {
				*ymtRot = 0;
			}
		}
	}
	if ((*DistShift + Motors.distance2MT) > (springRelax + radiusNC) || (*DistShift + Motors.distance2MT) < radiusNC) {
		*DistShift = 0;
	}

}
void SingleMotorAlign(struct MotorData* Motors, int whoattach) {
	double SinglemotorYrotate, SinglemotorXrotate;

	SinglemotorYrotate = atan(-Motors->motors[whoattach][1] / Motors->motors[whoattach][3]);
	SinglemotorXrotate = atan(Motors->motors[whoattach][2] / Motors->motors[whoattach][3]);
	rotateMotor(Motors, SinglemotorXrotate, SinglemotorYrotate, 0);

	if (Motors->motors[whoattach][3] > 0) {
		// Rotate 180° about the X-axis to flip the sphere.
		rotateMotor(Motors, PI, 0, 0);
	}
}

struct BalanceMovment GradientDescent(struct MotorData Motors) {
	double xRot = 0, yRot = 0, zRot = 0, xmtShift = 0, ymtRot = 0, disShift = 0;  // Rotation angles
	double gradx = 0, grady = 0, gradz = 0, gradxmt = 0, gradymt = 0, graddis = 0;      // Gradients for each rotation axis
	double prevEnergy, currentEnergy, momentum_x = 0, momentum_y = 0, momentum_z = 0, momentum_xmt = 0, momentum_ymt = 0, momentum_dis = 0; // Energy and momentum
	double beta = 0.9;  // Momentum decay factor
	double DeltaX = 0, DeltaY = 0;
	struct BalanceMovment BeadBalance;

	if (attMotors(Motors) == 0) {
		DeltaX = 0;
		DeltaY = 0;
		currentEnergy = dF;

	}
	else if (attMotors(Motors) == 1) {
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
		currentEnergy = energy(Motors);  // Initial energy

		for (int iter = 0; iter < MAX_GRADIENT_ITER; iter++) {
			// Calculate gradients
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


			// Update rotation angles
			xRot = -Angle_LEARNING_RATE * momentum_x;
			yRot = -Angle_LEARNING_RATE * momentum_y;
			zRot = -Angle_LEARNING_RATE * momentum_z;
			xmtShift = -Step_Learning_Rate * momentum_xmt;
			ymtRot = -Angle_LEARNING_RATE * 1e-2 * momentum_ymt;
			disShift = -Step_Learning_Rate * 1e-2 * momentum_dis;


			isLimiteisExceeded(&ymtRot, &disShift, Motors, &xRot, &yRot); // Clamp angles to valid ranges to avoid instability


			rotateMotor(&Motors, xRot, yRot, zRot);

			for (int i = 0; i < numMotors; i++) {
				Motors.motors[i][4] += xmtShift;
				Motors.motors[i][5] = radiusMT * sin(ymtRot + asin(Motors.motors[i][5] / radiusMT));
			}
			Motors.distance2MT += disShift;

			// Recalculate energy
			prevEnergy = currentEnergy;
			currentEnergy = energy(Motors);

			DeltaX -= xmtShift;
			DeltaY -= ymtRot;

			// Check convergence
			if (fabs(prevEnergy - currentEnergy) < GRADIENT_THRESHOLD) {
				break;
			}

		}

	}
	Motors.IterationDeltaX = DeltaX;
	Motors.IterationDeltaY = DeltaY; // Add total delta for attach calculation
	Motors.TotalDeltaX += DeltaX;
	Motors.TotalDeltaY += DeltaY;
	BeadBalance.NewLocations = Motors;
	BeadBalance.Energy = currentEnergy;

	return BeadBalance;
}





// Function to calculate the energy gradient
double CalculateEnergyGradient(struct MotorData Motors, int axis) {
	double energyPlus, energyMinus = energy(Motors), h = 1e-6;

	// Store original state
	struct MotorData Motorscopy;

	copyMotorDataStruct(&Motorscopy, &Motors);
	Motorscopy.distance2MT = Motors.distance2MT;

	if (axis == 1) {
		rotateMotor(&Motorscopy, h, 0, 0);
	}
	else if (axis == 2) {
		rotateMotor(&Motorscopy, 0, h, 0);

	}
	else if (axis == 3) {
		rotateMotor(&Motorscopy, 0, 0, h);

	}
	else if (axis == 4) {
		h = 1e-12;
		for (int i = 0; i < numMotors; i++) {
			Motorscopy.motors[i][4] += h;
		}
	}
	else if (axis == 5) {
		for (int i = 0; i < numMotors; i++) {
			Motorscopy.motors[i][5] = radiusMT * sin(h + asin(Motorscopy.motors[i][5] / radiusMT));

		}
	}
	else if (axis == 6) {
		h = 1e-12;
		Motorscopy.distance2MT += h;
	}


	energyPlus = energy(Motorscopy);

	// Calculate gradient
	return (energyPlus - energyMinus) / h;

}




/* This function checks excluded volume */
_Bool excludedVolume(int ChoosenMotor, double yNextLocation, double xNextLocation, double motors[numMotors][6])
{
	_Bool legalMove = true;
	{
		// notice! arrays start at 0. Thus, this function uses whoWalks minus 1 .
		for (int i = 0; i < numMotors; i++)
		{
			// If distance between motors is less the 
			if (i != ChoosenMotor) {
				if (motors[i][0] && (sqrt
				(
					/*X*/ pow((xNextLocation - motors[i][4]), 2) +
					/*Y*/ pow((yNextLocation - motors[i][5]), 2) +
					/*Z*/ pow(((sqrt(pow((radiusMT), 2) - pow((yNextLocation), 2))) - (sqrt(pow((radiusMT), 2) - pow((motors[i][5]), 2)))), 2)
				)

					<= dyenin_diameter
					)) {
					legalMove = false;
				}
			}
		}
	}
	return legalMove;
}

/* This function checks if any motor is in the stepping path */
_Bool isPathClear(int steppingMotor, double xCurr, double yCurr,
	double xNext, double yNext, double motors[numMotors][6])
{
	_Bool clearPath = true;
	int numSteps = 100;

	for (int i = 0; i < numMotors; i++) {
		if (i != steppingMotor && motors[i][0]) {
			// Motor position
			double xM = motors[i][4];
			double yM = motors[i][5];
			double zM = sqrt(pow(radiusMT, 2) - pow(yM, 2));

			for (int step = 0; step <= numSteps; step++) {
				double t = (double)step / numSteps;
				// Interpolated position along the curved path
				double xStep = xCurr + t * (xNext - xCurr);
				double yStep = yCurr + t * (yNext - yCurr);
				double zStep = sqrt(pow(radiusMT, 2) - pow(yStep, 2));

				// Distance from current step position to the motor position
				double distance = sqrt(pow(xM - xStep, 2) + pow(yM - yStep, 2) + pow(zM - zStep, 2));

				if (distance <= dyenin_diameter_path) {
					clearPath = false;
					break;
				}
			}
			if (!clearPath) {
				break;
			}

		}
	}

	return clearPath;
}



void* calculateDetachRates(void* arg) {
	struct ThreadData* data = (struct ThreadData*)arg;
	struct MotorData Motors = data->Motors; // Pointer to avoid copying
	struct DeltaTime* dTime = data->dTime;   // Shared structure
	double i_Energy = data->i_Energy;


	double localProbability = 0; // Thread-local aggregation
	double localDetachRate[numMotors] = { 0 }; // Thread-local detach rates

	for (int i = data->start; i < data->end; i++) {
		if (Motors.motors[i][0]) { // Motor is attached
			Motors.motors[i][0] = 0; // Temporarily detach the motor
			struct BalanceMovment BeadBalance = GradientDescent(Motors);
			double detachRate = 1 / (q0 * (1 + exp((BeadBalance.Energy - i_Energy) / (kb * Temper))));

			localDetachRate[i] = detachRate;
			localProbability += detachRate;

			Motors.motors[i][0] = 1; // Restore motor state
		}
	}

	// Lock once to update shared data
	pthread_mutex_lock(&mutex);
	for (int i = data->start; i < data->end; i++) {
		dTime->detachRate[i] += localDetachRate[i];
	}
	dTime->Probability[1] += localProbability;
	pthread_mutex_unlock(&mutex);

	pthread_exit(NULL);
}

void* calculateAttachRates(void* arg) {
	struct ThreadData* data = (struct ThreadData*)arg;
	struct MotorData Motors = data->Motors; // Pointer to avoid copying
	struct DeltaTime* dTime = data->dTime;   // Shared structure
	double i_Energy = data->i_Energy;
	double localProbability = 0; // Thread-local aggregation
	double localAttachRate[numMotors][meshLimit_att] = { {0} }; // Thread-local attach rates
	double localAttachProb[numMotors] = { 0 }; // Thread-local attach probabilities

	for (int i = 0; i < numMotors; i++) {
		for (int k = data->start; k < data->end; k++) {
			if (Motors.motors[i][0] == 0 && Motors.motors[i][3] < 0) { // (Motors.motors[i][0] == 0 && Motors.motors[i][3] < 0) 
				double PhiMeshShift = asin(Ymesh_att[k] / radiusMT) - fmod(Motors.TotalDeltaY, MTBSphiShift);
				Motors.motors[i][0] = 1; // Temporarily attach the motor
				Motors.motors[i][4] = Xmesh_att[k] - fmod(Motors.TotalDeltaX, MTBSxShift);
				Motors.motors[i][5] = radiusMT * sin(PhiMeshShift);

				if (excludedVolume(i, Motors.motors[i][5], Motors.motors[i][4], Motors.motors) && abs(PhiMeshShift) <= PI / 2) {
					struct BalanceMovment BeadBalance = GradientDescent(Motors);
					double attachRate = 1 / (q0 * (1 + exp((BeadBalance.Energy - i_Energy) / (kb * Temper))));



					localAttachRate[i][k] = attachRate;
					localAttachProb[i] += attachRate;
					localProbability += attachRate;
				}

				Motors.motors[i][0] = 0; // Restore motor state
				Motors.motors[i][4] = 0;
				Motors.motors[i][5] = 0;
			}
		}
	}

	// Lock once to update shared data
	pthread_mutex_lock(&mutex);
	for (int i = 0; i < numMotors; i++) {
		dTime->AttachProb[i] += localAttachProb[i];
		for (int k = data->start; k < data->end; k++) {
			dTime->attachRate[i][k] += localAttachRate[i][k];
		}
	}
	dTime->Probability[2] += localProbability;
	pthread_mutex_unlock(&mutex);

	pthread_exit(NULL);
}
// Function to assign ranges for each thread and initialize threadData
void assignThreadRanges(int numItems, int numThreads, int* Thread_active, int* ActionsPerThread) {


	if (numItems < numThreads) {
		*Thread_active = numItems;
		*ActionsPerThread = 1;
	}
	else {
		*Thread_active = numThreads;
		*ActionsPerThread = numItems / numThreads;
	}
}

struct DeltaTime deltaTime(struct MotorData Motors) {

	struct DeltaTime dTime = { 0 };
	int NumMotorsAttach = attMotors(Motors); // Number of attached motors
	double i_Energy = (NumMotorsAttach == 1) ? dF : energy(Motors); // Calculate energy before any changes
	int Thread_active, ActionsPerThread, start, end;


	// STEPPING RATE
	dTime.Probability[0] = (double)NumMotorsAttach / tauStep;

	// Detachment Rate Calculation
	if (NumMotorsAttach == 1) {
		// Special case: Only 1 motor is attached
		for (int i = 0; i < numMotors; i++) {
			dTime.detachRate[i] = (Motors.motors[i][0]) ? 1 : 0;
			dTime.Probability[1] += dTime.detachRate[i];
		}
	}
	else {
		assignThreadRanges(numMotors, numThreads, &Thread_active, &ActionsPerThread);
		// Parallel detachment calculation using threads
		for (int t = 0; t < Thread_active; t++) {
			start = t * ActionsPerThread;
			end = (t == (Thread_active - 1)) ? numMotors : (t + 1) * ActionsPerThread;
			threadData[t] = (struct ThreadData){ .Motors = Motors, .dTime = &dTime, .i_Energy = i_Energy, .start = start, .end = end };
			pthread_create(&threads[t], NULL, calculateDetachRates, &threadData[t]);
		}

		// Wait for all threads to complete detachment calculations
		for (int t = 0; t < Thread_active; t++) {
			pthread_join(threads[t], NULL);
		}
	}

	assignThreadRanges(meshLimit_att, numThreads, &Thread_active, &ActionsPerThread);
	// Attachment Rate Calculation
	for (int t = 0; t < Thread_active; t++) {
		start = t * ActionsPerThread;
		end = (t == (Thread_active - 1)) ? meshLimit_att : (t + 1) * ActionsPerThread;
		threadData[t] = (struct ThreadData){ .Motors = Motors, .dTime = &dTime, .i_Energy = i_Energy, .q0 = q0, .start = start, .end = end };
		pthread_create(&threads[t], NULL, calculateAttachRates, &threadData[t]);
	}

	// Wait for all threads to complete attachment calculations
	for (int t = 0; t < Thread_active; t++) {
		pthread_join(threads[t], NULL);
	}

	return dTime;
}

// Function to select one of the options based on the values as proportional probabilities
int monteCarloSelect(double* values, int size, MTRand* rand) {
	// Calculate the total sum of the values
	double total = 0.0;
	for (int i = 0; i < size; i++) {
		total += values[i];
	}
	// Generate a random number between 0 and the total sum
	double randNum = genRand(rand) * total;



	// Determine which option to select based on the random number
	double cumulativeSum = 0.0;
	for (int i = 0; i < size; i++) {
		cumulativeSum += values[i];  // Use values directly
		if (randNum < cumulativeSum) {
			return i; // Select option i
		}
	}
	return size - 1; // In case of an error (should not happen if probabilities are valid)
}

// Function to randomly select one of the attached motors
int selectAttachedMotor(struct MotorData Motors) {
	int numAttached = attMotors(Motors);  // Get the number of attached motors
	if (numAttached == 0) {
		printf("No attached motors found.\n");
		return -1;  // No attached motors, return error code
	}
	// Dynamically allocate memory for attachedMotors array
	int* attachedMotors = (int*)malloc(numAttached * sizeof(int));
	if (!attachedMotors) {
		printf("Memory allocation failed.\n");
		return -1;  // Memory allocation failed
	}
	int attachedCount = 0;
	// Populate the array with indices of attached motors
	for (int i = 0; i < numMotors; i++) {
		if (Motors.motors[i][0]) {  // Assuming motors[i][0] means "attached"
			attachedMotors[attachedCount++] = i; // Store attached motor index
		}
	}
	if (attachedCount == 0) {
		// Avoid returning invalid index if somehow no motors were attached
		free(attachedMotors);
		return -1;
	}
	// Proper random selection from attached motors
	int randomIndex = rand() % attachedCount;
	if (randomIndex == attachedCount) {
		randomIndex--;
	}
	int selectedMotor = attachedMotors[randomIndex];
	free(attachedMotors);  // Free allocated memory
	return selectedMotor;  // Return the selected motor index
}

void* calculateWalkProb(void* arg) {
	struct ThreadData* data = (struct ThreadData*)arg; // Cast input to ThreadData
	struct MotorData Motors = data->Motors;           // Local copy of Motors
	double i_Energy = data->i_Energy;                 // Initial energy
	int start = data->start;                          // Start index for thread
	int end = data->end;                              // End index for thread
	int whoWalks = data->whoWalks;                    // Index of the motor that walks
	double* walkProb = data->walkProb;                // Pointer to the walkProb array
	double TempX = Motors.motors[whoWalks][4];
	double TempY = Motors.motors[whoWalks][5];

	// Temporary array to store probabilities
	double localWalkProb[meshLimit_walk];

	for (int k = start; k < end; k++) {
		Motors.motors[whoWalks][4] += Xmesh_walk[k]; // Walk in the x-axis
		Motors.motors[whoWalks][5] -= Ymesh_walk[k] * cos(asin(Motors.motors[whoWalks][5] / radiusMT)); // Update y-axis position

		if (Motors.motors[whoWalks][5] < -radiusMT || Motors.motors[whoWalks][5] > radiusMT) {
			localWalkProb[k] = 0; // Out-of-bound condition
		}
		else {

			struct BalanceMovment BeadBalance = GradientDescent(Motors); // Calculate new energy
			double prob = 1 / (1 + exp((-fs * Xmesh_walk[k] + BeadBalance.Energy - i_Energy) / (kb * effTemp))) * xy_WalkProb[k];
			localWalkProb[k] = prob; // Store in local array

		}

		// Reset motor position
		Motors.motors[whoWalks][4] = TempX;
		Motors.motors[whoWalks][5] = TempY;
	}

	// Update the shared walkProb array with a single lock
	pthread_mutex_lock(&mutex);
	for (int k = start; k < end; k++) {
		walkProb[k] = localWalkProb[k];
	}
	pthread_mutex_unlock(&mutex);



	pthread_exit(NULL);
}

void walkAlgoritm(int whoWalks, struct MotorData Motors, double walkProb[meshLimit_walk]) {
	int NumMotorsAttach = attMotors(Motors); // Number of attached motors
	double i_Energy = (NumMotorsAttach == 1) ? dF : energy(Motors); // Calculate energy before any changes
	int Thread_active, ActionsPerThread, start, end;
	assignThreadRanges(meshLimit_walk, numThreads, &Thread_active, &ActionsPerThread);
	// Parallel detachment calculation using threads
	for (int t = 0; t < Thread_active; t++) {
		start = t * ActionsPerThread;
		end = (t == (Thread_active - 1)) ? meshLimit_walk : (t + 1) * ActionsPerThread;
		threadData[t] = (struct ThreadData){
			.Motors = Motors,          // Pass by value
			.i_Energy = i_Energy,
			.whoWalks = whoWalks,
			.start = start,
			.end = end, // Handle edge case for last thread
			.walkProb = walkProb
		};
		pthread_create(&threads[t], NULL, calculateWalkProb, &threadData[t]);
	}

	for (int t = 0; t < Thread_active; t++) {
		pthread_join(threads[t], NULL);
	}
}

double calculateAngle(double Motor1[6], double Motor2[6]) {
	double String1[2], String2[2]; // String vectors [0]==x, [1]==y
	String1[0] = Motor1[4] - Motor1[1];
	String1[1] = Motor1[5] - Motor1[2];
	String2[0] = Motor2[4] - Motor2[1];
	String2[1] = Motor2[5] - Motor2[2];
	return atan2(String2[1], String2[0]) - atan2(String1[1], String1[0]);
}


double normalizeAngle(double angle) {
	if (isnan(angle)) {
		fprintf(stderr, "Error: Angle is NaN\n");
		return 0.0; // or an appropriate default
	}
	while (angle > PI) angle -= 2 * PI;
	while (angle < -PI) angle += 2 * PI;
	return angle;
}
MotorPair* createNode(int motor1, int motor2) {
	MotorPair* newNode = (MotorPair*)malloc(sizeof(MotorPair));
	if (newNode == NULL) {
		// Handle memory allocation failure
		fprintf(stderr, "Memory allocation failed for MotorPair.\n");
		return NULL;
	}

	newNode->motor1 = motor1;
	newNode->motor2 = motor2;
	newNode->angle_sum = 0.0;
	newNode->last_angle = -1.0;
	newNode->next = NULL;

	return newNode;
}
// Add a pair to the linked list
void addMotorPair(MotorPair** head, int motor1, int motor2) {
	MotorPair* newNode = createNode(motor1, motor2);
	newNode->next = *head;
	*head = newNode;
}

// Remove a pair from the linked list
void removeMotorPair(MotorPair** head, int motor1, int motor2, double* angle_mean, int* angle_count, int* coil_count) {
	MotorPair* temp = *head;
	MotorPair* prev = NULL;

	while (temp != NULL) {
		if ((temp->motor1 == motor1 && temp->motor2 == motor2) ||
			(temp->motor1 == motor2 && temp->motor2 == motor1)) {
			*angle_mean += fabs(temp->angle_sum);
			(*angle_count)++;
			if (fabs(temp->angle_sum) >= 2 * PI) {
				(*coil_count)++;
			}
			if (prev == NULL) {
				*head = temp->next;
			}
			else {
				prev->next = temp->next;
			}
			free(temp);
			return;
		}
		prev = temp;
		temp = temp->next;
	}
}



// Update angles in the linked list
void updateMotorAngles(MotorPair* head, double Motors[][6]) {
	MotorPair* temp = head;
	while (temp != NULL) {

		double angle = normalizeAngle(calculateAngle(Motors[temp->motor1], Motors[temp->motor2]));
		if (temp->last_angle != -1.0) {
			temp->angle_sum += normalizeAngle(angle - temp->last_angle);


		}
		temp->last_angle = angle;
		temp = temp->next;
	}
}

double compute_b_value(int N, double TimeIntervalsLength, double msd_values[]) {
	double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;

	// Compute sums for linear regression
	for (int i = 0; i < N; i++) {
		double x = log(TimeIntervalsLength*(i+1));
		double y = log(msd_values[i]);

		sum_x += x;
		sum_y += y;
		sum_xy += x * y;
		sum_x2 += x * x;
	}

	// Compute slope (b value)
	double b = (N * sum_xy - sum_x * sum_y) / (N * sum_x2 - sum_x * sum_x);
	return b;
}

// * MAIN FUNCTION * // 
void main(void)
{
	int	trial, whoWalks, whoAttach, whoDetach, AttachLocation, Step, IndCount, NumMotorsAttach, StepsCount;
	int XpositiveIntTo = 0, XnegativeIntTo = 0, YpositiveIntTo = 0, YnegativeIntTo = 0, motorcount = 0, PitchCountTotal = 0;
	int TimeintervalCount, YpositiveIntTrial, YnegativeIntTrial, PitchCountTiral;
	struct DeltaTime dTime;
	struct BalanceMovment beadMovement;
	struct MotorData Motors;
	struct SimulationResults Results;
	double walkProb[meshLimit_walk];
	double OneMotorWalkProb[meshLimit_walk];
	double PositiveTotal = 0, NegativeTotal = 0;
	struct Statistics Stat = { 0 };
	clock_t start, end;
	double cpu_time_used, StepTime, timeinterval = 0, XstepInterval = 0, YstepInterval = 0, Pitchinterval = 0, PitchLastX = 0, PhiMeshShift;
	int OneMotorFlag = 0;
	// Dynamically allocate angle data
	int angle_count = 0, coil_count = 0;
	_Bool legalmove;
	double XMSDMeanlist[dataIntervalSize] = { 0 }, XMSDlist[dataIntervalSize] = { 0 }, TimeIntervalsLength = 0.27, MSDLastlocation = 0, MSDLastTime = 0, MSDSlope, MSDConst;
	int MSDcount[dataIntervalSize] = { 0 }, MSDmaxTime = 0;
	// Calculate q0.
	Na = (PolymerLength) / (aPolymer);  // Compute Na at runtime
	R0_Globe = sqrt(Na) * aPolymer;  // Compute R0_Globe at runtime
	dF = -attEnergy - kb * Temper * (log(9 * springRelax * pow(aPolymer, 3) / (2 * PI * pow(R0_Globe, 4))) - 3 * (pow(springRelax, 2) / (2 * pow(R0_Globe, 2))));
	q0 = 1 / (1 + exp(-dF / (kb * Temper)));

	uint32_t seed = (uint32_t)time(NULL);
	MTRand rand = seedRand(seed);

#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	numThreads = sysinfo.dwNumberOfProcessors;

#elif defined(__linux__)
	long numCores = sysconf(_SC_NPROCESSORS_ONLN);
	numThreads = (int)numCores;
#endif

	start = clock();  // Record the start time
	initializeThreads();
	allocateSimulationResults(&Results);
	printf("The simulation of %d motors as started \n", numMotors);





	for (trial = 0; trial < numTrials; trial++)
	{

		//reset
		IndCount = 0;
		timeinterval = 0;
		Pitchinterval = 0;
		XstepInterval = 0;
		YstepInterval = 0;
		PitchLastX = 0;
		PitchCountTiral = 0;
		StepsCount = 0;
		TimeintervalCount = 0;
		YpositiveIntTrial = 0;
		YnegativeIntTrial = 0;
		//Spread motors 
		spreadMotors(&Motors, &rand);
		NumMotorsAttach = 1;

		printf("NP number %d\n", trial);

		while (NumMotorsAttach && IndCount < dataSize)
		{


			dTime = deltaTime(Motors);// the function is auto defined as int. Thats why it shows error at the moment
			//printf("the time is %e, %e, %e\n", OneMotorWalkProb[0], dTime.Probability[1], dTime.Probability[2]);

			StepTime = 1 / (dTime.Probability[0] + dTime.Probability[1] + dTime.Probability[2]);
			// Perform one of the following actions:  walk,detach, attach or do nothing

			legalmove = false;

			while (!legalmove) {

				legalmove = true;
				Results.TrailTimetotal[trial] += StepTime; //Add to the total movment time 
				timeinterval += StepTime; //Add to the interval counting time

				// Perform one of the following actions:  walk,detach, attach or do nothing
				switch (monteCarloSelect(dTime.Probability, 3, &rand))
				{
				case 0: //walk
					//pick random motor.
				{
					StepsCount++;
					whoWalks = selectAttachedMotor(Motors);

					if (NumMotorsAttach == 1) {
						if (OneMotorFlag == 0) {
							walkAlgoritm(whoWalks, Motors, OneMotorWalkProb);
							OneMotorFlag = 1;
						}
						Step = monteCarloSelect(OneMotorWalkProb, meshLimit_walk, &rand);

					}
					else {

						//Get a probability list of all the possible steps for the detected motor
						walkAlgoritm(whoWalks, Motors, walkProb);
						// Choose Step with MonteCarlo
						Step = monteCarloSelect(walkProb, meshLimit_walk, &rand);

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
					whoDetach = monteCarloSelect(dTime.detachRate, numMotors, &rand);
					//detach
					Motors.motors[whoDetach][0] = 0;
					Motors.motors[whoDetach][4] = 0;
					Motors.motors[whoDetach][5] = 0;

					for (int i = 0; i < numMotors; i++) {
						if (Motors.motors[i][0] && i != whoDetach) {
							removeMotorPair(&motorPairs, whoDetach, i, &Stat.angle.mean, &angle_count, &coil_count);
						}
					}

					break;
				}

				case 2: //attach

				{
					// Pick the attaching motor

					whoAttach = monteCarloSelect(dTime.AttachProb, numMotors, &rand);
					//Pick the attaching location 
					AttachLocation = monteCarloSelect(dTime.attachRate[whoAttach], meshLimit_att, &rand);
					//attach
					Motors.motors[whoAttach][0] = 1;

					PhiMeshShift = asin(Ymesh_att[AttachLocation] / radiusMT) - fmod(Motors.TotalDeltaY, MTBSphiShift);//The angle that the NP has shifted from the original mesh
					Motors.motors[whoAttach][0] = 1; // Temporaly attach the selected motor
					Motors.motors[whoAttach][4] = Xmesh_att[AttachLocation] - fmod(Motors.TotalDeltaX, MTBSxShift); // Mesh shift on X axis
					Motors.motors[whoAttach][5] = radiusMT * sin(PhiMeshShift);// Mesh shift on Y axis


					for (int i = 0; i < numMotors; i++) {
						if (Motors.motors[i][0] && i != whoAttach) {
							addMotorPair(&motorPairs, whoAttach, i);
						}
					}

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
						//MSDSlope = (beadMovement.NewLocations.TotalDeltaX - MSDLastlocation) / (timeinterval - MSDLastTime);
						//MSDConst = beadMovement.NewLocations.TotalDeltaX - MSDSlope * timeinterval;

						Results.Xvelocity[trial][TimeintervalCount] = XstepInterval / timeinterval;
						Results.Yvelocity[trial][TimeintervalCount] = YstepInterval / timeinterval;
						Results.XIntervalMSD[trial][TimeintervalCount] = beadMovement.NewLocations.TotalDeltaX;

						if (XstepInterval > 0) {
							Results.XvelocityPlus[trial][TimeintervalCount] = XstepInterval / timeinterval;
							XpositiveIntTo++;

						}
						else if (XstepInterval <= 0) {
							Results.XvelocityMinus[trial][TimeintervalCount] = XstepInterval / timeinterval;
							XnegativeIntTo++;


						}


						timeinterval = 0;
						XstepInterval = 0;

						YstepInterval = 0;
						TimeintervalCount++;
					}
					else {
						MSDLastlocation = beadMovement.NewLocations.TotalDeltaX;
						MSDLastTime = timeinterval;
					}
				}

			}

			NumMotorsAttach = attMotors(Motors);

			updateMotorAngles(motorPairs, Motors.motors);
			// from this part on all the code is dedicated to documentation of the step data 

			//Calculation for mean motors attached
			if (NumMotorsAttach) {

				Results.IterationMotorsCount[trial][IndCount] = (double)NumMotorsAttach; //Create a list of the amount of motors attached in every round
				motorcount++;
			}

			//Calculation for forward or backward steps
			if (beadMovement.NewLocations.IterationDeltaX > 0) {
				PositiveTotal++;
			}
			else {
				NegativeTotal++;
			}

			//Calculation for Y positive and negative velocity in MC intervals 
			if (beadMovement.NewLocations.IterationDeltaY > 0) {
				Results.YvelocityPlus[trial] += beadMovement.NewLocations.IterationDeltaY / StepTime;

				YpositiveIntTrial++;
			}
			else if (beadMovement.NewLocations.IterationDeltaY < 0) {
				Results.YvelocityMinus[trial] += beadMovement.NewLocations.IterationDeltaY / StepTime;
				YnegativeIntTrial++;

			}

			//Calculation for Helical pitch 
			if (Pitchinterval > (2 * PI)) {
				Results.HelicalPitch[trial][PitchCountTiral] = fabs(beadMovement.NewLocations.TotalDeltaX - PitchLastX);
				PitchLastX = beadMovement.NewLocations.TotalDeltaX;
				PitchCountTiral++;
				PitchCountTotal++;
				Pitchinterval = 0;
			}

			IndCount++;

		}

		Results.TrailXtotal[trial] = beadMovement.NewLocations.TotalDeltaX; //Create a list of the X length
		Results.TrailYtotal[trial] = beadMovement.NewLocations.TotalDeltaY; //Create a list of the Y length
		if (YnegativeIntTrial) Results.YvelocityMinus[trial] /= YnegativeIntTrial; // Divide the negative Y velocity by the amount of velocities that measured
		if (YpositiveIntTrial) Results.YvelocityPlus[trial] /= YpositiveIntTrial; // Divide the positive Y velocity by the amount of velocities that measured
		Results.StepsNum[trial] = StepsCount;



	}


	for (int i = 0; i < numTrials; i++) {
		Stat.X.mean += Results.TrailXtotal[i];
		Stat.Y.mean += Results.TrailYtotal[i];
		Stat.Time.mean += Results.TrailTimetotal[i];
		Stat.Steps.mean += (double)Results.StepsNum[i];
		Stat.VelocityYplus.mean += Results.YvelocityPlus[i];
		Stat.VelocityYminus.mean += Results.YvelocityMinus[i];

		for (int j = 0; j < dataSize; j++) {
			if (Results.IterationMotorsCount[i][j]) {
				Stat.Motors.mean += Results.IterationMotorsCount[i][j];
			}
		}

		for (int j = 0; j < dataIntervalSize; j++) {
			if (Results.HelicalPitch[i][j]) Stat.HelicalPitch.mean += Results.HelicalPitch[i][j];
			if (Results.Xvelocity[i][j]) Stat.VelocityX.mean += Results.Xvelocity[i][j];
			if (Results.Yvelocity[i][j]) Stat.VelocityY.mean += Results.Yvelocity[i][j];
			if (Results.XvelocityPlus[i][j]) Stat.VelocityXplus.mean += Results.XvelocityPlus[i][j];
			if (Results.XvelocityMinus[i][j]) Stat.VelocityXminus.mean += Results.XvelocityMinus[i][j];
		}
	}

	// Compute means separately
	Stat.X.mean /= numTrials;
	Stat.Y.mean /= numTrials;
	Stat.Time.mean /= numTrials;
	Stat.Steps.mean /= numTrials;
	Stat.VelocityYplus.mean /= numTrials;
	Stat.VelocityYminus.mean /= numTrials;
	Stat.Motors.mean /= motorcount;
	Stat.HelicalPitch.mean /= PitchCountTotal;
	Stat.VelocityX.mean /= (XpositiveIntTo + XnegativeIntTo);
	Stat.VelocityY.mean /= (XpositiveIntTo + XnegativeIntTo);
	Stat.VelocityXplus.mean /= XpositiveIntTo;
	Stat.VelocityXminus.mean /= XnegativeIntTo;

	Stat.angle.mean /= angle_count;

	if (numTrials != 0) {
		for (int i = 0; i < numTrials; i++) // Calculate the STD values of all the NP trails 
		{
			Stat.X.STD += pow(Results.TrailXtotal[i] - Stat.X.mean, 2) / (numTrials - 1);
			Stat.Y.STD += pow(Results.TrailYtotal[i] - Stat.Y.mean, 2) / (numTrials - 1);
			Stat.Steps.STD += pow(Results.StepsNum[i] - Stat.Steps.mean, 2) / (numTrials - 1);
			Stat.Time.STD += pow(Results.TrailTimetotal[i] - Stat.Time.mean, 2) / (numTrials - 1);
			Stat.VelocityYplus.STD += pow(Results.YvelocityPlus[i] - Stat.VelocityYplus.mean, 2) / (numTrials - 1);
			Stat.VelocityYminus.STD += pow(Results.YvelocityMinus[i] - Stat.VelocityYminus.mean, 2) / (numTrials - 1);

			for (int j = 0; j < dataSize; j++) {

				if (Results.IterationMotorsCount[i][j]) Stat.Motors.STD += pow(Results.IterationMotorsCount[i][j] - Stat.Motors.mean, 2) / (motorcount - 1);

			}
			for (int j = 0; j < dataIntervalSize; j++) {

				if (Results.HelicalPitch[i][j]) Stat.HelicalPitch.STD += pow(Results.HelicalPitch[i][j] - Stat.HelicalPitch.mean, 2) / ((PitchCountTotal)-1);
				if (Results.Xvelocity[i][j]) Stat.VelocityX.STD += pow(Results.Xvelocity[i][j] - Stat.VelocityX.mean, 2) / ((XpositiveIntTo + XnegativeIntTo) - 1);
				if (Results.Yvelocity[i][j]) Stat.VelocityY.STD += pow(Results.Yvelocity[i][j] - Stat.VelocityY.mean, 2) / ((XpositiveIntTo + XnegativeIntTo) - 1);
				if (Results.XvelocityPlus[i][j]) Stat.VelocityXplus.STD += pow(Results.XvelocityPlus[i][j] - Stat.VelocityXplus.mean, 2) / (XpositiveIntTo - 1);
				if (Results.XvelocityMinus[i][j]) Stat.VelocityXminus.STD += pow(Results.XvelocityMinus[i][j] - Stat.VelocityXminus.mean, 2) / (XnegativeIntTo - 1);

			}
		}
		Stat.Motors.STD = sqrt(Stat.Motors.STD);
		Stat.X.STD = sqrt(Stat.X.STD);
		Stat.Y.STD = sqrt(Stat.Y.STD);
		Stat.Time.STD = sqrt(Stat.Time.STD);
		Stat.HelicalPitch.STD = sqrt(Stat.VelocityX.STD);
		Stat.VelocityX.STD = sqrt(Stat.VelocityX.STD);
		Stat.VelocityY.STD = sqrt(Stat.VelocityY.STD);
		Stat.VelocityXplus.STD = sqrt(Stat.VelocityXplus.STD);
		Stat.VelocityXminus.STD = sqrt(Stat.VelocityXminus.STD);
		Stat.VelocityYplus.STD = sqrt(Stat.VelocityYplus.STD);
		Stat.VelocityYminus.STD = sqrt(Stat.VelocityYminus.STD);
		Stat.Steps.STD = sqrt(Stat.Steps.STD);
	}


	Stat.angle.STD = sqrt(Stat.angle.STD);
	Stat.angle.SEM = Stat.angle.STD / sqrt(angle_count);



	Stat.X.SEM = Stat.X.STD / sqrt(numTrials);
	Stat.Y.SEM = Stat.Y.STD / sqrt(numTrials);
	Stat.Time.SEM = Stat.Time.STD / sqrt(numTrials);
	Stat.Motors.SEM = Stat.Motors.STD / sqrt(motorcount);
	Stat.HelicalPitch.SEM = Stat.HelicalPitch.STD / sqrt(PitchCountTotal);
	Stat.VelocityX.SEM = Stat.VelocityX.STD / sqrt((XpositiveIntTo + XnegativeIntTo));
	Stat.VelocityY.SEM = Stat.VelocityY.STD / sqrt((XpositiveIntTo + XnegativeIntTo));
	Stat.VelocityXplus.SEM = Stat.VelocityXplus.STD / sqrt(XpositiveIntTo);
	Stat.VelocityXminus.SEM = Stat.VelocityXminus.STD / sqrt(XnegativeIntTo);
	Stat.VelocityYplus.SEM = Stat.VelocityYplus.STD / sqrt(numTrials);
	Stat.VelocityYminus.SEM = Stat.VelocityYminus.STD / sqrt(numTrials);
	Stat.Steps.SEM = Stat.Steps.STD / sqrt(numTrials);

	// Claculate the MSD
	MSDmaxTime = Stat.Time.mean / TimeIntervalsLength;

	// First pass : Compute mean positions and counts
	for (int i = 0; i < numTrials; i++) {
		for (int j = 0; j <= MSDmaxTime; j++) {
			if (j < dataIntervalSize && Results.XIntervalMSD[i][j] != 0.0) {
				XMSDMeanlist[j] += Results.XIntervalMSD[i][j];
				MSDcount[j]++;
			}
		}
	}

	// Compute means
	for (int j = 0; j <= MSDmaxTime; j++) {
		if (MSDcount[j] > 0) {
			XMSDMeanlist[j] /= MSDcount[j];
		}
	}

	// Second pass: Compute MSD
	for (int i = 0; i < numTrials; i++) {
		for (int j = 0; j <= MSDmaxTime; j++) {
			if (j < dataIntervalSize && Results.XIntervalMSD[i][j] != 0.0 && MSDcount[j] > 0) {
				XMSDlist[j] += pow(Results.XIntervalMSD[i][j] - XMSDMeanlist[j], 2);
			}
		}
	}

	// Normalize MSD by counts
	for (int j = 0; j <= MSDmaxTime; j++) {
		if (MSDcount[j] > 0) {
			XMSDlist[j] /= MSDcount[j];
		}
	}

	double b_value= compute_b_value(MSDmaxTime,TimeIntervalsLength, XMSDlist);



	//save data
	FILE* file = fopen("output.txt", "w");  // Open the file in write mode
	if (file == NULL) {
		printf("Error opening file!\n");
		return 1;
	}
	// End time
	end = clock();  // Record the end time

	// Calculate the elapsed time
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

	// Print the runtime
	printf("Time taken: %f seconds\n", cpu_time_used);


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
	fprintf(file, "%-25s %-20.2e %-20.2e %-20.2e\n", "Number of steps", Stat.Steps.mean, Stat.Steps.SEM, Stat.Steps.STD);
	fprintf(file, "%-25s %-20.2e \n", "Angle between motors", Stat.angle.mean);
	fprintf(file, "%-25s %15.2e | %15.2e\n", "Fractions of forward & backward steps",
		(double)PositiveTotal / (double)(PositiveTotal + NegativeTotal),
		(double)NegativeTotal / (double)(PositiveTotal + NegativeTotal));
	fprintf(file, "There was %e 2Pi windings in average each run \n\n", (double)coil_count / numTrials);
	fprintf(file, "This is the MSD LIST \n");
	for (int i = 0; i < dataIntervalSize; i++) {
		if (XMSDlist[i]) {
			fprintf(file, "%e,", XMSDlist[i]);
		}
	}
	fprintf(file, " \n This is the b value %e \n",b_value);


	fclose(file);  // Close the file

	freeSimulationResults(&Results);
	// Free the linked list
	while (motorPairs != NULL) {
		MotorPair* temp = motorPairs;
		motorPairs = motorPairs->next;
		free(temp);
	}
	cleanupThreads();
}