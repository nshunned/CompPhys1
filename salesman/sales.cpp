#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>


using namespace std;

double const R = 6371;   // earth's mean radius = 6,371km
// time zones
double const E1 = -85.5; // easter edge (in deg)
double const E2 = -102;  // central edge (in deg)
double const E3 = -114;  // mountain edge (in deg)
bool zone = false;

// simple structure to store city coordinates
typedef struct {
    double lon, lat;
} COORD;

// fill the array of city locations
int GetData(char* fname, COORD* cities) {
    FILE* fp = fopen(fname,"r");
    const int bufsiz=1000;
    char line[bufsiz+1];
    int ncity=0;
    while(1) {
        fgets(line,bufsiz,fp);
        if (line[0] == '#') continue;  // skip comments
        if (feof(fp)) break;
        // we only scan for two numbers at start of each line
        sscanf(line,"%lf %lf",&cities[ncity].lon,&cities[ncity].lat);    
        ncity++;
    }
    fclose(fp);
    return ncity;
}

// print cities to the screen
void PrintCities(int ncity, COORD* cities) {
    printf("Longitude  Latitude\n");
    for (int i=0; i<ncity; i++)
        printf("%lf %lf\n",	cities[i].lon,cities[i].lat);
}

// get the distance between two cities
double GetDistance(const COORD city1, const COORD city2) {
    double lat1 = city1.lat * M_PI/180;
    double lat2 = city2.lat * M_PI/180;
    double lon1 = city1.lon * M_PI/180;
    double lon2 = city2.lon * M_PI/180;
    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;
    double a = sin(dLat/2)*sin(dLat/2) + cos(lat1)*cos(lat2)*sin(dLon/2)*sin(dLon/2);
    double c = 2*atan2(sqrt(a),sqrt(1-a));
    double l = R*c;
    if (zone) {
        lon1 = city1.lon;
        lon2 = city2.lon;
        if((lon1>E1 && lon2<E1) || (lon1<E1 && lon2>E1)) l += 100;
        if((lon1>E2 && lon2<E2) || (lon1<E2 && lon2>E2)) l += 100;
        if((lon1>E3 && lon2<E3) || (lon1<E3 && lon2>E3)) l += 100;
    }
    return l;
}

// get the total distance travelled through a list of cities
double GetTotalDistance(int ncity, const COORD* cities) {
    double distance=0;
    for (int i=0; i<ncity-1; i++) {
        distance += GetDistance(cities[i],cities[i+1]);
    }
    distance += GetDistance(cities[ncity-1],cities[0]);
    return distance;
}

double GetDL(int ncity, COORD* cities, int i, int j) {
    int il = (i==0 ? ncity-1: i-1), ir = (i==ncity-1 ? 0 : i+1);
    int jl = (j==0 ? ncity-1: j-1), jr = (j==ncity-1 ? 0 : j+1);
    double di = GetDistance(cities[il],cities[i]) + GetDistance(cities[i],cities[ir]) + GetDistance(cities[jl],cities[j]) + GetDistance(cities[j],cities[jr]);
    double df = GetDistance(cities[il],cities[j]) + GetDistance(cities[j],cities[ir]) + GetDistance(cities[jl],cities[i]) + GetDistance(cities[i],cities[jr]);
    return df-di;
}

// calculates the inital temperature T
double GetT(int ncity, COORD* cities) {
    double dL=0, dLMax=0;
    for (int i=0; i<ncity*10000; i++) {
        dL = GetDL(ncity,cities,(int)(drand48()*ncity),(int)(drand48()*ncity));
        if (dL > dLMax) dLMax = dL;
    }
    return dLMax *= 100;
}

// update the cities using Metropolis algorithm
// returns 1 if it's a successful reconfiguration
// returns 0 if it's not
int UpdateCities(int ncity, COORD* cities, double T) {
    int i = (int)(drand48()*ncity);
    int j = (int)(drand48()*ncity);
    double dL = GetDL(ncity,cities,i,j);
    if(dL <= 0 || drand48() < exp(-dL/T)) {
        COORD temp = cities[i];
        cities[i] = cities[j];
        cities[j] = temp;
        return 1;
    }
    return 0;
}

int main(int argc, char* argv[]) {
    // seed the random generator
    srand48((long)time(0));
    // max number of cities to read in
    const int NMAX=2500;
    COORD cities[NMAX];


    // read in citites
    if (argc < 2) {
        printf("Please include the name of the data file\n");
        exit (EXIT_FAILURE);
    }
    int ncity = GetData(argv[1],cities);
    printf("Read %d cities from data file\n",ncity);
    // consider time zones
    if (argc >= 3) {
        zone = true;
        cout << "Time zone penalty is turned on" << endl;
    }


    // initial distance
    const double L0 = GetTotalDistance(ncity,cities);
    cout << "The initial distance is: " << L0 << endl;
    // initial temperature
    double T = GetT(ncity,cities);
    cout << "The initial temperature is: " << T << endl;


    // Metropolis
    int success=0;              // # successful reconfigurations per trial
    double l,lMin=L0;           // distances
    COORD newCities[NMAX];      // optimization solution
    // start
    double sec;
    clock_t timer = clock();
    for (int i=0; i<1000000; i++) {
        for (int j=0; j<100*ncity; j++) {
            success += UpdateCities(ncity,cities,T);
            if (success==3) break;
        }
        success=0;
        T *= 0.999;

        l = GetTotalDistance(ncity,cities);
        if (l < lMin) {
            sec = (double)(clock()-timer)/CLOCKS_PER_SEC;
            lMin = l;
            copy(cities,cities+ncity,newCities);
            printf("The distance after optimization is: %-7f km (i = %-6i T = %-7.3f time (s) = %-6.2f)\n", l,i,T,sec);
        }
    }
    cout << (double)(clock()-timer)/CLOCKS_PER_SEC << " Fin." << endl;


    // name of output file
    string fileName(argv[1]); 
    fileName = fileName.substr(0,fileName.size()-4) + "_optimized";
    if(zone) fileName += "_zoned";
    fileName += ".dat";
    
    // print to file
    ofstream outFile(fileName.c_str());
    outFile << "#lontitude \t latitude" << endl;
    for (int i=0; i<ncity; i++) {
        outFile << newCities[i].lon << "\t \t \t" << newCities[i].lat << endl;
    }
    outFile.close();

    return 0;
}