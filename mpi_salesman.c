
// Brute force and heuristic to the Salesman problem
// Calculate the cylle from one city (Karlsruhe) to all others and back

/*--------------------------- INCLUDES ---------------------------*/

#include "math.h"
#include "time.h"
#include "stdio.h"
#include "mpi.h"

/*---------------------------- DEFINES ---------------------------*/

#define NUM_CITIES 8
#define KILL 1

/*-------------------------- STRUCTURES --------------------------*/

typedef struct {
                int x;
                int y;
               } COORD;

typedef struct {
                char name[15];
                COORD coord;
                } CITY;
                
typedef int SOLUTION[NUM_CITIES];
                
/*--------------------------- VARIABLES --------------------------*/

CITY cities[NUM_CITIES] = {
                           {"Karlsruhe", {   0,   0}},
                           {"Koeln"    , { -92, 183}},
                           {"Frankfurt", {  17, 106}},
		                   {"Muenchen" , { 200, -86}},
                           {"Freiburg" , { -36,-108}},
                           {"Dresden"  , { 331, 197}},
                           {"Berlin"   , { 311, 344}},
                           {"Hannover" , {  81, 331}}
                           /*{"Hamburg"  , {  95, 450}},
                           {"Stuttgart", {  47, -22}},
                           {"Ulm",       { 100, -58}},
                           {"Nurenberg", { 164,  42}},  // 12 - 25.00 seconds to calculate solutions
                           {"Bremmen",   {  25, 403}},  // 13 - 303.00 seconds to calculate solutions (5 minutes)
                           {"Dortmund",  { -58, 244}}   // 14 - 4258.00 seconds to calculate solutions (70 minutes)*/
                          }; // horizontal and vertical distances in Km
   
                          
SOLUTION final_solution;    /* Final Order of the Cities  */ 
double best_result = 100000L; /* Final cost (Global Minima) */                    
                         
/*--------------------------- PROTOTYPES -------------------------*/

double calc_dist( COORD point1, COORD point2 );
double calc_total_dist( SOLUTION solution ); 
void print_indexes ( SOLUTION solution );
void print_names ( SOLUTION solution );
unsigned int fat ( int number ); 
void perm_solution ( SOLUTION solution, int first, int last );
                            

double calc_dist( COORD point1, COORD point2 ) {
    double distx, disty, result;
  
    distx = point1.x - point2.x;
    disty = point1.y - point2.y;
    result = sqrt (pow(distx,2) + pow(disty,2));
    return( result );
} /* calc_dist */


double calc_total_dist( SOLUTION solution ) {
    int i;
    double total = 0L;
 
    for ( i = 0 ; i < NUM_CITIES-1; i++ )   // sum distance from each city from this solution to the next one
        total += calc_dist(cities[solution[i]].coord, cities[solution[i+1]].coord);
    total += calc_dist(cities[solution[i]].coord, cities[solution[0]].coord); // sum distance from last city to the first
    return(total);
} /* calc_total_dist */


void print_indexes ( SOLUTION solution ) {
    int i;
 
    for ( i = 0 ; i < NUM_CITIES ; i++ )
        printf("[%d] ", solution[i]);
    printf("\n");

} /* print_indexes */

void print_names ( SOLUTION solution ) {
    int i;
 
    for ( i = 0 ; i < NUM_CITIES ; i++ )
        printf("%s -> ",cities[solution[i]].name);
    printf("%s\n",cities[solution[0]].name);

} /* print_names */

void perm_solution ( SOLUTION solution, int first, int last ) {
    int i, j,temp; 
    double custo;

    if (first == last) {
        //print_indexes(solution);
        custo = calc_total_dist(solution);
        //printf("\nCusto %3.2f\n", custo );
        if ( custo < best_result ) {
            best_result = custo;
            for ( j=0 ; j < NUM_CITIES ; j++ )
                final_solution[j] = solution[j];
            }
        }
    else { 
        for (i = first; i <= last ; i++) { 
            // swap((a + l), (a + i)) 
            temp = solution[first];  
            solution[first] = solution[i];
            solution[i] = temp;
            perm_solution (solution, first+1, last); // permute(a, l + 1, r); 
            // swap((a + l), (a + i))
            temp = solution[first]; 
            solution[first] = solution[i];
            solution[i] = temp;
        } 
    } 
} 
  

/*----------------------------------------------------------------*
|                              fat                                |
*-----------------------------------------------------------------*/

unsigned int fat ( int number ) {
    unsigned int soma;
    int i;

    soma = number;
    for ( i = number-1 ; i > 1 ; i-- )
        soma *= i;
    return ( soma );
}  /* fat */

/*----------------------------------------------------------------*
|                               main                              |
*-----------------------------------------------------------------*/
                                
int main(int argc, char **argv) {

    int i, j, k, l, m, w, x, proc_index;
    int my_rank, proc_n;
    double t1, t2;
    SOLUTION received;

    MPI_Status status;

    MPI_Init (&argc , & argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);

    if(my_rank == 0)
    {
        t1 = MPI_Wtime();
        SOLUTION initial_solution;
        for(i = 0; i < NUM_CITIES; i++)
        {
            initial_solution[i] = i;
        }
        double test_result = 0L;
        SOLUTION workbag[NUM_CITIES*(NUM_CITIES-1)];
        SOLUTION best;
        double best_cost = 100000L;
        double returned_cost;
        w = 0;
        // montando o workbag
        for(i = 0; i < NUM_CITIES; i++)
        {
            for(j = 0; j < NUM_CITIES; j++)
            {
                if(i == j) continue;
                workbag[w][0] = initial_solution[i];
                workbag[w][1] = initial_solution[j];
                l = 2;
                m = 0;
                while(l < NUM_CITIES)
                {
                    if(m == i || m == j)
                    {
                        m++;
                        continue;
                    }
                    workbag[w][l] = m;
                    l++;
                    m++;
                }
                w++;
            }
        }

        // enviando primeira remessa
        w = 0;
        for(i = 1; i < proc_n; i++)
        {
            if(w == NUM_CITIES*(NUM_CITIES-1)) break;
            MPI_Send(&workbag[w], NUM_CITIES, MPI_INT, i, 0, MPI_COMM_WORLD);
            w++;
        }
        x = 0;
        proc_index = 1;
        while(1)
        {
            MPI_Recv(&received, NUM_CITIES, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            x++;
            returned_cost = calc_total_dist(received);
            if(returned_cost < best_cost)
            {
                for(i = 0; i < NUM_CITIES; i++)
                {
                    best[i] = received[i];
                }
                best_cost = returned_cost;
            }
            if(w == NUM_CITIES*(NUM_CITIES-1) && x == w)
            {
                int empty_message = 0;
                for(i = 1; i < proc_n; i++)
                {
                    MPI_Send(&empty_message, 1, MPI_INT, i, KILL, MPI_COMM_WORLD);
                }
                break;
            }
            if(w < NUM_CITIES*(NUM_CITIES-1))
            {
                MPI_Send(&workbag[w], NUM_CITIES, MPI_INT, proc_index, 0, MPI_COMM_WORLD);
                proc_index++;
                if(proc_index == proc_n) proc_index = 1;
                w++;
            }
        }
        t2 = MPI_Wtime();
        printf("\nMinimal for brute force solution: %3.2f Km roundtrip.\n\n", best_cost );
        //print_indexes (final_solution);   
        print_names (best);
        printf("\nTime: %lf seconds\n\n", t2-t1);
        
    }
    else
    {
        while(1)
        {
            MPI_Recv(&received, NUM_CITIES, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(status.MPI_TAG == KILL)
            {
                break;
            }
            perm_solution(received, 2, NUM_CITIES-1);
            MPI_Send(&final_solution, NUM_CITIES, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    
        
    MPI_Finalize();
}                                
                                
                        