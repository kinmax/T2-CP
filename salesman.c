
// Brute force and heuristic to the Salesman problem
// Calculate the cylle from one city (Karlsruhe) to all others and back

/*--------------------------- INCLUDES ---------------------------*/

#include "math.h"
#include "time.h"
#include "stdio.h"

/*---------------------------- DEFINES ---------------------------*/

#define NUM_CITIES 12

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
                           {"Hannover" , {  81, 331}},
                           {"Hamburg"  , {  95, 450}},
                           {"Stuttgart", {  47, -22}},
                           {"Ulm",       { 100, -58}},
                           {"Nurenberg", { 164,  42}}  // 12 - 25.00 seconds to calculate solutions
                           /*{"Bremmen",   {  25, 403}},  // 13 - 303.00 seconds to calculate solutions (5 minutes)
                           {"Dortmund",  { -58, 244}}   // 14 - 4258.00 seconds to calculate solutions (70 minutes)*/
                          }; // horizontal and vertical distances in Km
   
                          
SOLUTION final_solution;    /* Final Order of the Cities  */ 
double best_result = 1000L; /* Final cost (Global Minima) */                    
                         
/*--------------------------- PROTOTYPES -------------------------*/

double calc_dist( COORD point1, COORD point2 );
double calc_total_dist( SOLUTION solution ); 
void print_indexes ( SOLUTION solution );
void print_names ( SOLUTION solution );
unsigned int fat ( int number ); 
void perm_solution ( SOLUTION solution, int first, int last );
void heu_mais_perto ( void );
                            

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
        printf("%s -> ",cities[final_solution[i]].name);
    printf("%s\n",cities[final_solution[0]].name);

} /* print_names */

void perm_solution ( SOLUTION solution, int first, int last ) {
    int i, j,temp; 
    double custo;

    if (first == last) {
        //print_indexes(solution);
        custo = calc_total_dist(solution);
        // printf("\nCusto %3.2f\n", custo );
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

void heu_mais_perto (void) {
    /* rotina heuristica para a solucao aproximada do problema   */
    /* a partir da cidade inicial vou visitando sempre a cidade  */
    /* mais proxima ate esgotar as possibilidades; entao retorno */
  
    int    proxima_cidade, pos_solucao, j;
    double distancia, menor_distancia, custo_total;
    COORD  atual;
    SOLUTION visitadas;
  
    /* inicializo */
  
    pos_solucao = 1;  /* 0 ja eh Karlsruhe */
    final_solution[0] = 0;
  
    visitadas[0] = 1;
    for ( j=1 ; j < NUM_CITIES ; j++ )
        visitadas[j] = 0;
      
    /* posicao atual eh Karlsruhe */
  
    atual.x = cities[0].coord.x;
    atual.y = cities[0].coord.y;
  
    /* acho a cidade mais perto para cada posicao do vetor solucao */
  
    for ( pos_solucao = 1 ; pos_solucao < NUM_CITIES ; pos_solucao++ ) {
        /* procuro cidade mais perto do ponto atual,   */
        /* que ainda nao tenha sido visitada           */
      
        menor_distancia = 1000;
        for ( j = 1 ; j < NUM_CITIES ; j++ ) {
            /* se cidade nao foi visitada ainda */
          
            if ( visitadas[j] == 0 ) {
                /* calculo distancia para posicao atual */
                    
                distancia = calc_dist( atual, cities[j].coord);
                
                /* se distancia menor, marco */
                    
                if ( distancia < menor_distancia ) {
                    menor_distancia = distancia;
                    proxima_cidade = j;
                }
            }
        }  
    /* atualizo atual, vetor solucao e vetor visitadas */
      
    atual.x = cities[proxima_cidade].coord.x;
    atual.y = cities[proxima_cidade].coord.y;
    final_solution[pos_solucao] = proxima_cidade; 
    visitadas[proxima_cidade] = 1;
    }
      
    /* calculo custo da solucao encontrada */
  
    custo_total = calc_total_dist(final_solution);
  
    /* atualizo variavel global de resultado (custo) */
  
    best_result = custo_total;         
} /* heu_mais_perto */

/*----------------------------------------------------------------*
|                               main                              |
*-----------------------------------------------------------------*/
                                
void main ( void ) {
    SOLUTION initial_solution = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    double test_result = 0L;
    int i;
    time_t t1, t2;
  
    printf("\nVisiting %d cities from %s (brute force %u permutations):\n", NUM_CITIES-1, cities[0].name, fat(NUM_CITIES-1));

    time(&t1); /* disparo relogio */

    heu_mais_perto();
    printf("\nMinimal for heuristic solution: %3.2f Km for roundtrip.\n\n", best_result );
    //print_indexes (final_solution);
    print_names (final_solution);
  
    perm_solution(initial_solution, 1, NUM_CITIES-1);		  
    printf("\nMinimal for brute force solution: %3.2f Km roundtrip.\n\n", best_result );
    //print_indexes (final_solution);   
    print_names (final_solution);

    time(&t2); /* paro relogio */
  
    printf("\n%.2lf seconds to calculate solutions.\n\n", difftime(t2, t1));
}                                
                                
                        