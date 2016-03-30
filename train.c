#include <stdio.h>
#include "hmm.h"

inline void forward( HMM *hmm, char o[], int ob_len , double** a)
{   
    int i, j, t;
    
    for (i = 0; i < hmm->observ_num; i++)    // Alpha 1 Initial
        a[0][i] = hmm->initial[i] * hmm->observation[(int)(o[0]-'A')][i];

    for (t = 1; t < ob_len; t++) // Push forward
    {   for (j = 0; j < hmm->state_num; j++)
        {   for (i = 0; i < hmm->state_num; i++)
                a[t][j] += a[t-1][i] * hmm->transition[i][j];
            a[t-1][j] *= hmm->observation[(int)(o[t]-'A')][j];
        }
    }

}


int main(int argc, char* argv[])
{
    if (argc < 5) 
    {   fprintf( stderr, "Error: Arguments needed\n" );
        fprintf( stderr, "Usage: ./train [iteration] [Initial Model file] [Input sequence file] [output file]\n" );
        return 0;
    }

    HMM model;
    loadHMM( &model, argv[2] ); // Load the hmm

    FILE* f = fopen(argv[3], "r");
    char a[MAX_SEQ];

    double alpha[][model.state_num];
    
    
    return 0;
}
