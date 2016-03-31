#include <stdio.h> 
#include <string.h>
#include "hmm.h"

static inline void forward( HMM *hmm, char o[], int ob_len , double a[][MAX_STATE])
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

static inline void backward(HMM *hmm, char o[], int ob_len, double b[][MAX_STATE])
{   
    for (int i = 0; i < hmm->state_num; ++i) b[ob_len-1][i] = 1.;

    for (int t = ob_len-2; t >= 0; --t)
        for (int i = 0; i < hmm->state_num; ++i)
            for (int j = 0; j < hmm->state_num; ++j)
                b[t][i] += hmm->transition[i][j]*hmm->observation[(int)(o[t]-'A')][j]*b[t+1][j];
    

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
    double alpha[MAX_SEQ][MAX_STATE], beta[MAX_SEQ][MAX_STATE];
    double gamma[MAX_SEQ][MAX_STATE], e[MAX_SEQ][MAX_STATE][MAX_STATE];
    int sample_num = 0, len;

    while (fgets(a, MAX_SEQ, f) != NULL)
    {   ++sample_num;
        len = strlen(a);
        forward(&model,  a,  len, alpha);
        backward(&model, a,  len, beta);

        // Calculatin Gamma
        for (int t = 0; t < len; ++t)        
        {   double P = 0;
            for (int i = 0; i < model.state_num; ++i) 
                 P += alpha[t][i]*beta[t][i];
            
            for (int i = 0; i < model.state_num; ++i)
               gamma[t][i] += alpha[t][i]*beta[t][i]/P;
        }

        // Calculatin Epsilon
        for (int t = 0; t < len-1; ++t)
        {   double P = 0;
            for (int i = 0; i < model.state_num; ++i)
            {    for (int j = 0; j < model.state_num; ++j)

                    P += model.transition[i][j]*model.observation[(int)a[t+1]-'A'][j]*beta[t+1][j];
                    P *= alpha[t][i];
            } 

            for (int i = 0; i < model.state_num; ++i)
                for (int j = 0; j < model.state_num; ++j)
                    e[t][i][j] = alpha[t][i]*model.transition[i][j]*model.observation[(int)a[t+1]-'A'][j]*beta[t+1][j];
        }

    }
    
    for (int i = 0; i < model.state_num; ++i)
    {    model.initial[i] = gamma[0][i]/sample_num;
         
        double B = 0;
        for (int t = 0; t < len; ++t) B += gamma[t][i];

         for (int j = 0; j < model.state_num; ++j)
         {  double T = 0, B = 0;
             for (int t = 0; t < len-1; ++t)
                   T += e[t][i][j];
            model.transition[i][j] = T / B;
         }


    }

    printf("%s", a);

    return 0;
}
