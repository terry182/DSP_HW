#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "hmm.h"

int main(int argc, char* argv[])
{
    if (argc < 5) 
    {   fprintf( stderr, "Error: Arguments needed\n" );
        fprintf( stderr, "Usage: ./train [iteration] [Initial Model file] [Input sequence file] [output file]\n" );
        return 0;
    }

    HMM model;
    loadHMM( &model, argv[2] ); // Load the hmm
    
    int iter = atoi(argv[1]);
    fprintf( stderr, "iter = %d\n", iter);
    
    double accuGamma[MAX_SEQ][MAX_STATE], accuEpsilon[MAX_SEQ][MAX_STATE][MAX_STATE];
    double observGamma[MAX_OBSERV][MAX_STATE];

    // Initialize
    for (int i = 0; i < MAX_SEQ; ++i)
       for (int j = 0; j < MAX_STATE; ++j)
        {    accuGamma[i][j] = 0.0;
             for (int k = 0; k < MAX_STATE; ++k)
                 accuEpsilon[i][j][k] = 0.0;
        }

    while (iter--)
    {   FILE* f = open_or_die(argv[3], "r");
        char o[MAX_SEQ];
        int N = 0;
        int T = 0;
        while (fgets(o, MAX_SEQ, f) != NULL)
        {
            ++N;  // # of sample + 1
            T = strlen(o)-1; // Exclude the end-of-line

            // Alpha       
            double alpha[MAX_SEQ][MAX_STATE];
            for (int i = 0; i < model.state_num; ++i) 
                alpha[0][i] = model.initial[i] * model.observation[o[0]-'A'][i];

            for (int t = 0; t < T-1; ++t)
            {   for (int j = 0; j < model.state_num; ++j)
                {   alpha[t+1][j] = 0.0;
                    for (int i = 0; i < model.state_num; ++i)
                        alpha[t+1][j] += alpha[t][i] * model.transition[i][j];
                    alpha[t+1][j] *= model.observation[o[t+1]-'A'][j];
            //        fprintf(stdout, "%e ", alpha[t+1][j]);
                }
             //   fprintf(stdout, "\n");
            }

            // Beta
            double beta[MAX_SEQ][MAX_STATE];
            for (int i = 0; i < model.state_num; ++i)
                beta[T-1][i] = 1.0;
            
            for (int t = T-2; t >= 0; --t)
            {   
                for (int i = 0; i < model.state_num; ++i)
                {   beta[t][i] = 0.0;
                    for (int j = 0; j < model.state_num; ++j)
                        beta[t][i] += model.transition[i][j] * beta[t+1][j] * model.observation[o[t+1]-'A'][j];
                 //   fprintf(stdout, "%e ", beta[t][i]);
                }
                //fprintf(stdout, "\n");
            }

            // Gamma
            double gamma[MAX_SEQ][MAX_STATE];
            for (int t = 0; t < T; ++t)
            {   double all = 0.0;
                for (int i = 0; i < model.state_num; ++i)
                    all += alpha[t][i] * beta[t][i];
                for (int i = 0; i < model.state_num; ++i)
                    gamma[t][i] = alpha[t][i] * beta[t][i] / all;
            }

            // Epsilon
            double epsilon[MAX_SEQ][MAX_STATE][MAX_STATE];
            for (int t = 0; t < T-1; ++t)
            {   double P = 0.0; 
                for (int i = 0; i < model.state_num; ++i)
                    for (int j = 0; j < model.state_num; ++j)
                        P += alpha[t][i] * model.transition[i][j] * model.observation[o[t+1]-'A'][j] * beta[t+1][j];

                for (int i = 0; i < model.state_num; ++i)
                    for (int j = 0; j < model.state_num; ++j)
                      epsilon[t][i][j] = alpha[t][i] * model.transition[i][j] * model.observation[o[t+1]-'A'][j] * beta[t+1][j] / P;
            }
                
            // Acumulate
            for (int t = 0; t < T; ++t)
                for (int i = 0; i < model.state_num; ++i)
                {   accuGamma[t][i] += gamma[t][i];
                    observGamma[o[t]-'A'][i] += gamma[t][i];
                }
            
            for (int t = 0; t < T-1; ++t)
                for (int i = 0; i < model.state_num; ++i)
                    for (int j = 0; j < model.state_num; ++j)
                        accuEpsilon[t][i][j] += epsilon[t][i][j];

        }
        fclose(f);
        
        // Re-estimate
        for (int i = 0; i < model.state_num; ++i)
        {   model.initial[i] = accuGamma[0][i] / N;

            double B = 0.0;
            for (int t = 0; t < T-1; ++t) B += accuGamma[t][i];

            for (int j = 0; j < model.state_num; ++j)
            {   double E = 0.0;
                for (int t = 0; t < T-1; ++t) E += accuEpsilon[t][i][j];

                model.transition[i][j] = E/B;
            }    

            B += accuGamma[T-1][i]; // Shift one

            for (int k = 0; k < model.observ_num; ++k)
                model.observation[k][i] = observGamma[k][i] / B;

        }

    }
    dumpHMM(stderr, &model);

    return 0;
}
