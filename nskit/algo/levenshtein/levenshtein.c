# include <stdint.h>



double tmin(double a, double b, double c){
    double s;
    if (a<b){s = a;} else {s = b;}
    if (c<s){return c;} else {return s;}
}


double weighted_levenshtein(const char* a, const char* b, 
                              int aN, int bN, 
                              double ins, 
                              double rm, 
                              double sub
                             ){
    aN++;
    bN++;
    
    double result;
    double diagonal_cost = 0.;
    double* A = malloc(sizeof(double)*aN*bN);
    for (int i=1; i<aN; i++){A[i*bN] = (double)i * rm;}
    for (int j=0; j<bN; j++){A[j] = (double)j * ins;}
    
    for (int i=1; i<aN; i++){
        for (int j=1; j<bN; j++){
            
            if ( (*(a+i-1)) == (*(b+j-1)) ){
                diagonal_cost = 0.;
            } else {
                diagonal_cost = sub;
            }
            
            A[i*bN + j] = tmin((A[(i-1)*bN + (j-1)] + diagonal_cost), 
                               (A[i*bN + (j-1)] + ins), 
                               (A[(i-1)*bN + j] + rm)
                              );
        }
    }
    result = A[aN*bN - 1];
    free(A);
    
    return result;
}



