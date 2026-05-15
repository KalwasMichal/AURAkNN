#define R_NO_REMAP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <omp.h>

#define MIN(a,b) ((a) < (b)? (a):(b))
typedef struct {
    double distance;
    int row_index;
} Node;

void DownHeap(Node* root, int i, size_t n)
{
    int k= 2*i+1;
    Node v=root[i];
    while(k<n)
    {
        if(k+1<n)
        {
            if(root[k+1].distance>root[k].distance) k++;
        }
        if(root[k].distance >v.distance)
        {
            root[i]=root[k];
            i=k;
            k=2*i+1;
        }
        else break;
    }
    root[i]=v;

}
void build_heap(Node* root, size_t n)
{
    if(n <=1) return;
    for(int i=(n-1)/2; i>=0;i--)
    {
        DownHeap(root,i,n);
    }

}




//Różne metryki do obliczania sąsiadów 

Node* euclid_distance(double* matrix,int targetRow, int numRows, int numCols,int k, int* neigborsFound)
{
    double* distances = calloc(numRows,sizeof(double));
    if(distances==NULL) return NULL;
    int* validCols= calloc(numRows,sizeof(int));
    if(validCols==NULL)
    {
        free(distances);
        return NULL;
    }

    for(int i=0;i<numCols;i++)
    {
        size_t offset= (size_t)i * numRows;
        double targetVal=matrix[targetRow+offset];
        if(ISNA(targetVal)) continue;
        
        for(int j=0;j<targetRow;j++)
        {
           double curr= matrix[j+offset];
            if(!ISNA(curr))
            {
                double diff= targetVal-curr;
                distances[j]+=diff*diff;
                validCols[j]++;
            } 
        }
        for(int j=targetRow+1;j<numRows;j++)
        {
            double curr= matrix[j+offset];
            if(!ISNA(curr))
            {
                double diff= targetVal-curr;
                distances[j]+=diff*diff;
                validCols[j]++;
            }
        }
    }
    for(int i=0;i<numRows;i++) // skalujemy dystanse
    {
        if(i==targetRow || validCols[i]==0) continue;
        double scalingFactor = (double) numCols / (double) validCols[i];
        distances[i]= distances[i]*scalingFactor;
    }
    int Newk=MIN(k,numRows-1);
    Node* heap=malloc(Newk*sizeof(Node));
    if(heap==NULL)
    {
        free(distances);
        free(validCols);
        return NULL;
    }
    int rowIndex=0;
    int heapCount=0;
    while(rowIndex<numRows && heapCount<Newk)
    {
        if(rowIndex==targetRow || validCols[rowIndex]==0)
        {
            rowIndex++;
            continue;
        }
        heap[heapCount].distance=distances[rowIndex];
        heap[heapCount].row_index=rowIndex;   
        heapCount++;
        rowIndex++;
    }
    if(heapCount<Newk) Newk=heapCount;
    *neigborsFound = Newk;
    if(Newk == 0)
    {
        free(distances);
        free(validCols);
        free(heap);
        return NULL;
    }


    build_heap(heap,Newk);

    for(int i=rowIndex;i<numRows;i++)
    {
        if(i==targetRow || validCols[i]==0) continue;
        if(distances[i]<heap[0].distance)
        {
            heap[0].distance=distances[i];
            heap[0].row_index=i;
            DownHeap(heap,0,Newk);
        }
    }
    free(distances);
    free(validCols);
    return heap;
}

Node* manhattan_distance(double* matrix,int targetRow, int numRows, int numCols,int k, int* neigborsFound)
{
    double* distances = calloc(numRows,sizeof(double));
    if(distances==NULL) return NULL;
    int* validCols= calloc(numRows,sizeof(int));
    if(validCols==NULL)
    {
        free(distances);
        return NULL;
    }

    
    for(int i=0;i<numCols;i++)
    {
        size_t offset= (size_t)i * numRows;
        double targetVal=matrix[targetRow+offset];
        if(ISNA(targetVal)) continue;
        
        for(int j=0;j<targetRow;j++)
        {
           double curr= matrix[j+offset];
            if(!ISNA(curr))
            {
                double diff= targetVal-curr;
                distances[j]+=fabs(diff);
                validCols[j]++;
            } 
        }
        for(int j=targetRow+1;j<numRows;j++)
        {
            double curr= matrix[j+offset];
            if(!ISNA(curr))
            {
                double diff= targetVal-curr;
                distances[j]+=fabs(diff);
                validCols[j]++;
            }
        }
    }
    for(int i=0;i<numRows;i++) // skalujemy dystanse
    {
        if(i==targetRow || validCols[i]==0) continue;
        double scalingFactor = (double) numCols / (double) validCols[i];
        distances[i]= distances[i]*scalingFactor;
    }
    int Newk=MIN(k,numRows-1);
    Node* heap=malloc(Newk*sizeof(Node));
    if(heap==NULL)
    {
        free(distances);
        free(validCols);
        return NULL;
    }
    int rowIndex=0;
    int heapCount=0;
    while(rowIndex<numRows && heapCount<Newk)
    {
        if(rowIndex==targetRow || validCols[rowIndex]==0)
        {
            rowIndex++;
            continue;
        }
        heap[heapCount].distance=distances[rowIndex];
        heap[heapCount].row_index=rowIndex;   
        heapCount++;
        rowIndex++;
    }
    if(heapCount<Newk) Newk=heapCount;
    *neigborsFound = Newk;
    if(Newk == 0)
    {
        free(distances);
        free(validCols);
        free(heap);
        return NULL;
    }


    build_heap(heap,Newk);

    for(int i=rowIndex;i<numRows;i++)
    {
        if(i==targetRow || validCols[i]==0) continue;
        if(distances[i]<heap[0].distance)
        {
            heap[0].distance=distances[i];
            heap[0].row_index=i;
            DownHeap(heap,0,Newk);
        }
    }
    free(distances);
    free(validCols);
    return heap;
}


Node* gower_distance(double* matrix,int targetRow, int numRows, int numCols,int k, int* neigborsFound,int* colType, double* colRange)
{
    double* distances = calloc(numRows,sizeof(double));
    if(distances==NULL) return NULL;
    int* validCols= calloc(numRows,sizeof(int));
    if(validCols==NULL)
    {
        free(distances);
        return NULL;
    }
    for(int i=0;i<numCols;i++)
    {
        size_t offset=(size_t)i*numRows;
        double targetVal=matrix[targetRow+offset];
        if(ISNA(targetVal)) continue;

        if(colType[i])
        {
            
            for(int j=0;j<targetRow;j++)
            {
                double curr= matrix[j+offset];
                if(ISNA(curr)) continue;
                distances[j]+=(targetVal==curr)? 0.0:1.0;
                validCols[j]++;
            }
            for(int j=targetRow+1;j<numRows;j++)
            {
                double curr= matrix[j+offset];
                if(ISNA(curr)) continue;
                distances[j]+=(targetVal==curr)? 0.0:1.0;
                validCols[j]++;
            }
        }
        else
        {
            double range=colRange[i];
            for(int j=0;j<targetRow;j++)
            {
                double curr= matrix[j+offset];
                if(ISNA(curr)) continue;
                distances[j]+=fabs(targetVal-curr)/range;
                validCols[j]++;
            }
            for(int j=targetRow+1;j<numRows;j++)
            {
                double curr= matrix[j+offset];
                if(ISNA(curr)) continue;
                distances[j]+=fabs(targetVal-curr)/range;
                validCols[j]++;
            }
        }

    }
    for(int i=0;i<numRows;i++) // skalujemy dystanse
    {
        if(i==targetRow || validCols[i]==0) continue;
        double scalingFactor = (double) numCols / (double) validCols[i];
        distances[i]= distances[i]*scalingFactor;
    }
    int Newk=MIN(k,numRows-1);
    Node* heap=malloc(Newk*sizeof(Node));
    if(heap==NULL)
    {
        free(distances);
        free(validCols);
        return NULL;
    }
    int rowIndex=0;
    int heapCount=0;
    while(rowIndex<numRows && heapCount<Newk)
    {
        if(rowIndex==targetRow || validCols[rowIndex]==0)
        {
            rowIndex++;
            continue;
        }
        heap[heapCount].distance=distances[rowIndex];
        heap[heapCount].row_index=rowIndex;   
        heapCount++;
        rowIndex++;
    }
    if(heapCount<Newk) Newk=heapCount;
    *neigborsFound = Newk;
    if(Newk == 0)
    {
        free(distances);
        free(validCols);
        free(heap);
        return NULL;
    }


    build_heap(heap,Newk);

    for(int i=rowIndex;i<numRows;i++)
    {
        if(i==targetRow || validCols[i]==0) continue;
        if(distances[i]<heap[0].distance)
        {
            heap[0].distance=distances[i];
            heap[0].row_index=i;
            DownHeap(heap,0,Newk);
        }
    }
    free(distances);
    free(validCols);
    return heap;

}

SEXP knn_imputeC(SEXP rMatrix, SEXP rK, SEXP rMetricFlag,SEXP rThreads,SEXP rColType,SEXP rColRange)
{
    if(!(Rf_isReal(rMatrix)) || !Rf_isMatrix(rMatrix)) Rf_error("Matrix is not numeric");
    if(!(Rf_isInteger(rK))) Rf_error("K is not numeric");
    if(Rf_asInteger(rK)<=0) Rf_error("K must be positive");
    if(!(Rf_isInteger(rThreads))) Rf_error("rThreads is not numeric");
    if(Rf_asInteger(rThreads)<=0) Rf_error("rThreads must be positive");

    if(XLENGTH(rK)!=1) Rf_error("k must be scalar");

    if(XLENGTH(rMatrix)==0) Rf_error("There is nothing to impute");


    SEXP outputMatrix = Rf_duplicate(rMatrix);
    PROTECT(outputMatrix);
    double* matrixOut= REAL(outputMatrix);
    double* matrixIn=REAL(rMatrix);
    int numRows= Rf_nrows(outputMatrix);
    int numCols= Rf_ncols(outputMatrix);

    int k= Rf_asInteger(rK);
    int metricFlag= Rf_asInteger(rMetricFlag);
    int numThreads= Rf_asInteger(rThreads);


    int* colType= (rColType != R_NilValue)? INTEGER(rColType) : NULL;
    double* colRange= (rColRange != R_NilValue)? REAL(rColRange) : NULL;
    // tworzymy mape wierszy w których występują braki do uzupełnienia
    int* rowsToFIx = calloc(numRows,sizeof(int));
    if(rowsToFIx==NULL)
    {
        UNPROTECT(1);
        Rf_error("Allocation error");
    }

    for(int j=0; j<numCols;j++)
    {
        for(int i=0;i<numRows;i++)
        {
            if(ISNA(matrixIn[i+j*numRows])) rowsToFIx[i]=1;
        }
    }

    // Znajdujemy k sąsiadów dla wierszy które wymagają uzupełnienia 
    #pragma omp parallel for num_threads(numThreads)
    for(int target=0; target<numRows;target++)
    {
        if(rowsToFIx[target]==0) continue;

        Node* neighbors= NULL;
        int neighborsFound=0;
        switch (metricFlag)
        {
        case 1:
            neighbors=euclid_distance(matrixIn,target,numRows,numCols,k,&neighborsFound);
            break;
        
        case 2:
            neighbors=manhattan_distance(matrixIn,target,numRows,numCols,k,&neighborsFound);
            break;
        case 3:
            neighbors=gower_distance(matrixIn,target,numRows,numCols,k,&neighborsFound,colType,colRange);
            break;
        default:
        
            break;
        }
        if(neighbors==NULL || neighborsFound==0) continue;

        for(int j = 0; j < numCols; j++)
        {
            size_t offset= j*numRows;
            int targetIndex = target + offset;

            if(!ISNA(matrixIn[targetIndex])) continue;

            int validNeighbors = 0;

            if(colType != NULL && colType[j] == 1)
            {
                double* uniqueVal = malloc(neighborsFound*sizeof(double));
                if(uniqueVal==NULL)
                {
                    free(neighbors);
                    free(rowsToFIx);
                    UNPROTECT(1);
                    Rf_error("Allocation error");
                }
                int* counts = calloc(neighborsFound,sizeof(int));
                if(counts==NULL)
                {
                    free(uniqueVal);
                    free(neighbors);
                    free(rowsToFIx);
                    UNPROTECT(1);
                    Rf_error("Allocation error");
                }
                int uniqueCount=0;
                for(int n=0;n<neighborsFound;n++)
                {
                    double neighborVal= matrixIn[neighbors[n].row_index+offset];
                    if(ISNA(neighborVal)) continue;

                    validNeighbors++;
                    int found=0;

                    for(int u=0;u<uniqueCount;u++)
                    {
                        if(neighborVal==uniqueVal[u])
                        {
                            counts[u]++;
                            found=1;
                            break;
                        }
                    }
                    if(!found)
                    {
                        uniqueVal[uniqueCount]=neighborVal;
                        counts[uniqueCount]=1;
                        uniqueCount++;
                    }
                }
                if(validNeighbors>0)
                {
                    int maxVotes=-1;
                    double winner=0.0;

                    for(int u=0;u<uniqueCount;u++)
                    {
                        if(counts[u]>maxVotes)
                        {
                            maxVotes=counts[u];
                            winner=uniqueVal[u];
                        }
                    }
                    matrixOut[targetIndex]=winner;

                }
                free(uniqueVal);
                free(counts);

            }
            else 
            {
                double sum=0.0;
                for(int n=0;n<neighborsFound;n++) 
                {
                    double neighborVal= matrixIn[neighbors[n].row_index+offset];

                    if(!ISNA(neighborVal))
                    {
                    sum+=neighborVal;
                    validNeighbors++;
                    }
                }
                if(validNeighbors>0)
                {
                    matrixOut[targetIndex]= sum/(double)validNeighbors;
                }
            }
            }
            free(neighbors);
        }
    
    
    free(rowsToFIx);
    UNPROTECT(1);
    return outputMatrix;
}