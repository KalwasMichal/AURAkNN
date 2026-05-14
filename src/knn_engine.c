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
        double targetVal=matrix[targetRow+i*numRows];
        if(ISNA(targetVal)) continue;
        
        for(int j=0;j<numRows;j++)
        {
            double curr= matrix[j+i*numRows];
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
        double targetVal=matrix[targetRow+i*numRows];
        if(ISNA(targetVal)) continue;
        
        for(int j=0;j<numRows;j++)
        {
            double curr= matrix[j+i*numRows];
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