#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <math.h>
#define PNum 20 

using namespace std;

char *fichier, *fileOut;

int f, f_best, f_globe, f_local;
int *color, *best_color, *globe_color, *local_color;
int **gamma_;
int **edge;
int *adjactLen;
int **adjact_matrix;
int **tabuTenure;
int max_vtx, max_col, max_iter;
long long seed = 770;
int nbcfl;
double alpha = 0.6;

int *conflictVtx, *conflictIndex, conflictLen;

int *best_x, *best_v;
int *tabu_best_x, *tabu_best_v;

long long onceIter;

double total_iters;

/*-----------------------------------------------------------------------*/
int **freeColors, **freeColorsIndex, *freeColorsLen, **freeColorsFlag;
int *vertex, *vertexIndex, vertexLen;
int *vertexTmp, *vertexIndexTmp, vertexLenTmp;
bool *coloredVertex;
int *colorNum;
bool *colFlag;

/*-----------------------------------------------------------------------*/
int *onlyVtx;
int onlyVtxLen;
int *bestConfVtx, bestConfVtxLen;

int *unColVtx, *unColVtxEva, unColVtxLen;
long long *freqArray;

/*-----------------------------------------------------------------------*/
int **PSet, *PFbest;
int **distMatrix;
double *fitnessArray;
int **minDistance;
int *distanceTmp;
int minDist, maxDist;
int minF, maxF;

int **Concurrency, **dt, **fs ;
int *mc, *cc ; 
int cccnt;

void read()
{
     int x, y, i;
     
     ifstream FIC;
     string strLine;
     FIC.open(fichier);
     
     if ( FIC.fail() )
     {
           cout << "### Erreur ouverture, fichier " << fichier << endl;
           exit(0);
     }

     char StrReading[100];
     FIC >> StrReading;
     if ( FIC.eof() )
     {
           cout << "### Erreur ouverture, fichier " << fichier << endl;
           exit(0);
     }

     int nb_vtx=0, nb_edg=-1, max_edg=0;
     while ( ! FIC.eof() )
     {
           char bidon[50];
           if( strcmp(StrReading, "p" )==0 )
           {
                 FIC >> bidon >> max_vtx >> nb_edg;
                 cout << "Nombre de sommets = " << max_vtx << endl;
                 cout << "Nombre d'arretes = " << nb_edg << endl;
               
                 color=new int[max_vtx];
                 best_color=new int[max_vtx];
                 globe_color=new int[max_vtx];
                 local_color=new int[max_vtx];
                 edge=new int*[max_vtx];
                 adjact_matrix = new int *[max_vtx];
                 adjactLen = new int[max_vtx];
                 gamma_=new int*[max_vtx];   
                 tabuTenure=new int*[max_vtx];                               

                 conflictVtx=new int[max_vtx];
                 conflictIndex=new int[max_vtx];                        

                 best_x = new int[max_col * max_vtx]; 
                 best_v = new int[max_col * max_vtx];  
                 tabu_best_x = new int[max_col * max_vtx];  
                 tabu_best_v = new int[max_col * max_vtx];
                 
                 freeColors = new int *[max_vtx];
                 freeColorsIndex = new int *[max_vtx];
                 freeColorsLen = new int[max_vtx];
                 freeColorsFlag = new int *[max_vtx];               
                 vertex = new int[max_vtx];
                 vertexIndex = new int[max_vtx];
                 vertexTmp = new int[max_vtx];
                 vertexIndexTmp = new int[max_vtx];
                 onlyVtx = new int [max_vtx];
                 coloredVertex = new bool[max_vtx];
                 colorNum = new int [max_col];
                 colFlag = new bool [max_vtx];
                 
                 bestConfVtx = new int [max_vtx];
                 
                 unColVtx = new int [max_vtx];
                 unColVtxEva = new int [max_vtx];
                 
                 freqArray = new long long [max_vtx];
                 
                 PSet = new int *[PNum];
                 PFbest = new int[PNum];   
                 distMatrix = new int *[PNum]; 
                 fitnessArray = new double [PNum]; 
                 minDistance = new int *[max_vtx]; 
                 distanceTmp = new int[PNum]; 
                     
                 Concurrency = new int *[max_col];
                 dt = new int *[max_col];
                 fs = new int *[max_col];
                 mc = new int [max_col] ;
                 cc = new int [max_col] ; 				 
				  
                 for(x=0;x<max_vtx;x++) 
                 {
                           edge[x]=new int[max_vtx];
                           adjact_matrix[x] = new int [max_vtx];
                           gamma_[x]=new int[max_col];
                           tabuTenure[x]=new int[max_col];
                           freeColors[x] = new int [max_col];
                           freeColorsIndex[x] = new int [max_col];
                           freeColorsFlag[x] = new int [max_col]; 
						   minDistance[x] = new int [2];                      
                 }
                 
                 for(x=0; x<max_col; x++)
                 {
                           colorNum[x] = max_col;                                            
                 }  
               
                 for(x=0;x<max_vtx;x++)
                 {
                           color[x]=0;
                           adjactLen[x] = 0;
                           freeColorsLen[x] = 0;
                           coloredVertex[x] = false;
                           colFlag[x] = false;
                           freqArray[x] = 0;
                           for (y=0;y<max_vtx;y++)
                           {
                                edge[x][y]=0;
                                adjact_matrix[x][y] = 0;
                           }                        
                  }
                  
	               for(x=0; x<PNum; x++)
	               {
	                       PSet[x] = new int [max_vtx];
	                       distMatrix[x] = new int [PNum]; 
	                       PFbest[x] = 0;                
	               } 
	               
	               for(x=0; x<max_col; x++) 
	               {
	                        Concurrency[x] = new int[max_col];
	                        dt[x] = new int[max_col];
	                        fs[x] = new int[max_col];                                
	               }                         

           }

           int x1, x2;

           if ( strcmp(StrReading, "e")==0 )
           {
                 FIC >> x1 >> x2;
                 x1--; x2--;
                 if ( x1<0 || x2<0 || x1>=max_vtx || x2 >=max_vtx || x1 == x2 )
                 {
                       cout << "### Erreur de lecture : x1="<< x1 << ", x2=" << x2 << endl;
                       exit(0);
                 }

                 edge[x1][x2]=edge[x2][x1]=1;
                 max_edg++;
           }
           
            if ( strcmp(StrReading, "f")==0 )
            {
                       FIC >> x;
                       x--;
                       getline(FIC, strLine);

                       istringstream iss(strLine);
                       while(iss >> y)
                       {
                                 y--;
                                 freeColors[x][freeColorsLen[x]] = y;
                                 freeColorsIndex[x][y] = freeColorsLen[x];
                                 freeColorsLen[x]++;
                       }             
            }           
           
           FIC >> StrReading;
     }

     cout << "max_edge="<<max_edg<<endl;
     cout << "Denstite = " << (float) max_edg/max_vtx << endl;

     for(x=0; x<max_vtx; x++)
     {
             for(y=0; y<max_vtx; y++)
                     if(edge[x][y] != 0) 
                     {        
                                   adjact_matrix[x][adjactLen[x]++] = y; 
                     }
     } 
     
     vertexLen = 0;
     vertexLenTmp = 0;

     for(x=0; x<max_vtx; x++)
     {
            if(freeColorsLen[x] != 0 && freeColorsLen[x] != 1 && freeColorsLen[x] != max_col)
            {
                     cout<<"Error of the freeColorsLen!!! exit!"<<endl;
                     system("pause");
                     exit(-1);
            }
     }

     FIC.close();
     
     cout <<"Density = " << (float) 2 * max_edg/(max_vtx*(max_vtx-1)) << endl  ;
     cout << "File closed!" << endl ;
     cout << "The initializing is finished!" << endl; 
}

void q_sort_Descend(long long numbers[ ], int left, int right, int index[ ])
{
      int pivot, pivot_index, l_hold, r_hold;
    
      l_hold = left;
      r_hold = right;
      pivot = numbers[left];
      pivot_index = index[left] ;
      while (left < right)
      {
            while ((numbers[right] <= pivot) && (left < right))
              right--;
            if (left != right)
            {
              numbers[left] = numbers[right];
              index[left] = index[right];
              left++;
            }
            while ((numbers[left] >= pivot) && (left < right))
              left++;
            if (left != right)
            {
              numbers[right] = numbers[left];
              index[right] = index[left];
              right--;
            }
      }
      numbers[left] = pivot;
      index[left] = pivot_index;
      pivot = left;
      left = l_hold;
      right = r_hold;
      if (left < pivot)
        q_sort_Descend(numbers, left, pivot-1, index);
      if (right > pivot)
        q_sort_Descend(numbers, pivot+1, right, index);
}

void Quick_Sort_Descend(long long numbers[], int array_size, int index[])
{
      q_sort_Descend(numbers, 0, array_size - 1, index);
}

void q_sort_Ascend(int numbers[ ], int left, int right, int index[ ])
{
      int pivot, pivot_index, l_hold, r_hold;
    
      l_hold = left;
      r_hold = right;
      pivot = numbers[left];
      pivot_index = index[left] ;
      while (left < right)
      {
            while ((numbers[right] >= pivot) && (left < right))
              right--;
            if (left != right)
            {
              numbers[left] = numbers[right];
              index[left] = index[right];
              left++;
            }
            while ((numbers[left] <= pivot) && (left < right))
              left++;
            if (left != right)
            {
              numbers[right] = numbers[left];
              index[right] = index[left];
              right--;
            }
      }
      numbers[left] = pivot;
      index[left] = pivot_index;
      pivot = left;
      left = l_hold;
      right = r_hold;
      if (left < pivot)
        q_sort_Ascend(numbers, left, pivot-1, index);
      if (right > pivot)
        q_sort_Ascend(numbers, pivot+1, right, index);
}

void Quick_Sort_Ascend(int numbers[], int array_size, int index[])
{
      q_sort_Ascend(numbers, 0, array_size - 1, index);
}

void preReduction()
{
      int x, y, i, j;
      int adjV;
      int idx, col, colFin, val;
      int row, line;
      
      for(x=0; x<max_vtx; x++)
      {
          for(y=0; y<max_col; y++)
          {
                 freeColorsFlag[x][y] = false;
          }
      }
      
      vertexLen = 0;
      
      for(x=0; x<max_vtx; x++)
      {
               if(freeColorsLen[x] == 1)
               {
                     coloredVertex[x] = true;
                     colFlag[x] = true;
               }
      }    
      
      onlyVtxLen = 0;
      for(x=0; x<max_vtx; x++)
      {
               if(coloredVertex[x] == true)
               {
                     row = x / max_col;
                     line = x % max_col;
                     
                     col = freeColors[x][0];
                     best_color[x] = color[x] = col;
                     globe_color[x] = local_color[x] = col;
                     colorNum[col]--;
                     
                     for(y=0; y<adjactLen[x]; y++)
                     {
                             adjV = adjact_matrix[x][y];
                             
                             if(coloredVertex[adjV] == false)
                             {                                     
                                     if(freeColorsFlag[adjV][col] == false)
                                     {
                                                    freeColorsFlag[adjV][col] = true;
                                                    freeColorsLen[adjV]--;
                                                    idx = freeColorsIndex[adjV][col];
                                                    if(idx < 0 || idx > freeColorsLen[adjV])
                                                    {
                                                            cout<<"Something wrong!!!!!!!! idx = freeColorsIndex[adjV][col]"<<endl;
                                                            system("pause");
                                                            exit(-1);
                                                    }                                            
                                                    colFin = freeColors[adjV][freeColorsLen[adjV]];
                                                    freeColors[adjV][idx] = colFin;
                                                    freeColorsIndex[adjV][colFin] = idx;
                                                    
                                                    if(freeColorsLen[adjV] == 1)
                                                    {
                                                           onlyVtx[onlyVtxLen] = adjV;
                                                           onlyVtxLen++;                       
                                                    }                                            
                                     }
                             }
                     }
               }
               else
               {
                     vertex[vertexLen] = x;
                     vertexIndex[x] = vertexLen;
                     vertexLen++;
               }
               
      }

      if(vertexLen == max_vtx)
      {
                for(x=0; x<max_vtx; x++)
                {
                        for(y=0; y<max_col; y++)
                        {
                                 freeColors[x][freeColorsLen[x]] = y;
                                 freeColorsIndex[x][y] = freeColorsLen[x];
                                 freeColorsLen[x]++;
                        } 
                } 
                
      }
      
      while(onlyVtxLen > 0)
      {
                 x = onlyVtx[0];
                 coloredVertex[x] = true;
                 colFlag[x] = true;
                 
                 vertexLen--;
                 idx = vertexIndex[x];
                 if(idx < 0 || idx > vertexLen)
                 {
                        cout<<"Something wrong!!!!!!!! idx = vertexIndex[x]"<<endl;
                        system("pause");
                        exit(-1);
                 }                    
                 val = vertex[vertexLen];
                 vertex[idx] = val;
                 vertexIndex[val] = idx;                 

                 if(freeColorsLen[x] != 1)
                 {
                        cout<<"Something wrong!!!!!!!! freeColorsLen[x] should be 1"<<endl;
                        system("pause");
                        exit(-1);                                     
                 }
                 
                 row = x / max_col;
                 line = x % max_col;
                 
                 col = freeColors[x][0];
                 best_color[x] = color[x] = col;
                 globe_color[x] = local_color[x] = col;
                 colorNum[col]--;                           
                 
                 for(y=0; y<adjactLen[x]; y++)
                 {
                         adjV = adjact_matrix[x][y];                                                                              
                         
                         if(coloredVertex[adjV] == false)
                         {                                 
                                 if(freeColorsFlag[adjV][col] == false)
                                 {
                                                freeColorsFlag[adjV][col] = true;
                                                freeColorsLen[adjV]--;
                                                idx = freeColorsIndex[adjV][col];
                                                if(idx < 0 || idx > freeColorsLen[adjV])
                                                {
                                                            cout<<"Something wrong!!!!!!!! idx = freeColorsIndex[adjV][col]"<<endl;
                                                            system("pause");
                                                            exit(-1);
                                                }                                                  
                                                colFin = freeColors[adjV][freeColorsLen[adjV]];
                                                freeColors[adjV][idx] = colFin;
                                                freeColorsIndex[adjV][colFin] = idx;
                                                
                                                if(freeColorsLen[adjV] == 1)
                                                {
                                                       onlyVtx[onlyVtxLen] = adjV;
                                                       onlyVtxLen++;                       
                                                }
                                 }
                         }
                 }                  
                 
                 onlyVtxLen--;
                 val = onlyVtx[onlyVtxLen];
                 onlyVtx[0] = val;               
                       
      }
      
      cout<<"colorNum[x]:"<<endl; 
      for(x=0; x<max_col; x++)
           cout<<colorNum[x]<<" ";
      cout<<endl;            
      
      cout<<"vertexLen = "<<vertexLen<<endl;
         
      for(x=0; x<max_vtx; x++) 
      {
               adjactLen[x] = 0;
               vertexIndexTmp[x] = vertexIndex[x];
      }

      vertexLenTmp = vertexLen;
      for(i=0; i<vertexLen; i++)
      {  
         x = vertex[i];
         vertexTmp[i] = x;
         for(j=0; j<vertexLen; j++) 
         {    
                     y = vertex[j];
                     if(edge[x][y] != 0) 
                     {        
                                   adjact_matrix[x][adjactLen[x]] = y; 
                                   adjactLen[x]++;
                     }
         }
      }

}


void clearGamma()
{
     int x, v;
     
     for (x=0;x<max_vtx;x++)
           for (v=0;v<max_col;v++)
           {
                 gamma_[x][v]=0;
                 tabuTenure[x][v]=0;
           }               
     f=0;            
     conflictLen = 0;
}

void reduceConflictX(int x)
{
       int id, xValue;
       
       if(gamma_[x][color[x]] > 0)
       {
                conflictLen--;
                id = conflictIndex[x];
                xValue = conflictVtx[conflictLen];  
                conflictVtx[id] = xValue;
                conflictIndex[xValue] = id;                        
       } 
          
}

void addConflictX(int x)
{
       if(gamma_[x][color[x]] > 0)
       {
               conflictVtx[conflictLen] = x;
               conflictIndex[x] = conflictLen;
               conflictLen++;                           
       }      
}


void buildGamma(int *array, int arrayLen)
{
     int i, j, x, y;     
     clearGamma();

     for (i=0;i<arrayLen;i++)
     {
           x = array[i];
           for(j=i+1; j<arrayLen; j++)
           {
                 y = array[j];
                 
                 if(edge[x][y] == 1)
                 {
                       gamma_[x][color[y]] += 1;
                       gamma_[y][color[x]] += 1;
                       if(color[x] == color[y])
                                f += 1;        
                 }
           }
           if(gamma_[x][color[x]] > 0)
                  addConflictX(x);           
     } 	   
}

int delta(int x, int v)
{ 
    return gamma_[x][v]-gamma_[x][color[x]]; 
}


void updateGammaWhenMove(int x, int v0, int v1, int *array, int arrayLen)
{
     int j, y;
     
     for(j=0; j<arrayLen; j++)
     {
          y = array[j];
          
          if(edge[x][y] == 1 && y != x)
          {
                  reduceConflictX(y);
                  
                  gamma_[y][v0] -= 1;
                  gamma_[y][v1] += 1;
                  
                  addConflictX(y);          
          }         
     }
}


int isVertexInConflict(int x)
{  
    return gamma_[x][color[x]]; 
}

int isMoveTabu(int x, int v, int iter)
{  
    return (tabuTenure[x][v]>iter); 
}

void setTabuTenure(int x, int v, int iter)
{
     int res;
     
     res=rand() % 10 + (int)(alpha*nbcfl);
     tabuTenure[x][v]=iter+res;
}

int tabuSearch(int max_iters, int *array, int arrayLen)
{        
     int ix, x, j, v, iv;    
     int num_best, num_tabu_best;   
     int trouve=1;    
     int best_delta, tabu_best_delta, delt;     
     int v_x, old_value;
     int select_index;
     int select_x, select_v; 
     bool selFlag = false;
     int col;
     int lastId, lastX;
     int distNum;
	 
     buildGamma(array, arrayLen);
     
     f_best = 999999;
     
     if(f < f_best)
     {
             f_best = f;
             for(ix=0; ix<arrayLen; ix++)
             {
                     x = array[ix];
                     best_color[x] = color[x]; 
             }
             bestConfVtxLen = conflictLen;
             for(ix=0; ix<conflictLen; ix++)
             {
                   bestConfVtx[ix] = conflictVtx[ix]; 
				   freqArray[conflictVtx[ix]]++;
			}
     } 
     
     if(f_best == 0)
             return f_best;    
     

     for(int iters=1; iters<max_iters && f>0 ; iters++)
     {
           trouve = 1;
           tabu_best_delta = 9999999;
           best_delta = 9999999;
           num_tabu_best = 0;
           num_best = 0;
           nbcfl = 0;
           selFlag = false;

           for(ix=0; ix<conflictLen; ix++)
           {
                       x = conflictVtx[ix];
                       nbcfl++;
                       v_x = color[x];
                       for(j=0; j<freeColorsLen[x]; j++)
                       {
                            v = freeColors[x][j];
                            if(v != v_x)
                            {
                                  delt = delta(x,v);
                                  if (!isMoveTabu(x,v,iters))
                                  {
                                        if(delt < best_delta)
                                        {
                                               num_best = 0;  
                                               best_x[num_best] = x;
                                               best_v[num_best] = v;
                                               best_delta = delt;
                                               num_best++;
                                        }
                                        else if(delt == best_delta)
                                        {
                                                best_x[num_best] = x;
                                                best_v[num_best] = v;
                                                num_best++;                                                                                   
                                        }
                                  }
                                  else
                                  {
                                          if(delt < tabu_best_delta)
                                          {
                                                num_tabu_best = 0;
                                                tabu_best_x[num_tabu_best] = x;
                                                tabu_best_v[num_tabu_best] = v;
                                                tabu_best_delta = delt;
                                                num_tabu_best++;                                                   
                                          }
                                          else if(delt == tabu_best_delta)
                                          {
                                                    tabu_best_x[num_tabu_best] = x;
                                                    tabu_best_v[num_tabu_best] = v;
                                                    num_tabu_best++;                                                                                                      
                                          }                                         
                                  }
                            }
                       }
                 }
                 if(num_best == 0 && num_tabu_best == 0)
                 {
                        cout<<"num_best == 0 && num_tabu_best == 0"<<endl;
                        return f_best;
                 }

                 if( ( num_tabu_best > 0 && tabu_best_delta < best_delta && ( f + tabu_best_delta < f_best ) ) || (num_best == 0 && num_tabu_best > 0))
                 {
                        f += tabu_best_delta;
                        select_index = rand() % num_tabu_best;
                        select_x = tabu_best_x[select_index];
                        select_v = tabu_best_v[select_index]; 
                        selFlag = true;                    
                 }
                 else 
                 {
                        f += best_delta;
                        select_index = rand() % num_best;
                        
                        select_x = best_x[select_index];
                        select_v = best_v[select_index]; 
                        selFlag = true;                    
                 }

                 if(selFlag == true)
                 {
                         reduceConflictX(select_x);                
                         old_value = color[select_x];
                         updateGammaWhenMove(select_x, old_value, select_v, array, arrayLen);
                         color[select_x] = select_v;                
                         addConflictX(select_x);
        
                         setTabuTenure(select_x, old_value, iters);
        
                         if(f < f_best)
                         {
                                   f_best = f;
                                   for(ix=0; ix<arrayLen; ix++)
                                   {
                                           x = array[ix];  
                                           best_color[x] = color[x];   
                                   }
                                   if(f_best == 0)
                                           return f_best;
                                   bestConfVtxLen = conflictLen;
                                   for(ix=0; ix<conflictLen; ix++)
                                   {
                                           bestConfVtx[ix] = conflictVtx[ix];
                                           freqArray[conflictVtx[ix]]++;
                                   }
                                                                  
                         }
                         else if(f == f_best)
                         {
                                   for(ix=0; ix<arrayLen; ix++)
                                   {
                                           x = array[ix];  
                                           best_color[x] = color[x];   
                                   } 
                                   bestConfVtxLen = conflictLen;
                                   for(ix=0; ix<conflictLen; ix++)
                                   {
                                           bestConfVtx[ix] = conflictVtx[ix]; 
										   freqArray[conflictVtx[ix]]++;
								   }
                         }
                         
                 }
                 onceIter++;
     }
      
     return f_best;

}

int verifyGlobeSolution()
{
    int x,y;
    int err = 0;
    for(x=0; x<max_vtx; x++)
    {
            for(y=x+1; y<max_vtx; y++)
            {
                if(x!=y && edge[x][y] != 0)                       
                { 
                        if(globe_color[x] == globe_color[y])
                                 err += 1;                                                            
                }  
            }           
    }
    return err;
}



void clearFun()
{
     for (int x=0;x<max_vtx;x++) 
               color[x]=-1;         
}

void initSol()
{
     int x;
     
     for(x=0; x<max_vtx; x++)
     {
          if(coloredVertex[x] == true)
               color[x] = globe_color[x]; 
          else
          {
               color[x] = freeColors[x][rand() % freeColorsLen[x]];
               colFlag[x] = true;
          }
         
     }    
}

void reduceVtxFunc(int bestLen)
{
     int selV, selId;
     int idx, val;

     do
     {
             selId = rand() % bestConfVtxLen;
             selV = bestConfVtx[selId];
             
             vertexLen--;
             idx = vertexIndex[selV];
             if(idx < 0 || idx > vertexLen)
             {
                    cout<<"Something wrong!!!!!!!! idx = vertexIndex[x]"<<endl;
                    system("pause");
                    exit(-1);
             }                    
             val = vertex[vertexLen];
             vertex[idx] = val;
             vertexIndex[val] = idx; 
    
			 unColVtx[unColVtxLen] = selV;
		     unColVtxEva[unColVtxLen] = freeColorsLen[selV]; 
		     unColVtxLen++;
	         
	         colFlag[selV] = false;
             bestConfVtxLen--;
             bestConfVtx[selId] = bestConfVtx[bestConfVtxLen];   
               
     }while(bestConfVtxLen > bestLen);       
     
}

void copyColor()
{
     int i, x;
   
     for(i=0; i<vertexLenTmp; i++)
     {
            x = vertexTmp[i];
            if(colFlag[x] == true)
                  color[x] = best_color[x];  
		    else
			{
				  color[x] = best_color[x];
				  colFlag[x] = true;
			}     
     }
}


void recoverVtxArray()
{
     int i, x;
     
     vertexLen = vertexLenTmp;
     unColVtxLen = 0;
     for(i=0; i<vertexLenTmp; i++)
     {
            x = vertexTmp[i];
			vertex[i] =  x;

			freqArray[x] = 0; 
			if(colFlag[x] == false)
			{
			     cout<<"Something is wrong!!!!!!colFlag[x] == false"<<endl;
			     system("pause");
			     exit(-1);
			}   
     }
     for(i=0; i<max_vtx; i++)
            vertexIndex[i] = vertexIndexTmp[i];     
}

void initConfVtx()
{
     int i, x;
     
	 for(i=0; i<vertexLen; i++)
     {
            x = vertex[i];
            if(freqArray[x] != 0)
            {
            	 color[x] = freeColors[x][rand() % freeColorsLen[x]];
	    	}
	 }
}

bool localSearch(int maxLSIters)
{
     int LSIter;
     int i, x;
     int unImprovNum;
     
     LSIter = 0; 
     unImprovNum = 0;
     f_local = 999999;
  
     tabuSearch(max_iter, vertex, vertexLen);     
     if(f_best < f_local)
     {
              f_local = f_best;
              for(i=0; i<vertexLenTmp; i++)
              {
                     x = vertexTmp[i];
                     local_color[x] = best_color[x];         
              }
              unColVtxLen = 0;
              unImprovNum = 0;
              if(f_local == 0)
              {
              	    cout<<"1 f_local == 0"<<endl;
                     return true; 
			  }
     }  
   
     bool partFlag = false;
     int removeNum = 0; 
     do
     {                                   
			copyColor(); 
			if(bestConfVtxLen % 2 == 0)
			    removeNum = bestConfVtxLen/2;
			else
			    removeNum = bestConfVtxLen/2+1;    
			
			if(partFlag == false)					              
                 reduceVtxFunc(removeNum);
 
            tabuSearch(max_iter, vertex, vertexLen); /* Partial coloring */
            if(f_best  == 0)
	   		{
	            copyColor(); 
				recoverVtxArray();
	            tabuSearch(max_iter, vertex, vertexLen); /* Completed coloring */
	            partFlag = false;             
	            if(f_best < f_local)
	            {
	                  f_local = f_best;
	                  for(i=0; i<vertexLenTmp; i++)
	                  {
	                         x = vertexTmp[i];
	                         local_color[x] = best_color[x];         
	                  }
	                  unImprovNum = 0;
	                  if(f_local == 0)
	                  {
	                  	     cout<<"1 f_local == 0"<<endl;
	                         return true; }	                                    
	            }
	            else
	                   unImprovNum++;
            }
            else
            {
                unImprovNum++;
                copyColor();
                recoverVtxArray();
                partFlag = true;
            }

            LSIter++;                                                                
     }while(LSIter < maxLSIters);  
     return false;       
}

int augCycle( int ** cost)
{
         int i, j, k ;
         for( i = 0 ; i < max_col ; i ++ )
           for( j = 0 ; j < max_col ; j ++ )
           {
                dt[ i ][ j ] = cost[ i ][ mc[ i ] ] - cost[ i ][ mc[ j ] ] ;
                fs[ i ][ j ] = j ;    
           }

          for( k = 0 ; k < max_col ; k ++ )   
             for( i = 0 ; i < max_col ; i ++ )
                for( j = 0 ; j < max_col ; j ++ )
                  if( dt[ i ][ k ] + dt[ k ][ j ] < dt[ i ][ j ] )
                  {
                      dt[ i ][ j ] = dt[ i ][ k ] + dt[ k ][ j ];
                      fs[ i ][ j ] = fs[ i ][ k ] ;
                      if( i == j )      
                      {
                            cccnt = 0 ; 
                            do
                            {
                                 cc[cccnt] = i ; cccnt ++ ;
                                 i = fs[ i ][ j ] ; 
                            }while( i != j );
                            
                            return 1;
                      }
                  }
           return 0;
} 

int Bipartite_Matching( int ** cost )
{
     int i, j ;
     int ret = 0 ; 
     
     while(augCycle(cost))
     {
        j = mc[cc[ 0 ]];
        for( i = 0 ; i < cccnt - 1 ; i ++ )
           mc[ cc[ i ] ] = mc[ cc[ i + 1 ] ] ;
        mc[ cc[ i ] ] = j;
     }
     for( i = 0 ; i < max_col; i ++ )
         ret += cost[ i ][ mc[ i ] ] ;
     return ret;
} 

int Calculating_Distance(int *Color1, int *Color2)
{
           int i, j, h, b_i, b_j ;
           int c1, c2, cost, cost2, dist ; 
           int side_count ; 
           int * flag1, * flag2 ; 
           flag1 = new int[ max_col ] ;
           flag2 = new int[ max_col ] ;
           for( i = 0 ; i < max_col ;  i++ ) 
           {
                flag1[ i ] = 0 ;
                flag2[ i ] = 0 ;
           } 
           for( i = 0 ; i < max_col ; i ++ )   
             for( j = 0 ; j < max_col ; j ++ )
                  Concurrency[ i ][ j ] = 0 ;
           for( i = 0 ; i < max_vtx ; i ++ )      
           {
               c1 = Color1[ i ] ;
               c2 = Color2[ i ] ;
               Concurrency[ c1 ][ c2 ] ++ ;
           }
           for( h = 0 ; h < max_col ; h ++ )
           {
                cost = -1 ;
                for( i = 0 ; i < max_col ; i ++ )   
                   for( j = 0 ; j < max_col ; j ++ )
                      if(flag1[i] == 0 && flag2[j] == 0 && (cost < Concurrency[ i ][ j ] || cost == Concurrency[ i ][ j ] && rand() % 3 == 0 ) ) 
                      {
                          cost = Concurrency[ i ][ j ] ;
                          b_i = i ;
                          b_j = j ;
                      }
                mc[ b_i ] = b_j ;                
                flag1[ b_i ] = 1 ;   
                flag2[ b_j ] = 1 ; 
           }
           cost = 0 ; 
           for( i = 0 ; i<max_col ; i ++ )   
                cost += Concurrency[i][mc[i]];
           cost2 = Bipartite_Matching(Concurrency);
           dist = max_vtx - cost2 ;
           
           delete []flag1;
           delete []flag2;
           return dist;
}
        

// Memetic Algorithm Initializing
bool MemeticInitialize()
{
      int i, j, x;
      int dist = 0;
      int minD, minIndex;
      int pi, pj;
      bool sameFlag = false;
      int *sola = new int [max_vtx];
      int *solb = new int [max_vtx];
      
      f_globe = 999999;
      
      //Initialize the population
      for(pi=0; pi<PNum; pi++)
      {
                sameFlag = false;
				initSol();
                f_local = 999999;
                localSearch(10);            
                
                if(pi>0)
                {
	                for(pj=0;pj<pi;pj++)
	                {
						dist = Calculating_Distance(local_color, PSet[pj]);
						if(dist == 0)
						{
							sameFlag = true;
							pi--;
							break;					 	
						}   
					}
				}
				if(sameFlag == false)
				{
	                cout<<pi<<" ";
					PFbest[pi] = f_local;
	                for(x=0; x<max_vtx; x++)
	                         PSet[pi][x] = local_color[x];
	                if(f_local < f_globe)
	                {
	                         f_globe = f_local;
	                         for(x=0; x<max_vtx; x++)
	                                globe_color[x] = local_color[x];          
	                }
	                if(f_globe == 0)
	                {
	                	    cout<<"1 f_globe == 0"<<endl;
					        return true;
					}    	                         	                         
	                clearFun();    					
				}                 
      }
      cout<<endl;
      
      cout<<"PFbest:";
      for(i=0; i<PNum; i++)
              cout<<PFbest[i]<<" ";
      cout<<endl;       
      
      //initialize the distance matrix
      minDist = max_vtx;
      maxDist = 0;
      minF = 999999; 
      maxF = -1;
      
      for(i=0; i<PNum; i++)
      {
               if(PFbest[i] < minF)
                       minF = PFbest[i];
               if(PFbest[i] > maxF)
                       maxF = PFbest[i];
               
               for(x=0; x<max_vtx; x++)
                       sola[x] = PSet[i][x];
               
               for(j=i+1; j<PNum; j++)
               {
                         for(x=0; x<max_vtx; x++)
                                solb[x] = PSet[j][x];
                         dist = Calculating_Distance(sola, solb);
						 cout<<dist<<" ";                                                                     
                         distMatrix[i][j] = dist; 
                         distMatrix[j][i] = distMatrix[i][j];
                         if(distMatrix[i][j] < minDist)
                                  minDist = distMatrix[i][j];
                         if(distMatrix[i][j] > maxDist)
                                  maxDist = distMatrix[i][j];        
               }         
      } 
      cout<<endl;
     
      for(i=0; i<PNum; i++)
      {
               minD = max_vtx + 1;
               for(j=0; j<PNum; j++)
               {
                       if(j!=i && distMatrix[i][j] < minD)
                       {
                                  minD = distMatrix[i][j];
                                  minIndex = j;        
                       }         
               }  
               minDistance[i][0] = minD;
               minDistance[i][1] = minIndex;
               
               fitnessArray[i] = alpha * (PFbest[i] - minF) / (maxF - minF + 0.0001)
                                  + (1-alpha)*(maxDist - minDistance[i][0]) / (maxDist - minDist + 0.0001);     
      }
      cout<<"minD = "<<minD<<endl;

      delete []sola;
      delete []solb;
      return false;     
}

/*Crossover Operator*/
void crossOverFunc2()
{
     int selPNum = 2;
     int i, j, k, x;
     int maxP, maxK, maxNum, numTmp;
     int **countPC = new int *[selPNum];
     bool *pFlag = new bool [selPNum];
     int *pn = new int [selPNum]; 
     bool *vtxFlag = new bool [max_vtx];     
     int *tmp = new int [PNum];
     bool *boolCol = new bool [max_col];
     int uncolNum = 0;
     
     for(i=0; i<PNum; i++)
            tmp[i] = i;
            
     for(i=0; i<selPNum; i++)
             countPC[i] = new int [max_col];
             
     for(i=0; i<selPNum; i++)
     {
            pFlag[i] = false;  
            for(j=0; j<max_col; j++)
                    countPC[i][j] = 0;
     } 
     
     int idP = rand() % PNum;
	 pn[0] = tmp[idP];
     tmp[idP] = tmp[PNum-1];
	 idP = rand() % (PNum-1);
	 pn[1] = tmp[idP];
     
     for(i=0; i<max_vtx; i++)
     {
            if(coloredVertex[i] == true)
                      vtxFlag[i] = true;
            else
            {
		              for(j=0; j<selPNum; j++)
		                        countPC[j][PSet[pn[j]][i]]++;
                      vtxFlag[i] = false;
            }
     }
     
     for(i=0; i<max_col; i++)
            boolCol[i] = false;

     for(i=0; i<max_col; i++)
     {
            maxNum = 0;
            for(j=0; j<selPNum; j++)
            {
						   for(k=0; k<max_col; k++)
                           {
                                  if(boolCol[k] == false)
                                  {
									  numTmp = countPC[j][k];
									  if(numTmp <= colorNum[k])
		                              {   
									      if(maxNum < numTmp)
		                                  {
		                                        maxNum = numTmp;
		                                        maxP = j;
		                                        maxK = k;          
		                                  }  
		                                  else if(maxNum == numTmp && rand() % 2)
		                                  {
		                                        maxP = j;
		                                        maxK = k;
		                                  } 
									   }
								  }
                           }
            } 

            for(x=0; x<max_vtx; x++)
            {
                   if(vtxFlag[x] == false && PSet[pn[maxP]][x] == maxK)
                   {
                             color[x] = maxK;
                             vtxFlag[x] = true;
                             for(j=0; j<selPNum; j++)
                                     if(j != maxP)
                                           countPC[j][PSet[pn[j]][x]]--;                             
                   }  
            }	
            countPC[maxP][maxK] = 0;
            pFlag[maxP] = true;
            pFlag[1-maxP] = false;   			       
     }
     
     uncolNum = 0;
     for(x=0; x<max_vtx; x++)
     {
            if(vtxFlag[x] == false)
            {
                      uncolNum++;
					  if(PSet[pn[0]][x] == PSet[pn[1]][x])
                      {
                      	    color[x] = PSet[pn[0]][x];
					  }                              
					  else
					        color[x] = freeColors[x][rand()%freeColorsLen[x]];  
            }         
     }

     for(i=0; i<selPNum; i++)
     {
            delete []countPC[i];
     }
            
     delete []countPC;
     delete []pFlag;
     delete []vtxFlag;
     delete []tmp;
     delete []boolCol;
     delete []pn;
}


int updatePSet()
{
     int i, j, x;
     double fit, maxFit = 0.0;
     int maxFitIndex;
     int dist;
     int minD, minIndex;
     int accept = 0;
     int *sol = new int [max_vtx];
     
     for(i=0; i<PNum; i++)
     {
             if(fitnessArray[i] > maxFit)
             {
                   maxFit = fitnessArray[i];
                   maxFitIndex = i;                   
             }         
     }
     
     minD = max_vtx + 1;
     for(i=0; i<PNum; i++)
     {
             if(i != maxFitIndex)
             {
                         for(x=0; x<max_vtx; x++)
                                sol[x] = PSet[i][x];
                         dist = Calculating_Distance(best_color, sol);
                         distanceTmp[i] = dist;
                         if(dist < minD) 
                         {
                                  minD = dist;
                                  minIndex = i;        
                         }                          
             }         
     }
     cout<<"UpdatePSet minD = "<<minD<<endl;
     
     fit = 0.0;
     if(f_best >= maxF)
            fit *= alpha * 1.0;
     else if(f_best > maxF)
            fit += alpha * ((f_best - minF) / (maxF - minF + 0.0001));
     
     if(minD <= minDist)
            fit *= (1 - alpha) * 1.0;
     else if(minD < maxDist)
            fit += (1 - alpha) * ((maxDist - minD) / (maxDist - minDist + 0.0001));

     if(fit <= maxFit && minD > 0)
     {
             accept = 1;
             
             PFbest[maxFitIndex] = f_best;
             for(x=0; x<max_vtx; x++)
             {
                      PSet[maxFitIndex][x] = best_color[x];      
             }  
             
             minF = 999999;
             maxF = -1;
             for(i=0; i<PNum; i++)
             {
                   if(PFbest[i] < minF)
                             minF = PFbest[i];
                   if(PFbest[i] > maxF)
                             maxF = PFbest[i];         
             }
             
             for(i=0; i<PNum; i++)
             {
                    if(i != maxFitIndex)
                    {
                            distMatrix[i][maxFitIndex] = distanceTmp[i];
                            distMatrix[maxFitIndex][i] = distanceTmp[i];     
                    }         
             }
             
             minDistance[maxFitIndex][0] = minD;
             minDistance[maxFitIndex][1] = minIndex; 
             
             for(i=0; i<PNum; i++)
             {
                     if(i != maxFitIndex && minDistance[i][1] == maxFitIndex)
                     {
                            minD = max_vtx + 1;
                            for(j=0; j<PNum; j++)
                            {
                                     if(j != i)
                                            if(distMatrix[i][j] < minD)
                                            {
                                                    minD = distMatrix[i][j];
                                                    minIndex = j;                    
                                            }         
                            }     
                            minDistance[i][0] = minD;
                            minDistance[i][1] = minIndex; 
                     }         
             }
             
             minDist = max_vtx + 1;
             maxDist = 0;
             for(i=0; i<PNum; i++)
             {
                     if(minDistance[i][0] > maxDist)
                              maxDist = minDistance[i][0];
                     if(minDistance[i][0] < minDist)
                              minDist = minDistance[i][0];         
             }
             
             for(i=0; i<PNum; i++)
             {
                     fitnessArray[i] = alpha * ((PFbest[i] - minF) / (maxF - minF + 0.0001))
                                       + (1-alpha)*((maxDist - minDistance[i][0]) / (maxDist - minDist + 0.0001));                 
             }
                     
     } 
     else
        cout<<"Cannot accept the update!!!!!!!!!!!"<<endl;

     delete []sol;     
     return accept;
}


bool QWPmemetic()
{
     int Iters = 1, i, x;
     int totalIters = 100;
     int accept;
     int nonNum = 0;
     
     if(MemeticInitialize()== true)
             return true;
     
     while(Iters < totalIters)
     {
             clearFun();
             crossOverFunc2();

             localSearch(100);
             
             accept = updatePSet(); 
             if(accept == 1)
                       nonNum = 0;
             else
                       nonNum++;
                       
             if(f_local < f_globe)
             {
                       f_globe = f_local;
                       for(x=0; x<max_vtx; x++)
                                globe_color[x] = local_color[x];          
             }
             
             if(f_globe == 0)
             {
             	    cout<<"1 f_globe == 0"<<endl;
                    return true;
             }

		     cout<<"nonNum:"<<nonNum<<endl;
		 
		     cout<<"PFbest:";
		     for(i=0; i<PNum; i++)
		              cout<<PFbest[i]<<" ";
		     cout<<endl;   
                           
             Iters++;
     }
     
     return false;
}


int main(int argc, char **argv)
{
     if (argc == 4)
     {
           fichier = argv[1];          
           fileOut = argv[2];
           max_col = atoi(argv[3]);
           max_iter = 100000;
           seed = 770;
     }
     else if (argc == 5) 
     {
           fichier = argv[1];
           fileOut = argv[2];
           max_col = atoi(argv[3]);
           max_iter = atoi(argv[4]);
     }
     else 
     {
           cout << endl << "Syntaxes :" << endl;
           cout << "tabucol nom_fichier k max_iter seed" << endl;
           cout << "tabucol nom_fichier k max_iter seed alpha lgtl" << endl << endl;
           exit(0);
     }

     cout << endl << "Nom du fichier = " << fichier << endl;
     cout << "Nombre de couleurs = " << max_col << endl;
     cout << "Nombre max d'iterations = " << max_iter << endl;
     cout << "seed = " << seed << endl;
     cout << "alpha = " << alpha << endl;

     read();
     srand(time(0));

     FILE *fp; 

     int suc_iters = 0;
     double ave_times = 0.0;     
     double start_time, end_time, run_time;
     


     preReduction();

     FILE *fs_p = NULL;
     unsigned int seed = 0;

     //fs_p = fopen ("/dev/urandom", "r");
     //if (NULL == fs_p)
     //{
     //              printf("Can not open /dev/urandom\n");
     //              return -1;
     //}
     //else
     //{
     //              fread(&seed, sizeof(int), 1, fs_p);  //obtain one unsigned int data
     //             fclose(fs_p);
     //}
     srand(time(0));

     //cout<<"seed="<<seed<<endl;

     fp = fopen(fileOut, "w+");
     start_time = (double)clock();
     run_time = 0.0;
	//cout << "QWQ";
     if(QWPmemetic() == true)
     {
             end_time = (double)clock();
             run_time = (end_time - start_time) / CLOCKS_PER_SEC;

             total_iters += onceIter;

             fprintf(fp, "\n Success!!! \n max_col = %d: , run_time = %lf , once iters = %ld \n ", max_col, run_time, onceIter);
             cout<<"Success!!!"<<endl;
             //for(int x=0; x<max_vtx; x++)
                   //cout<<globe_color[x]<<" ";
             //cout<<endl;

             //for(int x = 0 ; x < max_vtx ; x ++ )
                       //fprintf(fp, " %d ", globe_color[ x ]);
             int qwq = sqrt(max_vtx);
             for(int ro=1;ro<=qwq;ro++)
             {
             	for(int co=1;co<=qwq;co++)
             	{
             		int posc=qwq*(ro-1)+co;
             		cout<<globe_color[posc-1]+1<<' ';
             		fprintf(fp, "%d ", globe_color[posc-1]+1);
				}
				cout<<endl;
				fprintf(fp,"\n");
        	 }
             //fprintf(fp, " \n ");
             suc_iters++;

             if(verifyGlobeSolution() == 0)
             {
                    cout<<"################proper##############"<<endl<<endl;
                    fprintf(fp, "################proper############## \n");
             }
             else
             {
                    cout<<"&&&&&&&&&&&&&&&&not proper&&&&&&&&&&&"<<endl<<endl;
                    fprintf(fp, "&&&&&&&&&&&&&&&&not proper&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n");
                    exit(-1);
             }
     }
     else
     {
             fprintf(fp, "\n Failed!!! \n max_col = %d: , run_time = %lf , iters = %d \n ", max_col, ((double)clock() - start_time) / CLOCKS_PER_SEC, onceIter);
             cout<<endl<<"Failed !!!"<<endl;
      }
      //fclose(fp);
      ave_times += run_time;

      if(suc_iters != 0)
      {
          ave_times = ave_times / suc_iters;
          total_iters = total_iters / suc_iters;
      }
      else
      {
          ave_times = 0.0;
          total_iters = 0.0;
      }
      cout<<"ave_times = "<<ave_times<<endl;
      cout<<"total_iters = "<<total_iters<<endl;
      //fp = fopen(fileOut, "w+");
      fprintf(fp, "\n average_time = %lf \n", ave_times);
      fprintf(fp, "\n success runs = %d \n", suc_iters);
      fprintf(fp, "\n total average iters = %lf \n", total_iters);
      fclose(fp);      
      
     
      //system("pause");
      return 1;

} 
