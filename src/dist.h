//Updated at Mar 11, 2020
//Updated by X. Su
//Added Bray-Curits dist

#ifndef DIST_H
#define DIST_H

#include <iostream>
#include <math.h>

using namespace std;

float Calc_Dist_Cos(float * abd_m, float * abd_n, int dim){
      
      float f = 0;
      
      float f_m_sum = 0;
      float f_n_sum = 0;
            
      for (int i = 0; i < dim; i ++){
          
          float fm = abd_m[i];
          float fn = abd_n[i];
          
          f += fm * fn;
          
          f_m_sum += fm * fm;
          f_n_sum += fn * fn;
                   
          }
      
      if (f_m_sum * f_n_sum == 0) return 1;
      
      f = f / (sqrt(f_m_sum * f_n_sum));
      
      return (1 - f);
      }

float Calc_Dist_E(float * abd_m, float * abd_n, int dim){
      
      float abd_m_norm[dim];
      float abd_n_norm[dim];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < dim; i ++){
          
          abd_m_norm[i] = abd_m[i];
          abd_n_norm[i] = abd_n[i];
          
          sum_m += abd_m_norm[i];
          sum_n += abd_n_norm[i];          
          }
      
      if (sum_m <= 0) return 1;
      if (sum_n <= 0) return 1;
      
      for (int i = 0; i < dim; i ++){
          abd_m_norm[i] /= sum_m;
          abd_n_norm[i] /= sum_n;
          }
      
      //calc
      float f = 0;
      
      for (int i = 0; i < dim; i ++)
          f += pow(abd_m_norm[i] - abd_n_norm[i], 2);
      
      return sqrt(f);
      }

float Calc_Dist_JSD(float * abd_m, float * abd_n, int dim){
      
      float abd_m_norm[dim];
      float abd_n_norm[dim];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < dim; i ++){
          
          abd_m_norm[i] = abd_m[i];
          abd_n_norm[i] = abd_n[i];
          
          sum_m += abd_m_norm[i];
          sum_n += abd_n_norm[i];          
          }
      
      for (int i = 0; i < dim; i ++){
          abd_m_norm[i] /= sum_m;
          abd_n_norm[i] /= sum_n;
          }
      
      //calc
      float dkl_m = 0;
      float dkl_n = 0;
      
      for (int i = 0; i < dim; i ++){
          
          if ((abd_m_norm[i] == 0) && (abd_n_norm[i] == 0)) continue;
          
          float abd_q = (abd_m_norm[i] +  abd_n_norm[i]) / 2;
          //debug
          //cout << abd_m_norm[i] << "\t" << abd_n_norm[i] << "\t" << abd_q << endl;
          
          if (abd_m_norm[i] > 0)
             dkl_m += abd_m_norm[i] * log(abd_m_norm[i] / abd_q);
          
          if (abd_n_norm[i] > 0)
             dkl_n += abd_n_norm[i] * log(abd_n_norm[i] / abd_q);
          
          //debug
          //cout  << dkl_m << "\t" << dkl_n << endl;
          }
          
      return sqrt((dkl_m + dkl_n)/2.0);
      }

float Calc_Dist_Bray_Curtis(float * abd_m, float * abd_n, int dim){
      
      float abd_m_norm[dim];
      float abd_n_norm[dim];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < dim; i ++){
          
          abd_m_norm[i] = abd_m[i];
          abd_n_norm[i] = abd_n[i];
          
          sum_m += abd_m_norm[i];
          sum_n += abd_n_norm[i];          
          }
      
      if (sum_m <= 0) return 1;
      if (sum_n <= 0) return 1;
      
      for (int i = 0; i < dim; i ++){
          abd_m_norm[i] /= sum_m;
          abd_n_norm[i] /= sum_n;
          }

      float sum = 0;
      float diff = 0;
            
      for (int i = 0; i < dim; i ++){
          sum += (abd_m_norm[i] + abd_n_norm[i]);
          float a_diff = abd_m_norm[i] - abd_n_norm[i];
          if (a_diff < 0) a_diff = a_diff * (-1.0);
          diff += a_diff;
          }
      if (sum <= 0) return 1;
      return diff / sum;
      }

#endif
