// Updated at Nov 26, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format

#ifndef TABLE_FORMAT_CORR_H
#define TABLE_FORMAT_CORR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>

#include "table_format.h"

using namespace std;

vector<string> selcoll(string selotu);

class _Table_Format_Meta{
    
      public:
      void Calc_Corr_Meta_Matrix(const char * infilename, const char * metafilename, const char * outfilename, int metrics, int coren);
	  void Calc_Corr_Meta_Matrix_S(const char * infilename, const char * metafilename, const char * outfilename, string selmeta, int metrics, int coren);
      float Calc_Corr_Meta_P(int sam_m, int sam_n);
      float Calc_Corr_Meta_S(int sam_m, int sam_n);
       
      private:
      _Table_Format DataAbdMeta;                   
      _Table_Format DataAbd;

};

void _Table_Format_Meta::Calc_Corr_Meta_Matrix(const char * infilename, const char * metafilename, const char * outfilename, int metrics, int coren){
       
     DataAbd.Load_Table_Rev(infilename);

     DataAbdMeta.Load_Table_Rev(metafilename);
    
 
     if (DataAbd.Get_Feature_Size() != DataAbdMeta.Get_Feature_Size()) {
        cerr << "Errors: Different size of two Samples, check your input files!" <<endl;
        return ;
        }
     
     vector<string> Features_Abd = DataAbd.Get_Sample_Names();
     vector<string> Features_Meta_Abd = DataAbdMeta.Get_Sample_Names();

     int abdsize = DataAbd.Get_Sample_Size();
     int abdmetasize = DataAbdMeta.Get_Sample_Size(); 

     int * order_m = new int [abdsize * abdmetasize ];                                                                                                      
     int * order_n = new int [abdsize * abdmetasize ];                                                                                                      
     int iter = 0;                                                                                                                                                           
                                                                                                                                                                             
     for (int i = 0; i < abdsize; i ++)                                                                                                                              
        for (int j = 0; j < abdmetasize; j ++){                                                                                                                     
            order_m[iter] = i;                                                                                                                                               
            order_n[iter] = j;                                                                                                                                               
            iter ++;                                                                                                                                                         
            }                                                                                                                                                                
                                                                                                                                                                             
    ofstream outfile(outfilename, ofstream::out);                                                                                                                            
    if (!outfile){                                                                                                                                                           
       cerr << "Error: Cannot open output file: " << outfilename << endl;                                                                                                           
       return;                                                                                                                                                               
       }                                                                                                                                                                     

      float * corr_matrix = new float [abdsize * abdmetasize];                                                                                              
      memset(corr_matrix, 0, abdsize * abdmetasize * sizeof(float));                                                                                        
                                                                                                                                                   
      omp_set_num_threads(coren);    
      
      #pragma omp parallel for schedule(dynamic, 1)                                                                                                                                
      for (int i = 0; i < iter; i++){                                                                                                                                        
          int m = order_m[i];                                                                                                                                                
          int n = order_n[i];         
                                                                                                                               
          if (metrics == 0)                                                                                                                                                  
               corr_matrix[m * abdmetasize + n] = Calc_Corr_Meta_S(m, n);       
          else corr_matrix[m * abdmetasize + n] = Calc_Corr_Meta_P(m, n);                                                                                               
        }                                                                                                                                                                    
  
      outfile << "Feature";
                                                                                                                                                         
      for (int i = 0; i < abdmetasize; i++)
          outfile <<"\t" <<Features_Meta_Abd[i];       
      outfile <<endl;                                                                                                                                                        
       for (int i = 0; i < abdsize; i++){   
           outfile <<Features_Abd[i];       

          for (int j = 0; j < abdmetasize; j++){          
              outfile <<"\t"<<corr_matrix[i*abdmetasize+j];                                                                                                         
              }                                                                                                                                                              
          outfile <<endl;                                                                                                                                                    
          }                                                                                                                                                                  
             
      outfile.close();   
      outfile.clear();                                         
      
     }

float _Table_Format_Meta::Calc_Corr_Meta_P(int sam_m, int sam_n){
          
     return DataAbdMeta.Calc_Corr_P(DataAbd.Get_Abd(sam_m), DataAbdMeta.Get_Abd(sam_n));

     }
     
float _Table_Format_Meta::Calc_Corr_Meta_S(int sam_m, int sam_n){

     return DataAbdMeta.Calc_Corr_S(DataAbd.Get_Abd(sam_m), DataAbdMeta.Get_Abd(sam_n));

     }

void _Table_Format_Meta::Calc_Corr_Meta_Matrix_S(const char * infilename, const char * metafilename, const char * outfilename, string selmeta, int metrics, int coren){
	
	DataAbd.Load_Table_Rev(infilename);

	DataAbdMeta.Load_Table_Rev(metafilename);


	if (DataAbd.Get_Feature_Size() != DataAbdMeta.Get_Feature_Size()) {
		cerr << "Errors: Different size of two Samples, check your input files!" <<endl;
		return ;
	}

	vector<string> Features_Abd = DataAbd.Get_Sample_Names();
	vector<string> Features_Meta_Abd = DataAbdMeta.Get_Sample_Names();

	int abdsize = DataAbd.Get_Sample_Size();
	int abdmetasize = DataAbdMeta.Get_Sample_Size();
 
	vector<string> svec=selcoll(selmeta);
	vector<int> Pos_sel;
	for(int i = 0; i < svec.size(); i++){
        
        bool is_contain = false;
            
		for (int j=0; j< Features_Meta_Abd.size();j++){

			if (Features_Meta_Abd[j]==svec[i])
			{
				Pos_sel.push_back(j);
				is_contain = true;
			}
		
		if (!is_contain)
		   cerr << "Warning: Feature " << svec[i] << " is not contained in the Correlation file (-m)" << endl;    	
		}
	}
	int * order_m = new int [abdsize * Pos_sel.size() ];                                                                                                      
	int * order_n = new int [abdsize * Pos_sel.size() ];                                                                                                      
	int iter = 0;                                                                                                                                                           

	for (int i = 0; i < abdsize; i ++)                                                                                                                              
		for (int j = 0; j < Pos_sel.size(); j ++){                                                                                                                     
			order_m[iter] = i;                                                                                                                                               
			order_n[iter] = Pos_sel[j];                                                                                                                                               
			iter ++;                                                                                                                                                         
		}                                                                                                                                                                

		ofstream outfile(outfilename, ofstream::out);                                                                                                                            
		if (!outfile){                                                                                                                                                           
			cerr << "Error: Error: Cannot open output file: " << outfilename << endl;                                                                                                           
			return;                                                                                                                                                               
		}                                                                                                                                                                     

		float * corr_matrix = new float [abdsize * Pos_sel.size()];                                                                                              
		memset(corr_matrix, 0, abdsize * Pos_sel.size() * sizeof(float));                                                                                        

		omp_set_num_threads(coren);    

        #pragma omp parallel for schedule(dynamic, 1)                                                                                                                                
		for (int i = 0; i < iter; i++){                                                                                                                                        
			int m = order_m[i];                                                                                                                                                
			int n = order_n[i];         

			if (metrics == 0)                                                                                                                                                  
				corr_matrix[m * Pos_sel.size() + n] = Calc_Corr_Meta_S(m, n);       
			else corr_matrix[m * Pos_sel.size() + n] = Calc_Corr_Meta_P(m, n);                                                                                               
		}                                                                                                                                                                    
        
        outfile << "Feature";
                                                                                                                                                    
		for (int i = 0; i < Pos_sel.size(); i++)
			outfile <<"\t" <<Features_Meta_Abd[Pos_sel[i]];       
		outfile <<endl;                                                                                                                                                        
		for (int i = 0; i < abdsize; i++){   
			outfile <<Features_Abd[i];       

			for (int j = 0; j < Pos_sel.size(); j++){          
				outfile <<"\t"<<corr_matrix[i*Pos_sel.size()+Pos_sel[j]];                                                                                                         
			}                                                                                                                                                              
			outfile <<endl;                                                                                                                                                    
		}                                                                                                                                                                  

		outfile.close();   
		outfile.clear();                                         
}

vector<string> selcoll(string selotu)

{
	vector<string> svec;

	stringstream str(selotu);

	string sub_str;

	while(getline(str,sub_str,','))

	{

		svec.push_back(sub_str);

	}

	return svec;

}
#endif
