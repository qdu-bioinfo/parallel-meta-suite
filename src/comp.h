// Updated at April 2, 2024
// Updated by Guosen Hou
// version 3.7 - 3.7.2 
// Change REG_SIZE to 120, To fit the Greengenes2 database of phylogenetic trees with more nodes

// Updated at Dec 23, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
// version 3.4.1 or above with OTU_Parser
// Added unweighted based on 3.4.1

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "hash.h"
#include "utility.h"
#include "version.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"
#include "dist.h"

#ifndef COMP_H
#define COMP_H

#define REG_SIZE 120//updated by Guosen Hou

#define MIN_DIST 0.00001

using namespace std;

class _Comp_Tree{
      
      public:
             _Comp_Tree(){                           
                           LeafN = 0;
                           OrderN = 0;                 
                           Init();
                          }
    
            _Comp_Tree(char db){
                        Database.Set_DB(db);                          
                        LeafN = 0;
                        OrderN = 0;
                        Init();
                        }
    
             int Load_abd(const char * infilename, float * Abd, bool is_cp_correct);
             int Load_abd(const char * infilename, float * Abd);
             int Load_abd(_Table_Format * table, float * Abd, int sample, bool is_cp_correct); //Load by table_format
             int Load_abd(_Table_Format * table, float * Abd, int sample); //Load by table_format

             float Calc_sim(float * Abd_1, float * Abd_2, bool w);
             float Calc_sim(float * Abd_1, float * Abd_2);
             float Calc_sim_unweight(float * Abd_1, float * Abd_2);
             float Calc_sim(float * Abd_1, float * Abd_2, int mode); //0: MS; 1: Unweighted MS; 2: Cos; 3: Eu; 4: JSD 
             
             int Get_LeafN(){
                 return LeafN;
                 }
                 
             string Get_Id(int n){
                    if ((n < 0) || (n >= LeafN)) return "NULL";
                    return Id[n];
                    }
             
      private:
              _PMDB Database;
    
              vector <float> Dist_1;
              vector <float> Dist_2;
              
              vector <int> Order_1;
              vector <int> Order_2;
              vector <int> Order_d;
              
              vector <string> Id; 
              
              hash_map <string, float, std_string_hash> Cp_number;
              hash_map <string, int, std_string_hash> Id_hash;
    
              void Init();
              int Load_id();
              int Load_order();   
              
              int LeafN;
              int OrderN;                         
              };

void _Comp_Tree::Init(){    
                    
     LeafN = Load_id();
     //load id hash
     for (int i = 0; i < LeafN; i ++)
        Id_hash[Id[i]] = i;
        
     OrderN = 0;
     //load tree  
     if (Database.Get_Is_Tree())      
        OrderN = Load_order();
    
     //load cp number         
     if (Database.Get_Is_Cp())
        Database.Load_Copy_Number(Cp_number);
     }

int _Comp_Tree::Load_id(){
     
     ifstream infile(Database.Get_Tree_Id().c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Id() << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           Id.push_back(buffer);
                           count ++;
                           }
     
     infile.close();
     infile.clear();
     return count;
     }


int _Comp_Tree::Load_order(){
    ifstream infile(Database.Get_Tree_Order().c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Order() << endl;
                 return 0;
                 }
    
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){                           
                           if(buffer.size() == 0) continue;
                           stringstream strin(buffer);
                           int order_1 = 0;
                           int order_2 = 0;
                           int order_d = 0;
                           float dist_1 = 0;
                           float dist_2 = 0;
                           strin >> order_1 >> dist_1 >> order_2 >> dist_2 >> order_d;
                           Order_1.push_back(order_1);
                           Order_2.push_back(order_2);
                           Order_d.push_back(order_d);
                           Dist_1.push_back(dist_1);
                           Dist_2.push_back(dist_2);
                           count ++;                           
                           }
       
    infile.close();
    infile.clear(); 
    
    return count;
    }

int _Comp_Tree::Load_abd(const char * infilename, float * Abd, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));
    
    hash_map<string, int, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
     
    _OTU_Parser otu_parser;
    
    otu_parser.Load_file_to_hash(infilename, otu_count);
    
    float total = 0;
        
    //cp_number_correct
    
    for (hash_map<string, int, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
        
        float cp_no = 1.0;
        
        if (is_cp_correct)
           if (Cp_number.count(miter->first) != 0)
                                          cp_no = Cp_number[miter->first];
        otu_abd[miter->first] = (float) miter->second / cp_no;
        total += otu_abd[miter->first];
        }
        
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
        miter->second /= total;
        //debug
        //cout << miter->first << endl;
        
        if (Id_hash.count(miter->first) != 0){
            Abd[Id_hash[miter->first]] = miter->second;
            mapped_otu_count ++;
            }
        }
        
    return mapped_otu_count;
    }

int _Comp_Tree::Load_abd(const char * infilename, float * Abd){
    
    return Load_abd(infilename, Abd, true);
    }

int _Comp_Tree::Load_abd(_Table_Format * table, float *Abd, int sample, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));

    vector <string> otus = table->Get_Feature_Names();
    vector <float> abds = table->Get_Abd(sample);
    
    hash_map<string, float, std_string_hash> otu_abd;
    
    float total = 0;
    
    for (int i = 0; i < otus.size(); i ++)
        if (abds[i] > 0){
            
            string a_otu = Check_OTU(otus[i]);
            float cp_no = 1.0;
            if (is_cp_correct)
               if (Cp_number.count(a_otu) != 0)
                                          cp_no = Cp_number[a_otu];
                                           
            if (Id_hash.count(a_otu) != 0){
                otu_abd[a_otu] = abds[i] / cp_no;
                total += otu_abd[a_otu];
                }                
            }    
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
                          miter->second /= total;
                          if (Id_hash.count(miter->first) != 0){
                             Abd[Id_hash[miter->first]] = miter->second;
                             mapped_otu_count ++;
                             }
                          }
    return mapped_otu_count;
    }

int _Comp_Tree::Load_abd(_Table_Format * table, float *Abd, int sample){
                                       return Load_abd(table, Abd, sample, true);
                                       }

float _Comp_Tree::Calc_sim(float * Abd_1, float * Abd_2, bool w){
      
      if (!Database.Get_Is_Tree()) return 1.0 - Calc_Dist_JSD(Abd_1, Abd_2, LeafN);
      
      if (w) return Calc_sim(Abd_1, Abd_2);
      else return Calc_sim_unweight(Abd_1, Abd_2);
      
      }

float _Comp_Tree::Calc_sim(float * Abd_1, float * Abd_2, int mode){ //0: MS; 1: Unweighted MS; 2: Cos; 3: Eu; 4: JSD; 5: Bray Curtis
      
      if (!Database.Get_Is_Tree()){
                                   if (mode < 2) mode = 4; //for no tree, default is 4: JSD
                                   }
      
      switch (mode){
             case 0: return Calc_sim(Abd_1, Abd_2); break;
             case 1: return Calc_sim_unweight(Abd_1, Abd_2); break;
             case 2: return 1.0 - Calc_Dist_Cos(Abd_1, Abd_2, LeafN); break;
             case 3: return 1.0 - Calc_Dist_E(Abd_1, Abd_2, LeafN); break;
             case 4: return 1.0 - Calc_Dist_JSD(Abd_1, Abd_2, LeafN); break;
	     case 5: return 1.0 - Calc_Dist_Bray_Curtis(Abd_1, Abd_2, LeafN); break;
             default: return 1.0 - Calc_Dist_JSD(Abd_1, Abd_2, LeafN); break;
             }
      return 0;
      }

float _Comp_Tree::Calc_sim(float * Abd_1, float * Abd_2){
      
      float Reg_1[REG_SIZE];
      float Reg_2[REG_SIZE];
      
      float total = 0;
      
      for(int i = 0; i < OrderN; i++){
              
              int order_1 = Order_1[i];
              int order_2 = Order_2[i];
              int order_d = Order_d[i] + REG_SIZE;
              
              float dist_1 = 1- Dist_1[i];
              float dist_2 = 1- Dist_2[i];
              
              dist_1 = (dist_1 < 0) ? MIN_DIST : dist_1;
              dist_2 = (dist_2 < 0) ? MIN_DIST : dist_2;
              
              float c1_1;
              float c1_2;
              
              float c2_1;
              float c2_2;
                                
              if (order_1 >= 0){
                          
                          c1_1 = Abd_1[order_1];
                          c1_2 = Abd_2[order_1];
                          
                          }
              else {
                   c1_1 = Reg_1[order_1 + REG_SIZE];
                   c1_2 = Reg_2[order_1 + REG_SIZE];
                   }
              
              if (order_2 >= 0){
                          
                          c2_1 = Abd_1[order_2];
                          c2_2 = Abd_2[order_2];
                          
                          }
              else {
                   c2_1 = Reg_1[order_2 + REG_SIZE];
                   c2_2 = Reg_2[order_2 + REG_SIZE];
                   }
              //min
              float min_1 = (c1_1 < c1_2)?c1_1:c1_2;
              float min_2 = (c2_1 < c2_2)?c2_1:c2_2;
              
              total += min_1;
              total += min_2;
              
              //reduce
              Reg_1[order_d] = (c1_1 - min_1) * dist_1 + (c2_1 - min_2) * dist_2;
              Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
              
              }
      
        total /= 100.0; //scale 0-1
        total = (total > 1.0) ? 1 : total;
        total = (total < 0.0) ? 0 : total;
        return total;
      }

float _Comp_Tree::Calc_sim_unweight(float * Abd_1, float * Abd_2){
      
      float Reg_1[REG_SIZE];
      float Reg_2[REG_SIZE];
      
      float total = 0;
      float total_count = 0;
      
      
      for(int i = 0; i < OrderN; i++){
              
              int order_1 = Order_1[i];
              int order_2 = Order_2[i];
              int order_d = Order_d[i] + REG_SIZE;
              
              float dist_1 = 1- Dist_1[i];
              float dist_2 = 1- Dist_2[i];
              
              dist_1 = (dist_1 < 0) ? MIN_DIST : dist_1;
              dist_2 = (dist_2 < 0) ? MIN_DIST : dist_2;
                            
              float c1_1 = 0;
              float c1_2 = 0;
              
              float c2_1 = 0;
              float c2_2 = 0;
                                
              if (order_1 >= 0){//leaf node
                          
                          c1_1 = (Abd_1[order_1] > 0)? 1 : 0;
                          c1_2 = (Abd_2[order_1] > 0)? 1 : 0;                                                    
                          
                          }
              else {//internal node
                   c1_1 = Reg_1[order_1 + REG_SIZE];
                   c1_2 = Reg_2[order_1 + REG_SIZE];
                   }
              
              total += (c1_1 > 0 && c1_2 > 0) ? dist_1 : 0;
              total_count += (c1_1 > 0 || c1_2 > 0) ? dist_1 : 0;
              
              if (order_2 >= 0){
                          
                          c2_1 = (Abd_1[order_2] > 0)? 1 : 0;
                          c2_2 = (Abd_2[order_2] > 0)? 1 : 0;
                                                    
                          }
              else {
                   c2_1 = Reg_1[order_2 + REG_SIZE];
                   c2_2 = Reg_2[order_2 + REG_SIZE];
                   }
              
              total += (c2_1 > 0 && c2_2 > 0) ? dist_2 : 0;              
              total_count += (c2_1 > 0 || c2_2 > 0) ? dist_2 : 0;              
              
              //reduce
              Reg_1[order_d] = (c1_1 > 0 || c2_1 > 0) ? 1: 0;
              Reg_2[order_d] = (c1_2 > 0 || c2_2 > 0) ? 1: 0;
                            
              }
      
        total /= total_count; //scale 0-1
        total = (total > 1.0) ? 1 : total;
        total = (total < 0.0) ? 0 : total;
        return total;
        
        /*
        for(int i = 0; i < LeafN; i ++){
                if ((Abd_1[i] > 0) && (Abd_2[i] > 0)) total ++;
                if ((Abd_1[i] > 0) || (Abd_2[i] > 0)) total_count ++;
                
                } 
       return total / total_count;
       */
      }
#endif
