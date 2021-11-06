// Updated at Aug 13, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include "hash.h"
#include "utility.h"
#include "db.h"

#include "table_format.h"
#include "dist.h"

#ifndef _COMP_FUNC_H
#define _COMP_FUNC_H

using namespace std;

class _Comp_Tree_Func{
      
      public:
             _Comp_Tree_Func(){                                                              
                               GeneN = Load_Id();  
                               Load_Pw();                   
                               }
             _Comp_Tree_Func(char db){
                               Database.Set_DB(db);
                               if (!Database.Get_Is_Func()){
                                                            cout << "The " << Database.Get_Description() << " domain is not supported yet" << endl;
                                                            exit(0);
                                                            }                                                             
                               GeneN = Load_Id();
                               Load_Pw();                               
                               }
    
             int Load_Gene_Count(const char * infilename, float * abd);
             int Load_Gene_Count(const char * table_file, vector<float * > & abds, vector<string> & samples);//load gene count by abd table
	     int Load_Gene_Count_Reverse(const char * table_file, vector<float * > & abds, vector<string> & samples);//load gene count by reverse abd table 

             int Load_Gene_Count(_Table_Format * table, float * abd, int sample);
             float Calc_sim(float * abd_m, float * abd_n, int mode);
             float Calc_sim_pw(float * abd_m, float * abd_n);
             
             int Get_GeneN(){
                 return GeneN;
                 }
             
			 string Get_Gene_Name(int n){
			 	if (n < 0) return ("NA");
			 	if (n >= Gene.size()) return ("NA");
			 	return Gene[n];
			 }                
             void Debug_output(const char * outfilename);
             
      private:
              int Load_Id();  
              int Load_Pw();            
              vector <string> Gene;
              
              vector <vector <string> > Pathway;
              vector <vector <set <int> > > Pw_arch;
              
              _PMDB Database;
              int GeneN; 
              int Pw_level;             
      };

int _Comp_Tree_Func::Load_Gene_Count(const char * infilename, float * abd){//load gene count by single abd file
    
    memset(abd, 0, GeneN * sizeof(float));
    
    hash_map<string, float, std_string_hash> table;
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open input file : " << infilename << endl;
                 exit(0);
                 }
    
    int count = 0;
    string buffer;
    getline(infile, buffer);//title
    while(getline(infile, buffer)){                          
                          stringstream strin(buffer);  
                          string gene;
                          string temp;
                          float gene_count;
        
                          strin >> gene >> gene_count;
                          
                          if (table.count(gene) == 0)
                                                table[gene] = gene_count;
                          else table[gene] += gene_count; 
                          count ++;
                          }
    
    for (int i = 0; i < GeneN; i ++){
        abd[i] = 0;
        if (table.count(Gene[i]) != 0)
           abd[i] = table[Gene[i]];
        }
                 
    infile.close();
    infile.clear();
    
    return count;
    }

int _Comp_Tree_Func::Load_Gene_Count(const char * table_file, vector<float * > & abds, vector<string> & samples) {//load gene count by abd table
   
    map<string, int> gene_id;
    for(int i = 0; i < GeneN; i++) {
	gene_id[Gene[i]] = i;
    }

    vector<string> genes;

    ifstream infile(table_file, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open input file : " << table_file << endl;
                 exit(0);
                 }

    int count = 0;
    string buffer, gene, sample;
    float gene_count;

    getline(infile, buffer);//title
    stringstream strin(buffer);

    strin >> gene;
    while(strin >> gene)
    	genes.push_back(gene);
 
    while(getline(infile, buffer)){
        float * abd = (float *) malloc(sizeof(float) * GeneN);
	memset(abd, 0, sizeof(float) * GeneN);
        stringstream strin(buffer);
        strin >> sample;
	samples.push_back(sample);
	int i = 0;
	while(strin >> gene_count) {
	    abd[gene_id[genes[i]]] = gene_count;
	    i++;
	}
	abds.push_back(abd);
    }

    return samples.size();

}

int _Comp_Tree_Func::Load_Gene_Count_Reverse(const char * table_file, vector<float * > & abds, vector<string> & samples) {//load gene count by reverse abd table 
    
    map<string, int> gene_id; 
    for(int i = 0; i < GeneN; i++) { 
        gene_id[Gene[i]] = i; 
    } 
 
    ifstream infile(table_file, ifstream::in); 
    if (!infile){ 
                 cerr << "Error: Cannot open input file : " << table_file << endl; 
                 exit(0); 
                 } 
 
    int count = 0; 
    string buffer, gene, sample; 
    float gene_count; 
 
    getline(infile, buffer);//title 
    stringstream strin(buffer); 
 
    strin >> sample; 
    while(strin >> sample) 
        samples.push_back(sample); 
    
    for(int i = 0; i < (int) samples.size(); i++) {
	float * abd = (float *) malloc(sizeof(float) * GeneN);
        memset(abd, 0, sizeof(float) * GeneN);  
	abds.push_back(abd);
    }
  
    while(getline(infile, buffer)){ 
        stringstream strin(buffer); 
        strin >> gene; 
        int i = 0; 
        while(strin >> gene_count) { 
            abds[i][gene_id[gene]] = gene_count; 
            i++; 
        } 
    }
    return samples.size(); 
 
} 


int _Comp_Tree_Func::Load_Gene_Count(_Table_Format *table, float * abd, int sample){
    
    memset(abd, 0, GeneN * sizeof(float));
    
    int count = 0;
    
    for (int i = 0; i < GeneN; i ++){
        abd[i] = table->Get_Abd_By_Feature(sample, Gene[i]);
        if (abd[i] > 0) count ++;
        }
    return count;
    }

float _Comp_Tree_Func::Calc_sim(float * abd_m, float * abd_n, int mode){
      
      switch (mode) {
             
             case 1: return 1.0 - Calc_Dist_Cos(abd_m, abd_n, GeneN); break;
             case 2: return 1.0 - Calc_Dist_E(abd_m, abd_n, GeneN); break;
             case 3: return 1.0 - Calc_Dist_JSD(abd_m, abd_n, GeneN); break;
             case 4: return 1.0 - Calc_Dist_Bray_Curtis(abd_m, abd_n, GeneN); break; 
             default:
             case 0: return Calc_sim_pw(abd_m, abd_n); break;
             }
      return 0;
      }

int _Comp_Tree_Func::Load_Id(){
    
    ifstream in_idfile(Database.Get_Func_Id().c_str(), ifstream::in);
    if (!in_idfile){
                    cerr << "Error: Open KO ID file error : " << Database.Get_Func_Id() << endl;
                    return 0;
                    }
    string buffer;
    int count = 0;
    while(getline(in_idfile, buffer)){
                          if (buffer.size() == 0) continue;
                          Gene.push_back(buffer);
                          count ++;
                          }
    in_idfile.close();
    in_idfile.clear();
    
    return count;
    }

int _Comp_Tree_Func::Load_Pw(){
    //
    ifstream in_pwfile(Database.Get_Func_Pw().c_str(), ifstream::in);
    if (!in_pwfile){
                    cerr << "Error: Open Pathway file error : " << Database.Get_Func_Pw() << endl;
                    return 0;
                    }
    
    string buffer;
    int count = 0;
    
    vector <vector <vector <string> > > pw_arch_vector;
    //read file
    while(getline(in_pwfile, buffer)){
                          
                          //make pw arch   
                          vector <vector <string> > pw_arch;
                          string a_pw_term;
                          vector <string> a_pw_arch;
                          
                          for (int i = 0; i < buffer.size(); i ++){
                              
                              if ((buffer[i] == '|') || (i == buffer.size() - 1)){
                                             
                                             a_pw_arch.push_back(a_pw_term);
                                             pw_arch.push_back(a_pw_arch);
                                             
                                             a_pw_term = "";
                                             //a_pw_arch.swap(vector <string>());
                                             vector<string>().swap(a_pw_arch);
                                             
                                             }
                              else if (buffer[i] == ';'){
                                   
                                   a_pw_arch.push_back(a_pw_term);
                                   a_pw_term = "";
                                   
                                   }
                              
                              else a_pw_term += buffer[i];
                              
                              }
                          
                          pw_arch_vector.push_back(pw_arch);    
                          
                          //debug display
                          /*
                          cout << Gene[count] << "\t";
                          for (int i = 0; i < pw_arch.size(); i ++){                              
                              if (i > 0) cout << "|";
                              
                              for (int j = 0; j < pw_arch[i].size(); j ++){
                                  if (j > 0) cout << ";";
                                  cout << pw_arch[i][j];
                                  }                                                                    
                              }
                          cout << endl;
                          */
                          count ++;                      
                          }
    
    Pw_level = pw_arch_vector[0][0].size();
    
    //debug
    //cout << Pw_level << " levels" << endl;
    
    hash_map <string, int, std_string_hash> pw_id_hash [Pw_level];
    vector <string> pw_id_vector [Pw_level];
    //
    vector < set <int> > pw_id_arch [Pw_level]; //KO is level 0
    
    //parse pw arch
    for (int i = 0; i <  pw_arch_vector.size(); i ++){
        
        pw_id_arch[0].push_back(set<int>());
        
        for (int j = 0; j < pw_arch_vector[i].size(); j ++){
            
            int last_level_index = i;            
            
            for (int k = Pw_level - 1; k >=0 ; k--){ //reversed order in pw_des.tab
                
                string this_pathway = pw_arch_vector[i][j][k];
                
                int level = Pw_level -k - 1;
                
                if (pw_id_hash[level].count(this_pathway) == 0){                                                                 
                                                                 pw_id_hash[level][this_pathway] = pw_id_vector[level].size();
                                                                 pw_id_vector[level].push_back(this_pathway);  
                                                                 
                                                                 if (k > 0)
                                                                    pw_id_arch[level + 1].push_back(set<int>());                                                               
                                                                 }
                                                                 
                int this_level_index = pw_id_hash[level][this_pathway];
                
                //debug
                //cout << this_level_index << " > " << last_level_index << endl;
                
                pw_id_arch[level][last_level_index].insert(this_level_index);
                
                last_level_index = this_level_index;
                }
            
            }    
    }
    
    for(int i = 0; i < Pw_level; i ++){
            
            Pathway.push_back(pw_id_vector[i]);
            Pw_arch.push_back(pw_id_arch[i]);
            
            }
    
    in_pwfile.close();
    in_pwfile.clear();
    
    
    //debug    
    /*
    for(int i = 0; i < Pw_level; i ++)
            cout << "Level " << i << " : " << Pathway[i].size() << " pathways" << "(last is " << Pw_arch[i].size() << ")" << endl;
    */
    return Pw_level;
    }

float _Comp_Tree_Func::Calc_sim_pw(float * abd_m, float * abd_n){
      
      float weight = 0;
      
      //diff = 2 * ( 1 - sim )
      //sim = 1 - diff / 2            
      vector <float> this_abd_m, this_abd_n, last_abd_m, last_abd_n;
      //normlize
      float sum_m = 0;
      float sum_n = 0;
      
      for (int i = 0; i < GeneN; i ++){
          last_abd_m.push_back(abd_m[i]);
          last_abd_n.push_back(abd_n[i]);
          
          sum_m += abd_m[i];
          sum_n += abd_n[i];
          }
      
      if (sum_m <= 0) return 1;
      if (sum_n <= 0) return 1;
      
      for (int i = 0; i < GeneN; i ++){
          last_abd_m[i] /= sum_m;
          last_abd_n[i] /= sum_n;
          }
          
      float sim = 1 - Calc_Dist_Bray_Curtis(&last_abd_m[0], &last_abd_n[0], GeneN);   
      
      //debug
      //cout << "KO sim: " << sim << endl;
      
      weight = (float) (Pw_level + 1);
      sim *= weight;
      
      //debug count
      //int debug_reduce_count = 0;
      
      for (int i = 0; i < Pw_level; i ++){
          
          //reduce
          this_abd_m = vector <float> (Pathway[i].size(), 0);
          this_abd_n = vector <float> (Pathway[i].size(), 0);
         
          for (int j = 0; j < Pw_arch[i].size(); j ++) { 
              
              //debug_reduce_count ++;            
                         
              if ((last_abd_m[j] > 0) || (last_abd_n[j] > 0)){
              
              //debug
              //cout << "Pw level " << i << " : " << j << " ";
              //cout << last_abd_m[j] << " vs " << last_abd_n[j] << " > ";
              //cout << " reduce #: " << Pw_arch[i][j].size() << " ";
              
              //float min = (last_abd_m[j] < last_abd_n[j]) ? last_abd_m[j] : last_abd_n[j];
              
              for (set<int>::iterator siter = Pw_arch[i][j].begin(); siter != Pw_arch[i][j].end(); siter ++){
                                      //debug
                                      //cout << *siter << "|";                                      
                                      this_abd_m[*siter] += last_abd_m[j];// - min;
                                      this_abd_n[*siter] += last_abd_n[j];// - min;         
                                                                                                        
                                      }
              //debug
              //cout << endl;
              }
              }
          //calc dist
          float res_sim = 1 - Calc_Dist_Bray_Curtis(&this_abd_m[0], &this_abd_n[0], Pathway[i].size());
          //debug
          //cout << "Pw level " << i + 1 << " sim: " << res_sim << endl;
          
          //add to total
          //sim += (1.0 - sim) * res_sim;
          sim += res_sim * (float) (Pw_level - i);
          weight += (float) (Pw_level - i);
                    
          last_abd_m = this_abd_m;
          last_abd_n = this_abd_n;
          }
      
      //debug
      //cout << "Reduce count: " << debug_reduce_count << endl;
                  
      return sim / weight;
      }


void _Comp_Tree_Func::Debug_output(const char * outfilename){
     
     ofstream outfile(outfilename, ios::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file: " << outfilename << endl;
                   return;
                   }
     
     for (int i = 0; i < Pw_arch[0].size(); i ++){
         
         outfile << Gene[i] << "\t";
         
         //for (int j = 0; j < Pw_arch[0][i].size(); j ++)
         for (set <int> ::iterator iter_j = Pw_arch[0][i].begin(); iter_j != Pw_arch[0][i].end(); iter_j ++)
             //for (int k = 0; k < Pw_arch[1][j].size(); k ++)
             for (set <int> ::iterator iter_k = Pw_arch[1][*iter_j].begin(); iter_k != Pw_arch[1][*iter_j].end(); iter_k ++)           
                 //for (int l = 0; l < Pw_arch[2][k].size(); l ++)
                 for (set <int> ::iterator iter_l = Pw_arch[2][*iter_k].begin(); iter_l != Pw_arch[2][*iter_k].end(); iter_l ++)           
                     outfile << Pathway[0][*iter_j] << ";" << Pathway[1][*iter_k] << ";" << Pathway[2][*iter_l] << "|" ;
                      
         outfile << endl;
         }
     
     outfile.close();
     outfile.clear();
     
     cout << "Test finished" << endl;
     }

#endif
