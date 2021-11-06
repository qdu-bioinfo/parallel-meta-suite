// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Version 3.4.2 or above with _OTU_Parser


#ifndef otu_parser_h
#define otu_parser_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "hash.h"
#include "db.h"

#define BUFF 10000000
#define STRLEN 1000
#define TAXA_LEVEL 8 //K, P, C, O, F, G, S, OTU

using namespace std;

class _OTU_Parser{
      
      public:
             
             _OTU_Parser(){};
             
             _OTU_Parser(_PMDB db){
                              Database = db;
                              Database.Read_Taxonomy(OTU_taxa_table);
                              Database.Load_Copy_Number(Cp_number);
                              };
             int Output_hash_to_table(const char * tablefilename, hash_map<string, int, std_string_hash> otu_count, bool is_cp);
             int Update_class_taxa(const char * infilename, const char * outfilename);
             int Load_file_to_hash(const char * classfilename, hash_map<string, int, std_string_hash> & otu_count);
             
             string Get_taxa_by_OTU(string otu);
             string Get_taxa_by_OTU(string otu, int level);
             void Get_taxa_by_OTU(string otu, string * taxa, int level);
             
             float Get_cp_by_OTU(string otu);                           
                          
      private:
              _PMDB Database;
              hash_map<string, string, std_string_hash> OTU_taxa_table;
              hash_map <string, float, std_string_hash> Cp_number;
      };

int _OTU_Parser::Output_hash_to_table(const char * tablefilename, hash_map<string, int, std_string_hash> otu_count, bool is_cp){
                 
                 hash_map <string, float, std_string_hash> otu_abd;
                 
                 //calc_abd
                 float seq_count = 0;
                 float otu_count_sum = 0;
                 for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){                     
                                  float cp = 1;
                                  if ((is_cp) && (Cp_number.count(miter->first) != 0))
                                            cp = Cp_number[miter->first];
                                  otu_count_sum += ((float) miter->second / cp );
                                  seq_count += miter->second;
                                  }
                 
                  for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                     
                     float cp = 1; 
                     if ((is_cp) && (Cp_number.count(miter->first) != 0))
                                 cp = Cp_number[miter->first];
                     
                     otu_abd[miter->first] = ((float) miter->second / cp ) / otu_count_sum;
                     }
                     
                 ofstream outfile(tablefilename, ios::out);
                 if (!outfile){
                               cerr << "Error: Cannot open output file: " << tablefilename << endl;
                               return 0;
                               }
                 
                 outfile << "#Database_OTU\tCount\tAbundance\tTaxonomy" << endl;
                 
                 for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                     
                     string taxa = "Unclassified; otu_" + miter->first;
                     if (OTU_taxa_table.count(miter->first) != 0) taxa = OTU_taxa_table[miter->first];
                     
                     outfile << miter->first << "\t" << miter->second << "\t" << otu_abd[miter->first] << "\t" << taxa << endl;
                     }
                                                   
                 outfile.close();
                 outfile.clear(); 
                 
                 return  seq_count;            
                 }
                                           
int _OTU_Parser::Update_class_taxa(const char * infilename, const char * outfilename){
                 hash_map<string, int, std_string_hash> otu_count;
                 int count = Load_file_to_hash(infilename, otu_count);
                 Output_hash_to_table(outfilename, otu_count, true);
                 
                 return count;   
                 }   
            
int _OTU_Parser::Load_file_to_hash(const char * classfilename, hash_map<string, int, std_string_hash> & otu_count){
                    //fgets
                    FILE * fptr = fopen(classfilename, "r");
                    if (fptr == NULL){
                       cerr << "Error: Cannot open input file: " << classfilename << endl;
                       return -1;
                       }
                 
                  char * str_buffer = new char [BUFF];
                  vector <string> file_buffer;
                  
                   while(fgets(str_buffer, BUFF-1, fptr) != NULL){
                                        
                                        string a_line = str_buffer;
                                        file_buffer.push_back(a_line);
                                        
                                        }
                   delete [] str_buffer;
                   fclose(fptr);
                 
                 int mode = 1; //default is new format
                 
                 string title = file_buffer[0];
                 if (title[1] == 'S') mode = 0;
                 
                 unsigned int count = 0;                              
                 
                 if (mode == 0){ //detail mode, old format                          
                   for(int i = 1; i < file_buffer.size(); i ++){
                   					   char str_seq_id [STRLEN];
                   					   char str_id [STRLEN];
                                       sscanf(file_buffer[i].c_str(), "%s%s", str_seq_id, str_id);
                                       string id = str_id;
                                       
                                       if (otu_count.count(id) == 0)
                                             otu_count[id] = 1;
                                       else otu_count[id] ++;
                                       
                                       count ++;
                                       } 
                     }
                 
                 else{//simple mode, new format 
                       for(int i = 1; i < file_buffer.size(); i ++){
                                       char str_id [STRLEN];
                                       int id_count = 0;                                       
                                       sscanf(file_buffer[i].c_str(), "%s%d", str_id, &id_count);									   
                                       string id = str_id;
                                       otu_count[id] = id_count;
                                       count += id_count;
                                       }                                             
                      } 
                 return count;
                 }
                                         
string _OTU_Parser::Get_taxa_by_OTU(string otu){
                    if (OTU_taxa_table.count(otu) != 0)
                       return OTU_taxa_table[otu];
                    return "Unclassified";
                    }
            
string _OTU_Parser::Get_taxa_by_OTU(string otu, int level){
                   string taxon = "Unclassified";
                   
                   if ((level > TAXA_LEVEL) || (level < 0)) return taxon;
                   
                   if (OTU_taxa_table.count(otu) == 0)
                   
                      return taxon;
                      
                   string otu_taxa = OTU_taxa_table[otu];
                   vector <string> taxa_buffer;
                   stringstream strin(otu_taxa);
                   
                   string temp;
    
                   if (level == 8){ //otu
                             taxon = "OTU_";
                             taxon += otu;
                   }
    
                   else for (int i = 0; i < level; i++){
                            strin >> taxon;
        
                            if (taxon.find("Unclassified") != string::npos)
                               return "Unclassified";
                            
                            if (taxon.find("otu_") != string::npos)
                               return "Unclassified";
                            
        
                            if(taxon[taxon.size()-1] != ';'){
                                                   strin >> temp;
                                                   taxon += "_";
                                                   taxon += temp;
                            }
                            taxa_buffer.push_back(taxon);
        					/*
                            //add genus info for species
                            if (i == 6){
                               taxon = taxa_buffer[5];
                               taxon = taxon.substr(0, taxon.size()-1);
                               taxon += "_";
                               taxon += taxa_buffer[6];
                            }
        
                            //add prefix for genus
                            if (i == 5){
                               string taxa_prefix = taxa_buffer[i-1];
                               taxa_prefix = taxa_prefix.substr(0, taxa_prefix.size() - 1);
                               if (taxa_prefix.size() > 3) taxa_prefix = taxa_prefix.substr(0, 3);
                               taxon = taxa_prefix + "_" + taxon;
                               }
                   			*/
                               }
    
                   return taxon;
                   }
             
void _OTU_Parser::Get_taxa_by_OTU(string otu, string * taxa, int level){

                 if ((level > TAXA_LEVEL)) level = TAXA_LEVEL;
                 
                 string otu_taxa = "Unclassified";
                 if (OTU_taxa_table.count(otu) != 0)
                      otu_taxa = OTU_taxa_table[otu];
    
                 vector <string> taxa_buffer;
                 stringstream strin(otu_taxa);
                 string taxon;
                 string temp;
                 bool is_taxa = true;

                 for (int i = 0; i < level; i++){        
        
                     if (!(strin >> taxon)){
                     taxon =  "Unclassified";
                     is_taxa = false;
                     }
        
               else if (taxon.find("Unclassified") != string::npos){
                    taxon =  "Unclassified";
                    is_taxa = false;
                    }
              else if (taxon.find("otu_") != string::npos){
                   taxon = "Unclassified";
                   is_taxa = false;
                   }
        
              else if(taxon[taxon.size()-1] != ';'){
                   strin >> temp;
                   taxon += "_";
                   taxon += temp;
              }
        
              taxa_buffer.push_back(taxon);
        		/*
              //add genus info for species
              if ((i == 6) && (is_taxa)) {
                 taxon = taxa_buffer[5];
                 taxon = taxon.substr(0, taxon.size()-1);
                 taxon += "_";
                 taxon += taxa_buffer[6];
                 }
        
              //add prefix for genus
              if ((i == 5) && (is_taxa)) {
                   string taxa_prefix = taxa_buffer[i-1];
                   taxa_prefix = taxa_prefix.substr(0, taxa_prefix.size() - 1);
                   if (taxa_prefix.size() > 3) taxa_prefix = taxa_prefix.substr(0, 3);
                   taxon = taxa_prefix + "_" + taxon;
                   }
                   
                */
               taxa[i] = taxon;
               }
            return;   
            }
                    
float _OTU_Parser::Get_cp_by_OTU(string otu){
                    if (Cp_number.count(otu) != 0)
                       return Cp_number[otu];
                    return 1;
                    } 
#endif
