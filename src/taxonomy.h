// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Added CIGAR
// _OTU_parser

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>

#include "hash.h"
#include "db.h"
#include "otu_parser.h"

#ifndef TAXONOMY_H
#define TAXONOMY_H

#define LEVEL 7 //(0: Kingdom, 5: Genus, 7: OTU)

using namespace std;

float Parse_CIGAR(string c){
         
         map <char, int> c_table;
         int value = 0;
         int total = 0;
         
         for (int i = 0; i < c.size(); i ++){
             
             if ((c[i] >= '0') && (c[i] <= '9'))
                value = value * 10 + c[i] -'0';
             else{
                  if (c_table.count(c[i]) == 0)
                     c_table[c[i]] = 0;
                  
                  c_table[c[i]] += value;
                  total += value;
                  value = 0;
                  }
             }
         
         if (c_table.count('M') == 0) return 0;
         else return (float) c_table['M'] / (float) total;
         }

    
void Parse_Taxonomy(const char * infilename, const char * outfilename, _OTU_Parser otu_parser, float Q, int * a_diver, bool is_paired, int & match_count, int & drop_count){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Open Mapping File error : " << infilename << endl;
                 exit(0);
                 }

    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                  cerr <<"Error: Open output file error : " << outfilename << endl;
                  exit(0);
                  }
    
    //outfile << "#Sequence_Id\tDatabase_OTU\tPOS\tCIGAR\tTaxonomy" << endl;
    outfile << "#Database_OTU\tCount\tTaxonomy" << endl;
    
    unsigned int total_count = 0;
    
    string buffer;
    // string last_seq_id = "";
    
    set <string> a_diver_table [LEVEL];
        
    while(getline(infile, buffer)){
                          if (buffer.size()==0) continue;
                          stringstream strin(buffer);
                          
                          //string seq_id, database_id;
                          //int flag, mapq, pos;
                          //string cigar;
                          //strin >> seq_id >> flag >> database_id >> pos >> mapq >> cigar;
                           
						  string database_id,count;      
						  strin >> database_id >> count;
						  
						  int m_count = atoi(count.c_str());
						  
						  if (database_id.find("#OTU") != string::npos){//table first line 
						  	continue;
						  }
						                     
                          if (is_paired)
                             getline(infile, buffer);
                          
                          if (database_id[0] == '*') continue; 
                          total_count ++;                               
                          /*
                          if (Parse_CIGAR(cigar) < Q) {
                                                 drop_count ++;
                                                 continue;
                                                 }                                                                                                        
                          */                                              
                          /*if (seq_id == last_seq_id) continue; //Duplications
                          else last_seq_id = seq_id;
                          */
                          //outfile << seq_id << "\t" << database_id << "\t" << pos << "\t" << cigar << "\t";
                          outfile << database_id << "\t" << count << "\t";
						      
                          //Output the taxonomy annotation
                          string taxa = otu_parser.Get_taxa_by_OTU(database_id);
                          outfile << taxa << endl;
                                                    
                          /*
                          stringstream strin_taxa(taxa);
                          for (int i = 0; i < LEVEL - 1; i ++){                              
                              string a_taxa;
                              if (!(strin_taxa >> a_taxa)) break;
                              if (a_taxa.find("otu_") != string::npos) break;
                              if (a_taxa.find("Unclassified") != string::npos) break;
                              a_diver_table[i].insert(a_taxa);
                              }
                          */
                          string taxa_vector[LEVEL - 1];
                          otu_parser.Get_taxa_by_OTU(database_id, taxa_vector, LEVEL - 1);
                          for (int i = 0; i < LEVEL - 1; i ++){                                                            
                              if (taxa_vector[i].find("otu_") != string::npos) break;
                              if (taxa_vector[i].find("Unclassified") != string::npos) break;
                              a_diver_table[i].insert(taxa_vector[i]);
                              }                                                                                        
                          
                          a_diver_table [LEVEL - 1].insert(database_id); //OTU
                          
                          match_count =match_count + m_count;                                                    
                          }
                          
    
    for (int i = 0; i < LEVEL; i ++)
        a_diver[i] = a_diver_table[i].size();
    
    infile.close();
    infile.clear();
    outfile.close();
    outfile.clear();
    
    cout << endl << match_count << " matches are parsed" << endl;
    }

unsigned int Out_Taxonomy(const char * infilename, string out_path, _PMDB db, int Q, int * a_diver, int is_paired, int & match_count, int & drop_count){
    
    _OTU_Parser otu_parser(db);
    
    string outfilename = out_path + "/classification.txt";
    string outfilename_detail = out_path + "/classification_detail.txt";
    
    //Parse_Taxonomy(infilename, outfilename_detail.c_str(), otu_parser, 0, a_diver, is_paired, match_count, drop_count);
    Parse_Taxonomy(infilename, outfilename.c_str(), otu_parser, 0, a_diver, is_paired, match_count, drop_count);
    
	//otu_parser.Update_class_taxa(outfilename_detail.c_str(), outfilename.c_str()); 
    otu_parser.Update_class_taxa(outfilename.c_str(), outfilename.c_str()); 
    
    cout << endl << match_count << " taxonomy annotations are parsed out" << endl << endl;
    
    return match_count;
    }

#endif
/*

int main(int argc, char * argv[]){
    
    string out_path = argv[3];
    string database = argv[2];
    
    Out_Taxonomy(argv[1], database, out_path, 0);

    
    return 0;
    
    }
*/
