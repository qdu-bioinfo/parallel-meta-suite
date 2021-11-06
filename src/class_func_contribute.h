// class_func_contribution
// Updated at Sept 20, 2018
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
// otu_parser

#include "class_func.h"
#include "utility.h"

#ifndef class_func_contribute_h
#define class_func_contribute_h

using namespace std;
      
class _KO_Contribute{
      
      public:
             friend class _KO_OTU_Table_Contribute;
             
             _KO_Contribute(){}
                                       
             void Set_Contribute(string otu, float ko_abd, float otu_abd, float ko_copy){
                                    KO_abd[otu] = ko_abd;
                                    OTU_abd[otu] = otu_abd;
                                    KO_copy[otu] = ko_copy;
                                    }
             
             void Add_Contribute(string otu, float ko_abd, float otu_abd, float ko_copy){
                                    if (KO_abd.count(otu) == 0){
                                                                    KO_abd[otu] = ko_abd;
                                                                    OTU_abd[otu] = otu_abd;
                                                                    }
                                    else{
                                        KO_abd[otu] += ko_abd;
                                        OTU_abd[otu] += otu_abd; 
                                        }
                                        
                                    KO_copy[otu] = ko_copy;     
                                    }
             float Get_Contribute(string otu){
                                    if (KO_abd.count(otu) == 0)
                                       return 0;
                                    else
                                        return KO_abd[otu];
                                    }
             float Get_OTU_Abd(string otu){
                                    if (OTU_abd.count(otu) == 0)
                                       return 0;
                                    else
                                        return OTU_abd[otu];
                                    }
             float Get_KO_Copy(string otu){
                                    if (KO_copy.count(otu) == 0)
                                       return 0;
                                    else
                                        return KO_copy[otu];
                                    }
                          
      private:
              hash_map <string, float, std_string_hash> KO_abd;
              hash_map <string, float, std_string_hash> OTU_abd;
              hash_map <string, float, std_string_hash> KO_copy;
      };

class _KO_OTU_Table_Contribute : public _KO_OTU_Table_All{
      
      public:
             _KO_OTU_Table_Contribute() : _KO_OTU_Table_All(){}
             _KO_OTU_Table_Contribute(_PMDB db, int count, int mode) : _KO_OTU_Table_All(db, count, mode){
                                            
                                            Sample_KO_Contribute = new hash_map <int, _KO_Contribute> [count];
                                            
                                            }
             
             int Load_Selected_KO(const char * infilename);
             int Output_Contribute(const char * outfilename);
             
      private:       
              hash_map <int, _KO_Contribute> * Sample_KO_Contribute;
              hash_set <int> Selected_KO;  
              
              int Load_Sample_By_OTU(hash_map <string, int, std_string_hash> * otu_seq_count, int sample); 
      };

int _KO_OTU_Table_Contribute::Load_Selected_KO(const char * infilename){
    
    vector <string> kos;
    Load_ID(infilename, kos);
    for (int i = 0; i < kos.size(); i ++){
        
        if (KO_Index.count(kos[i]) == 0){
                                   cerr << "Waring: Invalid input KO: " << kos[i] << endl;
                                   }
        else Selected_KO.insert(KO_Index[kos[i]]);
        }
    return Selected_KO.size();
    }

int _KO_OTU_Table_Contribute::Output_Contribute(const char * outfilename){
    
    //calc ko sum 
    hash_map <int, float> KO_abd_sum;
    for (int i = 0; i < KOs.size(); i ++){
        float this_ko_abd_sum = 0;
        for (int j = 0; j < Sample_count; j ++)
            this_ko_abd_sum += KO_Abd[j][KOs[i].Get_Index()];
        if (this_ko_abd_sum > 0) KO_abd_sum[KOs[i].Get_Index()] = this_ko_abd_sum;
        }
    
    ofstream outfile(outfilename, ofstream::out);
        if (!outfile){
            cerr << "Error: Cannot open output file : " << outfilename << endl;
            return 0;
        }
    
    int count = 0;
    
    outfile << "Gene\tSample\tOTU\tGeneCountPerGenome\tOTUAbundanceInSample\tCountContributebyOTU\tContributionPercentOfSample\tContributionPercentOfAllSamples\tTaxaonomy" << endl;
    
    for (hash_map <int, float>::iterator miter = KO_abd_sum.begin(); miter != KO_abd_sum.end(); miter ++){
        
        if ((Selected_KO.size() > 0) && (Selected_KO.count(miter->first) == 0)) continue;
                        
        count ++;
        
        int ko_index = miter->first;
        //debug
        //cout << ko_index << ": " << KOs[ko_index].Get_Name() << endl; 
                
        for (int i = 0; i < Sample_count; i ++){
            //debug
            //cout << Sample_KO_Contribute[i].size() << endl;
            
            if (Sample_KO_Contribute[i].count(ko_index) != 0){
                                                                                                                        
                                                            for (hash_map <string, float, std_string_hash> :: iterator miter_2 = Sample_KO_Contribute[i][ko_index].KO_abd.begin(); miter_2 != Sample_KO_Contribute[i][ko_index].KO_abd.end(); miter_2 ++){
                                                                string otu = miter_2->first;
                                                                float ko_abd = miter_2->second;
                                                                
                                                                outfile << KOs[ko_index].Get_Name() << "\t";
                                                                outfile << Sample_names[i] << "\t";
                                                                outfile << otu << "\t";
                                                                outfile << Sample_KO_Contribute[i][ko_index].Get_KO_Copy(otu) << "\t";
                                                                outfile << Sample_KO_Contribute[i][ko_index].Get_OTU_Abd(otu) << "\t";
                                                                outfile << ko_abd << "\t";
                                                                outfile << ko_abd / KO_Abd[i][ko_index] << "\t";
                                                                outfile << ko_abd / KO_abd_sum[ko_index] << "\t";
                                                                outfile << Otu_parser.Get_taxa_by_OTU(otu) << endl;
                                                                }
                                                            }
            }
                
        }
    
    outfile.close();
    outfile.clear();
    
    return count;
    }

int _KO_OTU_Table_Contribute::Load_Sample_By_OTU(hash_map <string, int, std_string_hash> * otu_seq_count, int sample){
    
    //debug
    //cout << "It starts" << endl; 
    //sum 16S count
    float otu_seq_abd_sum = 0;
    for (hash_map <string, int, std_string_hash> ::iterator miter = otu_seq_count->begin(); miter != otu_seq_count->end(); miter ++){
        float seq_abd = miter->second;
        float rna_cp = Otu_parser.Get_cp_by_OTU(miter->first);
        seq_abd /= rna_cp;
        otu_seq_abd_sum += seq_abd;
        }
    
    //hash_map <int, _KO_Contribute> this_sample_ko_contribute;
    
    for (hash_map <string, int, std_string_hash> ::iterator miter = otu_seq_count->begin(); miter != otu_seq_count->end(); miter ++){
            
            vector<_KO_Index_Copy> kos = OTU_KO_Index[miter->first];
            int seq_count = miter->second;
            
            for (int i = 0; i < kos.size(); i ++){
                
                int index = kos[i].Index;
                float copy = kos[i].Copy;
                float rna_cp = Otu_parser.Get_cp_by_OTU(miter->first);         
                float ko_abd = (float) seq_count * copy / rna_cp;
                KO_Abd[sample][index] += ko_abd;
                
                //contribute
                float otu_abd = ((float) seq_count / rna_cp ) / otu_seq_abd_sum;
                //this_sample_ko_contribute[index].Add_Contribute(miter->first, ko_abd, otu_abd, copy);
                Sample_KO_Contribute[sample][index].Add_Contribute(miter->first, ko_abd, otu_abd, copy);
                }
        }
        //Sample_KO_Contribute[sample] = this_sample_ko_contribute;
        return otu_seq_count->size();
    
    }
#endif
