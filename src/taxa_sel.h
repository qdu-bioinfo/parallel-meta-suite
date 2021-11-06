// Updated at Dec 20, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
// _Table_Format input support

#ifndef taxa_sel_h
#define taxa_sel_h

#include "table_format.h"
#include "hash.h"
#include "db.h"
using namespace std;

class _Table_Format_Seq{
    public:
    _Table_Format_Seq(){}
    _Table_Format_Seq(int sample_count){
        for (int i = 0; i < sample_count; i ++){
            vector <unsigned int> empty_vector;
            vector <float> empty_float;
            Seq_count.push_back(empty_vector);
            Abd_count.push_back(empty_float);
            }
        }
    
    _Table_Format_Seq(vector <string> sample_names){
        Abd_table.Samples = sample_names;
        for (int i = 0; i < Abd_table.Samples.size(); i ++){
            vector <unsigned int> empty_vector;
            vector <float> empty_float;
            Seq_count.push_back(empty_vector);
            Abd_count.push_back(empty_float);
            }
        }
    
    void Add_Feature(string feature, unsigned int count, int sample, float cp){
        
        if (Abd_table.Feature_hash.count(feature) == 0){
            Abd_table.Feature_hash[feature] = Abd_table.Features.size();
            Abd_table.Features.push_back(feature);
            
            for (int i = 0; i < Seq_count.size(); i ++){
                Seq_count[i].push_back(0);
                Abd_count[i].push_back(0);
                }
            }
        
        int feature_index = Abd_table.Feature_hash[feature];
        
        Seq_count[sample][feature_index] += count;
        Abd_count[sample][feature_index] += (float) count / cp;
        }
    
    void Add_Feature(string feature, int sample, float cp){
        
        Add_Feature(feature, 1, sample, cp);
        
        }

    
    void Filter_Seq_Count(unsigned int c){
        
        for (int i = 0; i < Seq_count.size(); i ++)
            for (int j = 0; j < Seq_count[i].size(); j ++)
                if (Seq_count[i][j] < c){
                    Seq_count[i][j] = 0;
                    Abd_count[i][j] = 0;
                    }
        }
    
    void Filter_Abd(float max, float min, float no_zero_rate, float ave_t){
         
         Abd_table.Init_Filter();
         
         Abd_table.Filter_Max(max);
         Abd_table.Filter_Min(min);
         Abd_table.Filter_Ave(ave_t);
         Abd_table.Filter_Zero(no_zero_rate);
         Abd_table.Filter_Empty();
         }
    
    //
    void Normalization(){                  
         
         for (int i = 0; i < Seq_count.size(); i ++){
             
             vector <float> abd;
             float seq_sum = 0;
             
             for (int j = 0; j < Abd_count[i].size(); j ++)
                 seq_sum += Abd_count[i][j];
                                  
             for (int j = 0; j < Abd_count[i].size(); j ++)
                 if (seq_sum > 0)
                    abd.push_back(Abd_count[i][j] / seq_sum);
                 else abd.push_back(0);
                 
             Abd_table.Abd.push_back(abd);             
             }
         }
         
    unsigned int Output_Abd(const char * outfilename){
         
         return Abd_table.Output_Table(outfilename);
          
         }
         
    unsigned int Output_Count(const char * outfilename){
         
         unsigned int out_count = 0;
         
         fstream outfile(outfilename, ofstream::out);
         if (!outfile){
                       cerr << "Error: Cannot open output file: outfilename" << endl;
                       return 0;
                       }
         
         outfile << "SampleID";
         for (int i = 0; i < Abd_table.Features.size(); i ++)
              if (!Abd_table.Check_Filter(Abd_table.Features[i])){
                 outfile << "\t" << Abd_table.Features[i];
                 out_count ++;
                 }
              outfile << endl;
    
         for (int i = 0; i < Abd_table.Samples.size(); i ++){
             outfile << Abd_table.Samples[i];
             for (int j = 0; j < Abd_table.Features.size(); j ++)
                 if (!Abd_table.Check_Filter(Abd_table.Features[j]))
                    outfile << "\t" << Seq_count[i][j];
                    outfile << endl;
                 }                  
         outfile.close();
         outfile.clear();
         
         return out_count;
         }             
    
    int Get_Taxa_Size(){
        return Abd_table.Get_Feature_Size();
        }
    
    private:
    _Table_Format Abd_table;
    vector < vector <unsigned int> > Seq_count;
    vector < vector <float> > Abd_count;
    
};

#endif /* table_format_seq_h */
