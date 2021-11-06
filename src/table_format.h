// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Add_feature and Add abd
#ifndef table_format_h
#define table_format_h

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

using namespace std;

class _Table_Format{
    
    public:
    
    friend class _Table_Format_Seq;
    
    _Table_Format(){}    
    _Table_Format(const char * infilename, bool order); //order: true:norm, false: rev
    _Table_Format(const char * infilename); //default is true
    _Table_Format(vector <string> features); //init by features
    
    int Load_Table(const char * infilename); //row: sample, column: feature
    int Load_Table_Rev(const char * infilename); //row: feature, column: sample
    void Add_Abd(vector <float> abd, string s); //row: sample, column: feature
    
    unsigned int Output_Table(ostream * outfile); //row: sample, column: feature
    unsigned int Output_Table(const char * outfilename);
    
     //Honglei, please implement the following function
    unsigned int Output_Table_Rev(ostream * outfile); //row: feature, column: sample
    unsigned int Output_Table_Rev(const char * outfilename); //row: feature, column: sample
    
    void Filter_Max(float max); // Filter features by the Maximum Value
    void Filter_Min(float min); // Filter features by the Minimum Value
    void Filter_Ave(float ave); // Filter features by the average Value
    void Filter_Zero(float zero); // Filter features by the proportion of non-zero values
    void Filter_Empty(); //Filter empty features
        
    vector <string> Get_Sample_Names();
    vector <string> Get_Feature_Names();
    int Get_Sample_Size();//Get the size of Samples
    int Get_Feature_Size();//Get the size of Features
    
    float Get_Abd_By_Order(unsigned int s, unsigned int i);
    float Get_Abd_By_Feature(unsigned int s, string a_feature);
    vector <float> Get_Abd(unsigned int s);
    
    float Calc_Dist_Cos(int sam_m, int sam_n); //Calculate the Cosin distance between two samples
    float Calc_Dist_E(int sam_m, int sam_n); //Calculate the Euclidean distance between two samples
    float Calc_Dist_JSD(int sam_m, int sam_n); //Calculate the Jeason-Shannon distance between two samples
    float Calc_Dist_Bray_Curtis(int sam_m, int sam_n); //Calculate the Bray-Curtis distance between two samples
    
    float Calc_Corr_S(int sam_m, int sam_n);//Calculate Spearman correlation coefficient
    float Calc_Corr_P(int sam_m, int sam_n);//Calculate Pearson correlation coefficient
    
    float Calc_Corr_S(vector <float> sam_m, vector <float> sam_n);
    float Calc_Corr_P(vector <float> sam_m, vector <float> sam_n);
    
    void Calc_Dist_Matrix(const char * outfilename, int metrics, int coren, bool is_sim); //0: cost 1: eu dist 2: jsd
    void Calc_Corr_Matrix(const char * outfilename, int metrics, int coren); //0: s-corr 1: p-corr
    
    protected:
    vector <string> Samples;
    vector <string> Features;
    vector < vector <float> > Abd;
    
    map <string, int> Feature_hash;
    
    bool * Max_filtered;
    bool * Min_filtered;
    bool * Ave_filtered;
    bool * Zero_filtered;
    bool * Empty_filtered;
    
    void Init_Filter(); //init features
    bool Check_Filter(string a_feature);//check features to determine whether to filter
    void Build_Feature_Hash();// Build the hash numbers of features.
    void BubbleSort(float *array1, int *rank1, int len);//Bubble sort code, was used by the function named 'Calc_Corr_S'

};

_Table_Format::_Table_Format(const char * infilename, bool order){ //order: true:norm, false: rev
    if (order) Load_Table(infilename);
    else Load_Table_Rev(infilename);
    
}

_Table_Format::_Table_Format(const char * infilename){
                                   Load_Table(infilename);
                                  
                                   } //default true

_Table_Format::_Table_Format(vector <string> features){
                                   Features = features;
                                   Build_Feature_Hash();
                                   Init_Filter();
                                   } //default true

int _Table_Format::Load_Table(const char * infilename){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
        cerr << "Error: Cannot open file : " << infilename << endl;
        return 0;
    }
    
    string buffer;
    
    //Samples
    getline(infile, buffer);
    stringstream strin(buffer);
    string feature_name;
    strin >> feature_name; // Title
    while(strin >> feature_name){
		//cout << feature_name<< "\t" ;
    	Features.push_back(feature_name);
	}
	//cout<<endl; 

    
    //Data
    while(getline(infile, buffer)){
        stringstream strin(buffer);
        string sample_id;
        strin >> sample_id;
        
        vector <float> abd;
        float a_abd;
        while(strin >> a_abd){
        	//cout << a_abd << "\t";
        	abd.push_back(a_abd);
		}
		//cout<<endl;

        //cout<< abd.size() <<endl;
        //cout<< Features.size() <<endl;
        //check feature number
        if (abd.size() != Features.size()){
            cerr << "Error: Sample: " << sample_id << " does not have " << Features.size() << " features" << endl;
            continue;
            }
        
        Samples.push_back(sample_id);
        Abd.push_back(abd);
        }
    
    infile.close();
    infile.clear();
    
    Build_Feature_Hash();
    Init_Filter();
    
    return Samples.size();
    }

int _Table_Format::Load_Table_Rev(const char * infilename){

    ifstream infile(infilename, ifstream::in);
    if (!infile){
        cerr << "Error: Cannot open file : " << infilename << endl;
        return 0;
    }

    Abd.clear();
    Samples.clear();
    Features.clear();

    vector< vector<float> > Abd_temp;

    string buffer;
    getline(infile, buffer);
    stringstream strin(buffer);
    string sample_name;
    strin >> sample_name; // Title
    while(strin >> sample_name)
        Samples.push_back(sample_name);

    while(getline(infile, buffer)){
        stringstream strin(buffer);
        string feature_id;
        strin >> feature_id;

        vector <float> abd;
        float a_abd;
        while(strin >> a_abd)
            abd.push_back(a_abd);


        if (abd.size() != Samples.size()){
            cerr << "Error: Features: " << feature_id << " does not have " << Samples.size() << " Samples" << endl;
            continue;
            }

        Features.push_back(feature_id);
        Abd_temp.push_back(abd);
        }


     for (int i=0; i<Samples.size(); i++){
         vector < float > abd1;
         for (int j=0; j<Features.size(); j++)
             abd1.push_back(Abd_temp[j][i]);

         Abd.push_back(abd1);
         }


    infile.close();
    infile.clear();
    Abd_temp.clear();
    
    Build_Feature_Hash();
    Init_Filter();

    return Samples.size();
        
    }

void _Table_Format::Add_Abd(vector <float> abd, string s){
     
     Samples.push_back(s);
     Abd.push_back(abd);
     
     } //row: sample, column: feature


unsigned int _Table_Format::Output_Table(ostream * outfile){
    
    unsigned int out_count = 0;
    
    *outfile << "SampleID";
    for (int i = 0; i < Features.size(); i ++)
        if (!Check_Filter(Features[i])){
           *outfile << "\t" << Features[i];
           out_count ++;
           }
    *outfile << endl;
    
    for (int i = 0; i < Samples.size(); i ++){
        *outfile << Samples[i];
        for (int j = 0; j < Features.size(); j ++)
            if (!Check_Filter(Features[j]))
               *outfile << "\t" << Abd[i][j];
        *outfile << endl;
        }
    return out_count;
    }

unsigned int _Table_Format::Output_Table(const char * outfilename){
            
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
        cerr << "Error: Cannot open output file: " << outfilename << endl;
        return 0;
        }
    
    unsigned int out_count = Output_Table(&outfile);

    outfile.close();
    outfile.clear();
    
    return out_count;
    }

// #############################################################################

vector <string> _Table_Format::Get_Sample_Names(){
    return Samples;
    }
vector <string> _Table_Format::Get_Feature_Names(){
       return Features;
       }
int _Table_Format::Get_Sample_Size(){
    return Samples.size();
    }
    
int _Table_Format::Get_Feature_Size(){
    return Features.size();
    }

unsigned int _Table_Format::Output_Table_Rev(ostream * outfile){
    
    unsigned int out_count = 0;
    
    *outfile << "SampleID";
    for (int i = 0; i < Samples.size(); i ++)
        *outfile << "\t" << Samples[i];
    *outfile << endl;

    for (int i = 0; i < Features.size(); i++){
        if (!Check_Filter(Features[i])){
           *outfile << Features[i];
           for (int j = 0; j < Samples.size(); j++)
               *outfile << "\t" << Abd[j][i];
           *outfile << endl;
           out_count ++;
           }
        }
     return out_count;
     }

unsigned int _Table_Format::Output_Table_Rev(const char * outfilename){

     ofstream outfile(outfilename, ofstream::out);
     if (!outfile){
        cerr << "Error: Cannot open output file: " << outfilename << endl;
        return 0;
        }

      unsigned int out_count = Output_Table_Rev(&outfile);

      outfile.close();
      outfile.clear();
      
      return out_count;
      }
// 
 
void _Table_Format::Filter_Max(float max){

      float max_temp;

      for (int i = 0; i < Features.size(); i ++){
          max_temp = 0.0;
          for (int j = 0; j < Samples.size(); j++)
               max_temp = max_temp > Abd[j][i] ? max_temp : Abd[j][i];

          if (max_temp < max)
             Max_filtered[i] = true;
           } 

      }
//
       
void _Table_Format::Filter_Min(float min){

      float min_temp;

      for (int i = 0; i < Features.size(); i ++){
          min_temp = 1.0;
          for (int j = 0; j < Samples.size(); j++)
               min_temp = min_temp < Abd[j][i] ? min_temp : Abd[j][i];
              
          if (min_temp < min) 
             Min_filtered[i] = true;
             
          } 

      }
// 

void _Table_Format::Filter_Ave(float ave){

      float samplesum, sampleave;      

      for (int i = 0; i < Features.size(); i ++){                                                                                                                               
          samplesum = 0;                                                                                                                                                      
          for (int j = 0; j < Samples.size(); j++)                                                                                                                          
              samplesum = samplesum + Abd[j][i];                                                                                                               
                                
              sampleave = samplesum/Samples.size();

          if (sampleave < ave)
             Ave_filtered[i] = true;
                                                
            }                                                                                                                                                                
                                                                                                                                                                             
      }
//
 
void _Table_Format::Filter_Zero(float zero){
     
      for (int i = 0; i < Features.size(); i ++){                                                                                                                               
          int samplesum = 0;                                                                                                                                         
          for (int j = 0; j < Samples.size(); j++)    
              if (Abd[j][i] > 0)
                 samplesum ++;                                                                                                                             
           
          if ((float) samplesum / (float) Samples.size() < zero)
             Zero_filtered[i] = true;
            }                                                                                                                                                                            
      }

void _Table_Format::Filter_Empty(){
     for (int i = 0; i < Features.size(); i ++){                                                                                                                               
          bool is_empty = true;                                                                                                                                         
          for (int j = 0; j < Samples.size(); j++)    
              if (Abd[j][i] > 0)
                 is_empty = false;                                                                                                                            
           
             Empty_filtered[i] = is_empty;
            }                                                                                                                                                                              
      }

float _Table_Format::Calc_Dist_E(int sam_m, int sam_n){

     float abd_m_norm[Features.size()], abd_n_norm[Features.size()];
     float sum_m = 0;
     float sum_n = 0;
     float     f = 0;
     

     for (int i = 0; i < Features.size(); i++){
         abd_m_norm[i] = Abd[sam_m][i];
         abd_n_norm[i] = Abd[sam_n][i];
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }

     
     for (int i = 0; i < Features.size(); i++)
         f += pow(abd_m_norm[i] - abd_n_norm[i], 2);

     return sqrt(f);
     }

float _Table_Format::Calc_Dist_Cos(int sam_m, int sam_n){

     float f_m_sum = 0;
     float f_n_sum = 0;
     float f = 0;

     for (int i = 0; i < Features.size(); i++){
         
         float fm = Abd[sam_m][i];
         float fn = Abd[sam_n][i];

         f += fm * fn;

         f_m_sum += fm * fm;
         f_n_sum += fn * fn;
         }
     
     float ff = sqrt(f_m_sum) * sqrt(f_n_sum);
     if (ff == 0) return 0;
     
     f = f / ff;
     
     return (1 - f);
     }

float _Table_Format::Calc_Dist_JSD(int sam_m, int sam_n){

      float abd_m_norm[Features.size()];
      float abd_n_norm[Features.size()];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < Features.size(); i++){
         if (Abd[sam_m][i] > 0) abd_m_norm[i] = Abd[sam_m][i];
         else abd_m_norm[i] = 0;
         
         if (Abd[sam_n][i] > 0) abd_n_norm[i] = Abd[sam_n][i];
         else abd_n_norm[i] = 0;
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }
      
      //calc
      float dkl_m = 0;
      float dkl_n = 0;
      
      for (int i = 0; i < Features.size(); i ++){
          
          if ((abd_m_norm[i] == 0) && (abd_n_norm[i] == 0)) continue;
          
          float abd_q = (abd_m_norm[i] +  abd_n_norm[i]) / 2;
          
          if (abd_m_norm[i] > 0)
             dkl_m += abd_m_norm[i] * log(abd_m_norm[i] / abd_q);
          
          if (abd_n_norm[i] > 0)
             dkl_n += abd_n_norm[i] * log(abd_n_norm[i] / abd_q);

          }
          
      return sqrt((dkl_m + dkl_n)/2.0);
     }

float _Table_Format::Calc_Dist_Bray_Curtis(int sam_m, int sam_n){

      float abd_m_norm[Features.size()];
      float abd_n_norm[Features.size()];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < Features.size(); i++){
         if (Abd[sam_m][i] > 0) abd_m_norm[i] = Abd[sam_m][i];
         else abd_m_norm[i] = 0;
         
         if (Abd[sam_n][i] > 0) abd_n_norm[i] = Abd[sam_n][i];
         else abd_n_norm[i] = 0;
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }
      
      //calc
      float sum = 0;
      float diff = 0;
            
      for (int i = 0; i < Features.size(); i ++){
          sum += (abd_m_norm[i] + abd_n_norm[i]);
          float a_diff = abd_m_norm[i] - abd_n_norm[i];
          if (a_diff < 0) a_diff = a_diff * (-1.0);
          diff += a_diff;
          }
      if (sum <= 0) return 1;
      return diff / sum;
     }

float _Table_Format::Calc_Corr_S(vector <float> sam_m, vector <float> sam_n){

      float sumValue;
      float partValue, rho;
      int sample_num1, sample_num2, rank_tempm, rank_tempn;
      int mline = sam_m.size();
      int nline = sam_n.size();
      
      if (mline!=nline){
          cerr << "Error: Different size of two vectors, check your input files! "<<endl;
          return 0;
                        }
      float abdm[mline], abdn[mline];
      int rankn[mline], rankm[mline];
      int rankn1[mline], rankm1[mline];
      float rankn2[mline], rankm2[mline];

      for (int i = 0; i < mline; i++){
          abdm[i] = sam_m[i];
          abdn[i] = sam_n[i];
          rankm[i] = i;
          rankn[i] = i;
          }
      
      BubbleSort(abdm, rankm, mline);
      BubbleSort(abdn, rankn, mline);
        
      for (int i = 0; i < mline; i++){
          for (int j = 0; j < mline; j++){
              if (sam_m[i] == abdm[j] && rankm[j] == i) rankm1[i]=j;
              if (sam_n[i] == abdn[j] && rankn[j] == i) rankn1[i]=j;
              }
          }
    
      for (int i = 0; i < mline; i++){

          rank_tempm = 0;
          rank_tempn = 0;
          sample_num1 = 0;
          sample_num2 = 0; 

          for (int j = 0; j < mline; j++){
              if (sam_m[i] == sam_m[j]) {
                 rank_tempm += rankm1[j];     
                 sample_num1 ++;
                 }
              if (sam_n[i] == sam_n[j]) {
                 rank_tempn += rankn1[j];
                 sample_num2 ++;
                 } 
              }
          rankm2[i] = (float) rank_tempm / sample_num1;
          rankn2[i] = (float) rank_tempn / sample_num2;
          
          }
      
      sumValue = 0.0; 

      for (int i = 0; i < mline; i++)
          sumValue += (rankm2[i] - rankn2[i])*(rankm2[i] - rankn2[i]); 
      
      partValue = mline*(mline*mline - 1);

      rho = 1 - ((6*sumValue)/partValue);
      
      return rho;

      }

float _Table_Format:: Calc_Corr_P(vector <float> sam_m, vector <float> sam_n){

      float ave_m = 0;
      float ave_n = 0;
      int mline = sam_m.size();
      int nline = sam_n.size();
      
      if (mline!=nline){
          cerr << "Error: Different size of two vectors, check your input files! "<<endl;
          return 0;
                        }
      
      for (int i = 0; i < mline; i++){
          
          ave_m += sam_m[i];
          ave_n += sam_n[i];

          }
      
      ave_m = ave_m / (float) mline;
      ave_n = ave_n / (float) nline;


      float deno_m = 0;
      float deno_n = 0;
 
      for (int i = 0; i < mline; i++){
           
          deno_m += pow(sam_m[i] - ave_m, 2);
          deno_n += pow(sam_n[i] - ave_n, 2);

          }

      float deno = sqrt(deno_m) * sqrt(deno_n);
      
      if (deno == 0) return 0;

      float nume = 0;

      for (int i = 0; i < mline; i++)
          nume += (sam_m[i] - ave_m)*(sam_n[i] - ave_n);


      return nume/deno;
      }

float _Table_Format::Calc_Corr_P(int sam_m, int sam_n){

      return Calc_Corr_P(Abd[sam_m], Abd[sam_n]);

      }

float _Table_Format::Calc_Corr_S(int sam_m, int sam_n){

           
      return Calc_Corr_S(Abd[sam_m], Abd[sam_n]);

      }


void _Table_Format::Calc_Dist_Matrix(const char * outfilename, int metrics, int coren, bool is_sim){
     
      //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < Samples.size() - 1; i ++)
        for (int j = i + 1; j < Samples.size(); j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
      
      ofstream outfile(outfilename, ofstream::out);                                 
      if (!outfile){                                                                
         cerr << "Error: Cannot open output file: " << outfilename << endl;               
         return;                                                                   
         } 
    
      //calc dist 
      vector <float>  dist_matrix;
      for (long i = 0; i < iter; i ++)
        dist_matrix.push_back(0);
      
      omp_set_num_threads(coren);
      
      #pragma omp parallel for schedule(dynamic, 1)
      for (long i = 0; i < iter; i ++){
        
         long m = order_m[i];
         long n = order_n[i];
         long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2;    
        
         switch (metrics){
                case 0: dist_matrix[p] = Calc_Dist_Cos(m, n); break;
                case 1: dist_matrix[p] = Calc_Dist_E(m, n); break;
                case 2: dist_matrix[p] = Calc_Dist_JSD(m, n); break;
                case 3: dist_matrix[p] = Calc_Dist_Bray_Curtis(m, n); break;
                default: dist_matrix[p] = Calc_Dist_Cos(m, n); break;
                }
        }
      
       for (int i = 0; i < Samples.size(); i++)
          outfile <<"\t" <<Samples[i];
       outfile <<endl;
      
       for (int i = 0; i < Samples.size(); i++){
          outfile <<Samples[i];
          for (int j = 0; j < Samples.size(); j++){
              long m = (i <= j) ? i : j;
              long n = (i > j) ? i : j;
              long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2; 
              if (is_sim) {
                          if (i == j) outfile << "\t" << 1;
                          else outfile << "\t" << 1 - dist_matrix[p];
                          }
              else {
                   if (i == j) outfile << "\t" << 0;
                   else outfile << "\t" << dist_matrix[p];
                   }
              }
          outfile <<endl;
          }

      outfile.close();
      outfile.clear();
      }

void _Table_Format::Calc_Corr_Matrix(const char * outfilename, int metrics, int coren){
     
      //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < Samples.size() - 1; i ++)
        for (int j = i + 1; j < Samples.size(); j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
      
      ofstream outfile(outfilename, ofstream::out);                                 
      if (!outfile){                                                                
         cerr << "Error: Cannot open output file: " << outfilename << endl;               
         return;                                                                   
         } 
    
      //calc dist 
      vector <float>  corr_matrix;
      for (long i = 0; i < iter; i ++)
        corr_matrix.push_back(0);
      
      omp_set_num_threads(coren);
      
      #pragma omp parallel for schedule(dynamic, 1)
      for (long i = 0; i < iter; i ++){
        
         long m = order_m[i];
         long n = order_n[i];
         long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2;    
        
         if (metrics == 0)
                    corr_matrix[p] = Calc_Corr_S(m, n);
         else corr_matrix[p] = Calc_Corr_P(m, n);

        }
      
       for (int i = 0; i < Samples.size(); i++)
          outfile <<"\t" <<Samples[i];
       outfile <<endl;
      
       for (int i = 0; i < Samples.size(); i++){
          outfile <<Samples[i];
          for (int j = 0; j < Samples.size(); j++){
              long m = (i <= j) ? i : j;
              long n = (i > j) ? i : j;
              long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2; 
              if (i == j) outfile << "\t" << 1;
                          else outfile << "\t" <<corr_matrix[p];
                          }             
          outfile <<endl;
          }

      outfile.close();
      outfile.clear();
      }

void _Table_Format::BubbleSort(float *array1, int *rank1, int len){

     float temp;
     int   temp_rank;

     for (int i = 0; i < len-1; i++){
         for (int j = 0; j < len-1-i; j++){
             if (array1[j] > array1[j+1]){
                temp        = array1[j];
                array1[j]   = array1[j+1];
                array1[j+1] = temp;
                temp_rank  = rank1[j];
                rank1[j]   = rank1[j+1];
                rank1[j+1] = temp_rank;
                }
             }
         }
     }
     
float _Table_Format::Get_Abd_By_Feature(unsigned int s, string a_feature){
            
      if (s >= Samples.size())
         return 0;
      if (Feature_hash.count(a_feature) == 0)
         return 0;
      if (Check_Filter(a_feature))
         return 0;
           
      return Abd[s][Feature_hash[a_feature]];
      }
      
float _Table_Format::Get_Abd_By_Order(unsigned int s, unsigned int i){
            
      if (s >= Samples.size() || (i >= Features.size()))
         return 0;
      return Abd[s][i];
      }      

vector <float> _Table_Format::Get_Abd(unsigned int s){
       
       if (s < Samples.size())
          return Abd[s];
       }

void _Table_Format::Init_Filter(){
     Max_filtered = new bool[Features.size()];
     Min_filtered = new bool[Features.size()];
     Ave_filtered = new bool[Features.size()];
     Zero_filtered = new bool[Features.size()];
     Empty_filtered = new bool[Features.size()];
     
     for (int i = 0; i < Features.size(); i ++){
         Max_filtered[i] = false;
         Min_filtered[i] = false;
         Ave_filtered[i] = false;
         Zero_filtered[i] = false;
         Empty_filtered[i] = false;
         }
     }

bool _Table_Format::Check_Filter(string a_feature){
     
     if (Feature_hash.count(a_feature) == 0) return true;
     int feature_index = Feature_hash[a_feature];    
     return (Max_filtered[feature_index] || Min_filtered[feature_index] || Ave_filtered[feature_index] || Zero_filtered[feature_index] || Empty_filtered[feature_index]);
     }
    
void _Table_Format::Build_Feature_Hash(){
         
     for (int i = 0; i < Features.size(); i ++)
         Feature_hash[Features[i]] = i;
         }
     
// #############################################################################
#endif /* table_format_hpp */
