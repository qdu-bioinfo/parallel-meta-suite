// Updated at July 2, 2024
// Updated by Xiaoquan Su
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// version 3.1 or above with _Table_Format
// version 3.4.1 or above wityh _OTU_Parser
// Added unweighted based on 3.4.1
// Notes: modify default db='G'

#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>

#include "comp.h"

using namespace std;

char Ref_db;

string Listfilename;
string Listprefix;
string Queryfile1;
string Queryfile2;

string Tablefilename;
string Outfilename;

int Coren = 0;

bool Is_cp_correct; //
bool Is_sim; //true: sim, false: dist;
bool Is_heatmap;
int Cluster = 2;

//bool Is_weight;
int Dist_metric = 0; //0: MS; 1: MS-uw; 2: cos 3: eu; 4: JSD; 5.Bray Curtis

int Mode = 0; //0: single, 1: multi_list, 2: multi_table
//bool Reversed_table = true;
 
int printhelp(){
    
    cout << "Comp-taxa version : " << Version << endl;
    cout << "\tCompute the Meta-Storms similarity/distance among samples" << endl;
    cout << "Usage: " << endl;
    cout << "PM-comp-taxa [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Two samples path for single sample comparison" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list for multi-sample comparison" << endl;
    cout << "\t  -p List files path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count) for multi-sample comparison" << endl;
    //cout << "\t  -R If the input table is reversed, T(rue) or F(alse), default is false [Optional for -T]" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file, default is to output on screen" << endl;
    cout << "\t  -d Output format, distance (T) or similarity (F), default is T" << endl;
    cout << "\t  -P (upper) Print heatmap and clusters, T(rue) or F(alse), default is F" << endl;
    
    cout << "\t[Other options]" << endl;
    //cout << "\t  -w weighted or unweighted, T(rue) or F(alse), default is T" << endl;
    cout << "\t  -M (upper) Distance Metric, 0: Meta-Storms; 1: Meta-Storms-unweighted; 2: Cosine; 3: Euclidean; 4: Jensen-Shannon; 5: Bray-Curtis, default is 0" << endl;
    cout << "\t  -r rRNA copy number correction, T(rue) or F(alse), default is T" << endl;
    cout << "\t  -c Cluster number, default is 2 [Optional for -P]" << endl;
    cout << "\t  -t Number of thread, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
void Parse_Para(int argc, char * argv[]){
    
    Ref_db = 'G';
    
    Coren = 0;
    Mode = 0; //default is single;
    
    Is_cp_correct = true;
    Is_sim = false;
    Is_heatmap = false;
    //Is_weight = true;
    Dist_metric = 0;
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D': Ref_db = argv[i+1][0]; break;
                 
                            case 'i': Queryfile1 = argv[i+1]; Queryfile2 = argv[i+2]; i++; Mode = 0; break;
                            
                            case 'l': Listfilename = argv[i+1]; Mode = 1; break;
                            case 'p': Listprefix = argv[i+1]; break;
                            
                            case 'T': Tablefilename = argv[i+1]; Mode = 2; break;
                            //case 'R': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Reversed_table = false; break;
                            
                            case 'o': Outfilename = argv[i+1]; break;
                            
                            //case 'w': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_weight = false; break;
                            case 'M': Dist_metric = atoi(argv[i+1]); break; 
                            case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;
                            case 'd': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_sim = true; break;
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_heatmap = true; break;
                            case 'c': Cluster = atoi(argv[i+1]); break;
                            
                            case 't': Coren = atoi(argv[i+1]); break;         
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = max_core_number;
                    } 
    if (Cluster <= 0){
                cerr << "Warning: cluster number must be larger than 0, change to default (2)" << endl;
                }
    }

void Output_Matrix(const char * outfilename, int n, vector <float> * sim_matrix, bool is_sim, vector <string> sam_name){
     
     ofstream outfile(outfilename, ofstream::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file : " << outfilename << endl;
                   return; 
                   }
     
     //Label
     for(int i = 0; i < n; i ++)
             outfile << "\t" << sam_name[i];
     outfile << endl;
     
     for(int i = 0; i < n; i ++){
             outfile << sam_name[i];
             
             for (int j = 0; j < n; j ++){                
                 
                 long ii = (i <= j) ? i : j;
                 long jj = (i >= j) ? i : j;
                 long p = ii * (long) n + jj - (1 + ii + 1) * (ii + 1) / 2;
                                  
                 if (is_sim){                 
                    if (ii == jj) outfile << "\t" << 1.0;
                    else outfile << "\t" << (*sim_matrix)[p];
                    }
                 else {
                      if (ii == jj) outfile << "\t" << 0.0;
                      else outfile << "\t" << 1.0 - (*sim_matrix)[p];
                      }
                 }
             outfile << endl;
             }
                 
     outfile.close();
     outfile.clear();
     }

void Single_Comp(){          
     
     _Comp_Tree comp_tree(Ref_db);
     
     float * Abd1 = new float [comp_tree.Get_LeafN()];
     float * Abd2 = new float [comp_tree.Get_LeafN()];
     
     cout << comp_tree.Load_abd(Queryfile1.c_str(), Abd1, Is_cp_correct) << " OTUs in file " << 1 << endl;
     cout << comp_tree.Load_abd(Queryfile2.c_str(), Abd2, Is_cp_correct) << " OTUs in file " << 2 << endl;
     
     float sim = comp_tree.Calc_sim(Abd1, Abd2, Dist_metric);
     
     if (Is_sim) 
                 cout << sim << endl;
     else 
                 cout << 1.0 - sim << endl;
     
     }

void Multi_Comp(){
    
    _Comp_Tree comp_tree(Ref_db);
         
     //load list
    vector <string> sam_name;
    vector <string> file_list;
    
    int file_count = 0;
    file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
            
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree.Get_LeafN()];
        cout << comp_tree.Load_abd(file_list[i].c_str(), Abd[i], Is_cp_correct) << " OTUs in file " << i + 1 << endl;
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
        
    vector <float>  sim_matrix;
    for (long i = 0; i < iter; i ++)
        sim_matrix.push_back(0);
    
    //openmp    
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < iter; i ++){
        
        long m = order_m[i];
        long n = order_n[i];
        long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;    
        sim_matrix[p] = comp_tree.Calc_sim(Abd[m], Abd[n], Dist_metric);
        }
    
    Output_Matrix(Outfilename.c_str(), file_count, &sim_matrix, Is_sim, sam_name);
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];
             
    if (Is_heatmap){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Heatmap.R -d %s -o %s", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".heatmap.pdf").c_str());
                    system(command);
                    sprintf(command, "Rscript %s/Rscript/PM_Hcluster.R -d %s -o %s -c %d", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".clusters.pdf").c_str(), Cluster); 
                    system(command);
                    }       
     };

void Multi_Comp_Table(_Table_Format abd_table){
    
    _Comp_Tree comp_tree(Ref_db);
    
    int file_count = abd_table.Get_Sample_Size();
         
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree.Get_LeafN()];
        cout << comp_tree.Load_abd(&abd_table, Abd[i], i, Is_cp_correct) << " OTUs in file " << i + 1 << endl;
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
        
    vector <float>  sim_matrix;
    for (long i = 0; i < iter; i ++)
        sim_matrix.push_back(0);
    
    //openmp
        
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < iter; i ++){
        
        long m = order_m[i];
        long n = order_n[i];
        long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;
        
        sim_matrix[p] = comp_tree.Calc_sim(Abd[m], Abd[n], Dist_metric);
        }
    
    Output_Matrix(Outfilename.c_str(), file_count, &sim_matrix, Is_sim, abd_table.Get_Sample_Names());
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];
             
    if (Is_heatmap){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Heatmap.R -d %s -o %s", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".heatmap.pdf").c_str());
                    system(command);
                    sprintf(command, "Rscript %s/Rscript/PM_Hcluster.R -d %s -o %s -c %d", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".clusters.pdf").c_str(), Cluster); 
                    system(command);
                    }       
     };

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);         
                  
    switch (Mode){
           case 0: Single_Comp(); break;
           case 1: Multi_Comp(); break;
           case 2:{
                   //_Table_Format table(Tablefilename.c_str(), Reversed_table);
                   _Table_Format table(Tablefilename.c_str());
                   Multi_Comp_Table(table); 
                   break;
                   }
           default: break;
           }
    return 0;
    }
