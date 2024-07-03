// Updated at July 2, 2024
// Updated by Xiaoquan Su, Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>

#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <ctime>

#include "utility.h"
#include "version.h"
#include "comp_sam_func.h"

using namespace std;

char Ref_db;

string Listfilename;
string Listprefix;

string Queryfile1;
string Queryfile2;

string Tablefilename;
string Outfilename;

int Coren = 0;

int Dist_matrix = 0; //0: cos 1: eu;

bool Is_sim; //true: sim, false: dist;

int Mode = 0; //0: single, 1: multi, 2: multi_table
bool Reversed_table = false;
bool Is_heatmap;
int Cluster = 2;

int printhelp(){
    
    cout << "PM-comp-func version : " << Version << endl;
    cout << "\tCompute the Meta-Storms functional distance/similarity among samples" << endl;
    cout << "Usage: " << endl;
    cout << "PM-comp-func [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Two samples path for single sample comparison" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list table for multi-sample comparison" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input KO count table (*.KO.Count) for multi-sample comparison" << endl;
    cout << "\t  -R If the input table is reversed, T(rue) or F(alse), default is false [Optional for -T]" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file, default is to output on screen" << endl;
    cout << "\t  -d Output format, distance (T) or similarity (F), default is T" << endl;
    cout << "\t  -P (upper) Print heatmap and clusters, T(rue) or F(alse), default is F" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -M (upper) Distance Metric, 0: Hierarchial Meta-Storms; 1: Cosine; 2: Euclidean; 3: Jensen-Shannon, 4: Bray-Curtis; default is 0" << endl;
    cout << "\t  -c Cluster number, default is 2 [Optional for -P]" << endl;
    cout << "\t  -t Number of thread, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
void Parse_Para(int argc, char * argv[]){
    
    Ref_db = 'G';
    
    Coren = 0;
    Mode = true; //default is single;
    
    Is_sim = false;
    Is_heatmap = false;
    Dist_matrix = 0;
    Reversed_table = false;
    
    int i = 1;
    if (argc == 1) 
	printhelp();

    while(i < argc){
     	if (argv[i][0] != '-') {
        	cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
         	exit(0);
        }
           
        switch(argv[i][1]){
        	case 'D': Ref_db = argv[i+1][0]; break;
	       	case 'i': Queryfile1 = argv[i+1]; Queryfile2 = argv[i+2]; i++; Mode = 0; break;
            	case 'l': Listfilename = argv[i+1]; Mode = 1; break;
            	case 'p': Listprefix = argv[i+1]; break;
           	case 'T': Tablefilename = argv[i+1]; Mode = 2; break;
            	case 'R': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Reversed_table = true; break;
            	case 'o': Outfilename = argv[i+1]; break;
                                        
            	case 'd': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_sim = true; break;
            	case 'M': Dist_matrix = atoi(argv[i+1]); break;
            	case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_heatmap = true; break;
                case 'c': Cluster = atoi(argv[i+1]); break;

            	case 't': Coren = atoi(argv[i+1]); break;         
            	case 'h': printhelp(); break;
            	default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
            	
	}

        i += 2;
    }
    if ((Dist_matrix > 4) || (Dist_matrix < 0)){
	    cerr << "Warning: Distance matrix (-M) must be 0-4, change to default (0)" << endl;
	    Dist_matrix = 0;
    }     

    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);    
    if ((Coren <= 0) || (Coren > max_core_number)){
	    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
    	    Coren = max_core_number;
    }
}

void Output_Matrix(const char * outfilename, unsigned long n, float * results, bool is_sim, vector <string> sam_name){
     
     FILE * outfile = fopen(outfilename, "w");
     if (outfile == NULL){
                   cerr << "Error: Cannot open output file : " << outfilename << endl;
                   return; 
                   }
     
     for(unsigned long i = 0; i < n; i ++)
             fprintf(outfile, "\t%s", sam_name[i].c_str());
     fprintf(outfile, "\n");
     
     for(unsigned long i = 0; i < n; i ++){
             fprintf(outfile, "%s", sam_name[i].c_str());
             
             for (unsigned long j = 0; j < n; j ++){                
                 
                 unsigned long ii = (i <= j) ? i : j;
                 unsigned long jj = (i >= j) ? i : j;
                 if (is_sim){                 
                    if (ii == jj) fprintf(outfile, "\t1.0");
                    else fprintf(outfile, "\t%f", results[ii *  n + jj - (1 + ii + 1) * (ii + 1) / 2]);
                    }
                 else {
                      if (ii == jj) fprintf(outfile, "\t0.0");
                      else fprintf(outfile, "\t%f", 1-results[ii *  n + jj - (1 + ii + 1) * (ii + 1) / 2]);
                      }
                 }
             fprintf(outfile, "\n");
             }
                 
     fclose(outfile);
     }
     
void Single_Comp(){
     
    //_Comp_Tree_Func comp_tree_func(Ref_db);
    _Comp_Tree_Func comp_tree_func;
    
    float * abd_1 = new float [comp_tree_func.Get_GeneN()];
    float * abd_2 = new float [comp_tree_func.Get_GeneN()];
    
    //cout << comp_tree_func.Load_Gene_Count(Queryfile1.c_str(), abd_1) << " loaded" << endl;
    //cout << comp_tree_func.Load_Gene_Count(Queryfile2.c_str(), abd_2) << " loaded" << endl;
    comp_tree_func.Load_Gene_Count(Queryfile1.c_str(), abd_1);
    comp_tree_func.Load_Gene_Count(Queryfile2.c_str(), abd_2);
    
    float sim = comp_tree_func.Calc_sim(abd_1, abd_2, Dist_matrix);
       
    if (Is_sim)
    	cout << sim << endl;
    else cout << 1.0 - sim << endl;
    
    }

void Multi_Comp(){
    
    //_Comp_Tree_Func comp_tree_func(Ref_db);
    _Comp_Tree_Func comp_tree_func;
         
     //load list
    vector <string> sam_name;
    vector <string> file_list;
    
    unsigned long file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
            
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree_func.Get_GeneN()];
        //cout << comp_tree_func.Load_Gene_Count(file_list[i].c_str(), Abd[i]) << " KOs in file " << i + 1 << endl;
        comp_tree_func.Load_Gene_Count(file_list[i].c_str(), Abd[i]);
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    unsigned long iter = (unsigned long) file_count * (file_count-1) / 2;
    float * results = (float *) malloc(sizeof(float) * iter); //sim_matrix
        
    //openmp    
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for
    for (unsigned long i = 0; i < file_count; i++) {
    	for(unsigned long j = i + 1; j < file_count; j++) {

       	 	results[(i * file_count + j - (i + 2) * (i + 1) / 2)] = comp_tree_func.Calc_sim(Abd[i], Abd[j], Dist_matrix);
        }
    }
    
    Output_Matrix(Outfilename.c_str(), file_count, results, Is_sim, sam_name);
    
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

void Multi_Comp_Table(const char * infilename, bool Reversed_table){
    
	//time_t t1 = time(0);

    //_Comp_Tree_Func comp_tree_func(Ref_db);
    _Comp_Tree_Func comp_tree_func;
            
    vector<string> samples;
    vector<float *> Abd; 
    unsigned long file_count;
    //load abd
    if(!Reversed_table){
    	file_count = (unsigned long) comp_tree_func.Load_Gene_Count(infilename, Abd, samples);
    }
    else {
	file_count = (unsigned long) comp_tree_func.Load_Gene_Count_Reverse(infilename, Abd, samples);
    }
    cout << file_count << " files loaded" << endl;

    unsigned long results_size = file_count * (file_count - 1) / 2;
    float * results = (float *) malloc(sizeof(float) * results_size);

    omp_set_num_threads(Coren);
    #pragma omp parallel for schedule(dynamic, 1)
    for(unsigned long i = 0; i < file_count; i++) {
	for(unsigned long j = i + 1; j < file_count; j++) {
		results[i*file_count+j-(2+i)*(i+1)/2] = comp_tree_func.Calc_sim(Abd[i], Abd[j], Dist_matrix);
	}
    }
    
	//time_t t2 = time(0);
	//cout << "time: " << endl << t2 - t1 << endl;

    Output_Matrix(Outfilename.c_str(), file_count, results, Is_sim, samples);
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];

    if (Is_heatmap){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Heatmap.R -d %s -o %s", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".heatmap.pdf").c_str());
                    system(command);
                    sprintf(command, "Rscript %s/Rscript/PM_Hcluster.R -d %s -o %s -c %d", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".clusters.pdf").c_str(), Cluster); 
                    system(command);
                    }  
}

int main(int argc, char * argv[]){
    
    //test();
    
    Parse_Para(argc, argv); 
    
    //debug
    //cout << "Method " << Dist_matrix << endl;
                  
    switch (Mode) {
       case 0: Single_Comp(); break;
       case 1: Multi_Comp(); break;
       case 2:{
            Multi_Comp_Table(Tablefilename.c_str(), Reversed_table); 
            break;
            }
            
       default: break;
       }

    return 0;
    }
