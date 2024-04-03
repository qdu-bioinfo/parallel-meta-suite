//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su

#ifndef _OPT_H
#define _OPT_H

#include "utility.h"

using namespace std;

string Dbseqfilename;
string Dbtaxafilename;
string Sampleseqfilename;

string Listfilename;
string Listprefix;

string Outfilename;

int Taxa_mode = 0; //0: HWL; 1: HWL + LCA;
int Run_mode = -1; //0: single; 1: Multi
int Strand_mode = 0;//0: Forward; 1: both forward and reverse;
bool Raw = false; //debug mode

float Sim = 0.98;
float Consistency = 0.0;
int K = 15;

int Coren = 0; //# of core

int printhelp(){
    
    cout << "PM profiler version : " << Version << endl;
    cout << "\tThe taxonomy profiler" << endl;
    cout << "Usage: " << endl;
    cout << "PM-profiler [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -d Reference database sequence file" << endl;
    cout << "\t  -m Reference database taxonomy file" << endl;
    
    cout << "\t  -i Input sequence file for single sample [conflicts with -l and -p]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input sequence file list for multiple samples [conflicts with -i]" << endl;
    cout << "\t  -p List files path prefix [Optional for -l]" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output name" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -M (upper) taxonomy mode, 0: HWL, 1: HWL + LCA, default is 0" << endl;
    cout << "\t  -K (upper) K-mer size (int), default is 15 , within 11-19" << endl;
    cout << "\t  -S (upper) search forward and reverse strand,default is F" << endl;
    cout << "\t  -s Similarity threshold (float value 0-1), default is 0.98" << endl;
    cout << "\t  -c HWL's Taxonomy Consistency threshold (float value 0-1), default is 0.0" << endl;
    cout << "\t  -t Number of thread, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };


void Parse_Para(int argc, char * argv[]){
    
    int i = 1;
      
      if (argc ==1)
        printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };
         switch(argv[i][1]){
                            case 'd': Dbseqfilename = argv[i+1]; break;
                            case 'm': Dbtaxafilename = argv[i+1];  break;
                            
                            case 'i': if (Run_mode == 1){
                                            cerr <<"Error: -i conflicts with -l"<< endl;
                                            exit(1);
                                       }
                                      Sampleseqfilename = argv[i+1];
                                      Run_mode = 0;
                                      break;
                            case 'l': if (Run_mode == 0){
                                            cerr <<"Error: -l conflicts with -i"<< endl;
                                            exit(1);
                                      }
                                      Listfilename = argv[i+1];
                                      Run_mode = 1;
                                      break;
                            case 'p': Listprefix = argv[i+1]; break;

                            case 'c': Consistency = atof(argv[i+1]); break;
                            
                            case 'o': Outfilename = argv[i+1]; break;
                            
                            case 'M': Taxa_mode = atoi(argv[i+1]); break;

                            case 's': Sim = atof(argv[i+1]); break;

                            case 'S': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Strand_mode = 1; break;

                            case 'K': if (atof(argv[i+1]) > 19 || atof(argv[i+1]) < 11){
                                            cerr <<"Error: The size of K-mer should be within 11-19"<< endl;
                                            exit(1);
                                      }
                                      K = atof(argv[i+1]); 
                                      break;
                            case 'R': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Raw = true; break;
                            
                            case 't': Coren = atoi(argv[i+1]); break;

                            case 'h': printhelp(); break;

                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
                            }
         i+=2;
         }
    
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    //int max_core_number = 160;
    
    if ((Coren <= 0) || (Coren > max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = max_core_number;
                    }
   
    }

#endif
