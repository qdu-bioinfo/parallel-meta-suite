// Updated at April 2, 2024
// Updated by Haobo Shi
// version 3.7 - 3.7.2 
// Update profiler,remove pair-end mode

// Updated at July 31, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <sys/dir.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include "utility.h"
#include "version.h"
#include "db.h"

using namespace std;

#ifndef INIT_H
#define INIT_H

#define LEVEL 7

//Parameters
class _Para{
      public:
             _Para(){
                     This_path = Check_Env();
                     //Align_exe_name = This_path + "/Aligner/bin/bowtie2-align-s";
                     Align_exe_name = "vsearch"; 
                     Length_filter = 0;
                     Core_number = 0;
                     Type = -1;
                     Out_path = "./Result";
                     profiler = 0;//Default profiler is vsearch
                 
             
                     Is_format_check = false;  
                     Is_paired = false; 
                     Is_func = true;
                     
                    //newadd parameters
                     Is_denoised = 'T';
                     Is_nonchimeras = 'T';
                    //Similarity when searching database, default is 0.97 
                     db_similarity = 0.99;
                     }
             
             string This_path;
             string Align_exe_name;
             string Infilename;
             string Infilename2;
             string Groupfilename;
             string Listfilename;
             string Out_path;
    
             int Type; //0: 16S 1: shotgun
             int Format; //0: fasta 1: fastq;
             int Length_filter;
             int Core_number;
             int profiler;//Profiler select
             bool Is_format_check;
             bool Is_paired;
             bool Is_func;
    		 
    		 //newadd
    		 char Is_denoised;
    		 char Is_nonchimeras;
    		 double db_similarity;
    		 
             _PMDB Database;
             };

int Print_Help(){
    
    cout << "Parallel-Meta Suite version: " << Version << endl;
    cout << "\tMicrobiome sample profiling" << endl;
    cout << "Usage:" << endl;
    cout << "PM-parallel-meta [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -m Input single sequence file (Shotgun) [Conflicts with -r]" << endl;
    cout << "\t  or" << endl;
	cout << "\t  -r Input single sequence file (rRNA targeted) [Conflicts with -m]" << endl;
	//cout << "\t  -R (upper) Input paired sequence file [Optional for -r, Conflicts with -m]" << endl;
    
    //newadd
    cout << "\t  -v ASV denoising, T(rue) or F(alse), default is T [optional for -i]" << endl;
    cout << "\t  -c Chimera removal, T(rue) or F(alse), default is T [optional for -i]" << endl;
    cout << "\t  -d Sequence alignment threshold (float value 0-1), default is 0.99 for ASV enabled and 0.97 for ASV disabled (-v F) [optional for -i]" << endl;
    cout << endl;
    
    cout << "\t[Output options]" << endl;
	cout << "\t  -o Output path, default is \"Result\"" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -k Sequence format check, T(rue) or F(alse), default is F" << endl;
    cout << "\t  -P (upper) Profiler selection. V (upper) or v for Vsearch, P (upper) or p for PM-profiler, default is Vsearch"  << endl;
	cout << "\t  -L (upper) rRNA length threshold of rRNA extraction. 0 is disabled, default is 0 [Optional for -m, Conflicts with -r]" << endl;
    cout << "\t  -f Functional analysis T(rue) or F(alse), default is T" << endl;
    cout << "\t  -t Number of thread, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    }

int Print_Config(_Para para){
    
    cout << "Reference database is ";
    cout << para.Database.Get_Description() << endl;
    
    cout << "Input sequence type is ";
    if (para.Type == 1) cout << "Metagenomic shotgun sequences" << endl;
    else cout << "rRNA targeted sequences" << endl;
    
    cout << "Input sequence is " << para.Infilename << endl;
    /* delete by Shi Haobo
    if (para.Is_paired)
       cout << "Input pair sequence is " << para.Infilename2 << endl;
    */
    cout << "Functional prediction is ";
    if (para.Is_func) cout << "On" << endl;
    else cout << "Off" << endl;
    
    if(para.Type != 1){
    	//newadd    
    	cout << "ASV denoising is ";
    	if (para.Is_denoised == 'T') cout << "On" << endl;
    	else cout << "Off" << endl;
    
    	cout << "Chimera removal is ";
    	if (para.Is_nonchimeras == 'T') cout << "On" << endl;
    	else cout << "Off" << endl;
	}   

    cout << "Profiler: ";
    cout << (para.profiler == 0 ? "Vsearch" : "PM-profiler") <<endl;

    cout << "Sequence alignment threshold is ";
    cout <<  para.db_similarity << endl;
    
    cout << "The thread number is " << para.Core_number << endl << endl;

    return 0;
    
    }

int Print_Report(_Para para, unsigned int seq_count, unsigned int rna_count, unsigned int match_rna_count, int * a_diver, int asv_count){
    
    ofstream outfile((para.Out_path + "/Analysis_Report.txt").c_str(), ofstream::out);
    
    if (!outfile){
                  cerr << "Error: Open Analysis Report file error : " << para.Out_path + "/Analysis_Report.txt" << endl;
                  return 0;
                  }
    
    string taxa_name [LEVEL] = {"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"};
    
    outfile << "Parallel-Meta Suite Analysis Report" << endl;
    outfile << "Version: " << Version << endl;
    
    outfile << "Reference database: ";
    
    outfile << para.Database.Get_Description() << endl;
    
    outfile << "Input type : ";
    if (para.Type == 1) outfile << "Metagenomic shotgun sequences" << endl;
    else outfile << "rRNA targeted sequences" << endl; 
    
    outfile << "Input file : " << para.Infilename << endl;
    /* // delete by Shi Haobo
    if (para.Is_paired){
       outfile << "Input paired file : " << para.Infilename2 << endl;
       }
    */
    if (para.Type == 1){
       outfile << "rRNA extraction length filter : ";
       if (para.Length_filter == 0) outfile << "Off" << endl;
       else outfile << para.Length_filter << endl;
       }
	if (para.Type !=1){
		//denoise chimera similarity
    	outfile << "ASV denoising:";
    	if (para.Is_denoised == 'T') outfile << " Yes" <<endl;
    	else outfile << " No" << endl;
     
    	outfile << "Chimera removal:";
    	if (para.Is_nonchimeras == 'T') outfile << " Yes" <<endl;
    	else outfile << " No" << endl;
	}

    outfile << "Profiler: ";
    outfile << (para.profiler == 0 ? "Vsearch" : "PM-profiler") <<endl;

    outfile << "Sequence alignment threshold: ";
    outfile <<  para.db_similarity <<endl;
	         
    if (para.Type == 1)
    outfile << "Metagenomic sequence number: " << seq_count << endl;
    
    outfile << "rRNA sequence number: " << rna_count << endl;
    
    outfile << "Mapped rRNA sequence number: " << match_rna_count << endl;
    
    if (para.Is_denoised == 'T' && para.Type!=1) outfile << "ASV number: " << asv_count << endl;
    else ;
    outfile << "Alpha diversity (# of taxa):" << endl;
    
    for (int i = 0; i < LEVEL; i ++)
        outfile << taxa_name[i] << "\t" << a_diver[i] << endl;
    
    outfile.close();
    outfile.clear();
    
    return 0;
    }

int Parse_Para(int argc, char * argv[], _Para &para){ //Parse Parameters

    if (argc ==1) 
		Print_Help();
    
    int i = 1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D' : para.Database.Set_DB(argv[i+1][0]);
                                       break; //Default is GG97
                 
                            case 'm' : if (para.Type != -1){
                                                     cerr << "Error: -m conflicts with -r" << endl;
                                                     exit(0);
                                                     }
                                       para.Infilename = argv[i+1];
                                       para.Type = 1;   
									   para.Is_denoised = 'F';
                                       para.Is_nonchimeras = 'F';
                                       para.db_similarity =0.97;                                                                                                       
                                       break;  
                            case 'P': if(argv[i+1][0] == 'P' || argv[i+1][0] == 'p') para.profiler = 1;//Profiler select,edit by Shi Haobo
                                        else if(argv[i+1][0] == 'V' || argv[i+1][0] == 'v') para.profiler = 0;
                                        else cerr << "Error: profiler only input P or V" <<endl;
                                        break;
                            case 'r' : if (para.Type != -1){
                                                     cerr << "Error: -r conflicts with -m" << endl;
                                                     exit(0);
                                                     }
                                       para.Infilename = argv[i+1];  
                                       para.Type = 0;
                                       break;  
                            
                            /*//delete by Shi Haobo 
                            case 'R' : para.Infilename2 = argv[i+1];                                        
                                       para.Is_paired = true;
                                       break;
                            */
                            case 'o' : para.Out_path = argv[i+1]; break; //Default is ./result                      

                            case 'k' : if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T' )) para.Is_format_check = true; 
                                       else if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F' )) para.Is_format_check = false;
                                       break;
                                        
                            case 'L' : para.Length_filter = atoi(argv[i+1]);
                                       if (para.Length_filter < 0){
                                                   cerr << "Error: Length Filter must be equal or larger than 0" << endl;
                                                   exit(0);      
                                                         }
                                       break; //Default is 0
                            case 't' : para.Core_number = atoi(argv[i+1]); break; //Default is Auto                                       
                            case 'f' : if ( (argv[i+1][0] == 'F') || (argv[i+1][0] == 'f') ) para.Is_func = false; break; //Default is On
                            case 'h' : Print_Help(); break;
							//newadd
							case 'v' : if(para.Type == 1){
											para.Is_denoised = 'F';
									   }else{
									   		if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) para.Is_denoised = 'T';  //Default is true
									   		else if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) para.Is_denoised = 'F';
									   		if(para.Is_denoised == 'F') para.db_similarity =0.97;
									   }
									   break; 
							case 'c' : if(para.Type == 1){
											para.Is_nonchimeras = 'F';
									   }else{
											if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) para.Is_nonchimeras = 'T';  //Default is true
									   		else if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) para.Is_nonchimeras = 'F';
									   }
									   break;	
							case 'd' : para.db_similarity = atof(argv[i+1]);
									   break; 	   
                            default : printf("Error: Unrec argument %s\n", argv[i]); Print_Help(); break; 
                            }
         i+=2;
         }
    //check func
    if (!para.Database.Get_Is_Func()) para.Is_func = false;
    
    //check input
    if (para.Infilename.size() == 0){
                                      cerr << "Error: Please input the sequence file by -m or -r ( and -R)" << endl;             
                                      exit(0);
                                      }
    //check search database similarity
    if (para.db_similarity <= 0 || para.db_similarity > 1){
    		cerr <<"Error: Please input right sequence alignment threshold, the value range is from 0 to 1" << endl;
    		exit(0);
	}
    //check pair
    /* //delete by Shi Haobo 
    if (para.Is_paired){
                        if (para.Type != 0){
                                      cerr << "Error: Pair-end sequences only support 16S rRNA" << endl;
                                      exit(0);
                                      } 
						if (Check_Format(para.Infilename.c_str()) == 0||Check_Format(para.Infilename2.c_str())==0){
							cerr << "Error: For pair-end sequences only support fastq format" << endl; 
							exit(0);
						}                    
                        }
    */
    Check_Path(para.Out_path.c_str(), 0); //Check output path 
    Check_Path((para.Out_path + "/tmp").c_str(), 0);    
    
    int max_Core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((para.Core_number <= 0) || (para.Core_number > max_Core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    para.Core_number = max_Core_number;
                    }
                         
    Print_Config(para);
    
    return para.Type;    
    }

#endif
