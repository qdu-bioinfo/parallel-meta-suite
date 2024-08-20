// Updated at Aug 16, 2024
// Updated by Haobo Shi,Xiaoquan Su
// version 3.7 - 3.7.2
// Update profiler,remove pair-end mode

// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>

#include "init.h"
#include "fastq.h"
#include "ExtractRNA.h"
#include "multialign.h"
#include "taxonomy.h"

using namespace std;

void Single_Run(_Para para){
	
	string command;
	
	int seq_count = 0;
	int rna_count = 0;
	int asv_count = 0;
	int match_rna_count = 0;
	int drop_rna_count = 0;

	int a_diver [LEVEL] = {0, 0, 0, 0, 0, 0, 0};
	
	string mergefile;
	
	string handlefile;
	
	//check format, infilename
	para.Format = Check_Format(para.Infilename.c_str());
	
	if(para.Format < 0) return;//Format error
	
	if(para.Is_format_check){//Check format
        cout << "Format Check Starts" << endl;
        command = "PM-format-seq -i " + para.Infilename;
        system(command.c_str());
        cout << endl;
    }
	
	//check format for pair ends, infilename2
    /* //delete
    if (para.Is_paired){
        para.Format = Check_Format(para.Infilename2.c_str()); 
        if (para.Format < 0) return;//Format error
                        
        if (para.Is_format_check){//Check format
            cout << "Format Check 2 Starts" << endl;
            command = "PM-format-seq -i " + para.Infilename2;
            system(command.c_str());
            cout << endl;
        }
    }  
	*/
	//Type, 0:16S  1:shotgun
	//Format, 0:fasta  1:fastq
	
	if(para.Type == 1){//if meta
		if (para.Format == 1){//If fastq
        	string tempfilename = para.Out_path + "/meta.fasta";
        	cout << "Pre-computation for Fastq Starts" << endl;
        	cout << endl << Fastq_2_Fasta(para.Infilename.c_str(), tempfilename.c_str()) << " sequences have been pre-computed" << endl << endl;
        	para.Infilename = tempfilename;
        }
		//Extract 16S r RNA
        seq_count = ExtractRNA(para.Database.Get_Domain(), para.Infilename, para.Out_path, para.Length_filter, para.This_path, para.Core_number);
        handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Out_path + "/meta.rna", para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, para.Format, para.Core_number);
		//search database
		Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt" , para.Database.Get_Path()+ "/database.fa" , para.db_similarity,'F', para.Core_number, para.profiler);
	}
    /* //delete by Shi Haobo 
    else if(para.Is_paired){//if paired
		if(para.Format == 0){//fasta cant merge in vsearch
			cerr << "Error: For pair ends you need to input fastq format file" << endl; 
    	 	return ;
		}else{//fastq
			//merge
			mergefile = Merge_Pairend(para.Align_exe_name.c_str(), para.Infilename, para.Infilename2, para.Out_path + "/tmp", para.Core_number);
			//dereplication, denoise, nonchimeras
			handlefile = Handle_seq(para.Align_exe_name.c_str(), mergefile , para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count,0, para.Core_number);
			//search database
			Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt" , para.Database.Get_Path()+ "/database.fa" , para.db_similarity,'F', para.Core_number, para.profiler);
			if (rna_count < 0) {                       
                       cerr << "Error: 2 ends contain different number of sequences" << endl;
                       return;
            }
		}	
		
	}
    */
    else{//single
		//dereplication, denoise, nonchimeras
		handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Infilename, para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, para.Format, para.Core_number);	
		
		//search database
		Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt" , para.Database.Get_Path()+ "/database.fa" , para.db_similarity,'F', para.Core_number, para.profiler);
	}
    if(para.profiler == 1)
        Otutab_Count(para.Out_path + "/tmp/PM.hwl.txt",para.Out_path + "/tmp/map_output.txt");
	//parse_taxonomy 
    Out_Taxonomy((para.Out_path + "/tmp/map_output.txt").c_str(), para.Out_path, para.Database, 0, a_diver, para.Is_paired, match_rna_count, drop_rna_count);
          
    //make_plot
    command = "PM-plot-taxa -D ";
    command += para.Database.Get_Id();
    command += " -i ";
    command += para.Out_path + "/classification.txt -o " + para.Out_path;
    system(command.c_str());
    
    //print report
    Print_Report(para, seq_count, rna_count, match_rna_count, a_diver,asv_count);
    
    //func_anno
    if (para.Is_func){
       //func
        command = "PM-predict-func -D ";
        command += para.Database.Get_Id();
        command += " -i ";
        command += para.Out_path + "/classification.txt -o " + para.Out_path;
        system(command.c_str());
        
        //nsti
        command = "PM-predict-func-nsti -D ";
        command += para.Database.Get_Id();
        command += " -i ";
        command += para.Out_path + "/classification.txt >> " + para.Out_path + "/Analysis_Report.txt";
        system(command.c_str());
       }
    //Remove the tmp 
    command = "rm -rf " + para.Out_path + "/tmp";
    system(command.c_str());	
}	
	

int main(int argc, char * argv[]){
    
    cout << endl << "Welcome to Parallel-Meta Suite version " << Version << endl << endl; 
    
    _Para para;         
    
    Parse_Para(argc, argv, para);
    
    Single_Run(para);
        
    cout << endl << "Parallel-Meta Suite Finished" << endl;
    cout << "Please check the analysis results and report at " << para.Out_path <<endl;
     
    return 0;
}
