// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Code by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// version 3.1 - 3.5.4 with Bowtie2


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>

#include <sys/types.h>
#include <sys/dir.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include "utility.h"

#ifndef MULTIALIGN_H
#define MULTIALIGN_H

#define Seq_Size 500

using namespace std;

string Merge_Pairend(string programname, string infilename_1, string infilename_2, string outpath, int coren){
	cout << "Merging starts" << endl;
        
    int seq_n_1 = Get_Count_fastq(infilename_1.c_str());  
    int seq_n_2 = Get_Count_fastq(infilename_2.c_str());
	  
    if (seq_n_1 != seq_n_2) return "";
    
    cout << "There are " << seq_n_1 << " paired sequences in total" << endl << endl;
    
    //mkdir((outpath + "/maptemp").c_str(), 0755);
	
	//command
    
    char command[BUFFER_SIZE];
    
	//sprintf(command,"%s --fastq_mergepairs %s --reverse %s --fastaout %s/merged",programname.c_str(),infilename_1.c_str(), infilename_2.c_str(),outpath.c_str()); 
	
	sprintf(command,"%s --fastq_mergepairs %s --reverse %s --fastq_minovlen 5 --fastaout %s/merged --threads %d",programname.c_str(),infilename_1.c_str(), infilename_2.c_str(),outpath.c_str(),coren); 	
	system(command);
	
	cout << "Merging Finished" << endl << endl;   
	string mergefile=outpath+"/merged" ;
    return  mergefile;
}

string Handle_seq(string programname, string infilename, string outpath, char Is_denoised, char Is_nonchimeras, int & rna_count, int & asv_count, int Format, int coren){
    if(Format == 0){
    	rna_count = Get_Count(infilename.c_str());
	}else{
		rna_count = Get_Count_fastq(infilename.c_str());
	}
    cout << "There are " << rna_count << " sequences in total" << endl << endl;
    
	cout << "Profiling starts" << endl << endl;
    //command  
    char command[BUFFER_SIZE];
    
    //vsearch 1.dereplication; 2.denoise; 3.nonchimeras; 4.db
    
	string fir_name="dereplication";
	string fir_path=outpath+'/'+fir_name;
	//v1:dereplication
		sprintf(command,"%s --derep_fulllength %s --sizeout --output %s --minuniquesize 1 -relabel sequence. --threads %d",programname.c_str(),infilename.c_str(),fir_path.c_str(),coren);
    	system(command); 	
    	//cout<< command << endl;
    	
    string sec_name;
    string sec_path;
    if(Is_denoised == 'T'){//=true,denoise 
    	sec_name="denoised";
    	sec_path=outpath+'/'+sec_name;
    	//v2:denoise
    	sprintf(command,"%s --cluster_unoise %s --sizein --sizeout --centroids %s --minsize 1 --threads %d",programname.c_str(),fir_path.c_str(),sec_path.c_str(),coren);
    	system(command);
    	//cout<< command << endl;
    	asv_count=Get_Count(sec_path.c_str());
	}else{
		sec_name=fir_name;
		sec_path=fir_path;
	}
	
	string thi_name;
	string thi_path;
	if(Is_nonchimeras == 'T'){//=true, remove chimeras
		thi_name="nonchimeras";
		thi_path=outpath+'/'+thi_name;
		//v3:nonchimeras
    	sprintf(command,"%s --uchime3_denovo %s --sizein --sizeout --nonchimeras %s --threads %d",programname.c_str(),sec_path.c_str(),thi_path.c_str(),coren);
    	system(command);
    	//cout<< command << endl;
	}else{
		thi_name=sec_name;
		thi_path=sec_path;
	}
	
	/*	
	string fou_name;
	string fou_path;
	fou_name="rereplicate";
	fou_path=outpath+'/'+fou_name;
	
	sprintf(command,"%s -rereplicate %s -output %s",programname.c_str(),thi_path.c_str(),fou_path.c_str());

	string handlefile = fou_path;
	*/
	string handlefile = thi_path;
	return handlefile; 
}
int Search_db(string programname, string infilename, string outpath, string database_path, double db_similarity, char Shotgun, int coren){
        
    int seq_n = Get_Count(infilename.c_str()); 
	
	//command  
    char command[BUFFER_SIZE];
    
    string fir_name="dereplication";
	string fir_path=outpath+'/'+fir_name;
	
    if(Shotgun == 'T'){
    	cout << "There are " << seq_n << " sequences in total" << endl << endl;
    	cout << "Profiling starts" << endl << endl;
    	
    	sprintf(command,"%s --derep_fulllength %s --sizeout --output %s --minuniquesize 1 -relabel sequence.",programname.c_str(),infilename.c_str(),fir_path.c_str());
    	system(command);
    	
    	infilename=fir_path;
	}

    //global search db, output format is sam 
    sprintf(command,"%s --id %.2f --db %s --usearch_global %s --otutabout %s/map_output.txt --threads %d",programname.c_str(), db_similarity, database_path.c_str(),infilename.c_str(),outpath.c_str(),coren); 
    system(command);
    //cout<< command << endl;
    
    cout << "Profiling finished" << endl ;          
    return seq_n;
    }
#endif
