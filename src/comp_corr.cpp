// Updated at July 2, 2024
// Updated by Xiaoquan Su
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#include "version.h"
#include "utility.h"
#include "table_format_corr.h"

using namespace std;
string Infilename;
string Infilenamecorr;
string Outfilenameself="self_corr_matrix.out";
string Outfilenamecorr="corr_matrix.out";
int Coeff = 0; //0: S 1: P
string Selotu;
int Flag=0;
int Tn=0;


bool Is_network = false;
float Network_t = 0.5;

void Print_Help(){

    cout << "Comp-corr version: " << Version << endl;
    cout << "\tCompute the correlation and network by feature table" << endl;
    cout << "Usage:" << endl;
    cout << "PM-comp-corr [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
	cout << "\t  -i or -T (upper) Input feature table file (*.Abd)" << endl;
	cout << "\t  -m Meta data name [Optional]" << endl;
    cout << "\t  -c Selected feature, separated by \",\" [Optional for -m]"<<endl;
    
    cout << "\t[Output options]" << endl;
	cout << "\t  -o Output prefix, default is \"corr_matrix\"" << endl;
	
    cout << "\t[Other options]" << endl;
	cout<<  "\t  -f 0:(Spearman) or 1:(Pearson) metrics,default is 0"<<endl;
	cout << "\t  -N (upper) Network based co-occurrence analysis, T(rue) or F(alse), default is F" << endl;
	cout << "\t  -G (upper) Netowrk analysis threshold, default is 0.5"<<endl;
	cout << "\t  -t Number of thread, default is auto"<<endl;
	cout << "\t  -h Help" << endl;

	exit(0);
}

void Para(int argc,char *argv[]){
	int i=1;
	if (argc==1)
		Print_Help();
	while (i<argc)
	{
		if (argv[i][0]!='-')
		{
			cout<<"Argument #"<<i<<" Error : Arguments must start with -\n"<<endl;
			exit(0);
		}
		switch (argv[i][1])
		{
        case 'i': 
        case 'T': Infilename=argv[i+1];break;
        case 'm': Infilenamecorr=argv[i+1];break;
		case 'o': {
            Outfilenameself=argv[i+1];
            Outfilenameself+=".self_matrix.out";
            Outfilenamecorr=argv[i+1];
            Outfilenamecorr+=".corr_matrix.out";
			break;
				  }
        case 'f': Coeff= atoi(argv[i+1]); break;
        case 'c': {Selotu=argv[i+1];Flag=1;break;}
        case 'N': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_network = true; break;
        case 'G': {Network_t = atof(argv[i+1]);break;}
        case 't': {Tn=atoi(argv[i+1]);break;
				  }
		default: cout<<"Error: Unrec arguments"<<argv[i]<<endl;Print_Help();break;
		}
		i+=2;
	}

     
int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Tn  <= 0) || (Tn  > Max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Tn  = Max_core_number;
                    }  
}



             
int main(int argc,char *argv[])
{

	Para(argc,argv);
       
	_Table_Format taxa(Infilename.c_str(),0);
	taxa.Calc_Corr_Matrix(Outfilenameself.c_str(),Coeff,Tn);
	
	
    if(Infilenamecorr.size()!=0)
    {
		_Table_Format_Meta Corr_With_Meta;
		if(Flag==0)
   	 Corr_With_Meta.Calc_Corr_Meta_Matrix(Infilename.c_str(),Infilenamecorr.c_str(),Outfilenamecorr.c_str(),Coeff,Tn);
		else
			 Corr_With_Meta.Calc_Corr_Meta_Matrix_S(Infilename.c_str(),Infilenamecorr.c_str(),Outfilenamecorr.c_str(),Selotu,Coeff,Tn);
    }
    
    if (Is_network){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Network.R -i %s -o %s -t %f", Check_Env().c_str(), Outfilenameself.c_str(), (Outfilenameself + ".network.pdf").c_str(), Network_t);
                    system(command);
                    }

	return 0;
}
