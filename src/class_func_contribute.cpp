// class_func_contribution
// Updated at Sep 20, 2018
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>

#include "class_func_contribute.h"

using namespace std;

// Parameters Def
_PMDB Database;

string Listfile;
string Listprefix;
string Tablefile;

vector <string> Infilename;
vector <string> Sam_name;

string SelectKOlist;
string Outfilename = "func.KO.contribute";
int Mode = -1; //0: single; 1: list; 2: table

int Coren = 0;
    
int printhelp(){
    
    cout << "Predict-func-contribute version: " << Version << endl;
    cout << "\tPredict the funtional contribution (16S only)" << endl;
    cout << "Usage:" << endl;
    cout << "PM-predict-func-contribute [Option] Value" << endl;
    cout << "Options: " << endl;

    cout << "\t-D (upper) ref database, " << _PMDB::Get_Func_Args() << endl;
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single file" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -L (upper) Selected KO list" << endl;
    cout << "\t  -o Output file name, default is \"func.KO.contribute\"" << endl;

    cout << "\t[Other options]" << endl;
    cout << "\t  -t Number of thread, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
    
    if (argc ==1) 
		printhelp();
    
    int i = 1;
    Mode = -1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D': Database.Set_DB(argv[i+1][0]); break;
                            case 'i':
                                 if (Mode != -1) {
                                             
                                             cerr << "Error: -i conflicts with -l and -T and -B" << endl;
                                             exit(0);
                                             
                                             }

                                 Infilename.push_back(argv[i+1]);
                                 Sam_name.push_back("Sample");
                                 Mode = 0;                                      
                                 break;
                                                      
                            case 'l':                                       
                                      if (Mode != -1){
                                                  
                                                  cerr << "Error: -l conflicts with -i and -T and -B" << endl;
                                                  exit(0);
                                                  
                                                  }
                                      Listfile = argv[i+1];
                                      Mode = 1;
                                      break;
                            case 'T':
                                      if (Mode != -1){
                     
                                                  cerr << "Error: -T conflicts with -i and -l and -B" << endl;
                                                  exit(0);
                     
                                                  }
                                      Tablefile = argv[i+1];
                                      Mode = 2;
                                      break;
                            case 'p': Listprefix = argv[i+1]; break;
                            case 'L': SelectKOlist = argv[i+1]; break;
                            case 'o': Outfilename = argv[i+1]; break;
                           
                            case 't': Coren = atoi(argv[i+1]); break;
                 
                            case 'h': printhelp(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    if (!Database.Get_Is_Func()) {
        cout << "The " << Database.Get_Description() << " domain is not supported yet" << endl;
        exit(0);
        }
    
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > max_core_number)){
        //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
        Coren = max_core_number;
    }
    
    return 0;
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    int sam_num = 0;
    
    _KO_OTU_Table_Contribute KOs;
    
    cout << endl << "Functional Annotation Starts" << endl;
    
    switch (Mode){
        case 0: {
		        sam_num = 1;
                KOs = _KO_OTU_Table_Contribute(Database, sam_num, 0);
                KOs.Load_Sample_Multi(Infilename, Sam_name, Coren);
                cout << "Total Sample Number is " << sam_num << endl;
                if (SelectKOlist.size() > 0)
                   cout << KOs.Load_Selected_KO(SelectKOlist.c_str()) << " KOs are selected" << endl;
                cout << KOs.Output_Contribute(Outfilename.c_str()) << " KOs have been parsed out" << endl; 	
                break;
            	}
        case 1: {
				sam_num = Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);
                KOs = _KO_OTU_Table_Contribute(Database, sam_num, 0);
                KOs.Load_Sample_Multi(Infilename, Sam_name, Coren);
                cout << "Total Sample Number is " << sam_num << endl;
                if (SelectKOlist.size() > 0)
                   cout << KOs.Load_Selected_KO(SelectKOlist.c_str()) << " KOs are selected" << endl;
                cout << KOs.Output_Contribute(Outfilename.c_str()) << " KOs have been parsed out" << endl; 	
                break;
            	}
        case 2: {
                _Table_Format table(Tablefile.c_str());
                sam_num = table.Get_Sample_Size();
                KOs = _KO_OTU_Table_Contribute(Database, sam_num, 0);
                KOs.Load_Sample_By_OTU_Table(&table, Coren);
                cout << "Total Sample Number is " << sam_num << endl;
                if (SelectKOlist.size() > 0)
                   cout << KOs.Load_Selected_KO(SelectKOlist.c_str()) << " KOs are selected" << endl;
                cout << KOs.Output_Contribute(Outfilename.c_str()) << " KOs have been parsed out" << endl; 	      
                break;
        	}
        default: cerr << "Error: Incorrect mode"; exit(0);
    }
    	    
    cout << endl << "Functional Annotation Finished"<< endl;

    return 0;
    }
