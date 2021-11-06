// Updated at July 31, 2021
// Updated by Yuzhu Chen
// Updated at Dec 23, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>

#include "class_func.h"

using namespace std;

// Parameters Def
_PMDB Database;

string Listfile;
string Listprefix;
string Tablefile;

vector <string> Infilename;
vector <string> Sam_name;
vector <string> Outfilename;

string Out_path = "Results_Func/";
string Outlistfile;

int Mode = -1; //0: single; 1: list; 2: table
bool Is_out_list = true;

int Coren = 0;
    
int printhelp(){
    
    cout << "Predict-func version: " << Version << endl;
    cout << "\tFunctional prediction" << endl;
    cout << "Usage:" << endl;
    cout << "PM-predict-func [Option] Value" << endl;
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
    cout << "\t  -o Output path, default is \"Results_Func\"" << endl;
    cout << "\t  -L (upper) If output list (at \"Output_path/func.list\"), T(rue) or F(alse), default is T [optional for -l]" << endl; 
    
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
                            case 'o': Out_path = argv[i+1]; break;
                            case 'L': if ((argv[i+1][0] == 'F') || (argv[i+1][0] == 'f')) Is_out_list = false; break;
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
        
    if (Is_out_list) {
	if(Out_path[Out_path.size()-1] != '/')
	        Outlistfile = Out_path + ".list";
 	else
        	Outlistfile = Out_path.substr(0, Out_path.size()-1) + ".list";
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
    
    _KO_OTU_Table_All KOs;
    
    cout << endl << "Functional Annotation Starts" << endl;
    
    switch (Mode){
        case 0: {
                Check_Path(Out_path.c_str(), 1);
		        sam_num = 1;
                Outfilename.push_back(Out_path + "/functions.txt");
                KOs = _KO_OTU_Table_All(Database, sam_num, 0);
                KOs.Load_Sample_Multi(Infilename, Sam_name, Coren);
                cout << "Total Sample Number is " << sam_num << endl;
                cout << KOs.Output_Multi(Outfilename) << " KOs have been parsed out" << endl; 	
                break;
                break;
            	}
        case 1: {
                Check_Path(Out_path.c_str(), 1);
				sam_num = Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);
                for (int i = 0; i < sam_num; i ++){
                    Check_Path((Out_path + "/" + Sam_name[i]).c_str(), 1);
                    Outfilename.push_back(Out_path + "/" + Sam_name[i] + "/functions.txt");
                    }
                KOs = _KO_OTU_Table_All(Database, sam_num, 0);
                KOs.Load_Sample_Multi(Infilename, Sam_name, Coren);
                if (Is_out_list){
                	//vector <string> sam_name_out;
    				//sam_name_out.push_back("samples"); // with id
                	//sam_name_out.insert(sam_name_out.end(), Sam_name.begin(), Sam_name.end());
                	//Make_list(Outlistfile.c_str(), Out_path.c_str(), sam_name_out, 1);
                	Make_list(Outlistfile.c_str(), Out_path.c_str(), Sam_name, 1);
					}	
                cout << "Total Sample Number is " << sam_num << endl;
                cout << KOs.Output_Multi(Outfilename) << " KOs have been parsed out" << endl; 	
                break;
            	}
        case 2: {
				_Table_Format table(Tablefile.c_str());
                sam_num = table.Get_Sample_Size();
                KOs = _KO_OTU_Table_All(Database, sam_num, 0);
                KOs.Load_Sample_By_OTU_Table(&table, Coren);
                cout << "Total Sample Number is " << sam_num << endl;
                cout << KOs.Output_By_Category((Out_path + ".KO").c_str(), 3, 0, 0) << " KOs have been parsed out" << endl;
                break;
        	}
        default: cerr << "Error: Incorrect mode"; exit(0);
    }
    	    
    cout << endl << "Functional Annotation Finished"<< endl;

    return 0;
    }
