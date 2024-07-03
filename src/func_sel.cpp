// Updated at July 2, 2024
// Updated by Xiaoquan Su
// Bioinformatics Group, Qingdao University

#include <iostream>
#include <fstream>
#include <sstream>

#include "class_func.h"

using namespace std;

// Parameters Def
#define FLevN 4

string Out_file = "functions_category";

_PMDB Database;

string Listfilename;
string Listprefix;
string Tablefilename;

vector <string> Infilename;
vector <string> Sam_name;

int Level = 2; //category level, default is 2

float Max_abd = 0;
float Min_abd = 0;

bool Is_print = false;

int Mode = -1; //0: list; 1: table

string Func_level[FLevN] = {"l1", "l2", "l3", "KO"};
string Func_level_display[FLevN] = {"Pathway level 1", "Pathway level 2", "Pathway level 3", "KO"};


int printhelp(){
    
    cout << "Select-func version " << Version << endl;
    cout << "\tSelect the functional features" << endl;
    cout << "Usage:" << endl;
    cout << "PM-select-func [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -l Input files list [Conflicts with -T and -B]" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input KO Absolute Count table (*.KO.Count) [Conflicts with -l and -B]" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file, default is \"functions_category\"" << endl;
    cout << "\t  -L (upper) KEGG Pathway level, Level 1, 2, 3 or 4 (KO number), default is 2" << endl;
    cout << "\t  -P (upper) Print distribution barchart, T(rue) or F(alse), default is F" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    return 0;
    }

void Parse_Para(int argc, char * argv[]){
    
    if (argc ==1) 
		printhelp();
    
    int i = 1;

    Level = 2;
    
    Mode = -1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){                            
                            case 'D': Database.Set_DB(argv[i+1][0]); break;
                 
                            case 'l': if (Mode != -1){
                                
                                                cerr << "Error: -l conflicts with -T and -B" << endl;
                                                exit(0);
                                
                                        }
                                        Listfilename = argv[i+1];
                                        Mode = 0;
                                        break;
                 
                            case 'p': Listprefix = argv[i+1]; break;
                 
                            case 'T': if (Mode != -1){
                 
                                                cerr << "Error: -T conflicts with -l and -B" << endl;
                                                exit(0);
                 
                                        }
                                        Tablefilename = argv[i+1];
                                        Mode = 1;
                                        break;
                 
                            case 'o': Out_file = argv[i+1]; break;
                            case 'L': Level = atoi(argv[i+1]); break; //by cate
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_print = true; break;                        
                            case 'h': printhelp(); break;

                            default : printf("Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    if (!Database.Get_Is_Func()){
        cout << "The " << Database.Get_Description() << " domain is not supported yet" << endl;
        exit(0);
        }
    
    //check
    if ((Level > 4) || (Level < 1)){               
               cerr << "Warning: Level must between 1 and 4, change to default (2)" << endl;
               Level = 2;
               }
    Out_file = Out_file + "." + Func_level[Level - 1];
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    _KO_OTU_Table_All KOs;
    _Table_Format table;
    int sam_num = 0;
    
    cout << endl << "Functional Selection Starts" << endl;
    cout << endl << "Level is " << Func_level_display[Level-1] << endl;
    
    switch (Mode){
        case 0: {
				sam_num = Load_List(Listfilename.c_str(), Infilename, Sam_name, Listprefix);
                KOs = _KO_OTU_Table_All(Database, sam_num, 1); // sim mode
                //load
                for (int i = 0; i< sam_num; i++){
                    cout << KOs.Load_Sample_By_Single_KO_Table(Infilename[i].c_str(), Sam_name[i], i) << " genes are loaded" << endl;
                    }
            
                break;
			}
        case 1: {
				table = _Table_Format(Tablefilename.c_str());
                sam_num = table.Get_Sample_Size();
                KOs = _KO_OTU_Table_All(Database, sam_num, 1); // sim mode
                //load
                KOs.Load_Sample_By_KO_Table(&table);
                break;
			}
        default: cerr << "Error: Incorrect mode"; exit(0);
        }
    
    cout << "Total Sample Number is " << sam_num << endl;

            
    if (Level > 0)
       cout << endl << KOs.Output_By_Category(Out_file.c_str(), Level - 1, Max_abd, Min_abd) << " Category have been parsed out" << endl;;
        
    if (Is_print){
                  char command[BUFFER_SIZE];
                  
                  sprintf(command, "Rscript %s/Rscript/PM_Distribution.R -i %s -o %s -p %s", Check_Env().c_str(), (Out_file + ".Abd").c_str(), (Out_file + ".Abd.distribution").c_str(), (Func_level[Level-1] + ".Abd").c_str());
                  system(command);
                  }
                  
    cout << endl << "Functional Selection Finished"<< endl;

    return 0;
    }
