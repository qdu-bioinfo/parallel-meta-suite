// Updated at Sep 4, 2018
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// _Table_Format supported

#include <iostream>
#include <stdlib.h>

#include <sys/dir.h>
#include <sys/stat.h>

#include <dirent.h>
#include <unistd.h>

#include "class_tax.h"

using namespace std;

// Parameters Def
string Out_path = "Result_Plot";
string Html_file = "taxonomy.html";
//string Txt_file = "taxonomy.txt";
//string Svg_file = "taxonomy.svg";

int Depth = 7;

string Listfile;
string Listprefix;
string Tablefile;

vector <string> Infilename;
vector <string> Sam_name;

int Mode = -1; //0: single; 1: list; 2: table

bool Is_cp_correct = true;

_PMDB Database;

int printhelp(){
    
    cout << "Plot-taxa version: " << Version << endl;
    cout << "\tTaxonomy profile visualization by Krona" << endl;
    cout << "Usage:" << endl;
    cout << "PM-plot-taxa [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single file" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output Path, default is \"Result_Plot\"" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -r rRNA copy number correction, T(rue) or F(alse), default is T" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){

    if (argc ==1) 
		printhelp();
    
    Mode = -1;
    Listprefix = "";
    
    Is_cp_correct = true;
    
    int i = 1;    
    int sam_num = 0;
    
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
                            
                            case 'p': Listprefix = argv[i+1]; break;
                 
                            case 'T':
                                    if (Mode != -1){
                     
                                            cerr << "Error: -T conflicts with -i and -l and -B" << endl;
                                            exit(0);
                     
                                    }
                 
                                    Tablefile = argv[i+1];
                                    Mode = 2;
                                    break;
                                    
                            
                            case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;
                            case 'o': Out_path = argv[i+1]; break;
                            case 'd': Depth = atoi(argv[i+1]); break;
                            case 'h': printhelp(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    Html_file = Out_path + "/" + Html_file;
    //Txt_file = Out_path + "/" + Txt_file;
    //Svg_file = Out_path + "/" + Svg_file;
    
    if (Depth <=0) Depth = 7;
    
    Check_Path(Out_path.c_str(), 1);

    return 0;
    }

int Copy_JS(const char * path){
    
    string js_copy;
    
    //copy js 
    js_copy ="cp "; 
    js_copy += Check_Env();
    js_copy += "/html/*js ";
    js_copy += path;
    
    system(js_copy.c_str());
    
    //copy png
    js_copy ="cp "; 
    js_copy += Check_Env();
    js_copy += "/html/*png ";
    js_copy += path;
    
    system(js_copy.c_str());
    
    return 0;
    
    }


int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    _Table_Format table;
    //load
    int sam_num = 0;
    
    cout << endl << "Taxonomic Classification Starts" << endl;
        
    switch (Mode){
            
        case 0: {
			sam_num = 1; 
			break;
			}
        case 1: {
			sam_num = Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix); 
			break;
      		}
	case 2: {
			table = _Table_Format(Tablefile.c_str()); 
			sam_num = table.Get_Sample_Size(); 
			break;
			}
        default: cerr << "Error: Incorrect mode" << endl; break;
        }
    
    TNode::Init(Database, sam_num);
    TNode root;
    
    switch (Mode){
            
        case 0:
        case 1:
            for (int i = 0; i < sam_num; i++)
                cout << root.Read_file(Infilename[i].c_str(), Sam_name[i], i, Is_cp_correct) << " sequences are loaded" << endl;
            break;
            
        case 2: // root.Read_table(&table, Is_cp_correct); break;
        case 3: root.Read_table(&table, Is_cp_correct); break;
        default: cerr << "Error: Incorrect mode" << endl; break;
        }
    
    root.Out_Tree_Html(Html_file.c_str(), Depth);
    Copy_JS(Out_path.c_str());
    
    cout << endl << "Taxonomic Classification Finished"<< endl;
    cout << "Plot Finished"<< endl;
    
    return 0;
    
    }

