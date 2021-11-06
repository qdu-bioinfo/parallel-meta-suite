// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
//  version 3.1 or above with _Table_Format
// _Table_Format input support
// _OTU_Parser 

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "taxa_sel.h"
#include "utility.h"
#include "version.h"
#include "taxonomy.h"
#include "otu_parser.h"

#define TLevN 8

using namespace std;

//_Taxa_Table Taxa_table;
_Table_Format_Seq Taxa_table;

string Listfilename;
string Listprefix;
string Tablefilename;

int Mode = -1; //0: list; 1: table

//string Listfilename;
string Outfilename = "taxaonomy_selection";

int Level = 5; ////(0: Kingdom, 5: Genus, 6: Species, 7: OTU)
float Max_abd = 0.001;
float Min_abd = 0;
float No_zero_rate = 0.1;
float Ave_t = 0.001;

int Min_seq = 2;

bool Is_cp_correct = true;
bool Is_print = false;

_PMDB Database;
_OTU_Parser otu_parser;

hash_map <string, float, std_string_hash> Cp_number;
hash_map <string, string, std_string_hash> OTU_taxa_hash;

string Taxa_level[TLevN] = {"kingdom", "phylum", "class", "order", "family", "genus", "species", "OTU"};

void Print_Help(){
    
    cout << "Select-taxa version:" << Version << endl;
    cout << "\tSelect the taxonomy features" << endl;
    cout << "Usage: " << endl;
    cout << "PM-select-taxa [Option] Value" << endl;
    cout << "Option: " << endl;
    
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    
    cout << "\t[Input options, requried]" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix for '-l' [Optional for -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file name, default is \"taxaonomy_selection\"" << endl;
    cout << "\t  -L (upper) Taxonomical level (1-6: Phylum - Species, 7: OTU). default is 5" << endl;
    cout << "\t  -P (upper) Print distribution barchart, T(rue) or F(alse), default is F" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -r rRNA copy number correction, T(rue) or F(alse), default is T" << endl;
    cout << "\t  -q Minimum sequence count threshold, default is 2" << endl;
    cout << "\t  -m Maximum abundance threshold, default is 0.001 (0.1%)" << endl;
    cout << "\t  -n Minimum abundance threshold, default is 0.0 (0%)" << endl;
    cout << "\t  -z Minimum No-Zero abundance threshold, default is 0.1 (10%)" << endl;
    cout << "\t  -v Minimum average abundance threshold, default is 0.001 (0.1%)" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    }

int Get_TLevel(string levels){
    switch(levels[0]){
                              case '1':
                              case 'P': 
                              case 'p': return 1; break;
                              case '2':
                              case 'C': 
                              case 'c': return 2; break;
                              case '3':
                              case 'O': 
                              case 'o': return 3; break;
                              case '4':
                              case 'F':
                              case 'f': return 4; break;
                              case '5':
                              case 'G':
                              case 'g': return 5; break;
                              case '6':
                              case 'S':
                              case 's': return 6; break;
                              case '7':
                              case 'U':
                              case 'u': return 7; break;
                              
                              default:
                                      cerr << "Warning: Unrec taxa level: " << levels[0] << endl; return 0; break; 
                              }
            
    return 0;
    }

void Parse_Para(int argc, char * argv[]){
      
      int i = 1;
      
      if (argc ==1)
		Print_Help();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D': Database.Set_DB(argv[i+1][0]); break;
                            case 'l': Listfilename = argv[i+1];
                                      if (Mode >=0){
                                          cerr << "Error: -l conflists with -T and -B" << endl;
                                          exit(0);
                                          }
                                      else Mode = 0; break;
                 
                            case 'p': Listprefix = argv[i+1]; break;
                 
                            case 'T': Tablefilename = argv[i+1];
                                      if (Mode >=0){
                                          cerr << "Error: -T conflists with -l and -B" << endl;
                                          exit(0);
                                          }
                                      else Mode = 1; break;
                                      
                            case 'o': Outfilename = argv[i+1]; break;
                            case 'L': Level = Get_TLevel(argv[i+1]); break;
                            case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;   
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_print = true; break;   
                            case 'q': Min_seq = atoi(argv[i+1]); break;   
                            case 'm': Max_abd = atof(argv[i+1]); break;
                            case 'n': Min_abd = atof(argv[i+1]); break;
                            case 'z': No_zero_rate = atof(argv[i+1]); break;
                            case 'v': Ave_t = atof(argv[i+1]); break;
                            case 'h': Print_Help(); break;
                            default : printf("Error: Unrec argument %s\n", argv[i]); Print_Help(); break; 
                            }
         i+=2;
         }
         
     //check
     if ((Level <= 0) || (Level > 7)){
               cerr << "Warning: Taxonomical level must be 1-7, change to default (5)" << endl;
               Level = 5;
               }
     
     if (Max_abd < Min_abd){
           cerr << "Warning: m must be larger than n, change to default (0, 0.001)" << endl;
           Max_abd = 0.001;
           Min_abd = 0.0;
        }
    
     Level ++; //1: Kingdom
     
     otu_parser = _OTU_Parser(Database);
     }

/*
string Get_Taxa(string otu, int level){
    
    string otu_taxa = otu_parser.Get_taxa_by_OTU(otu);
        
    vector <string> taxa_buffer;
    stringstream strin(otu_taxa);
    string taxa;
    string temp;
    
    if (level == 8){ //otu
        taxa = "OTU_";
        taxa += otu;
    }
    
    else for (int i = 0; i < level; i++){
        strin >> taxa;
        
        if (taxa.find("Unclassified") != string::npos){
            return "Unclassified";
        }
        if (taxa.find("otu_") != string::npos){
            return "Unclassified";
        }
        
        if(taxa[taxa.size()-1] != ';'){
            strin >> temp;
            taxa += "_";
            taxa += temp;
        }
        taxa_buffer.push_back(taxa);
        
        //add genus info for species
        if (i == 6){
            taxa = taxa_buffer[5];
            taxa = taxa.substr(0, taxa.size()-1);
            taxa += "_";
            taxa += taxa_buffer[6];
        }
        
        //add prefix for genus
        if (i == 5){
            string taxa_prefix = taxa_buffer[i-1];
            taxa_prefix = taxa_prefix.substr(0, taxa_prefix.size() - 1);
            if (taxa_prefix.size() > 3) taxa_prefix = taxa_prefix.substr(0, 3);
            taxa = taxa_prefix + "_" + taxa;
        }
    }
    
    return taxa;
}
*/
int Load_Taxa(string listfilename, string listprefix, int level){
    
    vector <string> infilename;
    vector <string> sam_name;
    
    //load
    Load_List(listfilename.c_str(), infilename, sam_name, listprefix);
    
    cout << "Total Sample Number is " << sam_name.size() << endl;
    
    Taxa_table = _Table_Format_Seq(sam_name);
        
    for (int i = 0; i < infilename.size(); i ++){
        
        hash_map <string, int, std_string_hash> otu_count;
        int seq_count = otu_parser.Load_file_to_hash(infilename[i].c_str(), otu_count);
        
        for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                              
                              //string taxa = Get_Taxa(miter->first, level);
                              string taxa = otu_parser.Get_taxa_by_OTU(miter->first, level);
                                                      
                              float cp_num = 1.0;
                              if (Is_cp_correct)
                                 cp_num = otu_parser.Get_cp_by_OTU(miter->first);
                                                                                                                                 
                              Taxa_table.Add_Feature(taxa, miter->second, i, cp_num);                                                                                       
                              }
        
        cout << seq_count << " sequences are loaded in file " << i + 1 << endl; 
        }
                 
    return 0;
    }

int Load_Taxa_Table(_Table_Format input_table, int level){
    
    cout << "Total Sample Number is " << input_table.Get_Sample_Size() << endl;
    
    Taxa_table = _Table_Format_Seq(input_table.Get_Sample_Names());
    
    
    vector <string> otus = input_table.Get_Feature_Names();
    
    vector <string> taxa;
    vector <float> cp_num;
    
    for (int i = 0; i < otus.size(); i ++){
        
        string a_otu = Check_OTU(otus[i]);
        //taxa.push_back(Get_Taxa(a_otu, level));
        taxa.push_back(otu_parser.Get_taxa_by_OTU(a_otu, level));
        
        float cp = 1.0;
        if (Is_cp_correct)
                cp = otu_parser.Get_cp_by_OTU(a_otu);
        cp_num.push_back(cp);
        }
    
    for (int i = 0; i < input_table.Get_Sample_Size(); i ++){
        
        unsigned int seq_sum = 0;
        
        //vector <float> seq_counts = input_table.Get_Abd(i);
        
        for (int j = 0; j < input_table.Get_Feature_Size(); j ++){
                    Taxa_table.Add_Feature(taxa[j], (int) input_table.Get_Abd_By_Order(i, j), i, cp_num[j]);
                    seq_sum += (int) input_table.Get_Abd_By_Order(i, j);
            }
        cout << seq_sum << " sequences are loaded in file " << i + 1 << endl;
        }
    return 0;
    }
    
int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);     
    
    cout << "Taxonomical Selection Starts" << endl << endl;
    cout << "Level is " << Taxa_level[Level-1] << endl;
    
        
    switch (Mode){
        case 0: Load_Taxa(Listfilename, Listprefix ,Level); break;
        case 1:{
                _Table_Format table(Tablefilename.c_str());
                Load_Taxa_Table(table, Level);
                break;
                }
        default: cerr << "Error: Incorrect mode"; exit(0);
        }
    
    Taxa_table.Filter_Seq_Count(Min_seq);
    
    Taxa_table.Normalization();
    
    Taxa_table.Filter_Abd(Max_abd, Min_abd, No_zero_rate, Ave_t);
    
    cout << "Total Taxa Number is " << Taxa_table.Get_Taxa_Size() << endl;
    //output abd
    string abd_file = Outfilename + "." + Taxa_level[Level-1] + ".Abd";
    string count_file = Outfilename + "." + Taxa_level[Level-1] + ".Count";
    cout << "Total Output Taxa Numer is " << Taxa_table.Output_Abd(abd_file.c_str()) << endl;
    Taxa_table.Output_Count(count_file.c_str());
    
    if (Is_print){
                  char command[BUFFER_SIZE];
                  sprintf(command, "Rscript %s/Rscript/PM_Distribution.R -i %s -o %s -p %s", Check_Env().c_str(), abd_file.c_str(), (abd_file + ".distribution").c_str(), (Taxa_level[Level-1] + ".Abd").c_str());
                  system(command);
                  }
    
    return 0;
    }
