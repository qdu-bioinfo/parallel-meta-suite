// Updated at Dec 23, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// _OTU_Parser
// _Table_Format

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <unistd.h>

#include "utility.h"
#include "version.h"
#include "db.h"
#include "otu_parser.h"
#include "table_format.h"

#define MAX_BOOT 1000
#define DEF_BOOT 200

using namespace std;

string Listfile;
string Listprefix;

vector <string> Infilename;
vector <string> Outfilename;
vector <string> Sam_name;

string Tablefilename;

_PMDB Database;

string Outpath = "Rare_Out";
string Outlistfile;

int Seq_depth = 0;
int Bootstrap = DEF_BOOT;

int Mode = 0; // 0: single; 1: list; 2: table
bool Is_out_list = true;

hash_map <string, string, std_string_hash> OTU_Taxa;

_OTU_Parser OTU_parser;

int printhelp(){
    
    cout << "Random-rare version: " << Version << endl;
    cout << "\tRarefy samples by taxonomy profiles" << endl;
    cout << "Usage: " << endl;
    cout << "PM-rand-rare [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single file name" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix [Optional or -l]" << endl;
    cout << "\tor" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output path (for -i and -l) or output table name (for -T), default is \"Rare_Out\"" << endl;
    cout << "\t  -L (upper) If output list, T(rue) or F(alse), default is T [optional for -l]" << endl; 
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -s Rarefaction depth [Required]" << endl;
    cout << "\t  -b Bootstrap for sequence number normalization, default is " << DEF_BOOT << ", maximum is " << MAX_BOOT << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
  
    Mode = 0;
    
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
             case 'D': Database.Set_DB(argv[i+1][0]); break;
             case 'i':
                 Infilename.push_back(argv[i+1]);
                 Sam_name.push_back("Sample");
                 Mode = 0;
                 break;
                 
             case 'l':
                 Listfile = argv[i+1];
                 Mode = 1;
                 break;
                                              
             case 'p': Listprefix = argv[i+1]; break;    
             case 'T': Tablefilename = argv[i+1]; Mode = 2; break; 
             case 'o': Outpath = argv[i+1]; break;
             case 'L': if ((argv[i+1][0] == 'F') || (argv[i+1][0] == 'f')) Is_out_list = false; break;
             case 's': Seq_depth = atoi(argv[i+1]); break;
             case 'b': Bootstrap = atoi(argv[i+1]); break;
                 
             case 'h': printhelp(); break;
             default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
             }
          
         i+=2;
         }
    
    if (Seq_depth <=0 ){
          cerr << "Please assign the rarefaction sequence depth by -s" << endl;
          exit(0);
          }
    if (Bootstrap <= 0){
          cerr << "Warning: The minimum bootstrap is 1, change to default " << DEF_BOOT << endl;
          Bootstrap = DEF_BOOT;
          }
    
    if (Bootstrap > MAX_BOOT){
        cerr << "Warning: The maximum bootstrap is " << MAX_BOOT << ", change to default " << DEF_BOOT << endl;
        Bootstrap = DEF_BOOT;
    }
    
    if ((Mode == 0) || (Mode == 1)){
       Check_Path(Outpath.c_str(), 1);
       if (Is_out_list) {
		if(Outpath[Outpath.size()-1] != '/')
	                Outlistfile = Outpath + ".list";
        	else
                	Outlistfile = Outpath.substr(0, Outpath.size()-1) + ".list";
	}
       if (Mode == 1) //list
       		Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);
       }
    return 0;
    }

/*
void Get_Random(int n, int * order){
    
    srand((int)time(NULL));
    
    int * order_table = new int [n];
    
    for(int i = 0; i< n; i++){
            order_table[i] = i;
            }
    for (int i = 0; i < n; i++ ){
         
         int r =(int)((float) (n-1-i)* rand()/(RAND_MAX+1.0)); 
         int temp = order_table[n-1-i]; //last 
         order_table[n-1-i] = order_table[r]; 
         order_table[r] = temp;
         }

    for (int i = 0; i < n; i++)
        order[i] = order_table[i];
    
    delete [] order_table;
    
    return;
    }
*/
void Get_Random(int n, int * order, int s, int loop){
    
    srand((int)time(NULL) + loop * 3000);
    
    int * order_table = new int [n];
    
    for(int i = 0; i< n; i++){
            order_table[i] = i;
            }
    for (int i = 0; i < s; i++ ){
         
         int r =(int)((float) (n-1-i)* rand()/(RAND_MAX+1.0)); 
         int temp = order_table[n-1-i]; //last 
         order_table[n-1-i] = order_table[r]; 
         order_table[r] = temp;
         }

    for (int i = 0; i < s; i++)
        order[i] = order_table[n-1-i];
    
    delete [] order_table;
    
    return;
    }

int Get_Int(float f){
    int n = (int) f;
    n += (f - n >= 0.5) ? 1 : 0;
    return n;
    }

unsigned int Load(const char * infilename, vector <string> & otus){
    
    hash_map <string, int, std_string_hash> otu_count;

    int count = OTU_parser.Load_file_to_hash(infilename, otu_count);
    for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++)
        for (int i = 0; i < miter->second; i ++)
            otus.push_back(miter->first);
    
    return otus.size();
    }

unsigned int Load(_Table_Format * table,  vector <string> & otus, int s){
         
         if (s >= table->Get_Sample_Size()) return 0;
         vector <string> otu_ids = table->Get_Feature_Names();
         vector <float> abds = table->Get_Abd(s);
         for (int i = 0; i < otu_ids.size(); i ++)
             for (int j = 0; j < abds[i]; j ++)
                 otus.push_back(otu_ids[i]);
         return otus.size();
         }

int Output(string outpath, string samname, int mode, hash_map <string, int, std_string_hash> count){
    
    string outfilename = outpath + "/";
    
    if (mode == 0)
            outfilename += "classification.rare.txt"; //Single mode
    else{
        Check_Path((outpath + "/" + samname).c_str(), 1);
        outfilename = outfilename + "/" + samname + "/" + "classification.rare.txt";
        }
    return OTU_parser.Output_hash_to_table(outfilename.c_str(), count, true);
    }

void Rare_Single_Boot(vector <string> otus, hash_map <string, int, std_string_hash> & otu_count, int s, int b){
    
    hash_map <string, float, std_string_hash> abd;
    
    if (s >= otus.size()){
          cerr << "Warning: Rarefaction sequence depth must be larger than input sequence number" << endl;
          return;
          } 
    
    for (int i = 0; i < b; i ++){
    
        int * order  = new int [s];
        Get_Random(otus.size(), order, s, i);
        for (int j = 0; j < s; j ++){
            if (abd.count(otus[order[j]]) == 0)
               abd[otus[order[j]]] == 0;
             abd[otus[order[j]]] ++;
             }
        }
    
    for (hash_map <string, float, std_string_hash> :: iterator miter = abd.begin(); miter != abd.end(); miter ++){
        miter->second = (float) Get_Int((float) miter->second / (float) b);
        }
    
    //correct abd
    float sum = 0;
    for (hash_map <string, float, std_string_hash> :: iterator miter = abd.begin(); miter != abd.end(); miter ++)
        sum += miter->second;
    
    float rate = (float) s / sum;
    //cout << rate << endl; //debug
    
    for (hash_map <string, float, std_string_hash> :: iterator miter = abd.begin(); miter != abd.end(); miter ++){
        miter->second *= rate;
        int a_count = Get_Int(miter->second);
        if (a_count > 0)
           otu_count[miter->first] = a_count;
        }
    }

int main(int argc, char * argv[]){        
    
    Parse_Para(argc, argv);
    
    Database.Read_Taxonomy(OTU_Taxa);
    
    OTU_parser = _OTU_Parser(Database);
    
    switch(Mode){
        case 0:
        case 1: {//single and list

          	int count_out = 0;
          	vector <string> sam_name_out;
          	//sam_name_out.push_back("samples"); // with id

          	for (int i = 0; i < Infilename.size(); i ++){
              	vector <string> otus;
              	hash_map <string, int, std_string_hash> otu_count;
              	int n = Load(Infilename[i].c_str(), otus);
              	cout << n << " input seqs loaded" << endl;
              	if (Seq_depth >= n){
                 	cerr << "Warning: Sample " << i + 1 << ": Rarefaction sequence depth must be larger than input sequence number, removed from the output" << endl;
                	continue;
                 	}
              	Rare_Single_Boot(otus, otu_count, Seq_depth, Bootstrap);
              	cout << Output(Outpath, Sam_name[i], Mode, otu_count) << " random seqs output" << endl;
              	sam_name_out.push_back(Sam_name[i]);
              	count_out ++;
              	}

          	if ((Is_out_list) && (Infilename.size() > 1))
          	Make_list(Outlistfile.c_str(), Outpath.c_str(), sam_name_out, 3);

          	cout << Infilename.size() << " Sample(s) input" << endl;
          	cout << count_out << " Sample(s) output" << endl;
          	}
          	break;
          
     	case 2:{//table
          
          	_Table_Format input_table(Tablefilename.c_str());
          	vector <string> otu_ids = input_table.Get_Feature_Names();
          	vector <string> sample_names = input_table.Get_Sample_Names();
          	_Table_Format output_table(otu_ids);
          
          	for (int i = 0; i < input_table.Get_Sample_Size(); i ++){
              	vector <string> otus;
              	hash_map <string, int, std_string_hash> otu_count;
              	int n = Load(&input_table, otus, i);
              	cout << n << " input seqs loaded" << endl;
              	if (Seq_depth >= n){
                    cerr << "Warning: Sample " << i + 1 << ": Rarefaction sequence depth must be larger than input sequence number, removed from the output" << endl;
                    continue;
                    }
              	Rare_Single_Boot(otus, otu_count, Seq_depth, Bootstrap);
              	vector <float> abds;
              	int output_n = 0;
              	for (int j = 0; j < otu_ids.size(); j ++)
                  	if (otu_count.count(otu_ids[j]) != 0) {
                                                output_n += otu_count[otu_ids[j]];
                                                abds.push_back(otu_count[otu_ids[j]]);
                                                }
                  	else abds.push_back(0);
                  
              	cout << output_n << " random seqs output" << endl;
              	output_table.Add_Abd(abds, sample_names[i]);
              	}
          
          	output_table.Filter_Empty();
          	output_table.Output_Table(Outpath.c_str());
          
         	cout << input_table.Get_Sample_Size() << " Sample(s) input" << endl;
          	cout << output_table.Get_Sample_Size() << " Sample(s) output" << endl;
          	}
    	break;
    } 
    return 0;
    }
