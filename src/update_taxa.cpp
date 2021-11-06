//Update Taxonomy of GG
//For PM 3 or above
//For GG_13-8
// Updated at July 31, 2021
// Updated by Yuzhu Chen
// Version 3.4.1 or above with _OTU_Parser
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include "utility.h"
#include "version.h"
#include "hash.h"
#include "db.h"
#include "otu_parser.h"

using namespace std;

_PMDB Database;

string Listfilename;
string Listprefix;
vector <string> Infilelist;

int Mode = -1; //0: single; 1: list

int printhelp(){

    cout << "Update-taxa version: " << Version << endl;
    cout << "\tUpdate the taxonomy annotation to the latest version" << endl;
    cout << "Usage:" << endl;
    cout << "PM-update-taxa [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) ref database, " << _PMDB::Get_Args() << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single file" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list" << endl;
    cout << "\t  -p List file path prefix [Optional for -l]" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;

    exit(0);

    return 0;

    };

void Parse_Para(int argc, char * argv[]){
    
    if(argc==1){
                printhelp();
                }

    int i=1;
    while(i<argc){
                 if(argv[i][0]!='-'){
                                     cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                                     exit(0);
                                     }
                 switch(argv[i][1]){
                                    case 'D': Database.Set_DB(argv[i+1][0]); break;
                                    case 'i' : if (Mode == 1){
                                                        cerr << "Error: -i conflists with -l" << endl;
                                                        exit(0);
                                                        } 
                                               else {                                                    
                                                    Infilelist.push_back(argv[i+1]);
                                                    //Infilename = argv[i+1]; 
                                                    Mode = 0;
                                                    }                                               
                                               break;
                                               
                                    case 'l' : if (Mode == 0){
                                                        cerr << "Error: -l conflists with -i" << endl;
                                                        exit(0);
                                                        } 
                                               else {
                                                    Listfilename = argv[i+1];
                                                    Mode = 1;
                                                    }
                                               break;
                                    case 'p' : Listprefix = argv[i+1]; break;
                                    case 'h' : printhelp(); break;           
                                    default : cerr << "Error: Unrec argument " << argv[i] << endl; break;
					      printhelp();
					      break;
                                    }
		i=i+2;
                  
                  }
    
    
    if (Mode == 1)
       Load_List(Listfilename.c_str(), Infilelist, Listprefix); 
    } 

    
int main(int argc, char * argv[]){
            
    Parse_Para(argc, argv);
    _OTU_Parser otu_parser(Database);
    
    for (int i = 0; i < Infilelist.size(); i ++){
         //copy the inputfile
         string command = "cp "+Infilelist[i] + " " + Infilelist[i] + ".bk";
	     system(command.c_str());
	     cout << otu_parser.Update_class_taxa((Infilelist[i] + ".bk").c_str(), Infilelist[i].c_str()) << " sequence output" << endl;
        }
    
    return 0;
    }
