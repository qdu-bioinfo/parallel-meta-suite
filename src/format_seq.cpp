// Updated at July 29, 2021
// Updated by Yuzhu Chen
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include "utility.h"
#include "version.h"

using namespace std;

string Listfilename;
string Listprefix;
vector <string> Infilelist;

int Mode = -1; //0: single; 1: list

int printhelp(){

    cout << "Format-seq version: " << Version << endl;
    cout << "\tSequence format check and reformat for Parallel-Meta Suite" << endl;
    cout << "Usage:" << endl;
    cout << "PM-format-seq [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -i Input single sequence file in FASTA or FASTQ format" << endl;
    cout << "\t  or" << endl;
    cout << "\t  -l Input sequence files list" << endl;
    cout << "\t  -p List file path prefix for '-l' [Optional for -l]" << endl;
    
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
    
bool Check_Id(string & id, int count){
     string id_new;
     string temp;
     bool p = true;
     
     stringstream strin(id);
     strin >> id_new;                                 
     while(strin >> temp){
                 id_new += '_';
                 id_new += temp;
                 p = false;
                 }
     if (id_new.size() == 1){
                       id_new += "seq_noname_";
                       stringstream strout;
                       string strout_id;
                       strout << count;
                       strout >> strout_id;
                       id_new += strout_id;
                       p = false;
                       }
     if (!p) id = id_new;
     return p;
     }

void Reformat(string infilename){
    ifstream infile(infilename.c_str(),ios::in);
    if(!infile){
                cerr << "Error: Cannot open infile : " << infilename << endl;
                return;
                }
    string outfilename = infilename+".tmp";
    ofstream outfile(outfilename.c_str(),ios::out);
    if(!outfile){
                 cerr << "Error: Cannot open infile : " << outfilename << endl;
                 return;
                 }
    
    string buffer;
    string id;
    string seq;
    unsigned int count = 0;
    unsigned int seq_count = 0;
    unsigned int seq_out_count = 0;
        
    bool is_ok=true;
    
   int format = Check_Format(infilename.c_str()); //0: fasta; 1: fastq
   if (format < 0) return;
   
   if (format == 0){ //fasta
      while(getline(infile, buffer)){
              if (buffer[0] == '>'){ // new seq
                            seq_count ++;
                            if (seq.size() > 0){
                                            outfile << id << endl;
                                            outfile << seq << endl;
                                            seq_out_count ++;                                            
                                            }
                            else if (seq_count != 1)
                                    is_ok = false;
                            
                            id = buffer;
                            is_ok = Check_Id(id, seq_count) && is_ok;
                            seq = "";
                            }
              else{
                   is_ok = (seq.size() == 0) && is_ok;
                   seq += buffer;
                   } 
              }
      if (seq.size() > 0){
                     outfile << id << endl;
                     outfile << seq << endl;
                     seq_out_count ++;                                            
                     }
      else is_ok = false;
      }
   
   else if (format == 1){ //fastq
        string qual_id;
        string qual;
        unsigned int seq_line_count = 0;
   
        while(getline(infile, buffer)){
                 if(buffer[0] == '@'){                              
                            seq_count ++;
                                                        
                            id = buffer;
                            is_ok = Check_Id(id, seq_count) && is_ok;
                            seq = "";
                            qual_id = "";
                            qual = "";  
                            seq_line_count = 0;
                            
                            while(getline(infile, buffer))
                                                  if (buffer[0] =='+'){
                                                                qual_id = buffer;
                                                                break;
                                                                }
                                                  else {
                                                       seq += buffer;
                                                       seq_line_count ++;
                                                       }
                            
                            if (seq_line_count > 1)
                               is_ok = false;
                            
                            //get qual
                            for (int i = 0; i < seq_line_count; i ++){
                                getline(infile, buffer);
                                qual += buffer;
                                }
                                                            
                            outfile << id << endl;
                            outfile << seq << endl;
                            outfile << qual_id << endl;
                            outfile << qual << endl;
                            seq_out_count ++;                                               
                            }             
                 
              else is_ok = false;
        }
   }
   cout << seq_count << " sequences in total" << endl;
   cout << seq_out_count << " sequences output" << endl;
   
   if(is_ok){
	string Opt="rm "+infilename+".tmp";
	system(Opt.c_str());
	cout<<infilename<<" is OK for Parallel-Meta Suite" << endl;
	}
	else{
         string Opt="cp "+infilename+" "+infilename+".bk";
	     system(Opt.c_str());
	     Opt="cp "+outfilename+" "+infilename;
	     system(Opt.c_str());
		 Opt="rm "+infilename+".tmp";
	     system(Opt.c_str());
         cout<<infilename<<" is formated for Parallel-Meta Suite" << endl;
         cout<<"The original file backup is "<<infilename<<".bk"<<endl; 
    }
    infile.close();
    outfile.close();        
    }
    
int main(int argc, char * argv[]){
    Parse_Para(argc, argv);
             
    for (int i = 0; i < Infilelist.size(); i ++)
                  Reformat(Infilelist[i]);         
    cout << Infilelist.size() << " Sample(s) in total" << endl;              
    return 0;
    }
