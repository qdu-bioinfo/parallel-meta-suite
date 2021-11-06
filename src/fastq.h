// Updated at Nov 25, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

int Fastq_2_Fasta(const char * infilename, const char * outfilename){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr <<"Error: Open infile error :" << infilename << endl;
                 //system("pause");
                 exit(0);
                 }
    
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                 cerr <<"Error: Open outfile error :" << outfilename << endl;
                 //system("pause");
                 exit(0);
                 }
    
    string buffer;
    unsigned int count = 0;
    
    while (getline(infile, buffer)){
          
          //count ++;
          if (buffer.size() == 0) continue;
          if (buffer[0]!='@'){
                              cerr << "Error: Input file Format Error: Line :" << count << endl;
                              //system("pause");
                              exit(0);
                              }
          buffer[0] = '>';
          outfile << buffer << endl;
          count ++;
          
          getline(infile, buffer);
          
          outfile << buffer << endl;
          
          getline(infile, buffer);
          
          if (buffer[0]!='+'){
                              cerr << "Error: Input file Format Error: Line :" << count << endl;
                              //system("pause");
                              exit(0);
                              }
          
          getline(infile, buffer);
                                                     
                              }
          
         infile.close();
         infile.clear();
         outfile.close();
         outfile.clear();
         
         return count; 
          
         }
/*
int main(int argc, char * argv[]){
    
    
    cout << Fastq_2_Fasta(argv[1], argv[2]) << endl;
    
    system("pause");
    
    return 0;
    
    
    }
*/
