// Updated at July 31, 2021
// Updated by Yuzhu Chen
// Code by: Yuzhu Chen, Honglei Wang, Gongchao Jing, Xiaoquan Su
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// _Table_Format,supported

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <omp.h>
#include <math.h>
#include <unistd.h>

#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>

#include "utility.h"
#include "version.h"
#include "table_format.h"

using namespace std;

vector< vector<int> > Otu_table;
vector<string> Sample;
int Bot=20;
int Skn=100;
int Threadnum=0;
int Mode = 0; // 0,1: OTU
string Outdirectory = "result";
string Resultname="out";
string Infilename;
bool Annotation=false;

int printhelp(){
    cout << "Rare-curv version : " << Version << endl;
    cout << "\tMake the rarefaction curve" << endl;
    cout << "Usage: " << endl;
    cout << "PM-rare-curv [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
	cout << "\t  -i or -T (upper) Input feature count table (*.OTU.Count)" << endl;
    
    cout << "\t[Output options]" << endl;
	cout << "\t  -o Output file directory, default is \"result\"" << endl;
    cout << "\t  -p Prefix name of output, default is \"out\""<<endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -b The bootstrap value, default is 20" << endl;
	cout << "\t  -s The rarefaction step, default is 100" << endl;
	cout << "\t  -l The rarefaction curve label, T is enable and F is disable, default is F"<<endl;
	cout << "\t  -t Number of thread, default is auto" <<endl;
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
		case 'i': Infilename = argv[i+1]; Mode = 0; break;
        case 'T': Infilename = argv[i+1]; Mode = 1; break;
        case 'B': Infilename = argv[i+1]; Mode = 2; break;
		case 'o': Outdirectory = argv[i+1];break;
        case 'p': Resultname = argv[i+1];break;
		case 't': Threadnum = atoi(argv[i+1]);break;
		case 's': Skn = atoi(argv[i+1]);break;
		case 'b': Bot = atoi(argv[i+1]);break;
		case 'h': printhelp();
		case 'l': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Annotation=true; break;
		default : cerr << "Error: Unrec argument " << argv[i] << endl;
			printhelp();
			break;
		}
		i=i+2;

	}
   int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
   if ((Threadnum <= 0) || (Threadnum > Max_core_number)){
                    Threadnum = Max_core_number;
                    }  
   
   if (Bot <= 0){
         cerr << "Warning: The bootstrap value should be larger than 0, change to default (20)";
         Bot = 20;
         }
}

void Load_Table(const char* infilename){
     ifstream infile(infilename, ifstream::in);
     if (!infile){
                  cout<<"Error: Cannot open file : " << infilename << endl;
                  return ;
                  }
     string buffer, Otu_id;
     vector<int> Otu_name;
     getline(infile,buffer);
	 stringstream strin(buffer);
	 strin>>Otu_id;
	 int Tri=0;
	 while(strin>>Otu_id){
          Otu_name.push_back(Tri);
          Tri++;
          }
     while(getline(infile,buffer)){
       vector<int> Otu_count;
       stringstream Str_Otu_count(buffer);
       string SampleID;
       Str_Otu_count>>SampleID;
       Sample.push_back(SampleID);
       string Taxa_count;
       int i=0;
       while(Str_Otu_count>>Taxa_count){
            if (atoi(Taxa_count.c_str())!=0){
                 for (int j=0;j<atoi(Taxa_count.c_str());j++){
                     Otu_count.push_back(Otu_name[i]);
                     }
                 } 
                 i++;
            }
       Otu_table.push_back(Otu_count);
       }
     }
     
void Load_Table(_Table_Format table){
     Sample= table.Get_Sample_Names();
            
     for(int i = 0; i < table.Get_Sample_Size(); i++) {
     	vector<float> Taxa_count = table.Get_Abd(i);
     	vector<int> Otu_count;
     	for(int i = 0; i < Taxa_count.size(); i++) {
     		if(Taxa_count[i] != 0) {
     			for(int j = 0; j < Taxa_count[i]; j++) {
     				Otu_count.push_back(i);
				 }
			 }
		 }
		 Otu_table.push_back(Otu_count);
	 }
}
     
void Calculater_Rarefraction(){
     
     string Outfilename_Shannon=Outdirectory+"/Shannon.txt";
     string Outfilename_Shannon_max=Outdirectory+"/Shannon_max.txt";
     string Outfilename_Observe_otu=Outdirectory+"/Observe_otu.txt";
     string Outfilename_Observe_otu_max=Outdirectory+"/Observe_otu_max.txt";
     ofstream Outfile_Shannon_max(Outfilename_Shannon_max.c_str(),ios::out);
     ofstream Outfile_Shannon(Outfilename_Shannon.c_str(),ios::out);
     ofstream Outfile_Observe_max(Outfilename_Observe_otu_max.c_str(),ios::out);
     ofstream Outfile_Observe(Outfilename_Observe_otu.c_str(),ios::out);
     Outfile_Shannon_max<<"SampleID"<<"\t"<<"xm"<<"\t"<<"ym"<<endl;
     Outfile_Shannon<<"SampleID"<<"\t"<<"x"<<"\t"<<"y"<<endl;
     Outfile_Observe_max<<"SampleID"<<"\t"<<"xm"<<"\t"<<"ym"<<endl;
     Outfile_Observe<<"SampleID"<<"\t"<<"x"<<"\t"<<"y"<<endl;
 
     int line_max=0;
	 for(int i=0;i<Otu_table.size();i++){
                if(Otu_table[i].size()>=line_max)
                line_max=Otu_table[i].size();
                }
                
     int **Observe_otu=new int*[Otu_table.size()];
	 double **Shannon_vec=new double*[Otu_table.size()];
	 
	 for(int i=0;i<Otu_table.size();i++){
             Shannon_vec[i]=new double[1+line_max/Skn];
             Observe_otu[i]=new int[1+line_max/Skn];
             }
    
     #pragma omp parallel for schedule(dynamic, 1)
     for(int i=0; i<Bot*Otu_table.size(); i++){
             int real_serials=i-floor((float)i/Otu_table.size())*Otu_table.size();
             int times=Otu_table[real_serials].size();
             set<int> Taxa;
             vector<int> line_taxa;
             vector<float> line_shannon;
             map<int,int> Map_shannon;
             
             srand((unsigned)time(0));
  			 for(int j=0;j<times;j++){
				     int ran_num=j+rand()%(times-j);
				     swap(Otu_table[real_serials][j],Otu_table[real_serials][ran_num]);
			         }
             for(int j=0;j<times;j+=Skn){
                  int kk;
                  if(j==0){kk=Skn;}else{kk=j;}
                  for(int k=j; k>kk-Skn; k--){
                             Taxa.insert(Otu_table[real_serials][k]);
                             if(Map_shannon.find(Otu_table[real_serials][k])!=Map_shannon.end()){
                                                                                                 Map_shannon[Otu_table[real_serials][k]]++;
                                                                                                 }
                                                                                                 else
                                                                                                 Map_shannon[Otu_table[real_serials][k]]=1;
                                   }
                     line_taxa.push_back(Taxa.size());
                     
                     float Shannon=0.0;
                     for(map<int,int>::iterator Map_shannon_it=Map_shannon.begin();Map_shannon_it!=Map_shannon.end();Map_shannon_it++){
                                                       Shannon+=(float(Map_shannon_it->second)/(float(j+1)))*(log((long double)(float(Map_shannon_it->second)/float(j+1))));
                                                       }
                     line_shannon.push_back(Shannon);
                     }//j
             if(i<Otu_table.size()){
                                    for(int j=0; j<line_shannon.size(); j++){
                                            Shannon_vec[i][j]=line_shannon[j];
                                            Observe_otu[i][j]=line_taxa[j];
                                            }
                                    }else{for(int j=0; j<line_shannon.size(); j++){
                                                  Shannon_vec[real_serials][j]+=line_shannon[j];
                                                  Observe_otu[real_serials][j]+=line_taxa[j];
                                                  }
                                          }
             }//i

     for(int i=0;i<Sample.size();i++){
             int ii=0;
             for(int j=0;j<Otu_table[i].size();j+=Skn){
                     Outfile_Observe<<Sample[i]<<"\t"<<j<<"\t"<<float(Observe_otu[i][ii])/float(Bot)<<endl;
                     Outfile_Shannon<<Sample[i]<<"\t"<<j<<"\t"<<-float(Shannon_vec[i][ii])/float(Bot)<<endl;
                     ii++;
                     }
             Outfile_Observe_max<<Sample[i]<<"\t"<<Otu_table[i].size()<<"\t"<<float(Observe_otu[i][ii-1])/float(Bot)<<endl;
             Outfile_Shannon_max<<Sample[i]<<"\t"<<Otu_table[i].size()<<"\t"<<-float(Shannon_vec[i][ii-1])/float(Bot)<<endl;
             }
      
   	Outfile_Observe.close();
	Outfile_Observe_max.close();
	Outfile_Shannon.close();
	Outfile_Shannon_max.close();
    }

int main(int argc, char *argv[]){
    Parse_Para(argc, argv);
    Check_Path(Outdirectory.c_str(), 1);
    
    switch(Mode) {
    	case 0: {
    		Load_Table(Infilename.c_str());
			break;
		}
		
    	case 1: {
    		_Table_Format table(Infilename.c_str());
			Load_Table(table);
			break;
		}
	}
	
    omp_set_num_threads(Threadnum);
    Calculater_Rarefraction();

	string Radd=Check_Env()+"/Rscript/PM_rarefaction.R ";
	string Opt;
	
	if(Annotation)
	 Opt="Rscript "+Radd+"-o "+Outdirectory+" -p "+Resultname+" -a T";
	else
	 Opt="Rscript "+Radd+"-o "+Outdirectory+" -p "+Resultname+" -a F";
	 
	system(Opt.c_str()); 
	return 0;
    }
     
