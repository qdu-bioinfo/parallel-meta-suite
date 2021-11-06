// Updated at Aug 21, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
// otu_parser

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string.h>
#include <omp.h>

#include "utility.h"
#include "version.h"
#include "hash.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"

#ifndef class_func_h
#define class_func_h

using namespace std;

#define BUFF 10000000
#define STRLEN 1000

class _KO_Index_Copy{
    
public:
    _KO_Index_Copy(){
        Index = 0;
        Copy = 0;
    }
    _KO_Index_Copy(int index, float copy){
        Index = index;
        Copy = copy;
    }
    
    int Index;
    float Copy;
};



class _KO{
    
public:
    _KO(){
        Name = "K_NA";
        Des = "None";
        Pathway = "None";
        Index = 0;
    }
    
    _KO(string name, string des, string pathway, int index){
        Name = name;
        Des = des;
        Pathway = pathway;
        Index = index;
    }
    
    string Get_Name(){
        return Name;
    }
    
    string Get_Des(){
        return Des;
    }
    
    string Get_Pathway(){
        return Pathway;
    }
    
    int Get_Index(){
        return Index;
    }
    
private:
    string Name;
    string Des;
    string Pathway;
    int Index;
};

class _KO_OTU_Table_All{

public:
    _KO_OTU_Table_All(){
        Sample_count = 0;
    }
    
    _KO_OTU_Table_All(_PMDB db, int count, int mode);
    int Load_KO_Id(const char * idfilename, const char * desfilename, const char * pathwayfilename);
    int Load_OTU_KO_Index(const char * abdfilename);
    int Load_Sample(const char * infilename, string sample_name, int sample);
    int Load_Sample_Multi(vector<string> infilename, vector <string> sample_name, int coren);
    int Load_Sample_By_OTU_Table(_Table_Format * table, int coren);
    int Load_Sample_By_Single_KO_Table(const char * infilename, string sample_name, int sample);
    int Load_Sample_By_KO_Table(_Table_Format * table);
    int Output(const char * outfilename);
    int Output_Multi(vector <string> outfilename);
    int Output_By_Category(const char * outfilename, int level, float max, float min);
    
protected:
    int Sample_count;
    string * Sample_names;
    vector <_KO> KOs; //KO information
    hash_map <string, int, std_string_hash> KO_Index; //KO number to index
    hash_map <string, vector<_KO_Index_Copy>, std_string_hash> OTU_KO_Index; //OTU to KO
    float ** KO_Abd; //Sample, KO
    _OTU_Parser Otu_parser;
    
    virtual int Load_Sample_By_OTU(hash_map <string, int, std_string_hash> * otu_seq_count, int sample);
    };

_KO_OTU_Table_All::_KO_OTU_Table_All(_PMDB db, int count, int mode){ //mode 0: norm; mode 1: sim
        
        string ko_id_file = db.Get_Func_Id();
        string ko_abd_file = db.Get_Func();
        string ko_des_file = db.Get_Func_Des();
        string ko_pw_file = db.Get_Func_Pw();
                
        Load_KO_Id(ko_id_file.c_str(), ko_des_file.c_str(), ko_pw_file.c_str());
        
        Sample_count = count;
        Sample_names = new string[count];
        
        if (mode ==  0)
            Load_OTU_KO_Index(ko_abd_file.c_str());
            
        
        Otu_parser = _OTU_Parser(db);
        
        //Init KO_Abd
        KO_Abd = new float * [count];
        for (int i = 0; i < count; i ++){
            KO_Abd[i] = new float [KOs.size()];
            for (int j = 0; j < KOs.size(); j ++)
                KO_Abd[i][j] = 0;
            }
        
    }
    
int _KO_OTU_Table_All::Load_KO_Id(const char * idfilename, const char * desfilename, const char * pathwayfilename){
        ifstream in_idfile(idfilename, ifstream::in);
        if (!in_idfile){
            cerr << "Error: Open KO ID file error : " << idfilename << endl;
            exit(0);
        }
        
        ifstream in_desfile(desfilename, ifstream::in);
        if (!in_desfile){
            cerr << "Error: Open KO description file error : " << desfilename << endl;
            exit(0);
        }
        
        ifstream in_pwfile(pathwayfilename, ifstream::in);
        if (!in_pwfile){
            cerr << "Error: Open KO pathway file error : " << pathwayfilename << endl;
            exit(0);
        }
        
        string buffer_id;
        string buffer_des;
        string buffer_pw;
        int i = 0;
        while(getline(in_idfile, buffer_id)){
            getline(in_desfile, buffer_des);
            getline(in_pwfile, buffer_pw);
            KOs.push_back(_KO(buffer_id, buffer_des, buffer_pw, i));
            KO_Index[buffer_id] = i;
            i ++;
        }
        
        in_idfile.close();
        in_idfile.clear();
        
        in_desfile.close();
        in_desfile.clear();
        
        in_pwfile.close();
        in_pwfile.clear();
        
        return KOs.size();
        
    }

   
int _KO_OTU_Table_All::Load_OTU_KO_Index(const char * abdfilename){
        ifstream infile(abdfilename, ifstream::in);
        if (!infile){
            cerr << "Error: Open KO Abduncance file error : " << abdfilename << endl;
            exit(0);
        }
        
        string buffer;
        while(getline(infile, buffer)){
            
            stringstream strin(buffer);
            string id;
            int index;
            float copy;
            
            strin >> id;
            
            while(strin >> index){
                strin >> copy;
                //OTU_KO_Index[id].push_back(_KO_Index_Copy(index, copy));
                //debug
                OTU_KO_Index[id].push_back(_KO_Index_Copy(index - 1, copy));
            }
            
        }
        
        infile.close();
        infile.clear();
        
        return OTU_KO_Index.size();
        
    }

/*
int _KO_OTU_Table_All::Load_OTU_KO_Index(const char * abdfilename){
        ifstream infile(abdfilename, ifstream::in);
        if (!infile){
            cerr << "Error: Open KO Abduncance file error : " << abdfilename << endl;
            exit(0);
        }
        
        string buffer;
        while(getline(infile, buffer)){
            
            //strtok
            char buf[buffer.size()-1];
            buffer.copy(buf,buffer.size()-1);
            //string del = "\t";
            char *p;
            
            string id;
            int index;
            float copy;
            
            //split for id
            p = strtok(buf, "\t");
            id = p;

            while(p){
                //split for index
                p = strtok(NULL, "\t");
                if(p != NULL)
                    index = atoi(p);
                //split for copy
                p = strtok(NULL, "\t");
                if(p != NULL)
                    copy = atof(p);
               OTU_KO_Index[id].push_back(_KO_Index_Copy(index, copy));
            }
            
        }
        
        infile.close();
        infile.clear();
        
        return OTU_KO_Index.size();
        
    }
 */   
int _KO_OTU_Table_All::Load_Sample(const char * infilename, string sample_name, int sample){
        
        Sample_names[sample] = sample_name;
        hash_map <string, int, std_string_hash> otu_seq_count;
        
        unsigned int count = Otu_parser.Load_file_to_hash(infilename, otu_seq_count);
        Load_Sample_By_OTU(&otu_seq_count, sample);
        
        cout << count << " seqs are loaded" << endl;
        return count;
        
    }
    
int _KO_OTU_Table_All::Load_Sample_Multi(vector<string> infilename, vector <string> sample_name, int coren){
        
        hash_map <string, int, std_string_hash> * otu_seq_count = new hash_map <string, int, std_string_hash> [infilename.size()];
        
        for(int i = 0; i < infilename.size(); i ++){
            
            Sample_names[i] = sample_name[i];
            
            unsigned int count = Otu_parser.Load_file_to_hash(infilename[i].c_str(), otu_seq_count[i]);
            
            cout << count << " seqs are loaded" << endl;
            
            }
        
        omp_set_num_threads(coren);
        #pragma omp parallel for schedule(dynamic, 1)
        for(int i = 0; i < infilename.size(); i ++)
            Load_Sample_By_OTU(&(otu_seq_count[i]), i);
        
        return infilename.size();
        }
    
int _KO_OTU_Table_All::Load_Sample_By_OTU_Table(_Table_Format * table, int coren){
        
        hash_map <string, int, std_string_hash> * otu_seq_count = new hash_map <string, int, std_string_hash> [table->Get_Sample_Size()];
        
        vector <string> otus = table->Get_Feature_Names();
        vector <string> sam_names = table->Get_Sample_Names();
        
        for (int i = 0; i < table->Get_Sample_Size(); i ++){
            
            Sample_names[i] = sam_names[i];
            //vector <float> abd= table->Get_Abd(i);
            
            for (int j = 0; j < table->Get_Feature_Size(); j ++)
                if (table->Get_Abd_By_Order(i, j) > 0){
                    string a_otu_j = Check_OTU(otus[j]);
                    otu_seq_count[i][a_otu_j] = (int) table->Get_Abd_By_Order(i, j);
                    }
            }
        
        omp_set_num_threads(coren);
        #pragma omp parallel for schedule(dynamic, 1)
        for(int i = 0; i < table->Get_Sample_Size(); i ++)
            Load_Sample_By_OTU(&(otu_seq_count[i]), i);
        
        return table->Get_Sample_Size();
        }
    
int _KO_OTU_Table_All::Load_Sample_By_Single_KO_Table(const char * infilename, string sample_name, int sample){
        
        ifstream infile(infilename, ifstream::in);
        if (!infile){
            cerr << "Error: Open input file error : " << infilename << endl;
            exit(0);
        }
        
        int count = 0;
        
        Sample_names[sample] = sample_name;
        string buffer;
        getline(infile, buffer);
        while(getline(infile, buffer)){
            
            stringstream strin(buffer);
            string ko;
            float ko_count;
            
            strin >> ko >> ko_count;
            
            if (KO_Index.count(ko) == 0) continue;
            int index = KO_Index[ko];
            
            KO_Abd[sample][index] = ko_count;
            
            count ++;
        }
        
        infile.close();
        infile.clear();
        
        return count;;
    }
    
int _KO_OTU_Table_All::Load_Sample_By_KO_Table(_Table_Format * table){
        
        vector <string> kos = table->Get_Feature_Names();
        vector <string> sam_names = table->Get_Sample_Names();
        
        for (int i = 0; i < table->Get_Sample_Size(); i ++){
            
            Sample_names[i] = sam_names[i];
            //vector <float> abd= table->Get_Abd(i);
            
            for (int j = 0; j < table->Get_Feature_Size(); j ++){
                if (table->Get_Abd_By_Order(i, j) <= 0) continue;
                if (KO_Index.count(kos[j]) == 0) continue;
                KO_Abd[i][KO_Index[kos[j]]] = table->Get_Abd_By_Order(i, j);
                }
            }

        return table->Get_Sample_Size();
        }
    
int _KO_OTU_Table_All::Output(const char * outfilename){
        
        ofstream outfile(outfilename, ofstream::out);
        if (!outfile){
            cerr << "Error: Cannot open output file : " << outfilename << endl;
            return 0;
        }
        
        int count = 0;
        
        outfile << "Gene";
        for (int i = 0; i < Sample_count; i ++)
            if (Sample_count == 1) outfile << "\tGene_Count";
            else outfile << "\t" << Sample_names[i] << "_Gene_Count";
        
        outfile << "\tGene_Description\tKEGG_Pathway" << endl;
                
                
        for (int i = 0; i < KOs.size(); i ++){
            
            bool No_Zero_Check = false;
            for (int j = 0; j < Sample_count; j ++)
                if (KO_Abd[j][i] > 0) No_Zero_Check = true;
            if (!No_Zero_Check) continue;
            count ++;    
            outfile << KOs[i].Get_Name();
                    
            for (int j = 0; j < Sample_count; j ++)
                        outfile << "\t" << KO_Abd[j][i];
                    
            outfile << "\t" << KOs[i].Get_Des() << "\t" << KOs[i].Get_Pathway() << endl;
            }
        
        
        outfile.close();
        outfile.clear();
        
        return count;
    }
    
int _KO_OTU_Table_All::Output_Multi(vector <string> outfilename){
        
        int count = 0;
        
        for (int i = 0; i < KOs.size(); i ++){
            
            bool No_Zero_Check = false;
            for (int j = 0; j < Sample_count; j ++)
                if (KO_Abd[j][i] > 0) No_Zero_Check = true;
            if (!No_Zero_Check) continue;
            count ++;
        }
        
        for (int i = 0; i < Sample_count; i ++){
        
            ofstream outfile(outfilename[i].c_str(), ofstream::out);
            if (!outfile){
            cerr << "Error: Cannot open output file : " << outfilename[i] << endl;
            return 0;
            }
            outfile << "Gene\tGene_Count\tGene_Description\tKEGG_Pathway" << endl;
        
            for (int j = 0; j < KOs.size(); j ++){
            
                if (KO_Abd[i][j] <= 0) continue;
            
                outfile << KOs[j].Get_Name() << "\t" << KO_Abd[i][j];
                outfile << "\t" << KOs[j].Get_Des() << "\t" << KOs[j].Get_Pathway() << endl;
                }
        
            outfile.close();
            outfile.clear();
            }
    
        return count;
        }

    
int _KO_OTU_Table_All::Output_By_Category(const char * outfilename, int level, float max, float min){ //2 output, abd and count
        
        
        hash_map <string, vector<string>, std_string_hash> pw_ko;
        hash_map <string, int, std_string_hash> pw_index;
        
        // open output
        //Abd file
        string outfilename_abd = outfilename;
        outfilename_abd += ".Abd";
        ofstream outfile_abd(outfilename_abd.c_str(), ofstream::out);
        if (!outfile_abd){
            cerr << "Error: Cannot open output file : " << outfilename_abd << endl;
            return 0;
        }
        //Count file
        string outfilename_count = outfilename;
        outfilename_count += ".Count";
        ofstream outfile_count(outfilename_count.c_str(), ofstream::out);
        if (!outfile_count){
            cerr << "Error: Cannot open output file : " << outfilename_count << endl;
            return 0;
        }
        
        //add_ko
        for (int i = 0; i < KOs.size(); i ++){
            
            if (level >= 3){
                
                pw_ko[KOs[i].Get_Name()].push_back(KOs[i].Get_Name());
                continue;
                
            }
            
            string pw = KOs[i].Get_Pathway();
            int l = 0;
            int begin = 0;
            int end = 0;
            while(end <= pw.size()){
                
                if ((pw[end] == ';') || (end == pw.size())){
                    if (l == level){
                        string pw_name = pw.substr(begin, end - begin);
                        pw_ko[pw_name].push_back(KOs[i].Get_Name());
                    }
                    begin = end + 1;
                    l ++;
                }
                else if (pw[end] == '|'){
                    l = 0;
                    begin = end + 1;
                }
                
                else if ((pw[end] == ' ') || (pw[end] == '\''))
                    pw[end] = '_';
                end ++;
            }
        }
        
        //make pw index
        int index = 0;
        for (hash_map <string, vector<string>, std_string_hash> ::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++){
            pw_index[miter->first] = index;
            index ++;
        }
        
        //calc the pathway sum
        float ** pw_count = new float * [Sample_count];
        float ** pw_abd = new float * [Sample_count];
        
        for (int i = 0; i < Sample_count; i ++){
            pw_count[i] = new float [pw_ko.size()];
            pw_abd[i] = new float [pw_ko.size()];
            
            for (int j = 0; j < pw_ko.size(); j ++){
                pw_count[i][j] = 0;
                pw_abd[i][j] = 0;
                }
                
            }
        
        for (int i = 0; i < Sample_count; i ++)
            for (hash_map <string, vector<string>, std_string_hash> ::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++){
                for (int j = 0; j < (miter->second).size(); j ++)
                    pw_count[i][pw_index[miter->first]] += KO_Abd[i][KO_Index[(miter->second)[j]]];
                }
       
        //norm for abd
        for (int i = 0; i < Sample_count; i ++){
            float sum = 0;
            for (hash_map <string, vector<string>, std_string_hash> ::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++)
                sum += pw_count[i][pw_index[miter->first]];
            for (hash_map <string, vector<string>, std_string_hash> ::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++)
                if (sum > 0)
                   pw_abd[i][pw_index[miter->first]] = pw_count[i][pw_index[miter->first]] / sum;
                else pw_abd[i][pw_index[miter->first]] = 0;
            }
        
        //check no zero
        bool * No_Zero_Check = new bool[pw_ko.size()];
        int count = 0;
        for (hash_map <string, vector<string>, std_string_hash> ::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++){
            index = pw_index[miter->first];
            No_Zero_Check[index] = false;
        
            for (int i = 0; i < Sample_count; i ++)
                if (pw_abd[i][index] > 0)
                    No_Zero_Check[index] = true;
            
            if (No_Zero_Check[index]) count ++;
            }
        
        //output
        outfile_abd << "Sample";
        outfile_count << "Sample";
        
        for (hash_map <string, vector<string>,std_string_hash>::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++){
            index = pw_index[miter->first];
            if (!No_Zero_Check[index]) continue;
            outfile_abd << "\t" << miter->first;
            outfile_count << "\t" << miter->first;
        }
        outfile_abd << endl;
        outfile_count << endl;

        for (int i = 0; i < Sample_count; i ++){
            
            outfile_abd << Sample_names[i];
            outfile_count << Sample_names[i];
            
            for (hash_map <string, vector<string>, std_string_hash>::iterator miter = pw_ko.begin(); miter != pw_ko.end(); miter ++){
                
                index = pw_index[miter->first];
                if (!No_Zero_Check[index]) continue;
                
                outfile_count << "\t" << pw_count[i][index];
                
                outfile_abd << "\t" << pw_abd[i][index];
                
            }
            
            outfile_abd << endl;
            outfile_count << endl;
        }
        
        outfile_abd.close();
        outfile_abd.clear();
        
        outfile_count.close();
        outfile_count.clear();
        
        return count;
    }

int _KO_OTU_Table_All::Load_Sample_By_OTU(hash_map <string, int, std_string_hash> * otu_seq_count, int sample){
        
        for (hash_map <string, int, std_string_hash> ::iterator miter = otu_seq_count->begin(); miter != otu_seq_count->end(); miter ++){
            
            vector<_KO_Index_Copy> kos = OTU_KO_Index[miter->first];
            int seq_count = miter->second;
            
            for (int i = 0; i < kos.size(); i ++){
                
                int index = kos[i].Index;
                float copy = kos[i].Copy;
                float rna_cp = Otu_parser.Get_cp_by_OTU(miter->first);
                KO_Abd[sample][index] += (float) seq_count * copy / rna_cp;
                }
            }
    
        return otu_seq_count->size();
        }
#endif /* class_func_h */
