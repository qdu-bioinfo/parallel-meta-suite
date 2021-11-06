// Updated at Sept 21, 2018
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// _Table_Format supported
// _OTU_Parser supported

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <stdlib.h>

#include "utility.h"
#include "hash.h"
#include "version.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"

using namespace std;

#ifndef _CLASS_TAX_H
#define _CLASS_TAX_H

#define Text_Length 120
#define Node_Width (Text_Length + Sam_num * 8 + 35)
#define X_Offset 54
#define Color_Count 10

#define TLevN 7

class TNode {
       public:
                 TNode(){
                         Count = new float [Sam_num];
                         for (int i = 0; i< Sam_num; i++)
                             Count[i] = 0;
                         }
                         
                 TNode(string _name){
                            Count = new float [Sam_num];
                            for (int i = 0; i< Sam_num; i++)
                                  Count[i] = 0;
                   
                            Name = _name;
                            }
                                                             
                 int Read_file(const char * infilename, string sam_name, int sam_n, bool is_cp_correct);
                 int Read_table(_Table_Format * table, bool is_cp_correct);
                 void Out_Tree_Html(const char * outfilename, int depth);
                 
                 static void Init(_PMDB db, int sam_num);      
                                  
       private:
               int Add_Node(string * branch, float count, string id, int sam_n, float copy_number);
               int Add_Node(string * branch, string id, int sam_n, float copy_number);
               
               string Name;
       
               map<string, TNode *> Table; //taxa to lower level

               float * Count;                                  
               
               static int Sam_num;
               static float * Sam_count;
               static string * Sam_name;               
               
               static _OTU_Parser OTU_parser;    
               
               static void Generate_Html(ostream & out);
               static void Show_Tree_Html(ostream & out, TNode * root, int depth, int i);
               //static void Get_Taxa(string otu, string * branch);                                                               
       };

int TNode::Sam_num = 0; 
float * TNode::Sam_count = NULL;
string * TNode::Sam_name = NULL;
_OTU_Parser TNode::OTU_parser = _OTU_Parser();

int TNode::Add_Node(string * branch, float count, string id, int sam_n, float copy_number){
    
    TNode * pt = this;
    
    for (int i = 0; i< 7; i++){
        
        //cin >> buffer;
        
        if (pt->Table.count(branch[i]) == 0){
            pt->Table[branch[i]] = new TNode(branch[i]);
            pt->Table[branch[i]]->Count[sam_n] = count / copy_number;
        }
        else pt->Table[branch[i]]->Count[sam_n] += (count / copy_number);
        
        pt = pt->Table[branch[i]];
        
        if ((branch[i] == "Unclassified")||(i == 6))
            break;
        
        }
    //pt = root;
    return 0;
}

int TNode::Add_Node(string * branch, string id, int sam_n, float copy_number){

    return Add_Node(branch, 1, id, sam_n, copy_number);
    
}

int TNode::Read_file(const char * infilename, string sam_name, int sam_n, bool is_cp_correct){
    
    hash_map <string, int, std_string_hash> otu_count;
    
    OTU_parser.Load_file_to_hash(infilename, otu_count);
    
    float count = 0;
    int line_count = 0;
    for (hash_map<string, int, std_string_hash> :: iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
                          
                          string branch[7];
                          //Get_Taxa(miter->first, branch);
                          OTU_parser.Get_taxa_by_OTU(miter->first, branch, TLevN);
                          
                          float copy_number = 1.0;
                          if (is_cp_correct)
                              copy_number = OTU_parser.Get_cp_by_OTU(miter->first);
                          
                          Add_Node(branch, (float) miter->second, miter->first, sam_n, copy_number);
                          line_count += miter->second;
                          count += miter->second / copy_number;
                          }
    
    Sam_count[sam_n] = count;
    Sam_name[sam_n] = sam_name;
    
    return line_count;
    }

int TNode::Read_table(_Table_Format *table, bool is_cp_correct){
    
    //file the content
    vector <string> otus = table->Get_Feature_Names();
    vector <string> sam_names = table->Get_Sample_Names();
    vector <string *> branch;
    vector <float> cp_num;
    
    for (int i = 0; i < otus.size(); i ++){
        
        string a_otu = Check_OTU(otus[i]);
        string * a_branch = new string [7];
        //Get_Taxa(a_otu, a_branch);
        OTU_parser.Get_taxa_by_OTU(a_otu, a_branch, TLevN);
        branch.push_back(a_branch);
        
        float cp = 1.0;
        if (is_cp_correct)
            cp = OTU_parser.Get_cp_by_OTU(a_otu);
        cp_num.push_back(cp);
        }
    
    for (int i = 0; i < table->Get_Sample_Size(); i ++){
        
        float count_sum = 0;
        
        //vector <float> seq_counts = table->Get_Abd(i);
        
        for (int j = 0; j < table->Get_Feature_Size(); j ++){
            Add_Node(branch[j], table->Get_Abd_By_Order(i, j), otus[j], i , cp_num[j]);
            count_sum += table->Get_Abd_By_Order(i, j) / cp_num[j];
            }
        
        Sam_count[i] = count_sum;
        Sam_name[i] = sam_names[i];
        
        //cout << seq_sum << " sequences are loaded" << endl;
        }

    return 0;
    }

void TNode::Init(_PMDB db, int sam_num){
     
     OTU_parser = _OTU_Parser(db);
     
     Sam_num = sam_num;
     Sam_count = new float [sam_num];
     for (int i = 0; i < sam_num; i ++)
         Sam_count[i] = 0;
     Sam_name = new string [sam_num];
     return;
     }

void TNode::Generate_Html(ostream & out){
    
    out << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" <<endl;
    out << "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"e\" lang=\"en\">" <<endl;
    out <<"\t"<< "<head>"<< endl;
    out <<"\t\t"<< "<meta charset=\"utf-8\"/>"<<endl;
    out <<"\t\t"<< "<style>"<<endl;
    out <<"\t\t\t"<< "body {margin:0;}"<<endl;
    out <<"\t\t"<< "</style>"<<endl;
    out <<"\t\t"<< "</head>"<<endl;
    out << "\t\t"<< "<body style=\"padding:0;position:relative\">" <<endl;
    out <<"\t\t\t"<< "<div id=\"options\" style=\"position:absolute;left:0;top:0;\">" << endl;
    out <<"\t\t\t"<< "</div>" << endl;
    out <<"\t\t\t"<< "<div id=\"details\" style=\"position:absolute;top:1%;right:2%;text-align:right;\">" << endl;
    out <<"\t\t\t"<< "</div>" << endl;
    out <<"\t\t\t"<< "<canvas id=\"canvas\" width=\"100%\" height=\"100%\">" << endl;
    out <<"\t\t\t"<< "This browser does not support HTML5." << endl;
    out <<"\t\t\t"<< "</canvas>" << endl;
    out <<"\t\t\t"<< "<img id=\"hiddenImage\" src=\"hidden.png\" visibility=\"hidden\"/>" << endl;
    out <<"\t\t\t"<< "<script name=\"meta_viewer\" src=\"meta_viewer.js\"></script>" << endl;
    out <<"\t\t"<< "</body>" << endl;
    out <<"\t\t"<< "<data>" << endl;
    out <<"\t\t\t"<< "<options collapse=\"false\" key=\"true\"></options>" << endl;
	out <<"\t\t\t"<< "<magnitude attribute=\"magnitude\"></magnitude>" << endl;
    out <<"\t\t\t"<< "<attributes magnitude=\"Total\"></attributes>" << endl;
    
    return;
    }

void TNode::Show_Tree_Html(ostream & out, TNode * root, int depth, int i){
    
    if (i >= depth) return;
           
    map<string, TNode *>::iterator map_iter = root->Table.begin();
    while (map_iter != root->Table.end()){
          
          out << "<node name=\"" <<map_iter->first << "\" magnitude=\""; 
          
          for (int j = 0; j< Sam_num; j++){
              out <<map_iter->second->Count[j] ;
              if (j != Sam_num-1) out << ",";
              }
              
          out <<"\">"<< endl;
          
          if (map_iter->first != "Unclassified") 
          
             Show_Tree_Html(out, map_iter->second, depth, i+1);
          
          out <<"</node>"<<endl;
          map_iter ++;
          
          }
    return;
    
    }


void TNode::Out_Tree_Html(const char * outfilename, int depth){
    
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                 cerr <<"Error: Open outfile error :" << outfilename << endl;
                 //system("pause");
                 exit(0);
                 }
                     
    Generate_Html(outfile);
    
    outfile <<"\t\t\t"<<"<datasets names=\"";
    for (int i = 0; i< Sam_num; i++){
        outfile <<Sam_name[i];
        if (i != Sam_num -1) outfile << ",";
        }
    outfile <<"\"></datasets>"<<endl;    
    
    outfile << "<node name=\"Root\" magnitude=\"";
    
    for (int i = 0; i< Sam_num; i++){
        outfile << Sam_count[i];
        if (i != Sam_num-1) outfile << ",";
        }
        
    outfile<<"\">"<< endl;
    
    Show_Tree_Html(outfile, this, depth, 0);
    
    outfile<<"</node>"<<endl;
    outfile <<"\t\t" << "</data>" << endl;
    outfile <<"</html>" << endl;
    
    outfile.close();
    outfile.clear();
    
    return;
    }

/*
void TNode::Get_Taxa(string otu, string * branch){
    
    string otu_taxa;
    int level = 7;
    
    otu_taxa = OTU_parser.Get_taxa_by_OTU(otu);
    
    vector <string> taxa_buffer;
    stringstream strin(otu_taxa);
    string taxa;
    string temp;
    bool is_taxa = true;

    for (int i = 0; i < level; i++){
        
        
        if (!(strin >> taxa)){
            taxa =  "Unclassified";
            is_taxa = false;
            }
        
        else if (taxa.find("Unclassified") != string::npos){
            taxa =  "Unclassified";
            is_taxa = false;
        }
        else if (taxa.find("otu_") != string::npos){
            taxa = "Unclassified";
            is_taxa = false;
        }
        
        else if(taxa[taxa.size()-1] != ';'){
            strin >> temp;
            taxa += "_";
            taxa += temp;
        }
        
        taxa_buffer.push_back(taxa);
        
        //add genus info for species
        if ((i == 6) && (is_taxa)) {
            taxa = taxa_buffer[5];
            taxa = taxa.substr(0, taxa.size()-1);
            taxa += "_";
            taxa += taxa_buffer[6];
        }
        
        //add prefix for genus
        if ((i == 5) && (is_taxa)) {
            string taxa_prefix = taxa_buffer[i-1];
            taxa_prefix = taxa_prefix.substr(0, taxa_prefix.size() - 1);
            if (taxa_prefix.size() > 3) taxa_prefix = taxa_prefix.substr(0, 3);
            taxa = taxa_prefix + "_" + taxa;
        }
        branch[i] = taxa;
    }
    return;
}
*/
#endif


