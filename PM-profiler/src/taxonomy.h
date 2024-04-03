//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su

#ifndef _TAXONOMY_H
#define _TAXONOMY_H

#include <unordered_map>
//#include <set>
#include <sstream>
#include "sequence.h"
#include "utility.h"

using namespace std;

#define TAXA_LEVEL 7 //k p c o f g s

vector <string> def_taxa = {"k__;", "p__;", "c__;", "o__;", "f__;", "g__;", "s__;"};
string unmapped_taxa = "unmapped";

class Taxonomy{

public:
	Taxonomy() {}
    Taxonomy(const char * infilename, vector <Sequence> & seqs);
	
	int Parse_taxa_LCA(vector<int> res, string & res_taxa);
	int Parse_taxa_HWL(vector<int> res, string & res_taxa , float consistency);
    
    int Parse_taxa_LCA(vector<int> * res, vector <Sequence> & seqs, const char * outfilename, int coren);
    int Parse_taxa_HWL(vector<int> * res, vector <Sequence> & seqs, const char * outfilename, int coren ,float consistency);
    
    int Parse_taxa_raw(vector <int> * res, vector<Sequence> & seqs, const char * outfilename);
	
private:
	
	vector<string> taxa_id;
	vector<string> * taxa_info; //taxonomy
};

Taxonomy::Taxonomy(const char * infilename, vector <Sequence> & seqs){
	
	ifstream infile(infilename, ios::in);
	if (!infile){
		cerr << "Error: Cannot open input file: " << infilename << endl;
		return;
	}
	
	unordered_map<string,vector<string> > tax_index;
	
    int line_count = 0;
	string buffer;    
    while(getline(infile, buffer)){
    	
        stringstream strin(buffer);
        
        string tax_name;
		string tax_tmp;
        strin >> tax_name;
        vector<string> tax_anno;
        
        for (int i = 0; i < TAXA_LEVEL; i ++){
        	strin >> tax_tmp;
            tax_anno.push_back(move(tax_tmp));
        }
        tax_index[tax_name] = move(tax_anno);
        line_count ++;
    }
	infile.close();
	infile.clear();
    
	taxa_info = new vector<string> [seqs.size()];
    
	for (int i = 0; i < seqs.size(); i ++){
        taxa_id.push_back(seqs[i].Get_seq_name());
        //taxa_info[i] = move(tax_index[taxa_id[i]]);
        
        if (tax_index.count(taxa_id[i]) == 0) cerr << "Error: Cannot find taxonomy for: " << taxa_id[i] << endl;
		else taxa_info[i] = move(tax_index[taxa_id[i]]);
        
    }
}

int Taxonomy::Parse_taxa_LCA(vector<int> res, string & res_taxa){
	
    if (res.size() == 0){ //no mapped res
        res_taxa = unmapped_taxa;
        return -1;
        }
	
	vector <string> res_taxa_list;
	
	for (int i = 0; i < TAXA_LEVEL; i++){
		
		string curr_taxa;
		bool lca_flag = true;
		
		for (int j = 0; j < res.size(); j ++){
            if (taxa_info[res[j]].size() <= i) continue;
			if (taxa_info[res[j]][i].size() <= 4) continue;
			if (curr_taxa.size() == 0) //init
				curr_taxa = taxa_info[res[j]][i]; 

			else if (curr_taxa != taxa_info[res[j]][i]){ //unmatched
				lca_flag = false;
				break;
			}
		}
		if (curr_taxa.size() ==0) lca_flag = false;
	
		if(lca_flag) 
			res_taxa_list.push_back(curr_taxa);
		else break;
	}
	
    if (res_taxa_list.size() == 0){
        res_taxa = def_taxa[0];
        return 0;
        }
	
	res_taxa = "";
	for (int i = 0; i < res_taxa_list.size(); i ++)//combine
		res_taxa += res_taxa_list[i] + "\t";
    
	return res_taxa_list.size();
}

/* old hwl_method
int Taxonomy::Parse_taxa_Retax(vector<int> res, string & res_taxa){
	
    if (res.size() == 0){ //no mapped res
        res_taxa = unmapped_taxa;
        return -1;
        }
	
	vector <int> curr_vote = res;
	int rep = curr_vote[0]; //default
	
	for (int i = 0 ; i < TAXA_LEVEL; i ++){
		
		unordered_map <string, vector <int> > taxa_vote;
		
        for (int j = 0; j < curr_vote.size(); j ++){
            if (taxa_info[curr_vote[j]].size() <= i) continue;
			if (taxa_info[curr_vote[j]][i].size() <= 4) continue;
			taxa_vote[taxa_info[curr_vote[j]][i]].push_back(curr_vote[j]);
		}
		
		int max_count = 0;
		string max_taxa;
		
		for (unordered_map <string, vector <int> >::iterator miter = taxa_vote.begin(); miter != taxa_vote.end(); miter ++)
				if ((miter->second).size() > max_count){
					
					max_count = (miter->second).size();
					max_taxa = miter->first;
				}
		
        //cout << max_taxa << endl; //debug
        
        if (max_count == 0) break;
        else{
            curr_vote = taxa_vote[max_taxa];
            rep = curr_vote[0];
        }
	}
    
	
    res_taxa = "";
	for (int i = 0; i < TAXA_LEVEL; i ++){//combine
        if(taxa_info[rep][i].size() > 4)
		    res_taxa += taxa_info[rep][i] + "\t";
		}
		
	return rep;
}
*/

int Taxonomy::Parse_taxa_HWL(vector<int> res, string & res_taxa , float consistency){
	
    if (res.size() == 0){ //no mapped res
        res_taxa = unmapped_taxa;
        return -1;
        }
	
	int rep = res[0]; //default
    int rep_count = 0,con_level = -1;
	
	for (int i = TAXA_LEVEL - 1 ; i >= 0; i --){
		
		unordered_map <string, vector <int> > taxa_vote;
		
        for (int j = 0; j < res.size(); j ++){
            if (taxa_info[res[j]].size() <= i) continue;
			if (taxa_info[res[j]][i].size() <= 4) continue;
			taxa_vote[taxa_info[res[j]][i]].push_back(res[j]);
		}
		
		int max_count = 0;
		string max_taxa;
		
		for (unordered_map <string, vector <int> >::iterator miter = taxa_vote.begin(); miter != taxa_vote.end(); miter ++)
				if ((miter->second).size() > max_count){
					
					max_count = (miter->second).size();
					max_taxa = miter->first;
				}
		
        //cout << max_taxa << endl; //debug
        
        if (max_count == 0 || (float)(max_count)/res.size() < consistency) continue;
        else{
            rep = taxa_vote[max_taxa][0];
            rep_count = max_count;
            con_level = i;
            break;
        }
	}
    
	
    res_taxa = "";
	for (int i = 0; i <= con_level; i ++){//combine
        if(taxa_info[rep][i].size() > 4)
		    res_taxa += taxa_info[rep][i] + " ";
		}
	res_taxa = res_taxa.substr(0,res_taxa.length()-1) + "\t" + to_string((float)(rep_count)/res.size());	
	return rep;
}

int Taxonomy::Parse_taxa_LCA(vector<int> * res, vector <Sequence> & seqs, const char * outfilename, int coren){
    
    string * res_taxa = new string [seqs.size()];
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 10)
    for(int i = 0; i < seqs.size(); i ++)
        Parse_taxa_LCA(res[i], res_taxa[i]);
    
    ofstream outfile(outfilename, ios::out);
    if (!outfile){
        cerr << "Error: Cannot open output file: " << outfilename << endl;
        delete [] res_taxa;
        return 0;
    }
    
    for (int i = 0; i < seqs.size(); i ++)
        outfile << seqs[i].Get_seq_name() << "\t" << res_taxa[i] << endl;
    
    outfile.close();
    outfile.clear();
    
    delete [] res_taxa;
    return seqs.size();
}

int Taxonomy::Parse_taxa_HWL(vector<int> * res, vector <Sequence> & seqs, const char * outfilename, int coren , float consistency){
    
    string * res_taxa = new string [seqs.size()];
    int * res_rep = new int [seqs.size()];
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 10)
    for(int i = 0; i < seqs.size(); i ++)
        res_rep[i] = Parse_taxa_HWL(res[i], res_taxa[i],consistency);
    
    ofstream outfile(outfilename, ios::out);
    if (!outfile){
        cerr << "Error: Cannot open output file: " << outfilename << endl;
        delete [] res_taxa;
        delete [] res_rep;
        return 0;
    }
    
    for (int i = 0; i < seqs.size(); i ++){
        outfile << seqs[i].Get_seq_name() << "\t";
        if (res_rep[i] >=0) outfile << taxa_id[res_rep[i]] << "\t" << res_taxa[i] << endl;    
        else outfile << unmapped_taxa << endl;
    }
    outfile.close();
    outfile.clear();
    
    delete [] res_taxa;
    delete [] res_rep;
    
    return seqs.size();
}

int Taxonomy::Parse_taxa_raw(vector <int> * res, vector<Sequence> & seqs, const char * outfilename){
    
    ofstream outfile(outfilename, ios::out);
    if (!outfile){
        
        cerr << "Error: Cannot open output file: " << outfilename << endl;
        return 0;
    }
    
    
    for(int i = 0; i < seqs.size(); i ++){
        
        outfile << seqs[i].Get_seq_name() << "\t" << res[i].size() << endl;
        unordered_map<string,int> taxa_list;

        for (int j = 0; j < res[i].size(); j++){
            string tmp_taxa = "";
            for(int k = 0;k<TAXA_LEVEL;k++){
                tmp_taxa += taxa_info[res[i][j]][k] + "\t";
            }
            taxa_list[tmp_taxa]++;
        }
        for (unordered_map <string, int >::iterator miter = taxa_list.begin(); miter != taxa_list.end(); miter ++){
            outfile << miter -> first << (float)(miter -> second)/res[i].size() << endl;
        }
    }
    
    /*
    for(int i = 0; i < sample_seqs.size(); i ++){
        
        outfile << seqs[i].Get_seq_name() << "\t";
        for (int j = 0; j < res[i].size(); j++)
            outfile << res[i][j] << "\t";
        outfile << endl;
    }
    */
    outfile.close();
    outfile.clear();
    
    return seqs.size();
}
#endif
