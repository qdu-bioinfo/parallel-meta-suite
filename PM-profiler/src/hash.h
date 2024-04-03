//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su

#ifndef _HASH_TABLE_H
#define _HASH_TABLE_H

#include "sequence.h"

using namespace std;

class Hash_Table;

class Compare_Info
{
public:
	Compare_Info() : match_count(0) , gap(0) {}
	Compare_Info(int m, short g) : match_count(m) , gap(g) {}
private:
	friend class Hash_Table;
	int match_count;
	short gap;
	bool operator < (Compare_Info y) const {
		if(this->match_count < y.match_count) return true;
		else if(this -> match_count > y.match_count) return false;
		else if(this -> gap > y.gap) return true;
		else return false;
	}
};

class Kmer_Info
{
public:
	Kmer_Info() : seq_index(0), pos(0) {}
	Kmer_Info(int s, short p): seq_index(s), pos(p){}
	
private:
	friend class Hash_Table;
	int seq_index; //position in the vector 
    short pos;
};


class Match_Info{
public:
	Match_Info() : match_count(0), gap(0) {}
	Match_Info(int m, short g) : match_count(m), gap(g) {}
	
	void Add_res(short pos){
		
		if (match_count != 0){	
			if ((pos - curr_pos > K + 6) || (pos <= curr_pos)) return;
			gap += pos - curr_pos -1;	
		}
		curr_pos = pos;
		
		match_count ++;
	}
	
private:
	friend class Hash_Table;
	int match_count; //number of match
	short gap;  //gap size
	short curr_pos; //tmp
};

class Hash_Table{
public:
    Hash_Table() : hash_index(NULL), kmer_hash(NULL), seq_count(0) {}
    Hash_Table(vector<Sequence> & seqs, int coren);
    Hash_Table(vector<Sequence> & seqs);
    
    int Hash_match(Sequence seq, vector <int> & res, float t, int strand_mode);
    int Hash_match(vector<Sequence> & seqs, vector <int> * res, float t, int coren, int strand_mode);
    int Hash_match(vector<Sequence> & seqs, vector <int> * res, float t, int strand_mode);
    
private:
	int * hash_index;
	Kmer_Info * kmer_hash;
	int seq_count;

};

Hash_Table::Hash_Table(vector<Sequence> & seqs, int coren){
			
	//make kmer
	omp_set_num_threads(coren);
	
	#pragma omp parallel for schedule(dynamic, 10)
	for (int i = 0; i < seqs.size(); i ++)
		seqs[i].Kmer_cal(0);
	
	unsigned int kmer_size = (pow(4,K)+1);
		
	//Calc hash index spaces
	hash_index = new int[kmer_size]();
	
	for (int i = 0; i < seqs.size(); i ++)
    	for (short j = 0; j < seqs[i].kmer_count; j ++)
    		hash_index[seqs[i].kmer_list[j]] ++;
    
    int total_kmer_count = 0;
    for (int i = 0; i < kmer_size; i ++){
    	total_kmer_count += hash_index[i];
		hash_index[i] = total_kmer_count;
	}
    
    //Make hash
    kmer_hash = new Kmer_Info[total_kmer_count];
    
    for (int i = 0; i < seqs.size(); i ++)
    	for (short j = 0; j < seqs[i].kmer_count; j ++)	
    		kmer_hash[--hash_index[seqs[i].kmer_list[j]]] = Kmer_Info(i, j);
    
    //cout << seqs[kmer_hash[237204916].seq_index].seq_name <<endl; //debug
    seq_count = seqs.size();	
}

Hash_Table::Hash_Table(vector<Sequence> & seqs){
	
	Hash_Table(seqs, 1);
	}

int Hash_Table::Hash_match(Sequence seq, vector <int> & res, float t,int strand_mode){
	seq.Kmer_cal(strand_mode);
	int l = seq.seq_line.length();
	Match_Info * res_list = new Match_Info[seq_count];
	Match_Info * rev_list;
	if(strand_mode) rev_list = new Match_Info[seq_count];
	
	int thd_miss = (1-t) * l + 0.5;
	
	for (int i = 0; i < seq.kmer_count; i ++){
		
		for(int j = hash_index[seq.kmer_list[i]]; j < hash_index[seq.kmer_list[i] + 1]; j ++){
	
			if ((i - res_list[kmer_hash[j].seq_index].match_count) / K > thd_miss) continue;
			res_list[kmer_hash[j].seq_index].Add_res(kmer_hash[j].pos); 
		}
	}
	
	if(strand_mode){
		for (int i = 0; i < seq.kmer_count; i ++){
			for(int j = hash_index[seq.rev_kmer_list[i]]; j < hash_index[seq.rev_kmer_list[i] + 1]; j ++){

				if ((i - rev_list[kmer_hash[j].seq_index].match_count) / K > thd_miss) continue;
					rev_list[kmer_hash[j].seq_index].Add_res(kmer_hash[j].pos); 
				}
		}
	}

	int max_match_count = 0;
	int min_gap = 10000;
	
	for (int i = 0; i < seq_count; i ++){
		
		//cout << res_list[i].match_count << "\t" << res_list[i].gap << endl; //debug

		if (((l - res_list[i].match_count) / K) + res_list[i].gap % K > thd_miss) continue;
		
		if (res_list[i].match_count == max_match_count){
			if (res_list[i].gap < min_gap) min_gap = res_list[i].gap;
		}
		

		if (res_list[i].match_count > max_match_count){			
			max_match_count = res_list[i].match_count;
			min_gap = res_list[i].gap;
		}
		
	}
	if(strand_mode){
		for (int i = 0; i < seq_count; i ++){
		
			//cout << res_list[i].match_count << "\t" << res_list[i].gap << endl; //debug

			if (((l - rev_list[i].match_count) / K) + rev_list[i].gap % K > thd_miss) continue;

			if (rev_list[i].match_count == max_match_count){
				if (rev_list[i].gap < min_gap) min_gap = rev_list[i].gap;
			}
			

			if (rev_list[i].match_count > max_match_count){			
				max_match_count = rev_list[i].match_count;
				min_gap = rev_list[i].gap;
			}	
			
		}
	}
		
    if (!max_match_count){
        delete [] res_list;
		if(strand_mode) delete [] rev_list;
        return 0;
    }


    
	for (int i = 0; i < seq_count; i ++){
		if ((res_list[i].match_count == max_match_count) && (res_list[i].gap == min_gap)) res.push_back(i);
		if ((strand_mode) && (rev_list[i].match_count == max_match_count) && (rev_list[i].gap == min_gap)) res.push_back(i);
	}
    delete [] res_list;
	if(strand_mode) delete [] rev_list;
	return res.size();
}

int Hash_Table::Hash_match(vector <Sequence> & seqs, vector <int> * res, float t, int coren,int strand_mode){
		
	omp_set_num_threads(coren);
	
	#pragma omp parallel for schedule(dynamic, 10)
	for (int i = 0; i < seqs.size(); i ++){
		
		Hash_match(seqs[i], res[i], t , strand_mode);
	}
	
	return 0;
}

int Hash_Table::Hash_match(vector <Sequence> & seqs, vector <int> * res, float t,int strand_mode){
	return Hash_match(seqs, res, t, 1 , strand_mode);
}
#endif
