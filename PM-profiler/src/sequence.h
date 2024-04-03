//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include "utility.h"

using namespace std;

extern int K;

unsigned int char2int[256] =
{
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

unsigned int char2int_rev[256] =
{
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  3221225472,  0,  2147483648,  0,  0,  0,  1073741824,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  2147483648,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  3221225472,  0,  2147483648,  0,  0,  0,  1073741824,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  2147483648,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

unsigned int char_mask[256] =
{
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  0,  1,  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1
};

class Sequence{
public:
    
    Sequence(string name, string line);
    ~Sequence();
    string Get_seq_name();
    string Get_seq_line();
    int Kmer_cal(int strand_mode);
    int Get_kmer(int i);
    int Get_kmer_count();
    
private:
	
	friend class Hash_Table;
	
	string seq_name;
	string seq_line;
	unsigned int * kmer_list;
	unsigned int * rev_kmer_list;
	int kmer_count;
};

Sequence::Sequence(string name, string line){
	seq_name = name;
	seq_line = line;
	kmer_list = NULL;
	rev_kmer_list = NULL;
	kmer_count = 0;
}

Sequence::~Sequence(){
    if (kmer_list != NULL)
        delete [] kmer_list;
	if (rev_kmer_list != NULL)
		delete [] rev_kmer_list;
}

int Sequence::Kmer_cal(int strand_mode){
	//calculate the k-mers for a sequence
	if (K > seq_line.size()) return 0;
	kmer_count = seq_line.size() - K + 1;
	kmer_list = new unsigned int[kmer_count];
	if(strand_mode) rev_kmer_list = new unsigned int[kmer_count];
	
	char * s =  & seq_line[0];
	unsigned int kmer_num = 0;
	unsigned int rev_kmer_num = 0;
    unsigned int err = 0;
    unsigned int kmer_mask = pow(4,K)-1;
	unsigned int rev_mask = 2*(32-K);
	
	for(int i = 0; i < K-1; i++){
        err <<= 2;
        err |= char_mask[(int)(seq_line[i])];

        kmer_num <<= 2;
        kmer_num |= char2int[(int)(seq_line[i])];

		if(strand_mode){
			rev_kmer_num >>=2;
			rev_kmer_num |= char2int_rev[(int)(seq_line[i])];
		}
    }

    for(int i = K-1; i < seq_line.size(); i++){
        err <<= 2;
        err |= char_mask[(int)(seq_line[i])];
        err &= kmer_mask;

        kmer_num <<= 2;
        kmer_num |= char2int[(int)(seq_line[i])];
        kmer_num &= kmer_mask;

		if(strand_mode){
			rev_kmer_num >>= 2;
			rev_kmer_num |= char2int_rev[(int)(seq_line[i])];
			if(!err) rev_kmer_list[seq_line.size() - i - 1] = (rev_kmer_num >> rev_mask);
			else rev_kmer_list[seq_line.size() - i - 1] = (pow(4,K));
		}

        if(!err){
            kmer_list[i-K+1] = kmer_num;
        }
        else{
            kmer_list[i-K+1] = (pow(4,K));
        }
    }
    
    return kmer_count;
}
string Sequence::Get_seq_name(){
	
	return seq_name;
}

string Sequence::Get_seq_line(){
	
	return seq_line;
}

int Sequence::Get_kmer(int i){
	
	if (i >= kmer_count) return -1;
	return kmer_list[i];
}

int Sequence::Get_kmer_count(){
	
	return kmer_count;
}

/*
int Load_Seq_1(const char * infilename){
	
	ifstream infile(infilename, ios::in);
	if (!infile){
		cerr << "Error: Cannot open input file: " << infilename << endl;
		return 0;
	}
	
	int line_count = 0;
	int seq_count = 0;
	string buffer;
	
	while(getline(infile, buffer)){
		
		if((buffer[0]!='>') && (buffer[0]!='@')) {
			cerr << "Error: Input format error: " << infilename << " at line: " << line_count + 1 << endl;
			infile.close();
			infile.clear();
			return seq_count;
			}
		getline(infile, buffer);
		line_count ++;
		
		if (buffer[0]=='@'){
			getline(infile, buffer);
			getline(infile, buffer);
			line_count += 2;
		}
		
		seq_count ++;
	}
	
	
	input_seq * seqs = new input_seq[seq_count];
	
	infile.clear();
	infile.seekg(0, ios::beg);
	
	seq_count = 0;
	
	while(getline(infile, buffer)){
		seqs[seq_count].seq_name = buffer;
		getline(infile, buffer);
		seqs[seq_count].seq_line = buffer;
		if (seqs[seq_count].seq_name[0] == '@'){
			getline(infile, buffer);
			getline(infile, buffer);
		}
		
	seq_count ++;	
	}
	
	infile.close();
	infile.clear();
	
	return seq_count;
}
*/

int Load_Seq_2(const char * infilename, vector <Sequence> & seqs){
	
	ifstream infile(infilename, ios::in);
	if (!infile){
		cerr << "Error: Cannot open input file: " << infilename << endl;
		return 0;
	}
	
	int line_count = 0;
	int seq_count = 0;
	string buffer;
	string seq_name,seq_line;
	bool fastq_flag = false;
		
	while(getline(infile, buffer)){
		line_count ++;
		
		if(seq_name.empty() && (buffer[0]!='>') && (buffer[0]!='@')) {
			cerr << "Error: Input format error: " << infilename << " at line: " << line_count  << endl;
			infile.close();
			infile.clear();
			return seq_count;
			}
		if(buffer[0] =='>' || buffer[0] =='@'){	
			if(!seq_name.empty()){
			seqs.push_back(Sequence(seq_name, seq_line));
			seq_line.clear();
			seq_count ++;
			}
			seq_name = buffer.substr(1);
			fastq_flag = (buffer[0] == '@');
		}
		else if(buffer[0] =='+') {
			getline(infile, buffer);
		}
		else{
			seq_line += buffer;
		}
		
		
	}
	if(!seq_name.empty()){
		seqs.push_back(Sequence(seq_name, seq_line));
		//cout << seq_name << endl << seq_line <<endl; //cout test
		seq_line.clear();
		seq_count ++;
	}
	
	infile.close();
	infile.clear();
	
	return seq_count;
}


#endif
