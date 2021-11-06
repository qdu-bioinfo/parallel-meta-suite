// Updated at Sept 8, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/mman.h>

#include <unistd.h>
#include <sys/types.h>

#include "hash.h"

#define MAX 100000000

using namespace std;

string Get_hash_label(string buffer, int p)
{

    stringstream strin(buffer);
    string res;
    for (int i = 0; i < p; i++)
    {
        strin >> res;
    }
    return res;
}

unsigned int Read_file_hash(const char *infilename, vector<string> &labels, vector<string> &seqs, hash_map<string, unsigned int, std_string_hash> &ID_seqs)
{
    ifstream infile(infilename, ifstream::in);
    if (!infile)
    {
        cerr << "Error: Open infile error : " << infilename << endl;
        exit(0);
    }
    unsigned int total_seq_count = 0;

    string buffer;
    string seq_temp = "";

    while (getline(infile, buffer))
    {
        //if (buffer[buffer.size()-1] == '\r') buffer = buffer.substr(0, buffer.size()-1);
        if (buffer[0] == '>')
        {

            labels.push_back(buffer);

            //(*ID_seqs)[Get_hash_label(buffer, 0)] = total_seq_count;
            ID_seqs[Get_hash_label(buffer, 1)] = total_seq_count;

            if (total_seq_count > 0)
            {

                seqs.push_back(seq_temp);
            }

            seq_temp = "";
            total_seq_count++;
        }
        else
        {

            seq_temp += buffer;
        };
    }

    seqs.push_back(seq_temp);
    infile.close();
    infile.clear();
    return total_seq_count;
}

unsigned int Hmm_parsing(const char *infilename, const char *outfilename, vector<string> *seqs, hash_map<string, unsigned int, std_string_hash> *ID_seqs, string label_pref, int length_filter)
{

    ifstream infile(infilename, ifstream::in);
    if (!infile)
    {
        cerr << "Error: Open Hmm file error : " << infilename << endl;
        exit(0);
    }

    ofstream outfile(outfilename, ofstream::out);
    if (!outfile)
    {
        cerr << "Error: Open Hmm outfile error :" << outfilename << endl;
        exit(0);
    }

    string buffer;
    unsigned int total_seq_count = 0;

    while (getline(infile, buffer))
    {

        if (buffer[0] == '>')
        {

            //get the hash_label

            //string hash_label = Get_hash_label(buffer, 3);
            string hash_label = Get_hash_label(buffer, 2);
            //cout << hash_label<< "\t"; //debug

            unsigned int ID = (*ID_seqs)[">" + hash_label];
            //cout << ID <<endl;

            getline(infile, buffer);
            getline(infile, buffer);
            getline(infile, buffer);

            while (buffer[5] == '!')
            {
                string skip;
                int begin;
                int end;

                istringstream strin(buffer);
                for (int i = 0; i < 9; i++)
                    strin >> skip;

                strin >> begin;
                strin >> end;

                if (end - begin > length_filter)
                { //Length_filter

                    //printf out the label
                    if (buffer[3] == '1')
                        outfile << ">" << total_seq_count + 1 << "|" << label_pref << hash_label << endl;
                    else
                        outfile << ">dup" << buffer[3] << label_pref << hash_label << endl;
                    //printf out the aligned seqs

                    string aligned_seq = ((*seqs)[ID]).substr(begin - 1, end - begin);
                    //cout << begin-1 << "\t" << end-begin << endl; //debug
                    outfile << aligned_seq << endl;

                    total_seq_count++;
                }
                buffer[5] = ' ';
                getline(infile, buffer);
            }

            //total_seq_count ++;
        }
    }

    infile.close();
    infile.clear();
    outfile.close();
    outfile.clear();
    return total_seq_count;
}

int Write_file(const char *outfilename, unsigned int seq_count, vector<string> *labels, vector<string> *seqs)
{

    ofstream outfile(outfilename, ofstream::out);
    if (!outfile)
    {
        cerr << "Error: Open outfile error :" << outfilename << endl;
        exit(0);
    }

    for (unsigned int i = 0; i < seq_count; i++)
    {
        outfile << (*labels)[i] << endl;
        outfile << (*seqs)[i] << endl;
    }

    outfile.close();
    outfile.clear();
    return 0;
}

int Init_comp_table(char *complement_table_A, char *complement_table_a)
{
    for (int i = 0; i < 26; i++)
    {
        complement_table_A[i] = 'N';
        complement_table_a[i] = 'n';
    }

    complement_table_A['A' - 'A'] = 'T';
    complement_table_A['C' - 'A'] = 'G';
    complement_table_A['G' - 'A'] = 'C';
    complement_table_A['T' - 'A'] = 'A';

    complement_table_a['a' - 'a'] = 't';
    complement_table_a['c' - 'a'] = 'g';
    complement_table_a['g' - 'a'] = 'c';
    complement_table_a['t' - 'a'] = 'a';

    return 0;
}

int Seq_reverse(unsigned int seq_count, vector<string> &seqs)
{

    char complement_table_A[26];
    char complement_table_a[26];

    Init_comp_table(complement_table_A, complement_table_a);

    string temp;

    for (unsigned int i = 0; i < seq_count; i++)
    {

        //For label
        //sprintf(temp, ">rev|%s", labels[i]+1);
        //strcpy(labels[i], temp);
        //
        int length = seqs[i].size();

        for (int j = 0; j < length / 2; j++)
        {
            char tempchar1 = seqs[i][j];
            char tempchar2 = seqs[i][length - j - 1];

            if ((tempchar1 >= 'A') && (tempchar1 <= 'Z'))
                seqs[i][length - j - 1] = complement_table_A[tempchar1 - 'A'];

            else if ((tempchar1 >= 'a') && (tempchar1 <= 'z'))
                seqs[i][length - j - 1] = complement_table_a[tempchar1 - 'a'];

            else
                seqs[i][length - j - 1] = 'X';

            if ((tempchar2 >= 'A') && (tempchar2 <= 'Z'))
                seqs[i][j] = complement_table_A[tempchar2 - 'A'];

            else if ((tempchar2 >= 'a') && (tempchar2 <= 'z'))
                seqs[i][j] = complement_table_a[tempchar2 - 'a'];

            else
                seqs[i][j] = 'X';
        }
        if (length / 2 * 2 < length)
        {
            if ((seqs[i][length / 2] >= 'A') && (seqs[i][length / 2] <= 'Z'))
                seqs[i][length / 2] = complement_table_A[seqs[i][length / 2] - 'A'];
            else if ((seqs[i][length / 2] >= 'a') && (seqs[i][length / 2] <= 'z'))
                seqs[i][length / 2] = complement_table_a[seqs[i][length / 2] - 'a'];
            else
                seqs[i][length / 2] = 'X';
        }
    }
    return 0;
}

unsigned int ExtractRNA(int domain, string infilename, string outpath, int length_filter, string this_path, int coren)
{

    pid_t cld_pid;
    int *pmap = (int *)mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);

    if ((cld_pid = fork()) == 0)
    {

        vector<string> labels;
        vector<string> seqs;

        hash_map<string, unsigned int, std_string_hash> ID_seqs;

        unsigned int rna_count1;
        unsigned int rna_count2;

        char command[3000];

        string exefilename = "hmmsearch";

        string modelfilename;
        if (domain == 0)
            modelfilename = this_path + "/models/bac_ssu.hmm";
        else
            modelfilename = this_path + "/models/euk_ssu.hmm";

        string hmmout1 = outpath + "/hmm.out1";
        string hmmout2 = outpath + "/hmm.out2";

        string fowardrna = outpath + "/tmp/foward.rna";
        string reverserna = outpath + "/tmp/reverse.rna";

        string reversefa = outpath + "/tmp/reverse.fa";

        string label_pref1 = "";
        string label_pref2 = "rev|";

        cout << "rRNA Extraction Starts" << endl;

        sprintf(command, "%s -E 1e-05 --cpu %d %s %s > %s", exefilename.c_str(), coren, modelfilename.c_str(), infilename.c_str(), hmmout1.c_str());
        system(command);

        //Extract Foward
        unsigned int seq_count = Read_file_hash(infilename.c_str(), labels, seqs, ID_seqs);
        cout << "There are " << seq_count << " sequences in total" << endl
             << endl;

        //cout <<"Start Parsing\n";
        rna_count1 = Hmm_parsing(hmmout1.c_str(), fowardrna.c_str(), &seqs, &ID_seqs, label_pref1, length_filter);
        //cout <<"Parsing Finished\n";

        //Extract Reverse
        //cout << "Reverse starts\n";
        Seq_reverse(seq_count, seqs);
        //cout <<"Reverse finished\n";

        Write_file(reversefa.c_str(), seq_count, &labels, &seqs);

        sprintf(command, "%s -E 1e-05 --cpu %d %s %s > %s", exefilename.c_str(), coren, modelfilename.c_str(), reversefa.c_str(), hmmout2.c_str());
        system(command);

        //cout <<"Start Parsing\n";
        rna_count2 = Hmm_parsing(hmmout2.c_str(), reverserna.c_str(), &seqs, &ID_seqs, label_pref2, length_filter);
        //cout <<"Parsing Finished\n";
        //Cat
        sprintf(command, "cat %s/tmp/*rna >> %s/meta.rna", outpath.c_str(), outpath.c_str());
        system(command);

        unsigned int rna_count = rna_count1 + rna_count2;

        cout << "rRNA Extraction Finished" << endl;

        cout << "There are " << rna_count << " rRNA sequences in total\n"
             << endl;
        ;

        sprintf(command, "rm %s/hmm.out*", outpath.c_str());
        system(command);

        ID_seqs.clear();

        pmap[0] = seq_count;

        exit(0);

        //system("pause");
    }
    while (wait(NULL) != -1)
        ;

    unsigned int seq_count = pmap[0];

    munmap(pmap, sizeof(int));

    return seq_count;
}
