//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su
#include "opt.h"
#include "hash.h"
#include "taxonomy.h"

int main(int argc, char * argv[]){
	//parse opts
	Parse_Para(argc, argv);
    
    //switch run_mode
    vector <string> sam_name;
    vector <string> infile_list;
    
    switch (Run_mode){
        case 0:
            sam_name.push_back("sample");
            infile_list.push_back(Sampleseqfilename);
            break;
        case 1:
            Load_List(Listfilename.c_str(), infile_list, sam_name, Listprefix);
            Check_Path(Outfilename.c_str(), 1);
            break;
        default:
            cerr << "Error: Please specify input sample file by -i or -l" << endl;
            exit(1);
        }
    
    cout << "PM profiler starts" << endl;
    
	cout << "Thread: " << Coren << endl;
	
	vector <Sequence> db_seqs;
	//load db seqs
	cout << Load_Seq_2(Dbseqfilename.c_str(), db_seqs) << " database sequence(s) loaded" << endl;
	
    //make kmer db
	cout << "Database construction starts" << endl;
	Hash_Table hash_table(db_seqs, Coren);
	cout << "Database construction finished" << endl;
	
    //load taxa
    Taxonomy taxa(Dbtaxafilename.c_str(), db_seqs);
    cout << "Database taxonomy loaded" << endl;
    
    for (int i = 0; i < infile_list.size(); i ++){
        //load sample seqs
        vector <Sequence> sample_seqs;
        cout << Load_Seq_2(infile_list[i].c_str(), sample_seqs) << " sequence(s) loaded" << endl;
        //cout <<"Input sample: " << infile_list[i] << endl;
        
        //search kmer db
        vector <int> * res = new vector <int> [sample_seqs.size()];
        cout << "Database search starts" << endl;
        hash_table.Hash_match(sample_seqs, res, Sim, Coren, Strand_mode);
        cout << "Database search finished" << endl;
        
        //make output filename
        string outfile;
        if (Run_mode == 0) outfile = Outfilename;
        else if (Run_mode == 1){
            outfile = Outfilename + "/" + sam_name[i];
            Check_Path(outfile.c_str(), 1);
            outfile += "/pm.out";
        }
        
        //parse taxa
        
        if (Raw)
            taxa.Parse_taxa_raw(res, sample_seqs, (outfile + ".raw.txt").c_str());
            //cout << "Output raw: " << outfile + ".raw.txt" << endl; //debug
        if (Taxa_mode == 1)
            taxa.Parse_taxa_LCA(res, sample_seqs, (outfile + ".lca.txt").c_str(), Coren);
            //cout << "Output LCA: " << outfile + ".lca.txt" << endl; //debug
        taxa.Parse_taxa_HWL(res, sample_seqs, (outfile + ".hwl.txt").c_str(), Coren , Consistency);
        //cout << "Output HWL: " << outfile + ".hwl.txt" << endl; //debug
         cout << "Taxonomy parse finished" << endl;
        
        delete [] res;
    }
    
    //debug taxonomy
    /*
    string a_test_taxa;
    int rep = 0;
    vector <int> res = {2385, 2386, 2387};
    
    taxa.Parse_taxa_LCA(res, a_test_taxa);
    cout << a_test_taxa << endl;
    rep = taxa.Parse_taxa_Retax(res, a_test_taxa);
    cout << rep << "\t" << a_test_taxa << endl;
	*/
	return 0;
	
}
