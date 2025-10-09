// Updated at July 6, 2025
// Updated by Yu Zhang and Jiaming Yu
// version 3.7.3 - 3.7.4
// update the issue of creating the database multiple times when calling PM-profiler,PM-parallel-meta support list input

// Updated at Aug 16, 2024
// Updated by Haobo Shi,Xiaoquan Su
// version 3.7 - 3.7.2
// Update profiler,remove pair-end mode

// Bioinformatics Group, College of Computer Science & Technology, Qingdao University

#include <iostream>
#include <tuple>

#include "init.h"
#include "fastq.h"
#include "ExtractRNA.h"
#include "multialign.h"
#include "taxonomy.h"

using namespace std;

struct SampleInfo {
	string name;
	string path;
};

void Single_Run(_Para para) {

	string command;

	int seq_count = 0;
	int rna_count = 0;
	int asv_count = 0;
	int match_rna_count = 0;
	int drop_rna_count = 0;

	int a_diver [LEVEL] = {0, 0, 0, 0, 0, 0, 0};

	string mergefile;

	string handlefile;

	//check format, infilename
	para.Format = Check_Format(para.Infilename.c_str());

	if(para.Format < 0) return;//Format error

	if(para.Is_format_check) { //Check format
		cout << "Format Check Starts" << endl;
		command = "PM-format-seq -i " + para.Infilename;
		system(command.c_str());
		cout << endl;
	}

	//check format for pair ends, infilename2
	/* //delete
	if (para.Is_paired){
	    para.Format = Check_Format(para.Infilename2.c_str());
	    if (para.Format < 0) return;//Format error

	    if (para.Is_format_check){//Check format
	        cout << "Format Check 2 Starts" << endl;
	        command = "PM-format-seq -i " + para.Infilename2;
	        system(command.c_str());
	        cout << endl;
	    }
	}
	*/
	//Type, 0:16S  1:shotgun
	//Format, 0:fasta  1:fastq

	if (para.Format == 1) { //If fastq
		string tempfilename = para.Out_path + "/meta.fasta";
		cout << "Pre-computation for Fastq Starts" << endl;
		cout << endl << Fastq_2_Fasta(para.Infilename.c_str(), tempfilename.c_str()) << " sequences have been pre-computed" << endl << endl;
		para.Infilename = tempfilename;
	}

	if(para.Type == 1) { //if meta
		//Extract 16S r RNA
		seq_count = ExtractRNA(para.Database.Get_Domain(), para.Infilename, para.Out_path, para.Length_filter, para.This_path, para.Core_number);
		handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Out_path + "/meta.rna", para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, 0, para.Core_number);
		//search database
		Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt", para.Database.Get_Path()+ "/database.fa", para.db_similarity,'F', para.Core_number, para.profiler);
	}
	/* //delete by Shi Haobo
	else if(para.Is_paired){//if paired
		if(para.Format == 0){//fasta cant merge in vsearch
			cerr << "Error: For pair ends you need to input fastq format file" << endl;
		 	return ;
		}else{//fastq
			//merge
			mergefile = Merge_Pairend(para.Align_exe_name.c_str(), para.Infilename, para.Infilename2, para.Out_path + "/tmp", para.Core_number);
			//dereplication, denoise, nonchimeras
			handlefile = Handle_seq(para.Align_exe_name.c_str(), mergefile , para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count,0, para.Core_number);
			//search database
			Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt" , para.Database.Get_Path()+ "/database.fa" , para.db_similarity,'F', para.Core_number, para.profiler);
			if (rna_count < 0) {
	                   cerr << "Error: 2 ends contain different number of sequences" << endl;
	                   return;
	        }
		}

	}
	*/
	else { //single
		//dereplication, denoise, nonchimeras
		handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Infilename, para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, 0, para.Core_number);

		//search database
		Search_db(para.Align_exe_name.c_str(), handlefile, para.Out_path + "/tmp", para.Database.Get_Path()+ "/taxonomy_annotation.txt", para.Database.Get_Path()+ "/database.fa", para.db_similarity,'F', para.Core_number, para.profiler);
	}
	if(para.profiler == 1)
		Otutab_Count(para.Out_path + "/tmp/PM.hwl.txt",para.Out_path + "/tmp/map_output.txt");
	//parse_taxonomy
	Out_Taxonomy((para.Out_path + "/tmp/map_output.txt").c_str(), para.Out_path, para.Database, 0, a_diver, para.Is_paired, match_rna_count, drop_rna_count);

	//make_plot
	command = "PM-plot-taxa -D ";
	command += para.Database.Get_Id();
	command += " -i ";
	command += para.Out_path + "/classification.txt -o " + para.Out_path;
	system(command.c_str());

	//print report
	Print_Report(para, seq_count, rna_count, match_rna_count, a_diver,asv_count);

	//func_anno
	if (para.Is_func) {
		//func
		command = "PM-predict-func -D ";
		command += para.Database.Get_Id();
		command += " -i ";
		command += para.Out_path + "/classification.txt -o " + para.Out_path;
		system(command.c_str());

		//nsti
		command = "PM-predict-func-nsti -D ";
		command += para.Database.Get_Id();
		command += " -i ";
		command += para.Out_path + "/classification.txt >> " + para.Out_path + "/Analysis_Report.txt";
		system(command.c_str());
	}
	//Remove the tmp
	command = "rm -rf " + para.Out_path + "/tmp";
	system(command.c_str());
}

void Preprocess(_Para para) {
	string command;

	int seq_count = 0;
	int rna_count = 0;
	int asv_count = 0;
	int match_rna_count = 0;
	int drop_rna_count = 0;

	int a_diver [LEVEL] = {0, 0, 0, 0, 0, 0, 0};

	string mergefile;

	string handlefile;

	//check format, infilename
	para.Format = Check_Format(para.Infilename.c_str());

	if(para.Format < 0) return;//Format error

	if(para.Is_format_check) { //Check format
		cout << "Format Check Starts" << endl;
		command = "PM-format-seq -i " + para.Infilename;
		system(command.c_str());
		cout << endl;
	}

	if (para.Format == 1) { //If fastq
		string tempfilename = para.Out_path + "/meta.fasta";
		cout << "Pre-computation for Fastq Starts" << endl;
		cout << endl << Fastq_2_Fasta(para.Infilename.c_str(), tempfilename.c_str()) << " sequences have been pre-computed" << endl << endl;
		para.Infilename = tempfilename;
	}

	if(para.Type == 1) { //if meta
		//Extract 16S r RNA
		seq_count = ExtractRNA(para.Database.Get_Domain(), para.Infilename, para.Out_path, para.Length_filter, para.This_path, para.Core_number);
		handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Out_path + "/meta.rna", para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, 0, para.Core_number);
	} else { //single
		//dereplication, denoise, nonchimeras
		handlefile = Handle_seq(para.Align_exe_name.c_str(), para.Infilename, para.Out_path + "/tmp", para.Is_denoised, para.Is_nonchimeras, rna_count, asv_count, 0, para.Core_number);
	}

	//print report
	Print_Report(para, seq_count, rna_count, match_rna_count, a_diver,asv_count);

	//write samplename and handlefile to list
	ofstream out_pm_sample((para.Singlesample_dir + "/pm_sample.list").c_str(), ofstream::app);
	out_pm_sample << para.Sample_name << '\t' << handlefile << std::endl;
	// close fstream
	out_pm_sample.close();

	// write the required parameters for the report to list
	ofstream report_list_out((para.Singlesample_dir + "/pm_report.list").c_str(), ofstream::app);
	report_list_out << para.Sample_name << '\t' << seq_count << '\t' << rna_count << '\t' << asv_count << endl;
	report_list_out.close();

}

void Run_Profiler(_Para para) {
	string command;
	// run PM-profiler
	//	sprintf(command,"PM-profiler -s %.2f -d %s -m %s -l %s -o %s -t %d -Q T", para.db_similarity, (para.Database.Get_Path()+ "/database.fa").c_str(), (para.Database.Get_Path()+ "/taxonomy_annotation.txt").c_str(),  (para.Singlesample_dir + "/pm_sample.list").c_str(), (para.Singlesample_dir).c_str(),para.Core_number);
	//	system(command.c_str());
	command = "PM-profiler";
	command += " -s " + std::to_string(para.db_similarity);
	command += " -d " + para.Database.Get_Path() + "/database.fa";
	command += " -m " + para.Database.Get_Path() + "/taxonomy_annotation.txt";
	command += " -l " + para.Singlesample_dir + "/pm_sample.list";
	command += " -o " + para.Singlesample_dir;
	command += " -t " + std::to_string(para.Core_number);
	command += " -Q T";

	system(command.c_str());
	//Run_With_Error(command, "PM-profiler", tmpError_file.c_str());
	// delete pm_sample.list
	string pm_sample_command;
	pm_sample_command = "rm ";
	pm_sample_command += para.Singlesample_dir + "/pm_sample.list";
	system(pm_sample_command.c_str());
}

void Parse_Plot_Func(_Para para,int seq_count, int rna_count, int asv_count) {
	string command;
	int a_diver [LEVEL] = {0, 0, 0, 0, 0, 0, 0};
	//parse_taxonomy
	int match_rna_count = 0;
	int drop_rna_count = 0;

	Otutab_Count(para.Out_path + "/tmp/PM.hwl.txt",para.Out_path + "/tmp/map_output.txt");
	Out_Taxonomy((para.Out_path + "/tmp/map_output.txt").c_str(), para.Out_path, para.Database, 0, a_diver, para.Is_paired, match_rna_count, drop_rna_count);

	//make_plot
	command = "PM-plot-taxa -D ";
	command += para.Database.Get_Id();
	command += " -i ";
	command += para.Out_path + "/classification.txt -o " + para.Out_path;
	system(command.c_str());

	//print report
	Print_Report(para, seq_count, rna_count, match_rna_count, a_diver,asv_count);

	//func_anno
	if (para.Is_func) {
		//func
		command = "PM-predict-func -D ";
		command += para.Database.Get_Id();
		command += " -i ";
		command += para.Out_path + "/classification.txt -o " + para.Out_path;
		system(command.c_str());

		//nsti
		command = "PM-predict-func-nsti -D ";
		command += para.Database.Get_Id();
		command += " -i ";
		command += para.Out_path + "/classification.txt >> " + para.Out_path + "/Analysis_Report.txt";
		system(command.c_str());
	}

	//Remove the tmp
	command = "rm -rf " + para.Out_path + "/tmp";
	system(command.c_str());
}

int main(int argc, char * argv[]) {

	cout << endl << "Welcome to Parallel-Meta Suite version " << Version << endl << endl;

	_Para para;
	int seq_count, rna_count, asv_count;

	Parse_Para(argc, argv, para);

	if (para.run_mode == 0) { // single
		cout << "Running in single-sample mode..." << endl;
		Single_Run(para);
	} else if (para.run_mode == 1) { // multiple
		cout << "Running in multiple-sample mode..." << endl;
		ifstream infile(para.Listfilename.c_str());
		if (!infile.is_open()) {
			cerr << "Error: Cannot open list file: " << para.Listfilename << endl;
			return 1;
		}

		vector<SampleInfo> samples;
		string line;

		while (getline(infile, line)) {
			if (line.empty()) continue;

			size_t tab_pos = line.find('\t');
			if (tab_pos == string::npos) {
				cerr << "Error: Invalid line format (no tab found): " << line << endl;
				continue;
			}

			SampleInfo sample;
			sample.name = line.substr(0, tab_pos);
			sample.path = line.substr(tab_pos + 1);

			sample.path.erase(0, sample.path.find_first_not_of(" \r\n"));
			sample.path.erase(sample.path.find_last_not_of(" \r\n") + 1);

			samples.push_back(sample);
		}

		infile.close();

		if(para.profiler == 1) { // pm-profiler
			// []

			// Preprocess
			for (size_t i = 0; i < samples.size(); ++i) {
				cout << "\nProcessing sample " << (i + 1) << ": " << samples[i].name
				     << "\t" << samples[i].path << endl;
				_Para single_para = para;
				single_para.Infilename = samples[i].path;
				single_para.Out_path = para.Out_path + "/" + samples[i].name;
				single_para.Sample_name = samples[i].name;
				string cmd = "mkdir -p " + single_para.Out_path + " && mkdir -p " + single_para.Out_path + "/tmp";
				system(cmd.c_str());

				Preprocess(single_para);
			}

			// Run PM-profiler
			Run_Profiler(para);

			map<string, tuple<int, int, int>> stat_map;
			ifstream stat_file((para.Singlesample_dir + "/pm_report.list").c_str());
			string stat_line;
			while (getline(stat_file, stat_line)) {
				istringstream iss(stat_line);
				string name;
				iss >> name >> seq_count >> rna_count >> asv_count;
				stat_map[name] = make_tuple(seq_count,rna_count,asv_count);
			}
			stat_file.close();

			// parse_taxonomy and make_plot and func
			for (size_t i = 0; i < samples.size(); ++i) {
				cout << "\nProcessing sample " << (i + 1) << ": " << samples[i].name
				     << "\t" << samples[i].path << endl;
				_Para single_para = para;
				single_para.Infilename = samples[i].path;
				single_para.Out_path = para.Out_path + "/" + samples[i].name;
				std::tie(seq_count,rna_count,asv_count) = stat_map[samples[i].name];

				Parse_Plot_Func(single_para, seq_count, rna_count, asv_count);
			}
			system(("rm -rf " + para.Out_path + "/tmp").c_str());
			system(("rm " + para.Singlesample_dir + "/pm_report.list").c_str());

		} else { // vsearch
			for (size_t i = 0; i < samples.size(); ++i) {
				cout << "\nProcessing sample " << (i + 1) << ": " << samples[i].name
				     << "\t" << samples[i].path << endl;
				_Para single_para = para;
				single_para.Infilename = samples[i].path;
				single_para.Out_path = para.Out_path + "/" + samples[i].name;
				string cmd = "mkdir -p " + single_para.Out_path + " && mkdir -p " + single_para.Out_path + "/tmp";
				system(cmd.c_str());
				Single_Run(single_para);
			}
			system(("rm -rf " + para.Out_path + "/tmp").c_str());
		}
	} else {
		cerr << "Error: Unknown run_mode: " << para.run_mode << endl;
		return 1;
	}

	cout << endl << "Parallel-Meta Suite Finished" << endl;
	cout << "Please check the analysis results and report at " << para.Out_path <<endl;

	return 0;
}
