// Updated at July 28, 2021
// Updated by Yuzhu Chen and Xiaoquan Su
// Code by Yuzhu Chen, Xiaoquan Su
// Bioinformatics Group, College of Computer Science and Technology, Qingdao University

#include "pipeline.h"
#include "utility.h"

using namespace std;

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    char command[BUFFER_SIZE];
    //Parallel for R
    vector <string> command_parallel_scripts; 
    vector <string> command_parallel_titles;
    
    vector <string> command_parallel_bdiversity_scripts; 
    vector <string> command_parallel_bdiversity_titles; 
	        
    Check_Path(Out_path.c_str(), 1);
        
    //open script file
    ofstream outscript((Out_path + "/scripts.sh").c_str(), ofstream::out);
    if (!outscript){
                    cerr << "Warning: Cannot open output file : " << Out_path << "/scripts.sh" << endl;                    
    }
    
    cout << "Parallel-Meta Suite Version: " << Version << endl;
    outscript << "#Parallel-Meta Suite Version: " << Version << endl;
    
    cout << "The reference sequence database is ";
    cout << Database.Get_Description() << endl;
    
    outscript << "#The reference sequence database is ";
    outscript << Database.Get_Description() << endl;
    
    //check metadata
    if (Load_ID(Meta_file.c_str(), Ids, 1) == 0){
        string error_info = "Error: Please check the Meta data file (-m): at least contains 1 columns of Sample ID";
        cerr << error_info << endl;
        Echo_Error(error_info.c_str(), Error_file.c_str());
        return 0;
    }
    
    if (!Check_Ids(Ids)){
        string error_info = "Warning: Sample ID starts by number";
        cerr << error_info << endl;
        Echo_Error(error_info.c_str(), Error_file.c_str());
        //return 0;
    }
    
    Input_sam_num = Ids.size();
	
    switch (Step){
        
    //Step 0: Parallel-META
    case 0:
               Check_Path(Singlesample_dir.c_str(), 1);
               Check_Path(Singlesamplelist_dir.c_str(), 1);
         	   Check_Path(Temp_dir.c_str(), 1);
               
               if (Load_List(Seq_list_file.c_str(), Seq_files, List_prefix) == 0){
                                              string error_info = "Error: Please check the sequence list file (-i) or the list path prefix (-p)";
                                              cerr << error_info << endl;
                                              Echo_Error(error_info.c_str(), Error_file.c_str());
                                              return 0;
                                              };
               
               //is_pair_end
               if (Seq_files.size() == Ids.size() * 2){
                                    Is_paired_seq = true;
                                    if (Seq_type != 'r'){                                                 
                                                 string error_info = "Error: Pair-end sequences only support 16S rRNA (-M F)";
                                                 cerr << error_info << endl;
                                                 Echo_Error(error_info.c_str(), Error_file.c_str());
                                                 return 0;
                                    }
                                    for (int i = 0; i < Ids.size(); i ++){
                                		if( Check_Format( Seq_files[i * 2].c_str() )== 0 || Check_Format(Seq_files[i * 2 +1].c_str() ) ==0 ){//fasta format
											cerr << "Error: For pair-end sequences only support fastq format" << endl;
											return 0;
										}
									}
                }
               else if (Seq_files.size() != Ids.size()){
                                    string error_info = "Error: Sequence files (pairs) and meta data should have the same sample number and order";
                                    cerr << error_info << endl;
                                    Echo_Error(error_info.c_str(), Error_file.c_str());
                                    return 0;
                }
            
               //format_seq
               /*
               if (Is_format_check){
               cout << endl << "Sequence format check" << endl;    
               outscript << endl << "#Sequence format check" << endl;                             
               
               if (List_prefix.size() > 0) // add prefix
                  sprintf(command, "%s/format-seq -l %s -p %s", Bin_path.c_str(), Seq_list_file.c_str(), List_prefix.c_str());          
               else
                  sprintf(command, "%s/format-seq -l %s", Bin_path.c_str(), Seq_list_file.c_str());                               
               Run_With_Error(command, "format-seq", Error_file.c_str());
               outscript << command << endl;
               }
               */
               
               //check search database similarity
    			if (db_similarity <= 0 || db_similarity > 1){
    						cerr <<"Error: Please input right similarity, the value range is from 0 to 1" << endl;
    						exit(0);
				}
				
               //profiling
               cout << endl << "Microbial Community profiling" << endl;    
               outscript << endl << "#Microbial Community profiling" << endl;   
               
               //taxa               
               for (int i = 0; i < Ids.size(); i ++){
                               
                               cout << endl << "Processing sample " << i + 1 << " of " << Ids.size() << endl;   
                               if (Is_paired_seq) //pair -end                            
                                  sprintf(command, "%s/PM-parallel-meta -r %s -R %s -o %s -t %d -f F -k %c -D %c -v %c -c %c -d %.2f", Bin_path.c_str(), Seq_files[i * 2].c_str(), Seq_files[i * 2 + 1].c_str(), (Singlesample_dir + "/" + Ids[i]).c_str(), Coren, Is_format_check, Ref_db, Is_denoised, Is_nonchimeras, db_similarity);
                               else
                                   sprintf(command, "%s/PM-parallel-meta -%c %s -o %s -t %d -f F -L %d -k %c -D %c -v %c -c %c -d %.2f", Bin_path.c_str(), Seq_type, Seq_files[i].c_str(), (Singlesample_dir + "/" + Ids[i]).c_str(), Coren, Length_t, Is_format_check, Ref_db, Is_denoised, Is_nonchimeras, db_similarity);
                               Run_With_Error(command, "PM-parallel-meta", tmpError_file.c_str());
                               
                               //system(command);
                               outscript << command << endl;
                               
                               }                              
              //taxa list
              Taxa_list_file = Singlesamplelist_dir + "/taxa.list";
              Make_list(Taxa_list_file.c_str(), Singlesample_dir.c_str(), Ids, 0);
              
              List_prefix = "";
              //Step 0 finished
            
    //Step 1: check path and make OTU table
    case 1: //check output dir
           Check_Path(Abd_dir.c_str(), 1);
           Check_Path(Dist_dir.c_str(), 1);
           Check_Path(Clust_dir.c_str(), 1);   
           Check_Path(Marker_dir.c_str(), 1);
           Check_Path(Network_dir.c_str(), 1);
           Check_Path(Alpha_dir.c_str(), 1);
           Check_Path(Beta_dir.c_str(), 1);
           Check_Path(Sampleview_dir.c_str(), 1);
           Check_Path(Temp_dir.c_str(), 1);
           
           //make OTU table
           switch (Mode){
                  case 0: //list
                         if (Taxa_list_file.size() == 0){
                              string error_info = "Error: Please check the taxa list (-l)";
                              cerr << error_info << endl;
                              Echo_Error(error_info.c_str(), Error_file.c_str());
                              return 0;
                              }
                                 
                         //If prefix
                         if (List_prefix.size() > 0){ // Add prefix to list
                                  Check_Path(Singlesamplelist_dir.c_str(), 1);
                                  Add_list_prefix(Taxa_list_file.c_str(), List_prefix.c_str(), (Singlesamplelist_dir + "/taxa.list").c_str());
                                  Taxa_list_file = Singlesamplelist_dir + "/taxa.list";
                                  }
                                                                                                             
                         sprintf(command, "%s/PM-select-taxa -l %s -o %s/taxa -L %d -m 0 -n 0 -z 0 -v 0 -q 0 -D %c -r %c", Bin_path.c_str(), Taxa_list_file.c_str(), Abd_dir.c_str(), 7, Ref_db, Is_cp);
                         Run_With_Error(command, "PM-select-taxa", Error_file.c_str());
                         outscript << command << endl;
                         break;
                         
                  case 1: //table
                         sprintf(command, "%s/PM-select-taxa -T %s -o %s/taxa -L %d -m 0 -n 0 -z 0 -v 0 -q 0 -D %c -r %c", Bin_path.c_str(), Table_file.c_str(), Abd_dir.c_str(), 7, Ref_db, Is_cp);
                         Run_With_Error(command, "PM-select-taxa", Error_file.c_str());
                         outscript << command << endl;
                         break;
                  }
                                                                                           
    //Step 2: Rarefaction & Normalization
    if (Is_rare){
            cout << endl << "Sequence Normalization" << endl;
            outscript << endl << "#Sequence Normalization" << endl;
            //rand-rare                 
            sprintf(command, "%s/PM-rand-rare -T %s/taxa.OTU.Count -o %s/taxa.OTU.Count -b %d -s %d -D %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), Bootstrap, Rare_depth, Ref_db);
            Run_With_Error(command, "PM-rand-rare", Error_file.c_str());
            outscript << command << endl;
            
            //make OTU table
            sprintf(command, "%s/PM-select-taxa -T %s/taxa.OTU.Count -o %s/taxa -L %d -m 0 -n 0 -z 0 -v 0 -q 0 -D %c -r %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), 7, Ref_db, Is_cp);
            Run_With_Error(command, "PM-select-taxa", Error_file.c_str());
            outscript << command << endl;
            }
            
    //Ids reload from table file
    Ids.clear();
    Load_ID((Abd_dir + "/taxa.OTU.Count").c_str(), Ids, 1); 
    //Update metadata
    if (!(Check_Metadata_By_Ids(Ids, Meta_file, Out_path + "/meta.txt")))
                return 0;
    
    Filter_sam_num = Ids.size();
    
    if (Filter_sam_num < MIN_SAM_NUM){
        
       char error_info [BUFFER_SIZE];
       sprintf(error_info, "Error: %d sample(s) passed the filtering", Filter_sam_num);
       cerr << error_info << endl;
       Echo_Error(error_info, Error_file.c_str());
       
       sprintf(error_info, "The required minimum sample number is %d", MIN_SAM_NUM); 
       cerr << error_info << endl;
       Echo_Error(error_info, Error_file.c_str()); 
       return 0;
       }
        
    //Step 3: Feature select
    cout << endl << "Feature Selection" << endl;
    outscript << endl << "#Feature Selection" << endl;
           //taxa-sel
           if (Is_taxa){                 
           for (int i = 0; i < TLevN - 1; i ++) // no OTU
               if (TLevel_Set[i]){             
               sprintf(command, "%s/PM-select-taxa -T %s/taxa.OTU.Count -o %s/taxa -L %d -m 0 -n 0 -D %c -r %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), i + 1, Ref_db, Is_cp);
               Run_With_Error(command, "PM-select-taxa", Error_file.c_str());
               outscript << command << endl;
               //plot
               sprintf(command, "Rscript %s/PM_Distribution.R -m %s -i %s/taxa.%s.Abd -o %s -p taxa.%s.Abd", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Abd_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Distribution.R");
               outscript << command << endl;
               }
           
            }
    //Step 4: Function Prediction
        if (Is_func){
        
        cout << endl << "Function Prediction" << endl;
        outscript << endl << "#Function Prediction" << endl;
                
        //Parse func
            //Check_Path(Singlesample_dir.c_str(), 1);
            sprintf(command, "%s/PM-predict-func -T %s/taxa.OTU.Count -o %s/func -t %d -D %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), Coren, Ref_db);
            Run_With_Error(command, "PM-predict-func", Error_file.c_str());
            outscript << command << endl;
        
        //calc func-nsti
            sprintf(command, "%s/PM-predict-func-nsti -T %s/taxa.OTU.Count -o %s/func.nsti -D %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), Ref_db);
            Run_With_Error(command, "PM-predict-func-nsti", Error_file.c_str());
            outscript << command << endl;
            
        //func list
            //Check_Path(Singlesamplelist_dir.c_str(), 1);
            //Func_list_file = Singlesamplelist_dir + "/func.list";
            //Make_list(Func_list_file.c_str(), Singlesample_dir.c_str(), Ids, 1);
            
        //func-sel
           // KO table
           //sprintf(command, "%s/PM-select-func -l %s -o %s/func -L %d -D %c", Bin_path.c_str(), Func_list_file.c_str(), Abd_dir.c_str(), 4, Ref_db);
           //Run_With_Error(command, "PM-select-func", Error_file.c_str());
           //outscript << command << endl;
           
           for (int i = 0; i < FLevN - 1; i ++) // no KO
               if (FLevel_Set[i]){
               sprintf(command, "%s/PM-select-func -T %s/func.KO.Count -o %s/func -L %d -D %c", Bin_path.c_str(), Abd_dir.c_str(), Abd_dir.c_str(), i + 1, Ref_db);
               Run_With_Error(command, "PM-select-func", Error_file.c_str());
               outscript << command << endl;
               //plot
               sprintf(command, "Rscript %s/PM_Distribution.R -m %s -i %s/func.%s.Abd -o %s -p func.%s.Abd", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Abd_dir.c_str(), FLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Distribution.R");
               outscript << command << endl;
               }
           }
    
    //Step 5: Sample visualization
    //if (Ids.size() <= 300){
    cout << endl << "Sample visualization" << endl;
    outscript << endl << "#Sample visualization" << endl;
            //class-tax
            if (Is_taxa){
                sprintf(command, "%s/PM-plot-taxa -T %s/taxa.OTU.Count -o %s/ -D %c", Bin_path.c_str(), Abd_dir.c_str(), Sampleview_dir.c_str(), Ref_db);
                Run_With_Error(command, "PM-plot-tax", Error_file.c_str());
                outscript << command << endl;
            }
    //    }

    //Step 6: Dist
    cout << endl << "Dist Calculation" << endl;
    outscript << endl << "#Dist Calculation" << endl;
           //calc taxa dist
           if (Is_taxa){
           
           if ((Taxa_dist_type == 0) || (Taxa_dist_type == 2)){
           sprintf(command, "%s/PM-comp-taxa -T %s/taxa.OTU.Count -o %s/taxa.dist -t %d -P T -c %d -D %c", Bin_path.c_str(), Abd_dir.c_str(), Dist_dir.c_str(), Coren, Cluster, Ref_db);
           Run_With_Error(command, "PM-comp-taxa", Error_file.c_str());
           outscript << command << endl;
           }
           
           if ((Taxa_dist_type == 1) || (Taxa_dist_type == 2)){
           sprintf(command, "%s/PM-comp-taxa -T %s/taxa.OTU.Count -o %s/taxa.unweighted.dist -t %d -P T -c %d -D %c -M 1", Bin_path.c_str(), Abd_dir.c_str(), Dist_dir.c_str(), Coren, Cluster, Ref_db);
           Run_With_Error(command, "PM-comp-taxa", Error_file.c_str());
           outscript << command << endl;
           }
           }
           
           if (Is_func){      
           //calc func dist
           sprintf(command, "%s/PM-comp-func -T %s/func.KO.Count -o %s/func.dist -t %d -P T -c %d -D %c", Bin_path.c_str(), Abd_dir.c_str(), Dist_dir.c_str(), Coren, Cluster, Ref_db);
           Run_With_Error(command, "PM-comp-func", Error_file.c_str());
           outscript << command << endl;
           }
           
    //Step 7: Co-occur network
    cout << endl << "Correlation Calculation" << endl;
    outscript << endl << "#Correlation Calculation" << endl;
           //taxa-corr                 
           if (Is_taxa){
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){               
               sprintf(command, "%s/PM-comp-corr -i %s/taxa.%s.Count -o %s/taxa.%s -f 1 -N T -G %f", Bin_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Network_dir.c_str(), TLevel[i].c_str(), Network_t);
               Run_With_Error(command, "PM-comp-corr", Error_file.c_str());
               outscript << command << endl;
               }
           }
           
    //Step 8: Rarefaction curve
    //might be slow
    if (Is_rare_curve){
    cout << endl << "Rarefaction Analysis" << endl;
    outscript << endl << "#Rarefaction Analysis" << endl;
           //OTU based
           if (Is_taxa){            
           sprintf(command, "%s/PM-rare-curv -i %s/taxa.OTU.Count -o %s -p taxa.OTU -b 10", Bin_path.c_str(), Abd_dir.c_str(), Alpha_dir.c_str());
           Run_With_Error(command, "PM-rare-curv", Error_file.c_str());
           outscript << command << endl;
           }
    }     
    
    //Parallel Rscripts
    cout << endl << "Statistical Analysis" << endl;  
    //Step 9: PCoA
    //cout << endl << "PCoA Calculation" << endl;
    outscript << endl << "#PCoA Calculation" << endl;
           //taxa PCoA
           if (Is_taxa){
           
           if ((Taxa_dist_type == 0) || (Taxa_dist_type == 2)){
           sprintf(command, "Rscript %s/PM_Pcoa.R -m %s -d %s/taxa.dist -o %s/taxa.pcoa.pdf", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Clust_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Pcoa.R");
           outscript << command << endl;
           }
           
           if ((Taxa_dist_type == 1) || (Taxa_dist_type == 2)){
           sprintf(command, "Rscript %s/PM_Pcoa.R -m %s -d %s/taxa.unweighted.dist -o %s/taxa.unweighted.pcoa.pdf", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Clust_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Pcoa.R");
           outscript << command << endl;
           }
           }
           
           //func PCoA
           if (Is_func){
           sprintf(command, "Rscript %s/PM_Pcoa.R -m %s -d %s/func.dist -o %s/func.pcoa.pdf", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Clust_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Pcoa.R");
           outscript << command << endl;
           }
           
    //Step 10: PCA
    //cout << endl << "PCA Calculation" << endl;
    outscript << endl << "#PCA Calculation" << endl;
            
           //taxa PCA
           if (Is_taxa){
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){           
               sprintf(command, "Rscript %s/PM_Pca.R -m %s -i %s/taxa.%s.Abd -o %s/taxa.%s.pca.pdf", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Clust_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Pca.R");
               outscript << command << endl;
               }
           }
           
           //func PCA
           if (Is_func){
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){                        
                sprintf(command, "Rscript %s/PM_Pca.R -m %s -i %s/func.%s.Abd -o %s/func.%s.pca.pdf", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Clust_dir.c_str(), FLevel[i].c_str());
                command_parallel_scripts.push_back(command);
                command_parallel_titles.push_back("PM_Pca.R");
                outscript << command << endl;
               }
           }
           
    //Step 11: Multivariate Statistical Analysis
    //cout << endl << "Multivariate Statistical Analysis" << endl;
    outscript << endl << "#Multivariate Statistical Analysis" << endl;           
    
           //Parallel  
           //Alpha: taxa         
           if (Is_taxa){                        
               for (int i = 0; i < TLevN; i ++) //contain OTU
               if ((TLevel_Set[i]) || (i == TLevN - 1)){ // OTU                  
               sprintf(command, "Rscript %s/PM_Adiversity.R -m %s -i %s/taxa.%s.Abd -o %s -p taxa.%s", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Alpha_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Adiversity.R");
               outscript << command << endl;
               }                      
            
           //Beta: taxa dist
           if ((Taxa_dist_type == 0) || (Taxa_dist_type == 2)){
           sprintf(command, "Rscript %s/PM_Bdiversity.R -m %s -d %s/taxa.dist -o %s -p taxa.dist -n %s -t %d", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Beta_dir.c_str(), Taxa_dist_name.c_str(), Coren);
           command_parallel_bdiversity_scripts.push_back(command);
           command_parallel_bdiversity_titles.push_back("PM_Bdiversity.R");
		   //command_parallel_scripts.push_back(command);
           //command_parallel_titles.push_back("PM_Bdiversity.R");
           outscript << command << endl;   
           }       
           
           if ((Taxa_dist_type == 1) || (Taxa_dist_type == 2)){
           sprintf(command, "Rscript %s/PM_Bdiversity.R -m %s -d %s/taxa.unweighted.dist -o %s -p taxa.unweighted.dist -n %s -t %d", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Beta_dir.c_str(), Taxa_dist_name_uw.c_str(), Coren);
           command_parallel_bdiversity_scripts.push_back(command);
           command_parallel_bdiversity_titles.push_back("PM_Bdiversity.R");
		   //command_parallel_scripts.push_back(command);
           //command_parallel_titles.push_back("PM_Bdiversity.R");
           outscript << command << endl;                    
           }
           }
           
           //Alpha: func
           if (Is_func){
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){               
               sprintf(command, "Rscript %s/PM_Adiversity.R -m %s -i %s/func.%s.Abd -o %s -p func.%s", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Alpha_dir.c_str(), FLevel[i].c_str()); 
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Adiversity.R");
               outscript << command << endl;           
               }                     
                       
           //Beta: func dist
           sprintf(command, "Rscript %s/PM_Bdiversity.R -m %s -d %s/func.dist -o %s -p func.dist -n %s -t %d", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Beta_dir.c_str(), Func_dist_name.c_str(), Coren);
           command_parallel_bdiversity_scripts.push_back(command);
           command_parallel_bdiversity_titles.push_back("PM_Bdiversity.R");
		   //command_parallel_scripts.push_back(command);
           //command_parallel_titles.push_back("PM_Bdiversity.R");
           outscript << command << endl;           
           }
                   
    //Step 12: Test
    //cout << endl << "Marker Analysis" << endl;
    outscript << endl << "#Marker Analysis" << endl;
    
           //Parallel
           //taxa-test
           if (Is_taxa){                 
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){ 
               //Test         
               sprintf(command, "Rscript %s/PM_Marker_Test.R -i %s/taxa.%s.Abd -m %s -o %s -p taxa.%s -P %c", R_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), TLevel[i].c_str(), Is_pair);
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_Test.R");
               outscript << command << endl;     
               
               //RF
               sprintf(command, "Rscript %s/PM_Marker_RFscore.R -i %s/taxa.%s.Abd -m %s -o %s -p taxa.%s", R_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_RFscore.R");
               outscript << command << endl;    
               
               //Corr
               sprintf(command, "Rscript %s/PM_Marker_Corr.R -i %s/taxa.%s.Abd -m %s -o %s -p taxa.%s", R_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_Corr.R");
               outscript << command << endl;         
               }                                              
           }
           
           //func-test
           if (Is_func){                     
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){     
               //Test                                     
               sprintf(command, "Rscript %s/PM_Marker_Test.R -i %s/func.%s.Abd -m %s -o %s -p func.%s -P %c", R_path.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), FLevel[i].c_str(), Is_pair);
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_Test.R");
               outscript << command << endl;
               
               //RF
               sprintf(command, "Rscript %s/PM_Marker_RFscore.R -i %s/func.%s.Abd -m %s -o %s -p func.%s", R_path.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), FLevel[i].c_str());               
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_RFscore.R");
               outscript << command << endl;
               
               //Corr
               sprintf(command, "Rscript %s/PM_Marker_Corr.R -i %s/func.%s.Abd -m %s -o %s -p func.%s", R_path.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), FLevel[i].c_str());               
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker_Corr.R");
               outscript << command << endl;
               }                        
           }
    //Parallel bdiversity R
    omp_set_num_threads(Min(command_parallel_bdiversity_scripts.size(), Coren));    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < command_parallel_bdiversity_scripts.size(); i ++){
        char error_bdiv_parallel[BUFFER_SIZE];
        sprintf(error_bdiv_parallel, "%s/%d.log", Temp_dir.c_str(), i);
        Run_With_Error(command_parallel_bdiversity_scripts[i].c_str(), command_parallel_bdiversity_titles[i].c_str(),error_bdiv_parallel);
    }
     
    //Combine error
    for (int i = 0; i < command_parallel_bdiversity_scripts.size(); i ++){
        sprintf(command, "cat %s/%d.log >> %s", Temp_dir.c_str(), i, Error_file.c_str());
        system(command);
    }
         
    //Parallel R
     omp_set_num_threads(Min(command_parallel_scripts.size(), Coren));    
     #pragma omp parallel for schedule(dynamic, 1)
     for (int i = 0; i < command_parallel_scripts.size(); i ++){
         char error_parallel[BUFFER_SIZE];
         sprintf(error_parallel, "%s/%d.log", Temp_dir.c_str(), i);
         Run_With_Error(command_parallel_scripts[i].c_str(), command_parallel_titles[i].c_str(),error_parallel);
         }
     
     //Combine error
     for (int i = 0; i < command_parallel_scripts.size(); i ++){
         sprintf(command, "cat %s/%d.log >> %s", Temp_dir.c_str(), i, Error_file.c_str());
         system(command);
         }
                                                                                    
    default: break;
    }
    
    //if (rmdir(Temp_dir.c_str()) < 0){
        sprintf(command, "rm -rf %s", Temp_dir.c_str());
        system(command);
   //     }
    
    if (outscript){
                   outscript.close();
                   outscript.clear();
                   }
    
    Print_Report(Report_file.c_str());
    Copy_Index(Out_path.c_str());
    
    cout << endl << "Parallel-Meta Suite Pipeline Finished" << endl;
    outscript << endl << "#Parallel-Meta Suite Pipeline Finished" << endl;
    cout << "Please check the analysis results and report at " << Out_path <<endl;
    
    return 0;
    }

