// Creat at Nov 1, 2021
// Creat by Jian Lee
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Notes: Only for Bioconda version, download runtime environment of Parallel-META

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include <sys/stat.h>

#include <unistd.h>
#include <sys/types.h>

#include "utility.h"
#include "hash.h"

void Check_Database()
{
    string env_Path = Check_Env();
    string gg_13_path = env_Path + "/databases/gg_13/";
    string gg_13_99_path = env_Path + "/databases/gg_13_99/";
    bool this_part_ok = true;

    const int gg_13_files = 8;
    string name[gg_13_files] = {"copy_number.txt", "database.fa", "otu_nsti.tab", "taxonomy_annotation.txt", "KO/ko.tab", "tree/dist.txt", "tree/id.txt", "tree/order.txt"};

    for (int i = 0; i < gg_13_files; i++)
    {
        if (!Check_File((gg_13_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    for (int i = 0; i < gg_13_files; i++)
    {
        if (!Check_File((gg_13_99_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/databases_gg_13.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/databases_gg_13.tar.gz");
        system("mkdir $ParallelMETA/databases");
        system("tar -xzvf $ParallelMETA/databases_gg_13.tar.gz -C $ParallelMETA/databases");
    }
}

void Check_Database2()
{
    string env_Path = Check_Env();
    string silva_16s_path = env_Path + "/databases/silva_16s";
    string silva_18s_path = env_Path + "/databases/silva_18s";
    bool this_part_ok = true;

    const int silva_files = 8;
    string name[silva_files] = {"copy_number.txt", "database.fa", "otu_nsti.tab", "taxonomy_annotation.txt", "KO/ko.tab", "tree/dist.txt", "tree/id.txt", "tree/order.txt"};

    for (int i = 0; i < silva_files; i++)
    {
        if (!Check_File((silva_16s_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }
    for (int i = 0; i < silva_files; i++)
    {
        if (!Check_File((silva_18s_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/databases_silva.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/databases_silva.tar.gz");
        system("mkdir $ParallelMETA/databases");
        system("tar -xzvf $ParallelMETA/databases_silva.tar.gz -C $ParallelMETA/databases");
    }
}

void Check_Database3()
{
    string env_Path = Check_Env();
    string its_path = env_Path + "/databases/its/";
    string KO_path = env_Path + "/databases/KO/";
    string oral_core_path = env_Path + "/databases/oral_core/";
    bool this_part_ok = true;

    const int its_files = 4;
    const int KO_files = 3;
    const int oral_core_files = 9;
    string name[its_files] = {"copy_number.txt", "database.fa", "taxonomy_annotation.txt", "tree/id.txt"};
    string name2[KO_files] = {"ko_des.tab", "ko_id.tab", "ko_pw.tab"};
    string name3[oral_core_files] = {"copy_number.txt", "database.fa", "otu_nsti.tab", "taxonomy_annotation.txt", "KO/ko.tab", "tree/core.tre", "tree/dist.txt", "tree/id.txt", "tree/order.txt"};

    for (int i = 0; i < its_files; i++)
    {
        if (!Check_File((its_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }
    for (int i = 0; i < KO_files; i++)
    {
        if (!Check_File((KO_path + name2[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }
    for (int i = 0; i < oral_core_files; i++)
    {
        if (!Check_File((oral_core_path + name3[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    this_part_ok = this_part_ok & Check_File((env_Path + "/databases/db.config").c_str());

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/databases_other.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/databases_other.tar.gz");
        system("mkdir $ParallelMETA/databases");
        system("tar -xzvf $ParallelMETA/databases_other.tar.gz -C $ParallelMETA/databases");
    }
}
void Check_Rscript()
{
    string env_Path = Check_Env();
    string rscript_path = env_Path + "/Rscript/";
    bool this_part_ok = true;

    const int r_files = 13;
    string name[r_files] = {"config.R", "PM_Adiversity.R", "PM_Bdiversity.R",
                            "PM_Distribution.R", "PM_Hcluster.R", "PM_Heatmap.R",
                            "PM_Marker_Corr.R", "PM_Marker_RFscore.R", "PM_Marker_Test.R",
                            "PM_Network.R", "PM_Pca.R", "PM_Pcoa.R", "rarefaction.R"};

    for (int i = 0; i < r_files; i++)
    {
        if (!Check_File((rscript_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/Rscript.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/Rscript.tar.gz");
        system("tar -xzvf $ParallelMETA/Rscript.tar.gz -C $ParallelMETA/");
        system("Rscript $ParallelMETA/Rscript/config.R");
    }
}

void Check_Models()
{
    string env_Path = Check_Env();
    string model_path = env_Path + "/models/";
    bool this_part_ok = true;

    const int model_files = 4;
    string name[model_files] = {"arc_ssu.hmm", "bac_ssu.hmm", "bac_ssu.hmm.bk2",
                                "euk_ssu.hmm"};

    for (int i = 0; i < model_files; i++)
    {
        if (!Check_File((model_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/models.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/models.tar.gz");
        system("tar -xzvf $ParallelMETA/models.tar.gz -C $ParallelMETA/");
    }
}

void Check_Html()
{
    string env_Path = Check_Env();
    bool this_part_ok = true;

    const int r_files = 2;
    string name[r_files] = {"/html/index.html", "/PMS-config/index.html"};

    for (int i = 0; i < r_files; i++)
    {
        if (!Check_File((env_Path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/html.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/html.tar.gz");
        system("tar -xzvf $ParallelMETA/html.tar.gz -C $ParallelMETA/");
    }
}

void Check_Example()
{
    string env_Path = Check_Env();
    string rscript_path = env_Path + "/example/";
    bool this_part_ok = true;

    const int r_files = 3;
    string name[r_files] = {"meta.txt", "seqs.list", "Readme"};

    for (int i = 0; i < r_files; i++)
    {
        if (!Check_File((rscript_path + name[i]).c_str()))
        {
            this_part_ok = false;
            break;
        }
    }

    if (!this_part_ok)
    {
        system("rm -rf $ParallelMETA/example.tar.gz");
        system("wget -P $ParallelMETA/ http://bioinfo.single-cell.cn/Released_Software/parallel-meta/conda/example.tar.gz");
        system("tar -xzvf $ParallelMETA/example.tar.gz -C $ParallelMETA/");
    }
}

int main(int argc, char *argv[])
{

    cout << "Parallel-META Suite Environment Check Start..." << endl;

    cout << "Checking Rscripts..." << endl;

    Check_Rscript();

    cout << "Checking Databases..." << endl;
    Check_Database();
    Check_Database2();
    Check_Database3();

    cout << "Checking HTML GUI..." << endl;

    Check_Html();

    cout << "Checking Examples..." << endl;

    Check_Example();

    cout << "Checking models..." << endl;

    Check_Models();

    return 0;
}