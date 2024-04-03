###Parallel-META installer
###Updated at July 31, 2021
###Updated by Shihao Bo at April 3, 2024. Added compilation related to PM-profiler
###Bioinformatics Group, College of Computer Science and Technology, Qingdao University
###Code by: Yuzhu Chen, Xiaoquan Su, Honglei Wang, Gongchao Jing
#!/bin/bash
##Users can change the default environment variables configuration file here
if [[ $SHELL = '/bin/zsh' ]];
then
        PATH_File=~/.zshrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.zsh_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi
else
        PATH_File=~/.bashrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.bash_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi

fi
PM_PATH=`pwd`
Sys_ver=`uname`
###Checking that environment variable of Parallel-META exists###
Check_old_pm=`grep "export ParallelMETA"  $PATH_File|awk -F '=' '{print $1}'`
Check_old_path=`grep "ParallelMETA/bin"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
Check_old_Rscript_path=`grep "ParallelMETA/Rscript"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
Add_Part="####DisabledbyParallelMetaSuite####"
echo "**Parallel-Meta Suite Installation**"
echo "**version 3.7.2**"

###Build source code for src package###
if [ -f "Makefile" ]
   then
       echo -e "\n**Parallel-Meta Suite src package**"
       make
       echo -e "\n**Build Complete**"
else
   echo -e "\n**Parallel-Meta Suite bin package**"
fi
chmod +x $PM_PATH/bin/vsearch
chmod +x $PM_PATH/bin/hmmsearch
###Build source code for PM-profiler###
cd PM-profiler
if [ -f "Makefile" ]
   then
       echo -e "\n**PM-profiler src package**"
       make
       echo -e "\n**Build Complete**"
else
   echo -e "\n**PM-profiler bin package**"
fi
cd ..
###Configure environment variables###

if [ "$Check_old_pm" != "" ]
   then
      Checking=`grep ^export\ ParallelMETA  $PATH_File|awk -F '=' '{print $2}'`
      if [ "$Checking" != "$PM_PATH" ]
         then
         if [ "$Sys_ver" = "Darwin" ]
            then
            sed -i "" "s/^export\ ParallelMETA/$Add_Part\ &/g" $PATH_File
            sed -i "" -e "`grep -n "$Add_Part" $PATH_File | cut -d ":" -f 1 | head -1` a\ 
            export\ ParallelMETA=$PM_PATH
            " $PATH_File
         else
             sed -i "s/^export\ ParallelMETA/$Add_Part\ &/g" $PATH_File
             sed -i "/$Add_Part\ export\ ParallelMETA/a export\ ParallelMETA=$PM_PATH" $PATH_File
         fi
     fi    
elif [ "$Check_old_pm" = "" ]
    then
      echo "export ParallelMETA="${PM_PATH} >> $PATH_File
fi
if [ "$Check_old_path" = "" ]
    then
      echo "export PATH=\$PATH:\$ParallelMETA/bin" >> $PATH_File
fi
if [ "$Check_old_Rscript_path" = "" ]
    then
	  echo "export PATH=\$PATH:\$ParallelMETA/Rscript" >> $PATH_File
fi
###Source the environment variable file###
source $PATH_File
echo -e "\n**Environment Variables Configuration Complete**"
###Configurate the R packages###
echo -e ""
Rscript $PM_PATH/Rscript/config.R
chmod +x $PM_PATH/Rscript/PM_*.R
###End
echo -e "\n**Parallel-Meta Suite Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"


