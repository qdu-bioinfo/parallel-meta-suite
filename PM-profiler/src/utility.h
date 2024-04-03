//PM profiler
//QDU Bioinfo Group
//Jan 9, 2024
//Haobo Shi, Xiaoquan Su
#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/dir.h>
#include <unistd.h>
#include <queue>
#include <set>

#include "version.h"

using namespace std;

int Check_Path(const char * path, int type){//0: remove and mkdir; 1: mkdir
    
    if (strlen(path) < 1) return 0;
    
    DIR *pDir = opendir(path);
            
    if(pDir!=NULL){
                  closedir(pDir);
                  if (type == 0){
                     string command = "rm -rf ";
                     command += path;
                     system(command.c_str());
                     mkdir(path, 0755);
                     }
                  }
    else
         mkdir(path, 0755);
    return 0;
    
    }

bool Check_Path(const char * path){
    
    if (strlen(path) < 1) return false;
    
    DIR *pDir = opendir(path);
    
    if (pDir != NULL){
             closedir(pDir);
             return true;
             }
    
    return false;
    }

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids, string prefix){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(prefix + str_temp[1]);
                                              
                                              }
                          
                          else{
                               list.push_back(prefix + buffer);
                               int name_begin = buffer.find_first_of('/');
                               int name_end = buffer.find_last_of('/');
                               if (name_begin == name_end)
                                              name_begin = 0;
                               else name_begin = buffer.rfind('/', name_end-1)+1;
                               ids.push_back(buffer.substr(name_begin, name_end - name_begin));
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();
    
    return count;
    }
#endif
