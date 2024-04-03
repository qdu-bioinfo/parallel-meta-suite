#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<math.h>
#include<string.h>
#include<limits.h>
#include<malloc.h>
#include<algorithm>
#include<time.h>
#include<set>
#include<unordered_set>
#include<cstring>
#include<map>
#include<unordered_map>

using namespace std;

int main(){
    clock_t f,e;
    f = clock();
    int m = 1,s = 0;
    for(int i = 0;i < 10000000;i++){
        if(m != 0) s++;
    }
    e = clock();
    cout << "The run time is: " <<(double)(e - f) / CLOCKS_PER_SEC << "s" << endl; 
}