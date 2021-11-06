// Updated at Nov 25, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <ext/hash_set>
#include <ext/hash_map>

#ifndef HASH_H
#define HASH_H

using namespace __gnu_cxx;

struct std_string_hash
	{                                                                                           
		size_t operator()( const std::string& x ) const                                           
		{                                                                                         
                                                                                                  
			return  __gnu_cxx::hash< const char* >()( x.c_str() );                                              
		}                                                                                         
	};

#endif
