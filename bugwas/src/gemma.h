/*
	Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011  Xiang Zhou and Matthew Stephens

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __GEMMA_H__                
#define __GEMMA_H__

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif

using namespace std;

class GEMMA {

public:			
	//parameters
	string version;
	string date;
	string year;
	
	//constructor
	GEMMA(void);
	
	//functions
	void PrintHeader (void);
	void PrintHelp (size_t option);
	void PrintLicense (void);
	void Assign (int argc, char **argv, PARAM &cPar);
	void BatchRun (PARAM &cPar);
	void WriteLog (int argc, char **argv, PARAM &cPar);
};


#endif

