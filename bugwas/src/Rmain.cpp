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
#include "Rmain.h"
#include <R.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "param.h"

#ifdef FORCE_FLOAT
#include "gemma_float.h"
#else
#include "gemma.h"
#endif

using namespace std;



void Rmain(int* argv_nrow, int* argv_ncol, char** argv, int* EXIT_STATUS)
{
	int argc = *argv_nrow;
	
	/* Debug: output arguments to cout, bookended by '|'
	cout << "Arguments:\n";
	for(int i=0;i<argc;i++) {
		cout << "|";
		for(int j=0;j<*argv_ncol;j++) {
			if(argv[i][j]=='\0') break;
			cout << argv[i][j];
		}
		cout << "|\n";
	}
	cout << "End of arguments\n"; */
	
	GEMMA cGemma;	
	PARAM cPar;

	if (argc <= 1) {
		cGemma.PrintHeader();
		*EXIT_STATUS = EXIT_SUCCESS;
		return;
	}
	if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		cGemma.PrintHelp(0);
		*EXIT_STATUS = EXIT_SUCCESS;
		return;
	}
	if (argc==3 && argv[1][0] == '-' && argv[1][1] == 'h') {
		string str;
		str.assign(argv[2]);
		cGemma.PrintHelp(atoi(str.c_str()));
		*EXIT_STATUS = EXIT_SUCCESS;
		return;
	}
	if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'l') {
		cGemma.PrintLicense();
		*EXIT_STATUS = EXIT_SUCCESS;
		return;
	}
	
	ifstream check_dir("output/");
	if (!check_dir) {
		mkdir("output", S_IRWXU);
	}	
	
	cGemma.Assign(argc, argv, cPar); 
	
	if (cPar.error==true) {
		*EXIT_STATUS = EXIT_FAILURE;
		return;
	}
	
	cPar.CheckParam();
	
	if (cPar.error==true) {
		*EXIT_STATUS = EXIT_FAILURE;
		return;
	}
	
	cGemma.BatchRun(cPar);
	
	if (cPar.error==true) {
		*EXIT_STATUS = EXIT_FAILURE;
		return;
	}
	
	cGemma.WriteLog(argc, argv, cPar);
	return;
}


 
