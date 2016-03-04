// =============================================================================
// CD-HIT
// http://cd-hit.org/
// http://bioinformatics.burnham-inst.org/cd-hi
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
//
// modified by:
//                    Thomas Lin Pedersen
//                    Center for Biological Sequencing (CBS), DTU
//                    2300 Kongens Lyngby, Denmark
//                    Email: thomasp85@gmail.com
//                    
// =============================================================================

#include "cdhit-common.h"
#include <Rcpp.h>
using namespace Rcpp;

Options options;
SequenceDB seq_db;

//[[Rcpp::export]]
IntegerVector cdhitC(List opts, CharacterVector name, bool showProgress) {
    CharacterVector arguments = opts.names();
    std::string argn, argv;
    std::vector<int> clusters;
    for (int i = 0; i < opts.size(); i++) {
        argn = "-" + as<std::string>(arguments[i]);
        argv = as<std::string>(opts[i]);
        options.SetOption(argn.c_str(), argv.c_str());
    }
    options.Validate();
    
    InitNAA( MAX_UAA );
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];
    
    seq_db.Read(options.input.c_str(), options);
    seq_db.SortDivide(options);
    seq_db.DoClustering(options, as<std::string>(name), showProgress);
    
    //seq_db.WriteClusters(options.input.c_str(), options.input.c_str(), options);
    clusters = seq_db.GetClusters(options);
    
    return wrap(clusters);
}

// ////////////////////////////////////  MAIN /////////////////////////////////////
// int main(int argc, char *argv[])
// {
// 	string db_in;
// 	string db_out;
// 
// 	float begin_time = current_time();
// 	float end_time;
// 
// 	// ***********************************    parse command line and open file
// 	if (argc < 5) print_usage(argv[0]);
// 	if (options.SetOptions( argc, argv ) == 0) print_usage(argv[0]);
// 	options.Validate();
// 
// 	db_in = options.input;
// 	db_out = options.output;
// 
// 	InitNAA( MAX_UAA );
// 	options.NAAN = NAAN_array[options.NAA];
// 	seq_db.NAAN = NAAN_array[options.NAA];
// 
// 	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );
// 
// 	seq_db.Read( db_in.c_str(), options );
// 	cout << "total seq: " << seq_db.sequences.size() << endl;
// 	seq_db.SortDivide( options );
// 
// 	seq_db.DoClustering( options );
// 
// 	printf( "writing new database\n" );
// 	seq_db.WriteClusters( db_in.c_str(), db_out.c_str(), options );
// 
// 	// write a backup clstr file in case next step crashes
// 	seq_db.WriteExtra1D( options );
// 	cout << "program completed !" << endl << endl;
// 	end_time = current_time();
// 	printf( "Total CPU time %.2f\n", end_time - begin_time );
// 	return 0;
// } // END int main
