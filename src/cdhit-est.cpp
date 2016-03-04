// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
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
// //over-write some defs in cd-hi.h for est version
// #undef MAX_UAA
// #define MAX_UAA 4
// 
//over-write some defs in cd-hi-init.h for est version

void setaa_to_na();
void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> & Comp_AAN_idx);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);

#include <Rcpp.h>
using namespace Rcpp;

#define MAX_UAA_EST 4

Options options_est;
SequenceDB seq_db_est;

//[[Rcpp::export]]
IntegerVector cdhitestC(List opts, CharacterVector name, bool showProgress) {
    options_est.cluster_thd = 0.95;
    options_est.NAA = 10;
    options_est.NAAN = NAA8;
    seq_db_est.NAAN = NAA8;
    options_est.NAA_top_limit = 12;
    setaa_to_na();
    mat.set_to_na(); //mat.set_gap(-6,-1);
    
    CharacterVector arguments = opts.names();
    std::string argn, argv;
    std::vector<int> clusters;
    for (int i = 0; i < opts.size(); i++) {
        argn = "-" + as<std::string>(arguments[i]);
        argv = as<std::string>(opts[i]);
        options_est.SetOption(argn.c_str(), argv.c_str());
    }
    options_est.Validate();
    
    InitNAA( MAX_UAA_EST );
    seq_db_est.NAAN = NAAN_array[options_est.NAA];
    
    seq_db_est.Read(options_est.input.c_str(), options_est);
    seq_db_est.SortDivide(options_est);
    seq_db_est.DoClustering(options_est, as<std::string>(name), showProgress);
    
    //seq_db.WriteClusters(options.input.c_str(), options.input.c_str(), options);
    clusters = seq_db_est.GetClusters(options_est);
    
    return wrap(clusters);
}

////////////////////////////////////  MAIN /////////////////////////////////////
// int main(int argc, char **argv) 
// {
// 	string db_in;
// 	string db_out;
// 
// 	options.cluster_thd = 0.95;
// 	options.NAA = 10;
// 	options.NAAN = NAA8;
// 	seq_db.NAAN = NAA8;
// 	options.NAA_top_limit = 12;
// 	setaa_to_na();
// 	mat.set_to_na(); //mat.set_gap(-6,-1);
// 
// 	float begin_time = current_time();
// 	float end_time;
// 
// 	// ***********************************    parse command line and open file
// 	if (argc < 5) print_usage_est(argv[0]);
// 	if (options.SetOptions( argc, argv, false, true ) == 0) print_usage_est(argv[0]);
// 	options.Validate();
// 
// 	db_in = options.input;
// 	db_out = options.output;
// 
// 	InitNAA( MAX_UAA );
// 	seq_db.NAAN = NAAN_array[options.NAA];
// 
// 	if ( options.option_r ) {
// 		Comp_AAN_idx.resize( seq_db.NAAN );
// 		make_comp_short_word_index(options.NAA, NAAN_array, Comp_AAN_idx);
// 	}
// 
// 	seq_db.Read( db_in.c_str(), options );
// 	cout << "total seq: " << seq_db.sequences.size() << endl;
// 	seq_db.SortDivide( options );
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
// }
