#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "parameters.h"
#include "susypar.h"
#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;

//////////////////////////////////////////////////////////////////
///
/// modification of input parameters
/// default file to read input from: input.cfg
///
//////////////////////////////////////////////////////////////////
// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}


int configure(int ac, char* av[], string configfile, bool runupdate = true)
{
    try {

        int order; //order for LO, NLO etc.
        int resummed; // option for resummation (1) or scet (2)
        string resum_order; // LL, NLL, NNLL
        string process; //DY, higgs, dihiggs, WW, ZZ
        double comTeV; // sqrt(S) in TeV
        string config_file = configfile;
		// Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "produce help message")
            ("config,c", po::value<string>(&config_file)->default_value(configfile), "name of a file of a configuration.")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("include-path,I",
                 po::value< vector<string> >()->composing(),
                 "include path")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-file", po::value< vector<string> >(), "input file")
            //output
            ("output", po::value<string>(&result_map)->default_value("/home/mbeekveld/DY_numerical_NE/Numerical_Code_try/results"), "output directory")
            //PDFstuf
            ("setname", po::value<string>(&setname)->default_value("PDF4LHC15_nnlo_100"), "Name of PDFset")
            ("usemember", po::value<int>(&use_member)->default_value(0),"member of PDFset to use")
            ("fitPDF", po::value<bool>(&fitPDF)->default_value(false),"use own fits to PDF grid")
            ("toy_pdfs", po::value<bool>(&toy_pdfs)->default_value(false),"use toy PDFs")
            ("chebPDF", po::value<bool>(&chebPDF)->default_value(false),"use cheby fits to PDF grid")
            ("Nfix", po::value<bool>(&Nfixed)->default_value(false),"Use a fixed value for N in the resummation functions")
            //process
            ("sqrtS", po::value<double>(&comTeV)->default_value(13.), "center-of-mass energy [TeV]")
            ("process", po::value<string>(&process)->default_value("DY"),"process")
            ("observable", po::value<string>(&observable)->default_value("Qinv"),"observable")
            //perturbative order
            ("fixed_order", po::value<int>(&order)->default_value(2), "order of calculation")
            //resummation
            ("resum", po::value<int>(&resummed)->default_value(1), "resummation")
            ("inceuler", po::value<double>(&INCEULER)->default_value(1), "resum euler constant in pQCD")
            ("resum_order", po::value<string>(&resum_order)->default_value("NNLL"), "order to resum to")
            ("NLP", po::value<double>(&ISNLP)->default_value(0.), "include NLP contributions")
            ("INCSQRTZ", po::value<bool>(&INCSQRTZ)->default_value(false), "include sqrt(z) in SCET")
            ("INCHARD", po::value<bool>(&INCHARD)->default_value(false), "include hardpart in NLP NNLL SCET")
            //scales
            ("muR", po::value<double>(&muR)->default_value(500.), "renormalization scale [GeV]")
            ("muF", po::value<double>(&muF)->default_value(500.), "factorization scale [GeV]")
            ("Q", po::value<double>(&Q)->default_value(500.), "hard scale [GeV]")
            ("muh", po::value<double>(&muh)->default_value(500.), "SCET hard scale [GeV]")
            ("mus", po::value<double>(&mus)->default_value(500.), "SCET soft scale [GeV]")
            ("setdym", po::value<bool>(&setdym)->default_value(true), "set soft scale / factorization scale dynamically")
            ("highscale", po::value<bool>(&highscale)->default_value(true), "set soft scale / factorization scale dynamically")
            //masses
            ("MZ", po::value<double>(&mZ)->default_value(91.1876), "Z boson mass [GeV]")
            ("MW", po::value<double>(&mW)->default_value(80.419002), "W boson mass [GeV]")
            ("mh0", po::value<double>(&mH)->default_value(125.), "Higgs boson mass [GeV]")
            ("mt", po::value<double>(&mt)->default_value(173.1), "Top quark mass [GeV]")
            ("mb", po::value<double>(&mb)->default_value(4.18), "Bottom quark mass [GeV]")
            //SUSY
            ("setSUSY", po::value<double>(&SUSYset)->default_value(false), "relevant for dihiggs, use SUSY?")
            ("tanb", po::value<double>(&tb)->default_value(15.), "tan(beta)")
            //integration
            ("phiMP", po::value<double>(&phiMP)->default_value(3./4.*M_PI), "phiMP for inverse Mellin transform")
            ("CMP", po::value<double>(&CMP)->default_value(2.1), "CMP for inverse Mellin transform")
            ;


        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        ifstream ifs(config_file.c_str());
        if (!ifs)
        {
            cout << "Cannot open config file: " << config_file << "\n";
            exit(0);
        }
        else
        {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }

        if (vm.count("help")) {
			  cout << endl
				   << "fill in the input.cfg file" << endl;
			  exit(0);
        }
        if (vm.count("sqrtS")){
            S = comTeV*1000.;
            S2 = pow(S,2);
        }
        if (vm.count("fixed_order")){
          if(order == -1){
  			   LO = false;
  			   NLO = false;
  			   NNLO = false;
  		   }
            if(order == 0){
			   LO = true;
			   NLO = false;
			   NNLO = false;
		   }
		   else if(order == 1){
			   LO = true;
			   NLO = true;
			   NNLO = false;
		   }
		   else if(order == 2){
			   LO = true;
			   NLO = true;
			   NNLO = true;
        }
		}
        if( vm.count("resum")){
		   if(resummed == 0){
			   RES = false;
			   SCET = false;
		   }
		   else if(resummed == 1){
			   RES = true;
			   SCET = false;
		   }
		   else if(resummed == 2){
			   RES = false;
			   SCET = true;
		   }
		   else if(resummed == 3){
			   RES = true;
			   SCET = true;
		   }
	   }
	    if( vm.count("resum_order")){
		   if(resum_order.compare("LL") == 0){
			   ISLL = 1;
			   ISNLL = 0;
			   ISNNLL = 0;
		   }
		   else if(resum_order.compare("NLL")==0){
			   ISLL = 1;
			   ISNLL = 1;
			   ISNNLL = 0;
		   }
		   else if(resum_order.compare("NNLL")==0){
			   ISLL = 1;
			   ISNLL = 1;
			   ISNNLL = 1;
		   }

	   }
	    if( vm.count("process")){
		   if(process.compare("DY") == 0){
			   DY = true;
			   higgs = false, WW = false, ZZ = false, hh = false;
		   }
		   else if(process.compare("higgs") == 0){
			   higgs = true;
			   DY = false, WW = false, ZZ = false, hh = false;
		   }
		   else if(process.compare("dihiggs") == 0){
			   hh = true;
			   DY = false, WW = false, ZZ = false, higgs = false;
		   }
		   else if(process.compare("wwfull") == 0){
			   WW = true;
			   full = true;
			   diff = false;
			   DY = false, hh = false, ZZ = false, higgs = false;
		   }
		   else if(process.compare("zzfull") == 0){
			   ZZ = true;
			   full = true;
			   diff = false;
			   DY = false, hh = false, WW = false, higgs = false;
		   }
		   else if(process.compare("wwdiff") == 0){
			   WW = true;
			   full = false;
			   diff = true;
			   DY = false, hh = false, ZZ = false, higgs = false;
		   }
		   else if(process.compare("zzdiff") == 0){
			   ZZ = true;
			   full = false;
			   diff = true;
			   DY = false, hh = false, WW = false, higgs = false;
		   }
       else if(process.compare("ttH") == 0){
         ttH = true;
  			 ZZ = false;
  			 full = false;
  			 diff = false;
  			 DY = false, hh = false, WW = false, higgs = false;
       }
		   else
			{
				cout << "Option not available, choose from: DY, higgs, dihiggs, zzfull, zzdiff, wwfull, wwdiff" << endl;
				exit(0);
			}
	   }
        if( vm.count("MZ")){mZ2 = pow(mZ,2);}
        if( vm.count("MW")){mW2 = pow(mW,2);}
        if( vm.count("mh0")){mH2 = pow(mH,2);}
        if( vm.count("mt")){mt2 = pow(mt,2);}
        if( vm.count("mb")){mb2 = pow(mb,2);}
        if( vm.count("setSUSY")){if(SUSYset==true) set_SUSY(sqrt(mA2), tb);}
        if (runupdate){update_defaults();}

    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
    return 0;
}
