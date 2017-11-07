/*****************************************************
 * Train and test the expression model
 * Input: sequence, expression, motif, factor expr, cooperativity rules,
 *   activator list, repression rules
 * Output: trained model and expression of training data
 * File formats:
 * (1) Sequence: Fasta format
 * (2) Expression: one line per sequence
 *       <seq_name expr_1 ... expr_m>
 * (3) Motif: Stubb format
 * (4) Factor expression: one line per factor
 *       <factor expr_1 ... expr_m>
 * (5) Cooperativity: list of cooperative factor pairs
 * (6) Factor roles: the role field is either 1 (yes) or 0 (no)
 *       <factor activator_role repressor_role>
 * (7) Repression: list of <R A> where R represses A
 * (8) Parameters: the format is:
 *       <factor binding activation repression>
 *       <factor1 factor2 coop>
 *       <basal_transcription = x>
 *     where repression is optional, and the coop. lines are optional.
 * Note that (5), (6), (7) and (8) may be empty
 ******************************************************/
#include <stdexcept>

#include "utils/gs_errors.h"
#include "IO.h"

#include "ExprModel.h"
#include "ExprPredictor.h"

#include "ObjFunc.h"

#include <stdexcept>

int main( int argc, char* argv[] )
{
    // command line processing
    string seqFile, r_seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, axis_wtFile;
    string outFile;                               // output file
    string dnase_file;
    string factor_thr_file;
    string par_out_file; // the learned parameters will get stored here
    ofstream par_out_stream; // Uninitialized at first.

    string train_weights_filename;
    bool train_weights_loaded = false;

    ModelType cmdline_modelOption = LOGISTIC;
    string cmdline_interaction_option_str = "BINARY";
    double coopDistThr = 50;
    double factorIntSigma = 50.0;                 // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 0;
    int maxContact = 1;
    bool read_factor_thresh_file = false;
    bool read_factor_thresh_eTF = false;
    double eTF = 0.60;
    unsigned long initialSeed = time(0);

    bool read_par_init_file = false;

    ObjType cmdline_obj_option = SSE;
    double l1 = 0.0;
    double l2 = 0.0;

    bool cmdline_one_qbtm_per_crm = false;
    bool cmdline_one_beta = false;
    bool cmdline_write_gt = true;

    string lower_bound_file; ExprPar lower_bound_par; bool lower_bound_par_read = false;
    string upper_bound_file; ExprPar upper_bound_par; bool upper_bound_par_read = false;
    string free_fix_indicator_filename;
    //ExprPar::one_qbtm_per_crm = false;
    //ExprFunc::one_qbtm_per_crm = false;

    // additional control parameters
    double gcContent = 0.5;
    FactorIntType intOption = BINARY;             // type of interaction function
    SearchType cmdline_search_option = CONSTRAINED;//TODO:: actually read this from the commandline

    int cmdline_n_alternations = 5;
    int cmdline_n_random_starts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    int cmdline_max_simplex_iterations = 400;
    int cmdline_max_gradient_iterations = 50;
    for ( int i = 1; i < argc; i++ )
    {
        if ( !strcmp( "-s", argv[ i ] ) )
            seqFile = argv[ ++i ];
        else if ( !strcmp( "-rs", argv[ i ] ) )
            r_seqFile = argv[ ++i ];
        else if ( !strcmp( "-a", argv[ i ] ) )
            annFile = argv[ ++i ];
        else if ( !strcmp( "-e", argv[ i ] ) )
            exprFile = argv[ ++i ];
        else if ( !strcmp( "-m", argv[ i ] ) )
            motifFile = argv[ ++i ];
        else if ( !strcmp( "-f", argv[ i ] ) )
            factorExprFile = argv[ ++i ];
        else if ( !strcmp( "-o", argv[ i ] ) )
            cmdline_modelOption = getModelOption( argv[++i] );
        else if ( !strcmp( "-c", argv[ i ] ) )
            coopFile = argv[ ++i ];
        else if ( !strcmp( "-i", argv[ i ] ) )
            factorInfoFile = argv[ ++i ];
        else if ( !strcmp( "-r", argv[ i ] ) )
            repressionFile = argv[ ++i ];
        else if ( !strcmp( "-oo", argv[ i ] ) )
            cmdline_obj_option = getObjOption( argv[++i] );
        else if ( !strcmp( "-mc", argv[i] ) )
            maxContact = atoi( argv[++i] );
        else if ( !strcmp( "-fo", argv[i] ) )
            outFile = argv[++i];
        else if ( !strcmp( "-p", argv[i] ) )
            parFile = argv[++i];
        else if ( !strcmp( "-wt", argv[ i ]) )
            axis_wtFile = argv[ ++ i ];
        else if ( !strcmp( "-ct", argv[i] ) )
            coopDistThr = atof( argv[++i] );
        else if ( !strcmp( "-sigma", argv[i] ) )
            factorIntSigma = atof( argv[++i] );
        else if ( !strcmp( "-rt", argv[i] ) )
            repressionDistThr = atof( argv[++i] );
        else if ( !strcmp( "-na", argv[i] ) )
            cmdline_n_alternations = atoi( argv[++i] );
        else if ( !strcmp( "-ff", argv[i] ) )
            free_fix_indicator_filename = argv[++i];
        else if ( !strcmp( "-oq", argv[i] ) )
        {
            cmdline_one_qbtm_per_crm = true;

            //ExprPar::one_qbtm_per_crm = true;
            //ExprFunc::one_qbtm_per_crm = true;
        }
        else if ( !strcmp( "-et", argv[i] ) ){
            eTF = atof( argv[ ++i ] );
	           read_factor_thresh_eTF = true;
        }
        else if ( !strcmp( "-df", argv[ i ]))
            dnase_file = argv[ ++i ];
        else if ( !strcmp( "-ft", argv[ i ]))
            factor_thr_file = argv[ ++i ];
	else if ( !strcmp( "--seed", argv[ i ]))
	    initialSeed = atol( argv[++i] );
	else if ( !strcmp("-po", argv[ i ]))
	    par_out_file = argv[ ++i ]; //output file for pars at the en
  else if ( !strcmp("-onebeta", argv[ i ]))
      cmdline_one_beta = true;
  else if ( !strcmp("-l1", argv[ i ]))
      l1 = atof(argv[ ++i ]);
  else if ( !strcmp("-l2", argv[ i ]))
      l2 = atof(argv[ ++i ]);
	else if ( !strcmp("-lower_bound", argv[ i ]))
	    lower_bound_file = argv[ ++i ];
  else if ( !strcmp("-upper_bound", argv[ i ]))
	    upper_bound_file = argv[ ++i ];
        else if( !strcmp("-no_gt_out", argv[ i ]))
            cmdline_write_gt = false;
        else if( !strcmp("-int", argv[ i ]))
            cmdline_interaction_option_str = argv[ ++i ];
    else if( !strcmp("-train_weights", argv[ i ]))
        train_weights_filename = argv[ ++i ];
    }

    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( cmdline_modelOption == QUENCHING || cmdline_modelOption == CHRMOD_UNLIMITED || cmdline_modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( cmdline_modelOption == QUENCHING && repressionFile.empty() ) )
    {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -na nAlternations -ct coopDistThr -sigma factorIntSigma --seed RNG_SEED]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }

    //     bool readSites = false;     // whether read sites (if true) or read sequences

    int rval;
    vector< vector< double > > data;              // buffer for reading matrix data
    vector< string > labels;                      // buffer for reading the labels of matrix data
    string factor1, factor2;                      // buffer for reading factor pairs

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    ASSERT_MESSAGE(rval != RET_ERROR, "Could not read the sequence file.");
    int nSeqs = seqs.size();

    /*
    //read the random sequences
    vector< Sequence > r_seqs;
    vector< string > r_seqNames;
    if( !r_seqFile.empty() ){
    	rval = readSequences( r_seqFile, r_seqs, r_seqNames );
    	ASSERT_MESSAGE( rval != RET_ERROR , "Coule not read the random sequences file.");
    }
    int r_nSeqs = r_seqs.size();
    */

    // read the expression data
    vector< string > condNames;
    rval = readMatrix( exprFile, labels, condNames, data );
    ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the expression data file.");
    ASSERT_MESSAGE( labels.size() == nSeqs , "Mismatch between number of labels and number of sequences");
    for ( int i = 0; i < nSeqs; i++ )
    {
        if( labels[ i ] != seqNames[ i ] ) cout << labels[i] << seqNames[i] << endl;
        ASSERT_MESSAGE( labels[i] == seqNames[i] , "A label and a sequence name did not agree.");
    }
    Matrix exprData( data );
    int nConds = exprData.nCols();
    vector < string > expr_condNames = condNames;

    //read the weights if provided
    Matrix *training_weights = NULL;
    if( ! train_weights_filename.empty() ){
        vector<vector< double> > weights_data(0);
        vector<string> weights_labels(0);
        vector<string> weights_condNames(0);
        rval = readMatrix( train_weights_filename, weights_labels, weights_condNames, weights_data );
        ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the weights data file.");
        ASSERT_MESSAGE( labels.size() == nSeqs , "Mismatch between number of labels and number of sequences");
        for ( int i = 0; i < nSeqs; i++ )
        {
            if( weights_labels[ i ] != seqNames[ i ] ) cout << weights_labels[i] << seqNames[i] << endl;
            ASSERT_MESSAGE( weights_labels[i] == seqNames[i] , "A label and a sequence name did not agree.");
        }
        training_weights = new Matrix( weights_data );
        train_weights_loaded = true;
    }


    // read the motifs
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    try{
        rval = readMotifs( motifFile, background, motifs, motifNames );
        ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the motifs.");
    }catch( std::runtime_error e){
        cerr << "Unable to read the PWM file because of error: " << e.what() << endl;
        exit(1);
    }
    int nFactors = motifs.size();

    // factor name to index mapping
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ )
    {
        factorIdxMap[motifNames[i]] = i;
    }

    // read the factor expression data
    labels.clear();
    data.clear();
    rval = readMatrix( factorExprFile, labels, condNames, data );
    ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the factor expression matrix");
    cout << labels.size() << " " << nFactors << " " << condNames.size() << " " << nConds << endl;
    ASSERT_MESSAGE( labels.size() == nFactors, "Number of labels and number of transcriptions factors differ.");
    ASSERT_MESSAGE( condNames.size() == nConds, "Number of condition-names and number of conditions differ.");
    for ( int i = 0; i < nFactors; i++ ){ ASSERT_MESSAGE( labels[i] == motifNames[i], "A label and a motif name disagree in the factor expression file." ); }
    Matrix factorExprData( data );
    ASSERT_MESSAGE( factorExprData.nCols() == nConds , "Number of columns in factor expression data differs from the number of conditions.");

    //Initialize the dataset that is actually provided
    DataSet training_dataset(factorExprData,exprData);

    //****** MODEL ********

    // read the roles of factors
    vector< bool > actIndicators( nFactors, true );
    vector< bool > repIndicators( nFactors, false );
    if ( !factorInfoFile.empty() )
    {
	int readRet = readFactorRoleFile(factorInfoFile, factorIdxMap, actIndicators, repIndicators);
        ASSERT_MESSAGE(0 == readRet, "Could not parse the factor information file.");
    }

    // read the repression matrix
    IntMatrix repressionMat( nFactors, nFactors, false );
    if ( !repressionFile.empty() )
    {
	int readRet = readEdgelistGraph(repressionFile, factorIdxMap, repressionMat, true);
	ASSERT_MESSAGE(0 == readRet, "Error reading the repression file.");
    }


    //Create a new ExprModel with all of the selected options.
    //TODO: Continue here
    if( cmdline_modelOption == DIRECT || cmdline_modelOption == LOGISTIC || cmdline_modelOption == MARKOV || cmdline_modelOption == RATES){
        //TODO: Kind of a hacky workaround, the models/DP implementations should know that they need to ignore this during their setup.
        repressionDistThr = 0;
    }
    ExprModel expr_model( cmdline_modelOption, cmdline_one_qbtm_per_crm, motifs, motifNames, maxContact, actIndicators, repIndicators, repressionMat, repressionDistThr);
    expr_model.shared_scaling = cmdline_one_beta;

    //********* SETUP COOPERTIVITIES ********

    //TODO: break this block into a separate function
    {
        FactorIntFunc* default_int_func;
        intOption = getIntOption(cmdline_interaction_option_str);
        if ( intOption == BINARY ) default_int_func = new FactorIntFuncBinary( coopDistThr );
        else if ( intOption == GAUSSIAN ) default_int_func = new FactorIntFuncGaussian( coopDistThr, factorIntSigma );
        else if ( intOption == HELICAL ) default_int_func = new FactorIntFuncHelical( coopDistThr );
        else
        {
            cerr << "Interaction Function is invalid " << endl; exit( 1 );
        }
        expr_model.coop_setup->set_default_interaction(default_int_func);
    }

    // read the cooperativity matrix
    if ( !coopFile.empty() )
    {
        expr_model.coop_setup->read_coop_file(coopFile, factorIdxMap);
    }
    //***** END COOPS *****

    //Deleted AXIS_WEIGHTS from here

    cerr << "Created the parameter factory...";
    //Setup a parameter factory
    ParFactory *param_factory = new ParFactory(expr_model, nSeqs);//This param_factory is used for loading/unloading, it should be unconstrained.
    cerr << "DONE." << endl;

    cerr << "Creating the initial parameters...";
    // read the initial parameter values
    ExprPar par_init = param_factory->create_expr_par(); //Currently, code further down expects par_init to be in PROB_SPACE.
    cerr << " ... " << par_init.my_pars;
    par_init = param_factory->changeSpace(par_init, PROB_SPACE); //This will cause the expected behaviour, but may hide underlying bugs.
                                                                //Code that needs par_init in a particular space should use an assertion, and do the space conversion itself.
    cerr << "DONE." << endl;

    cerr << "Created the parameter factory." << endl;

    if ( !parFile.empty() ){
        try{
          par_init = param_factory->load( parFile );
          read_par_init_file = true;
        }catch (int& e){
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        }
    }

    /******** FREE fix
    */
    //Load free_fix from the same format as parameter vectors!
    ExprPar param_ff = param_factory->createDefaultFreeFix();

    if( !free_fix_indicator_filename.empty() )
    {
        //ExprPar param_ff;
        try{
          cerr << "loading free fix" << endl;
          param_ff = param_factory->load( free_fix_indicator_filename );
          cerr << "loaded free fix" << endl;
        }catch (int& e){
          cerr << "Could not parse/read the free_fix file " << free_fix_indicator_filename << endl;
          exit(1);
        }
    }
        #ifndef REANNOTATE_EACH_PREDICTION
        //prevent optimization of annotation thresholds if that will be useless.
        //param_ff.energyThrFactors.assign(param_ff.energyThrFactors.size(),0.0);
        #endif

        vector <bool> indicator_bool;
        vector < double > tmp_ff;
        param_ff.getRawPars(tmp_ff);
        indicator_bool.clear();
        for(vector<double>::iterator iter = tmp_ff.begin();iter != tmp_ff.end();++iter){
          double one_val = *iter;
          if( -0.00000001 < one_val && 0.0000001 > one_val){ indicator_bool.push_back(false); }
          else if (0.9999999 < one_val && 1.0000001 > one_val){ indicator_bool.push_back(true);}
          else{ ASSERT_MESSAGE(false,"Illegal value in indicator_bool file");}
        }




    /******* END OF FREE free_fix
    */



    /*
    //Make sure that parameters use the energy thresholds that were specified at either the command-line or factor thresh file.
    if( read_factor_thresh_eTF ){
        par_init = param_factory->changeSpace(par_init, PROB_SPACE);
        ASSERT_MESSAGE(par_init.my_space == PROB_SPACE,"This should never happen: Preconditions not met for -et option. This is a programming error, and not the fault of the user. For now, you can try avoiding the -et commandline option, and contact the software maintainer.");
        //par_init.energyThrFactors = energyThrFactors;
        //I really don't want to have code dependent on the structure of the parameters here.
        //dammit.
    }
    */

    if ( !upper_bound_file.empty() ){
	try{
		upper_bound_par = param_factory->load( upper_bound_file );
		upper_bound_par = param_factory->changeSpace(upper_bound_par, ENERGY_SPACE);
		upper_bound_par_read = true;
	}catch (int& e){
		cerr << "Cannot read upper bounds from " << upper_bound_file << endl;
		exit( 1 );
	}
    }

    if ( !lower_bound_file.empty() ){
	try{
		lower_bound_par = param_factory->load( lower_bound_file );
		lower_bound_par = param_factory->changeSpace(lower_bound_par, ENERGY_SPACE);
		lower_bound_par_read = true;
	}catch (int& e){
		cerr << "Cannot read lower bounds from " << lower_bound_file << endl;
		exit( 1 );
	}
    }

    //Check AGAIN that the indicator_bool will be the right shape for the parameters that are read.
    vector < double > all_pars_for_test;
    par_init.getRawPars(all_pars_for_test );
    ASSERT_MESSAGE(all_pars_for_test.size() == indicator_bool.size(), "For some reason, the number of entries in free_fix did not match the number of free parameters.\n"
		  "Remember that whatever model, there are 3 parameters for every transcription factor\n");
    all_pars_for_test.clear();//Won't be used again.
    //It is possible that the user wants to write out to the same par file, doing this after reading the par file means we won't have overridden it before reading

    //Check that we can access and write to the par outfile now, so that we can warn the user before a lot of time was spent on the optimization
    if( !par_out_file.empty() ){
        par_out_stream.open( par_out_file.c_str() );
        if( ! par_out_stream ){
                cerr << "Cannot open the parameter output file " << par_out_file << " for writing." << endl;
                exit(1);
        }
    }

    //initialize the energy threshold factors
    vector < double > energyThrFactors(nFactors, eTF);
    //TODO: par_init
    if(read_par_init_file){
      assert(par_init.my_space == PROB_SPACE);
      //assert(energyThrFactors.size() == par_init.energyThrFactors.size());

      //energyThrFactors = par_init.energyThrFactors;
      for(int i = 0;i<((gsparams::DictList&)par_init.my_pars)["tfs"].size();i++){
          energyThrFactors[i] = ((gsparams::DictList&)par_init.my_pars)["tfs"][i]["annot_thresh"];
      }
    }

    //TODO: move this to after the reading of the factor_thr_file? Meh. We should never use factor_thr_file anyway.
    if(read_factor_thresh_eTF){
      energyThrFactors.assign(energyThrFactors.size(),eTF);
    }


    if( ! factor_thr_file.empty() )//TODO: Totally eliminate the factor_thr_file.
    {
      int readFactorRet = readFactorThresholdFile(factor_thr_file, energyThrFactors, nFactors);
      ASSERT_MESSAGE( 0==readFactorRet , "Difficulty opening the factor_thr_input file.");
      read_factor_thresh_file = true;
    }



    //assign that back to the initial par file so that it receives any changes made.
    assert(par_init.my_space == PROB_SPACE);
    //assert(energyThrFactors.size() == par_init.energyThrFactors.size());//TODO: restore

    //New way to do this
    //par_init.energyThrFactors = energyThrFactors;

    for(int i = 0;i<((gsparams::DictList&)par_init.my_pars)["tfs"].size();i++){
        par_init.my_pars["tfs"][i]["annot_thresh"] = energyThrFactors.at(i);
    }

    // site representation of the sequences
    // TODO: Should this code be removed? If we are using this code, and no command-line option was provided for energyThrFactors, but a .par file was provided, shouldn't it use the thresholds learned there? (So, shouldn't it happen after reading the par file?)
    // TODO: Relates to issue #19
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrFactors );
    if ( annFile.empty() )                        // construct site representation
    {

        if( dnase_file.empty() )
        {
            for ( int i = 0; i < nSeqs; i++ )
            {
                //cout << "Annotated sites for CRM: " << seqNames[i] << endl;
                ann.annot( seqs[ i ], seqSites[ i ] );
                seqLengths[i] = seqs[i].size();
            }
        }
        else
        {
            for ( int i = 0; i < nSeqs; i++ )
            {
                //cout << "Annotated sites for CRM: " << seqNames[i] << endl;
                ifstream dnase_input( dnase_file.c_str() );
                assert( dnase_input.is_open());

                string temp_s;
                string temp_gen;
                string chr;
                double temp_start, temp_end;
                vector < double > dnase_start;
                vector < double > dnase_end;
                vector < double > scores;

                while( dnase_input >> temp_s )
                {
                    dnase_input >> temp_gen >> temp_gen >> temp_gen >> temp_gen >> temp_gen;
                    dnase_input >> chr;
                    dnase_input >> temp_gen >> temp_start >> temp_end >> temp_gen;
                    if( temp_s == seqNames[ i ] )
                    {
                        //cout << "Processing for:\t" << temp_s << endl;
                        dnase_start.clear();
                        dnase_end.clear();
                        scores.clear();
                        ifstream chr_input( ("chr" + chr + ".bed").c_str() );
                        //cout << "File: " << "chr" + chr + ".bed" << "opened" << endl;
                        assert( chr_input.is_open() );
                        double chr_start, chr_end, chr_score;
                        //cout << "Starting location on chromosome: " << (long long int)temp_start << endl;
                        //cout << "Ending location on chromosome: " << (long long int)temp_end << endl;
                        while( chr_input >> temp_gen >> chr_start >> chr_end >> temp_gen >> chr_score )
                        {
                            if( ( chr_start < temp_start && chr_end < temp_start ) || ( chr_start > temp_end && chr_end > temp_end ) )
                            {
                                ;
                            }
                            else
                            {
                                dnase_start.push_back( chr_start );
                                dnase_end.push_back( chr_end );
                                scores.push_back( chr_score );
                                //cout << "Inserting: " << (long long int)chr_start << "\t" << (long long int)chr_end << "\t" << chr_score << endl;
                            }
                        }
                        chr_input.close();
                        break;
                    }
                }

                dnase_input.close();
                ann.annot( seqs[ i ], seqSites[ i ], dnase_start, dnase_end, scores, temp_start );
                seqLengths[i] = seqs[i].size();
            }
        }
    }                                             // read the site representation and compute the energy of sites
    else
    {
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ )
        {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }

    //TODO: R_SEQ Either remove this feature or un-comment it.
    /*
    //site representation of the random sequences
    vector< SiteVec > r_seqSites( r_nSeqs );
    vector< int > r_seqLengths( r_nSeqs );

    if( r_seqs.size() > 0){
    SeqAnnotator r_ann( motifs, energyThrFactors );
        for ( int i = 0; i < r_nSeqs; i++ ) {
        //cout << "Annotated sites for CRM: " << seqNames[i] << endl;
                r_ann.annot( r_seqs[ i ], r_seqSites[ i ] );
                r_seqLengths[i] = r_seqs[i].size();
        }

    }
    */


    // CHECK POINT
    //     cout << "Sequences:" << endl;
    //     for ( int i = 0; i < seqs.size(); i++ ) cout << seqNames[i] << endl << seqs[i] << endl;
    //     cout << "Expression: " << endl << exprData << endl;
    //     cout << "Factor motifs:" << endl;
    //     for ( int i = 0; i < motifs.size(); i++ ) cout << motifNames[i] << endl << motifs[i] << endl;
    //     cout << "Factor expression:" << endl << factorExprData << endl;
    //     cout << "Cooperativity matrix:" << endl << expr_model.coop_setup->coop_matrix << endl;
    //     cout << "Activators:" << endl << actIndicators << endl;
    //     cout << "Repressors:" << endl << repIndicators << endl;
    //     cout << "Repression matrix:" << endl << repressionMat << endl;
    //     cout << "Site representation of sequences:" << endl;
    //     for ( int i = 0; i < nSeqs; i++ ) {
    //         cout << ">" << seqNames[i] << endl;
    //         for ( int j = 0; j < seqSites[i].size(); j++ ) cout << seqSites[i][j] << endl;
    //     }

    // print the parameters for running the analysis
    cout << "Parameters for running the program: " << endl;
    cout << "Model = " << getModelOptionStr( cmdline_modelOption ) << endl;
    if ( cmdline_modelOption == QUENCHING || cmdline_modelOption == CHRMOD_LIMITED )
    {
        cout << "Maximum_Contact = " << maxContact << endl;
    }
    if ( cmdline_modelOption == QUENCHING || cmdline_modelOption == CHRMOD_LIMITED || cmdline_modelOption == CHRMOD_UNLIMITED )
    {
        cout << "Repression_Distance_Threshold = " << repressionDistThr << endl;
    }
    cout << "Objective_Function = " << getObjOptionStr( cmdline_obj_option ) << endl;
    if ( !coopFile.empty() )
    {
        cout << "Interaction_Model = " << getIntOptionStr( intOption ) << endl;
        cout << "Interaction_Distance_Threshold = " << coopDistThr << endl;
        if ( intOption == GAUSSIAN ) cout << "Sigma = " << factorIntSigma << endl;
    }
    //cout << "Search_Option = " << getSearchOptionStr( ExprPar::searchOption ) << endl; //TODO: restore




    // create the expression predictor
    ExprPredictor* predictor = new ExprPredictor( seqs, seqSites, seqLengths, training_dataset, motifs, expr_model, indicator_bool, motifNames );
    //And setup parameters from the commandline
    predictor->search_option = cmdline_search_option;
    predictor->set_objective_option(cmdline_obj_option);
    predictor->n_alternations = cmdline_n_alternations;
    predictor->n_random_starts = cmdline_n_random_starts;
    predictor->max_simplex_iterations = cmdline_max_simplex_iterations;
    predictor->max_gradient_iterations = cmdline_max_gradient_iterations;


    //Setup a weighted objective if that is appropriate
    if( predictor->objOption == WEIGHTED_SSE) {
        ASSERT_MESSAGE( train_weights_loaded , "User requested WEIGHTED_SSE objective, but provided no weights.");
        delete predictor->trainingObjective;
        Weighted_RMSEObjFunc *tmp_ptr = new Weighted_RMSEObjFunc();
        tmp_ptr->set_weights(training_weights);
        predictor->trainingObjective = tmp_ptr;
        //TODO: log a message about this.
    }


    //Setup regularization objective function
    if(0.0 != l1 || 0.0 != l2){
      cerr << "INFO: Regularization was turned on and will be used. l1 = " << l1 << " l2 = " << l2 << " ."<< endl;

      ExprPar tmp_centers = predictor->param_factory->create_expr_par();
      ExprPar tmp_l1 = predictor->param_factory->create_expr_par();
      ExprPar tmp_l2 = predictor->param_factory->create_expr_par();

      //TODO: add an option to read l1 and l2 values from a file.
      vector< double > tmp_l12_vector;
      tmp_l1.getRawPars(tmp_l12_vector);
      std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l1);
      tmp_l1 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);

      tmp_l2.getRawPars(tmp_l12_vector);
      std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l2);
      tmp_l2 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);


      RegularizedObjFunc *tmp_reg_obj_func = new RegularizedObjFunc(predictor->trainingObjective,
                                              tmp_centers,
                                              tmp_l1,
                                              tmp_l2
                                            );
      predictor->trainingObjective = tmp_reg_obj_func;
    }

    if(upper_bound_par_read){
    	predictor->param_factory->setMaximums(upper_bound_par);
    }
    if(lower_bound_par_read){
    	predictor->param_factory->setMinimums(lower_bound_par);
    }

    // random number generator
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;     // create rng type
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, initialSeed );                // set the seed equal to simulTime(0)

    // model fitting

    predictor->train( par_init, rng );

    gsl_rng_free( rng );
    // print the training results
    ExprPar par = predictor->getPar();
    par = predictor->param_factory->changeSpace(par, PROB_SPACE);
    if( par_out_stream){
        par.print( par_out_stream, motifNames, expr_model.coop_setup->coop_matrix );
        par_out_stream.close();
    }
    cout << "Estimated values of parameters:" << endl;
    par.print( cout, motifNames, expr_model.coop_setup->coop_matrix );
    cout << "Performance = " << setprecision( 5 ) << ( ( cmdline_obj_option == SSE || cmdline_obj_option == PGP ) ? predictor->getObj() : -predictor->getObj() ) << endl;

    // print the predictions
    writePredictions(outFile, *predictor, training_dataset.exprData, expr_condNames, cmdline_write_gt, true);

    //TODO: R_SEQ Either remove this feature or make it conditional.
    /*
        cout << "Max expressions of the random sequences:" << endl;

        for( int i = 0; i < nSeqs; i++ ){
            vector < double > targetExprs;
        predictor -> predict(r_seqSites[i], r_seqLengths[i], targetExprs, i);
        double max = 0;
        for( int j = 0; j < targetExprs.size(); j++ ){
            if( targetExprs[ j ] > max  ){
                max = targetExprs[ j ];
            }
    }
    cout << i + 1 << "\t" << max << endl;
    }
    */
    return 0;
}
