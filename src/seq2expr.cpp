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
#include "Utils.h"
#include "IO.h"

#include "ExprModel.h"
#include "ExprPredictor.h"

#include "ObjFunc.h"

int main( int argc, char* argv[] )
{
    // command line processing
    string seqFile, r_seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, axis_wtFile;
    string outFile;                               // output file
    string dnase_file;
    string factor_thr_file;
    string par_out_file; // the learned parameters will get stored here
    ofstream par_out_stream; // Uninitialized at first.

    ModelType cmdline_modelOption = LOGISTIC;
    double coopDistThr = 50;
    double factorIntSigma = 50.0;                 // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 250;
    int maxContact = 1;
    bool read_factor_thresh = false;
    double eTF = 0.60;
    unsigned long initialSeed = time(0);

    double l1 = 0.0;
    double l2 = 0.0;

    bool cmdline_one_qbtm_per_crm = false;
    bool cmdline_one_beta = false;

    string lower_bound_file; ExprPar lower_bound_par; bool lower_bound_par_read = false;
    string upper_bound_file; ExprPar upper_bound_par; bool upper_bound_par_read = false;
    string free_fix_indicator_filename;
    ExprPar::one_qbtm_per_crm = false;
    ExprFunc::one_qbtm_per_crm = false;

    ExprPredictor::nAlternations = 5;
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
            ExprPredictor::objOption = getObjOption( argv[++i] );
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
            ExprPredictor::nAlternations = atoi( argv[++i] );
        else if ( !strcmp( "-ff", argv[i] ) )
            free_fix_indicator_filename = argv[++i];
        else if ( !strcmp( "-oq", argv[i] ) )
        {
	    cmdline_one_qbtm_per_crm = true;

            ExprPar::one_qbtm_per_crm = true;
            ExprFunc::one_qbtm_per_crm = true;
        }
        else if ( !strcmp( "-et", argv[i] ) ){
            eTF = atof( argv[ ++i ] );
	    read_factor_thresh = true;
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
    }

    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( cmdline_modelOption == QUENCHING || cmdline_modelOption == CHRMOD_UNLIMITED || cmdline_modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( cmdline_modelOption == QUENCHING && repressionFile.empty() ) )
    {
        cerr << "Usage: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -na nAlternations -ct coopDistThr -sigma factorIntSigma --seed RNG_SEED]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        exit( 1 );
    }

    //     bool readSites = false;     // whether read sites (if true) or read sequences

    // additional control parameters
    double gcContent = 0.5;
    FactorIntType intOption = BINARY;             // type of interaction function
    ExprPar::searchOption = CONSTRAINED;          // search option: unconstrained; constrained.

    ExprPredictor::nRandStarts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 400;
    ExprPredictor::nGradientIters = 50;

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

    //read the random sequences
    vector< Sequence > r_seqs;
    vector< string > r_seqNames;
    if( !r_seqFile.empty() ){
    	rval = readSequences( r_seqFile, r_seqs, r_seqNames );
    	ASSERT_MESSAGE( rval != RET_ERROR , "Coule not read the random sequences file.");
    }
    int r_nSeqs = r_seqs.size();

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

    // read the motifs
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    rval = readMotifs( motifFile, background, motifs, motifNames );
    ASSERT_MESSAGE( rval != RET_ERROR , "Could not read the motifs.");
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

    //initialize the energy threshold factors
    vector < double > energyThrFactors(nFactors, eTF);

    if( ! factor_thr_file.empty() )
    {
	int readFactorRet = readFactorThresholdFile(factor_thr_file, energyThrFactors, nFactors);
	ASSERT_MESSAGE( 0==readFactorRet , "Difficulty opening the factor_thr_input file.");
	read_factor_thresh = true;
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
    //site representation of the random sequences
    vector< SiteVec > r_seqSites( r_nSeqs );
    vector< int > r_seqLengths( r_nSeqs );

    if( r_seqs.size() > 0){
    /*SeqAnnotator r_ann( motifs, energyThrFactors );
        for ( int i = 0; i < r_nSeqs; i++ ) {
        //cout << "Annotated sites for CRM: " << seqNames[i] << endl;
                r_ann.annot( r_seqs[ i ], r_seqSites[ i ] );
                r_seqLengths[i] = r_seqs[i].size();
        }
    */
    }

    // read the cooperativity matrix
    int num_of_coop_pairs = 0;
    IntMatrix coopMat( nFactors, nFactors, false );
    if ( !coopFile.empty() )
    {
	int readRet = readEdgelistGraph(coopFile, factorIdxMap, coopMat, false);
	ASSERT_MESSAGE(0 == readRet, "Error reading the cooperativity file");

	//Calculate the number of cooperative pairs.
	for(int i = 0;i < nFactors;i++)
		for(int j=i;j<nFactors;j++)
			num_of_coop_pairs += coopMat.getElement(i,j);
    }

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


    FactorIntFunc* intFunc;
    if ( intOption == BINARY ) intFunc = new FactorIntFuncBinary( coopDistThr );
    else if ( intOption == GAUSSIAN ) intFunc = new FactorIntFuncGaussian( coopDistThr, factorIntSigma );
    else
    {
        cerr << "Interaction Function is invalid " << endl; exit( 1 );
    }

    //Create a new ExprModel with all of the selected options.
    //TODO: Continue here
    ExprModel expr_model( cmdline_modelOption, cmdline_one_qbtm_per_crm, motifs, intFunc, maxContact, coopMat, actIndicators, repIndicators, repressionMat, repressionDistThr);
    expr_model.shared_scaling = cmdline_one_beta;

    // read the axis wt file
    vector < int > axis_start;
    vector < int > axis_end;
    vector < double > axis_wts;

    axis_start.clear();
    axis_end.clear();
    axis_wts.clear();

    if( !axis_wtFile.empty() )
    {
	    int readRet = readAxisWeights(axis_wtFile, axis_start, axis_end, axis_wts);
	    ASSERT_MESSAGE(0 == readRet, "Error reading the axis weights (-aw) file.");
    }
    else
    {//Alternative intialization.
        axis_start.push_back( 0 );
        axis_end.push_back( condNames.size() - 1 );
        axis_wts.push_back( 100 );
    }


    //Setup a parameter factory
    ParFactory *param_factory = new ParFactory(expr_model, nSeqs);

    // read the initial parameter values
    ExprPar par_init = param_factory->create_expr_par(); //Currently, code further down expects par_init to be in PROB_SPACE.
    par_init = param_factory->changeSpace(par_init, PROB_SPACE); //This will cause the expected behaviour, but may hide underlying bugs.
                                                                //Code that needs par_init in a particular space should use an assertion, and do the space conversion itself.
    if ( !parFile.empty() ){
        try{
          par_init = param_factory->load( parFile );
        }catch (int& e){
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        }
    }

    //Load free_fix from the same format as parameter vectors!
    vector< double > tmp_vector;
    par_init.getRawPars(tmp_vector);
    int num_indicators = tmp_vector.size();
    vector <bool> indicator_bool(num_indicators, true);
    if( !free_fix_indicator_filename.empty() )
    {
        ExprPar param_ff;
        try{
          param_ff = param_factory->load( free_fix_indicator_filename );
        }catch (int& e){
          cerr << "Could not parse/read the free_fix file " << free_fix_indicator_filename << endl;
          exit(1);
        }
        vector < double > tmp_ff;
        param_ff.getRawPars(tmp_ff);
        indicator_bool.clear();
        for(vector<double>::iterator iter = tmp_ff.begin();iter != tmp_ff.end();++iter){
          double one_val = *iter;
          if( -0.00000001 < one_val && 0.0000001 > one_val){ indicator_bool.push_back(false); }
          else if (0.9999999 < one_val && 1.0000001 > one_val){ indicator_bool.push_back(true);}
          else{ ASSERT_MESSAGE(false,"Illegal value in indicator_bool file");}
        }
    }



    //Make sure that parameters use the energy thresholds that were specified at either the command-line or factor thresh file.
    if( read_factor_thresh ){
        par_init = param_factory->changeSpace(par_init, PROB_SPACE);
        ASSERT_MESSAGE(par_init.my_space == PROB_SPACE,"This should never happen: Preconditions not met for -et option. This is a programming error, and not the fault of the user. For now, you can try avoiding the -et commandline option, and contact the software maintainer.");
        par_init.energyThrFactors = energyThrFactors;
    }

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
    // CHECK POINT
    //     cout << "Sequences:" << endl;
    //     for ( int i = 0; i < seqs.size(); i++ ) cout << seqNames[i] << endl << seqs[i] << endl;
    //     cout << "Expression: " << endl << exprData << endl;
    //     cout << "Factor motifs:" << endl;
    //     for ( int i = 0; i < motifs.size(); i++ ) cout << motifNames[i] << endl << motifs[i] << endl;
    //     cout << "Factor expression:" << endl << factorExprData << endl;
    //     cout << "Cooperativity matrix:" << endl << coopMat << endl;
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
    cout << "Objective_Function = " << getObjOptionStr( ExprPredictor::objOption ) << endl;
    if ( !coopFile.empty() )
    {
        cout << "Interaction_Model = " << getIntOptionStr( intOption ) << endl;
        cout << "Interaction_Distance_Threshold = " << coopDistThr << endl;
        if ( intOption == GAUSSIAN ) cout << "Sigma = " << factorIntSigma << endl;
    }
    cout << "Search_Option = " << getSearchOptionStr( ExprPar::searchOption ) << endl;




    // create the expression predictor
    ExprPredictor* predictor = new ExprPredictor( seqs, seqSites, r_seqSites, seqLengths, r_seqLengths, exprData, motifs, factorExprData, expr_model, indicator_bool, motifNames, axis_start, axis_end, axis_wts );

    //Setup regularization
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
    if( par_out_stream){
        par.print( par_out_stream, motifNames, coopMat );
        par_out_stream.close();
    }
    cout << "Estimated values of parameters:" << endl;
    par.print( cout, motifNames, coopMat );
    cout << "Performance = " << setprecision( 5 ) << ( ( ExprPredictor::objOption == SSE || ExprPredictor::objOption == PGP ) ? predictor->getObj() : -predictor->getObj() ) << endl;

    // print the predictions
    writePredictions(outFile, *predictor, exprData, expr_condNames, true);

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
