/*
 * PredictorTrainer.h
 *
 *  Created on: Jul 31, 2015
 *      Author: lunt
 */

#ifndef SRC_PREDICTORTRAINER_H_
#define SRC_PREDICTORTRAINER_H_

#include <string>

using namespace std;

enum ObjType
{
    SSE,                                          // sum of squared error
    CORR,                                         // Pearson correlation
    CROSS_CORR,                                   // cross correlation (maximum in a range of shifts)
    PGP,                                           // PGP score
    LOGISTIC_REGRESSION,                            // Logistic Regression
    PEAK_WEIGHTED,                                  // SSE with equal weight to peaks and non-peaks
    WEIGHTED_SSE                                  //User provides weights for sse.
};

ObjType getObjOption( const string& objOptionStr );
string getObjOptionStr( ObjType objOption );

enum SearchType
{
    UNCONSTRAINED,                                // unconstrained search
    CONSTRAINED                                   // constrained search
};

string getSearchOptionStr( SearchType searchOption );


class TrainingAware {
	protected:
		bool in_training;
		int epoch_number;
		int batch_number;

		virtual void training_updated(bool new_value, bool old_value){}
		virtual void epoch_begun(int new_epoch){}
		virtual void batch_begun(int new_batch){}
	public:
		TrainingAware() : in_training(false), epoch_number(0), batch_number(0) {}

		inline bool is_training() const { return in_training; }

		void set_training(bool do_training){bool old_train = in_training; this->in_training = do_training; this->training_updated(this->in_training, old_train);}
		void start_training(){this->set_training(true);}
		void end_training(){this->set_training(false);}
		void begin_epoch(int epoch){this->epoch_number = epoch; this->batch_number = 0; this->epoch_begun(this->epoch_number);};
		void begin_epoch(){this->begin_epoch(this->epoch_number+1);}
		void begin_batch(int batch){this->batch_number = batch; this->batch_begun(this->batch_number);};
		void begin_batch(){this->begin_batch(this->batch_number+1);}

};


#include "SeqAnnotator.h"


#endif /* SRC_PREDICTORTRAINER_H_ */
