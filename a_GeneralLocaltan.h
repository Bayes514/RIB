#ifndef A_GENERALLOCALTAN_H
#define A_GENERALLOCALTAN_H

#pragma once
#include <limits.h>

#include "incrementalLearner.h"
#include "distributionTree.h"
#include "xxyDist.h"
#include "xxxyDist.h"
#include "yDist.h"
#include <limits>
class a_GeneralLocaltan: public IncrementalLearner
{
    public:
        a_GeneralLocaltan();
        a_GeneralLocaltan(char* const *& argv, char* const * end);
        virtual ~a_GeneralLocaltan();
        void reset(InstanceStream &is);
        void initialisePass();
        void train(const instance &inst);
        void finalisePass();
        bool trainingIsFinished();
        void getCapabilities(capabilities &c);
        virtual void classify(const instance &inst, std::vector<double> &classDist);
        std::vector<CategoricalAttribute> parents_;
        std::vector<CategoricalAttribute> parentsLocal_;
        unsigned int noCatAtts_;
        unsigned int noClasses_;
    protected:

    private:
//        unsigned int noCatAtts_;
//        unsigned int noClasses_;
        double t_;
        InstanceStream* instanceStream_;
//        std::vector<CategoricalAttribute> parents_;
//        std::vector<CategoricalAttribute> parentsLocal_;
        xxyDist xxyDist_;

        bool trainingIsFinished_;
        const static CategoricalAttribute NOPARENT = 0xFFFFFFFFUL;
};

#endif // A_GENERALLOCALTAN_H
