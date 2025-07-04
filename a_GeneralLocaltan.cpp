#include "a_GeneralLocaltan.h"
#include "utils.h"
#include "correlationMeasures.h"
#include <assert.h>
#include <math.h>
#include <iostream>
#include <set>
#include <stdlib.h>
#include <string>
#include <cstdlib>


#include <algorithm>

#include "globals.h"
a_GeneralLocaltan::a_GeneralLocaltan() :
trainingIsFinished_(false)
{
    //ctor
}
a_GeneralLocaltan::a_GeneralLocaltan(char* const *& , char* const * ):
xxyDist_(), trainingIsFinished_(false)
{
    name_ = "A_GENERALLOCALTAN";
//    while (argv != end)
//    {
//        if (*argv[0] != '-')
//        {
//            break;
//        }
//        else if(argv[0][1] == 'a')
//        {
//            const char * str = argv[0] + 2;
//            char *ptr;
//            t_ = strtod(str, &ptr);
////            std::cout <<"********"<< str <<std::endl;
////            std::cout <<"&&&&&&&&"<< t_ <<std::endl;
//        }
//        else
//        {
//            break;
//        }
//
//        name_ += argv[0];
//
//        ++argv;
//    }
}
a_GeneralLocaltan::~a_GeneralLocaltan()
{
    //dtor
}
class miCmpClass
{
public:
    miCmpClass(std::vector<float> *m)
    {
        mi = m;
    }
    bool operator() (CategoricalAttribute a, CategoricalAttribute b)
    {
        return (*mi)[a] > (*mi)[b];
    }
private:
    std::vector<float> *mi;
};
void a_GeneralLocaltan::reset(InstanceStream  &is)
{
    instanceStream_ = &is;
    const unsigned int noCatAtts = is.getNoCatAtts();
    noCatAtts_ = noCatAtts;
    noClasses_ = is.getNoClasses();

    trainingIsFinished_ = false;

    parentsLocal_.resize(noCatAtts);
    for (CategoricalAttribute a = 0; a < noCatAtts_; a++)
    {
        parentsLocal_[a] = NOPARENT;
    }
    parents_.resize(noCatAtts);
    for (CategoricalAttribute a = 0; a < noCatAtts_; a++)
    {
        parents_[a] = NOPARENT;
    }
//    dist.reset(is);
    xxyDist_.reset(is);
}
void a_GeneralLocaltan::getCapabilities(capabilities &c)
{
    c.setCatAtts(true);
}
void a_GeneralLocaltan::initialisePass()
{
    assert(trainingIsFinished_ == false);
}
void a_GeneralLocaltan::train(const instance &inst)
{
    xxyDist_.update(inst);
}
void a_GeneralLocaltan::finalisePass()
{
    assert(trainingIsFinished_ == false);



    crosstab<float> cmi = crosstab<float>(noCatAtts_);
//    getCondMutualInf(xxyDist_, cmi);//原
//    getxxi1CondMutualInf(xxyDist_, cmi);//5-3
//    getxxi12CondMutualInf(xxyDist_, cmi);//5-4      000048
//    getxx63maxiCondMutualInf(xxyDist_, cmi);//6-3 000048
//    getxx64maxiCondMutualInf(xxyDist_, cmi);//6-4  000192
//    getxxI2maxi1CondMutualInf(xxyDist_, cmi);//6-2
    getxxI2maxiCondMutualInf(xxyDist_, cmi);//6-1
//    getxxi2CondMutualInf(xxyDist_, cmi);//5-2
//    getxxiCondMutualInf(xxyDist_, cmi);//5-1
    // find the maximum spanning tree
cmi.print();
    CategoricalAttribute firstAtt = 0;

    parents_[firstAtt] = NOPARENT;

    float *maxWeight;
    CategoricalAttribute *bestSoFar;
    CategoricalAttribute topCandidate = firstAtt;
    std::set<CategoricalAttribute> available;

    safeAlloc(maxWeight, noCatAtts_);
    safeAlloc(bestSoFar, noCatAtts_);

    maxWeight[firstAtt] = -std::numeric_limits<float>::max();

    for (CategoricalAttribute a = firstAtt + 1; a < noCatAtts_; a++)
    {
//        if(cmi[firstAtt][a] > t_)//不加弧
//        {
            maxWeight[a] = cmi[firstAtt][a];
            if (cmi[firstAtt][a] > maxWeight[topCandidate])
                topCandidate = a;
            bestSoFar[a] = firstAtt;
            available.insert(a);
//        }

    }

    while (!available.empty())
    {
        const CategoricalAttribute current = topCandidate;
        parents_[current] = bestSoFar[current];
        available.erase(current);

        if (!available.empty())
        {
            topCandidate = *available.begin();
            for (std::set<CategoricalAttribute>::const_iterator it =
                    available.begin(); it != available.end(); it++)
            {
                    if (maxWeight[*it] < cmi[current][*it])
                    {
                        maxWeight[*it] = cmi[current][*it];
                        bestSoFar[*it] = current;
                    }

                    if (maxWeight[*it] > maxWeight[topCandidate] && maxWeight[*it] > 0)
                        topCandidate = *it;

            }
        }
    }

    delete[] bestSoFar;
    delete[] maxWeight;

    trainingIsFinished_ = true;//训练完成
}

bool a_GeneralLocaltan::trainingIsFinished()
{
    return trainingIsFinished_;
}
void a_GeneralLocaltan::classify(const instance &inst, std::vector<double> &classDist)
{
//        assert(trainingIsFinished_ == false);

    crosstab<float> cmi = crosstab<float>(noCatAtts_);
    getxxLocal61CondMutualInf(xxyDist_, cmi,inst);
//    getlocalCondMutualInf(xxyDist_, cmi, inst);//原

//    cmi.print();
    CategoricalAttribute firstAtt = 0;

    parentsLocal_[firstAtt] = NOPARENT;

    float *maxWeight;
    CategoricalAttribute *bestSoFar;
    CategoricalAttribute topCandidate = firstAtt;
    std::set<CategoricalAttribute> available;

    safeAlloc(maxWeight, noCatAtts_);
    safeAlloc(bestSoFar, noCatAtts_);

    maxWeight[firstAtt] = -std::numeric_limits<float>::max();

    for (CategoricalAttribute a = firstAtt + 1; a < noCatAtts_; a++)
    {
//        if(cmi[firstAtt][a] > t_)//不加弧
//        {
            maxWeight[a] = cmi[firstAtt][a];
            if (cmi[firstAtt][a] > maxWeight[topCandidate])
                topCandidate = a;
            bestSoFar[a] = firstAtt;
            available.insert(a);
//        }

    }

    while (!available.empty())
    {
        const CategoricalAttribute current = topCandidate;
        parentsLocal_[current] = bestSoFar[current];
        available.erase(current);

        if (!available.empty())
        {
            topCandidate = *available.begin();
            for (std::set<CategoricalAttribute>::const_iterator it =
                    available.begin(); it != available.end(); it++)
            {
                    if (maxWeight[*it] < cmi[current][*it])
                    {
                        maxWeight[*it] = cmi[current][*it];
                        bestSoFar[*it] = current;
                    }

                    if (maxWeight[*it] > maxWeight[topCandidate] && maxWeight[*it] > 0)
                        topCandidate = *it;

            }
        }
    }

    delete[] bestSoFar;
    delete[] maxWeight;


for (CatValue y = 0; y < noClasses_; y++)
    {
        classDist[y] = xxyDist_.xyCounts.p(y);
    }
    for (unsigned int x1 = 0; x1 < noCatAtts_; x1++)
    {

        const CategoricalAttribute parent = parents_[x1];

        if (parent == NOPARENT)
        {
            for (CatValue y = 0; y < noClasses_; y++)
            {
                classDist[y] *= xxyDist_.xyCounts.p(x1, inst.getCatVal(x1), y);

            }
        } else
        {
            for (CatValue y = 0; y < noClasses_; y++)
            {
                classDist[y] *= xxyDist_.p(x1, inst.getCatVal(x1), parent,
                        inst.getCatVal(parent), y);
            }
        }

    }
    normalise(classDist);
std::vector<double> classDistLocal;
classDistLocal.assign(noClasses_, 0);
for (CatValue y = 0; y < noClasses_; y++)
    {
        classDistLocal[y] = xxyDist_.xyCounts.p(y);
    }
    for (unsigned int x1 = 0; x1 < noCatAtts_; x1++)
    {

        const CategoricalAttribute parent = parentsLocal_[x1];

        if (parent == NOPARENT)
        {
            for (CatValue y = 0; y < noClasses_; y++)
            {
                classDistLocal[y] *= xxyDist_.xyCounts.p(x1, inst.getCatVal(x1), y);

            }
        } else
        {
            for (CatValue y = 0; y < noClasses_; y++)
            {
                classDistLocal[y] *= xxyDist_.p(x1, inst.getCatVal(x1), parent,
                        inst.getCatVal(parent), y);
            }
        }

    }


    for(int i=0;i<noCatAtts_;i++)
            if(parentsLocal_[i]!=NOPARENT)
            std::cout<<i<<" --> "<<parentsLocal_[i]<<std::endl;
    normalise(classDistLocal);
    for (int classno = 0; classno < noClasses_; classno++)
        {
            classDist[classno] += classDistLocal[classno];
            classDist[classno] = classDist[classno] / 2;
        }
        for(int i=0;i<noClasses_;i++)
        std::cout<<i<<" -- "<<classDistLocal[i]<<std::endl;

}















