/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef OptParser_h
#define OptParser_h

#include <iostream>
#include <stdlib.h>

#include "TRegexp.h"
#include "Utils.h"
#include "tclap/CmdLine.h"

using namespace Utils;
using namespace TCLAP;

class OptParser
{
    public:
        OptParser();
        ~OptParser();

        void bookOption(TString opt);
        void bookAllOptions();
        void bookPlottingOptions();
        void bookPluginOptions();
        void bookProbOptions();
        void bookFlowcontrolOptions();
        void parseArguments(int argc, char* argv[]);
        bool isAction(TString s);
        bool isAsimovCombiner(int id);
        bool isQuickhack(int id);

        std::vector<TString> action;
        std::vector<int>     asimov;
        std::vector<TString> asimovfile;
        bool            cacheStartingValues;
        std::vector<float>   CL;
        std::vector<int>     cls;
        std::vector<int>     color;
        std::vector<int>     combid;
        std::vector<std::vector<int> >    combmodifications; // encodes requested modifications to the combiner ID through the -c 26:+12 syntax,format is [cmbid:[+pdf1,-pdf2,...]]
        bool            compare;
        bool            confirmsols;
        bool            controlplot;
        int             coverageCorrectionID;
        int             coverageCorrectionPoint;
        bool            debug;
        int             digits;
        bool            enforcePhysRange;
        std::vector<int>     fillstyle;
        std::vector<int>     fillcolor;
        std::vector<float>   filltransparency;
        std::vector<int>     linewidth;
        std::vector<int>     linecolor;
        std::vector<int>     linestyle;
        std::vector<std::string>  hexfillcolor;
        std::vector<std::string>  hexlinecolor;
        TString         filenamechange;
        TString         filenameaddition;
        std::vector<std::vector<FixPar> >   fixParameters;
        std::vector<std::vector<StartPar> > startVals;
        std::vector<std::vector<RangePar> > physRanges;
        std::vector<std::vector<TString> >  removeRanges;
        std::vector<std::vector<TString> >  randomizeToyVars;
        bool            grid;
        TString         group;
        TString         groupPos;
        TString         hfagLabel;
        TString         hfagLabelPos;
        int             id;
        bool            importance;
        bool            info;
        bool            interactive;
        std::vector<int>     jmax;
        std::vector<int>     jmin;
        TString         jobdir;
        bool            largest;
        bool            latex;
        std::vector<TString> loadParamsFile;
        bool            lightfiles;
        int             batchstartn;
        bool            batcheos;
        bool            batchsubmit;
        TString         batchout;
        TString         batchreqs;
        int             nbatchjobs;
        int             nBBpoints;
        int             ndiv;
        int             ndivy;
        bool            nosyst;
        int             npoints1d;
        int             npoints2dx;
        int             npoints2dy;
        int             npointstoy;
        int             ncoveragetoys;
        int             nrun;
        int             ntoys;
        int             nsmooth;
        TString         parsavefile;
        bool            parevol;
        std::vector<int>     pevid;
        std::vector<int>     plot2dcl;
        TString         plotdate;
        TString         plotext;
        int             plotid;
        bool            plotlog;
        bool            plotlegend;
        float           plotlegx;
        float           plotlegy;
        float           plotlegsizex;
        float           plotlegsizey;
        TString         plotlegstyle;
        int             plotlegcols;
        bool            plotlegbox;
        float           plotlegboxx;
        float           plotlegboxy;
        float           plotgroupx;
        float           plotgroupy;
        Double_t        plotHFAGLabelPosX;
        Double_t        plotHFAGLabelPosY;
        Double_t        plotHFAGLabelScale;
        bool            plotmagnetic;
        int             plotnsigmacont;
        std::map<int, std::vector<int> > contourlabels;
        bool            plotpluginonly;
        bool            plotpulls;
        bool            plotprelim;
        float           plotoriginx;
        float           plotoriginy;
        std::vector<int>     plotsolutions;
        std::vector<int>     plotsoln;
        bool            plotunoff;
        float           plotymin;
        float           plotymax;
        bool            intprob;
        float           pluginPlotRangeMin;
        float           pluginPlotRangeMax;
        bool            probforce;
        bool            probimprove;
        TString         probScanResult;
        bool            printcor;
        float           printSolX;
        float           printSolY;
        std::vector<int>     qh;
        /*TString         queue;*/
        std::vector<std::vector<TString> > readfromfile;
        std::vector<TString> relation;
        bool            runCLs;
        TString         save;
        bool            saveAtMin;
        std::vector<float>   savenuisances1d;
        std::vector<float>   savenuisances2dx;
        std::vector<float>   savenuisances2dy;
        bool            scanforce;
        float           scanrangeMin;
        float           scanrangeMax;
        float           scanrangeyMin;
        float           scanrangeyMax;
        float           scaleerr;
        float           scalestaterr;
        bool            smooth2d;
        bool            square;
        int             teststatistic;
        std::vector<TString> title;
        TString         xtitle;
        TString         ytitle;
        TString         toyFiles;
        int             updateFreq;
        bool            usage;
        std::vector<TString> var;
        bool            verbose;

        CmdLine cmd;

    private:
        int         convertToDigitWithCheck(TString parseMe, TString usage);
        int         convertToIntWithCheck(TString parseMe, TString usage);
        void        defineOptions();
        void        parsePosition(TString parseMe, float &x, float &y, TString usage);
        void        parsePositionAndScale(TString parseMe, Double_t &x, Double_t &y, Double_t &scale, TString usage);
        bool        parseRange(TString parseMe, float &min, float &max);
        bool        parseAssignment(TString parseMe, TString &name, TString &value);
        bool        parseAssignment(TString parseMe, TString &name, float &value);
        void        parseCombinerString(TString parseMe, int& resultCmbId, std::vector<int>& resultAddDelPdf);
        std::vector<TString> availableOptions;
        std::vector<TString> bookedOptions;
};

#endif
