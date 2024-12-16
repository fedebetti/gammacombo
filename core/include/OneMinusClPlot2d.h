/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef OneMinusClPlot2d_h
#define OneMinusClPlot2d_h

#include <TROOT.h>
#include <TMultiGraph.h>
#include <TSystem.h>

#include "OneMinusClPlotAbs.h"
#include "Utils.h"
#include "ColorBuilder.h"
#include "ConfidenceContours.h"

class OneMinusClPlot2d : public OneMinusClPlotAbs
{
    public:

        OneMinusClPlot2d(OptParser *arg, TString name="c1", TString title="c1");

        void        addScanner(MethodAbsScan* s, int CLsType=0);
        void        addFile(TString fName);
        void        Draw();
        void        DrawFull();
        void        drawCLcontent(bool isFull=false);
        void        drawMarker(float x, float y, int color=0, int style=3, float size=2.0);
        void        drawGroup();
        void        drawSolutions();
        inline int  getNumberOfDefinedColors(){return linecolor[0].size();}
        inline void setContoursOnly(){contoursOnly = true;};
        inline void setXaxisTitle(TString s){xTitle=s;};
        inline void setYaxisTitle(TString s){yTitle=s;};

    protected:

        std::vector<TH2F*> histos;
        TString       xTitle;
        TString       yTitle;
        bool          contoursOnly;
        std::vector<std::vector<int> > linecolor;          ///< defines colors of 1 sigma lines and solutions of different scanners
        std::vector<std::vector<int> > fillcolor;          ///< defines colors of 1 sigma areas of different scanners
        std::vector<std::vector<int> > linestyle;          ///< defines the line style of 1 sigma line of different scanners
        std::vector<std::vector<int> > fillstyle;          ///< defines the fill style of
        std::vector<std::vector<int> > linewidth;          ///< defines the line width
        std::vector<std::vector<float> > filltransparency; ///< defines the fill transparency
        std::vector<int>          markerstyle;        ///< defines marker styles of the solutions of different scanners
        std::vector<float>        markersize;

    private:

        void drawLegend();
        bool hasHistoType(Utils::histogramType t);
        void makeNewPlotStyle(TString htmlColor, int ROOTColor=-1);
        void makeOneColorPlotStyle(TString htmlColor, int ROOTColor=-1);

        std::vector<Utils::histogramType>       histosType;          ///< defines if histogram is interpreted as p-value or chi2
        std::vector<ConfidenceContours*> m_contours;          ///< holds the contours for each scanner
        std::vector<bool>                m_contours_computed; ///< true if the contours were computed for that scanner by computeContours()
        TLegend*                    m_legend;            ///< pointer to the plot legend. Filled by drawLegend().
};

#endif
