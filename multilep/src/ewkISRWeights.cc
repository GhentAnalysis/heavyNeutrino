#include "heavyNeutrino/multilep/interface/ewkISRWeights.h"
#include "heavyNeutrino/multilep/interface/RangedMap.h"


double ewkISRWeight2016( const double ptISR ){

    //2016 ISR weights taken from https://indico.cern.ch/event/616816/contributions/2489809/attachments/1418579/2174166/17-02-22_ana_isr_ewk.pdf
    static const RangedMap< double > ewkWeights2016( { { 0., 1. }, { 50., 1.052 }, { 100., 1.179 }, { 150., 1.150 }, { 200., 1.057 }, { 300., 1.000 }, { 400., 0.912 }, { 600., 0.783 } } );

    return ewkWeights2016[ ptISR ];
}
