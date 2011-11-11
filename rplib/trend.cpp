#include "rplib/trend.h"

Trend::Trend() {
}

Trend::~Trend() {

}


ConstantTrend::ConstantTrend(double trend) : trend_(trend) {

}

ConstantTrend::~ConstantTrend() {

}


Trend1D::Trend1D(const std::vector<double>& trend) : trend_(trend) {

}

Trend1D::~Trend1D() {

}

Trend2D::Trend2D(const std::vector<double>& trend, int ns1) : trend_(trend), ns1_(ns1) {

}

Trend2D::~Trend2D() {

}

