#include "GaugeField.h"
#include "EigenVector.h"

int main() {
  for (size_t t = 0; t < 24; ++t){
    EigenVector V_t(1,dim_row,48);
    V_t.read_eigen_vector("",0,0);
    // Read in timeslice of gauge
    GaugeField gauge = GaugeField(24,12,12,12,"",0,1,4);
    gauge.read_gauge_field(400,t,t);
    // Calculate Umu_V
    umu_v = gauge.Umu_times_V(v,t,0);
    // Calculate Umu_V_shifted
    umu_vshift = gauge.Umu_times_V(v,t,0);
    // Write them to disk
    EigenVector UV(1,dim_row,48);
    UV.set_V(umu_v;)
    EigenVector UVshift(1,dim_row,48);
    UVshift.set_V(umu_vshift);
  }
}
