#include "algo_top_parameters.h"
#include "algo_top.h"
#include "cordic.h"
#include <algorithm>
#include <utility>

#include "objects.h"
#include "ap_fixed.h"
using namespace std;

TowersInEta unpackInputLink(ap_uint<576> &link) {
#pragma HLS INLINE

TowersInEta tEta_;

  tEta_.towers[0]  = Tower(link( 31,   0));
  tEta_.towers[1]  = Tower(link( 63,  32));
  tEta_.towers[2]  = Tower(link( 95,  64));
  tEta_.towers[3]  = Tower(link(127,  96));
  tEta_.towers[4]  = Tower(link(159, 128));
  tEta_.towers[5]  = Tower(link(191, 160));
  tEta_.towers[6]  = Tower(link(223, 192));
  tEta_.towers[7]  = Tower(link(255, 224));
  tEta_.towers[8]  = Tower(link(287, 256));
  tEta_.towers[9]  = Tower(link(319, 288));
  tEta_.towers[10] = Tower(link(351, 320));
  tEta_.towers[11] = Tower(link(383, 352));
  tEta_.towers[12] = Tower(link(415, 384));
  tEta_.towers[13] = Tower(link(447, 416));
  tEta_.towers[14] = Tower(link(479, 448));
  tEta_.towers[15] = Tower(link(511, 480));
  tEta_.towers[16] = Tower(link(543, 512));

  return tEta_;
}

void OutputLinkPack(ap_uint<16> Exs[24], ap_uint<16> Eys[24], ap_uint<576> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=Exs complete dim=0
#pragma HLS ARRAY_PARTITION variable=Eys complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS LATENCY min=9

link_out[0] = 0;
size_t start = 0;
link0OutputLoop: for (size_t tower = 0; tower < 24; tower++) {
                  #pragma HLS UNROLL
                  size_t end = start + 15;
                  link_out[0].range(end, start) = Exs[tower];
                  start += 16;
}

link_out[1] = 0;
start = 0;
link1OutputLoop: for (size_t tower = 0; tower < 24; tower++) {
                  #pragma HLS UNROLL
                  size_t end = start + 15;
                  link_out[1].range(end, start) = Eys[tower];
                  start += 16;
}
}

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]) {
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return

// Step 1: Unpack links
// Input is 32 links carrying 32phix17eta towers
TowersInEta towersInEta[TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=towersInEta complete dim=0


for (size_t ilink = 0; ilink < N_INPUT_LINKS; ilink++) {
  #pragma LOOP UNROLL
  size_t iEta = ilink;
  towersInEta[ilink] = unpackInputLink(link_in[iEta]);

}

   // Step 2: MET Algo goes here
ap_uint<16> Exs[24];
ap_uint<16> Eys[24];
#pragma HLS ARRAY_PARTITION variable=Exs complete dim=0
#pragma HLS ARRAY_PARTITION variable=Eys complete dim=0

COS_SIN_TYPE sinphi[angle];
COS_SIN_TYPE cosphi[angle];
#pragma HLS ARRAY_PARTITION variable=sinphi complete dim=0
#pragma HLS ARRAY_PARTITION variable=cosphi complete dim=0

/*---sinphi and cosphi calculation through CORDIC---*/

for(ap_uint<8> degree=0; degree < angle; degree++){
  #pragma HLS unroll
  THETA_TYPE theta = (degree+2.5)*0.01745329;  /*---degree to radian conversion---*/
  cordic(theta, sinphi[degree], cosphi[degree]);
}

for(ap_uint<6> b = 0; b < 24; b++) {
  #pragma hls unroll
  ap_fixed<16,12> Ey;
  ap_fixed<16,12> Ex;
  ap_uint<16> j;
  ap_uint<8> phi;

 	       if (b<36)
 	  		 phi=180*(b+4)/36; /*calulating theta of the respective tower*/
 	       else
 	         phi= -180*(72-(b+4))/36;

  j = towersInEta[b+4].towers[0].tower_et()  + towersInEta[b+4].towers[1].tower_et()  + towersInEta[b+4].towers[2].tower_et()  + towersInEta[b+4].towers[3].tower_et()  +
      towersInEta[b+4].towers[4].tower_et()  + towersInEta[b+4].towers[5].tower_et()  + towersInEta[b+4].towers[6].tower_et()  + towersInEta[b+4].towers[7].tower_et()  +
      towersInEta[b+4].towers[8].tower_et()  + towersInEta[b+4].towers[9].tower_et()  + towersInEta[b+4].towers[10].tower_et() + towersInEta[b+4].towers[11].tower_et() +
      towersInEta[b+4].towers[12].tower_et() + towersInEta[b+4].towers[13].tower_et() + towersInEta[b+4].towers[14].tower_et() + towersInEta[b+4].towers[15].tower_et() +
      towersInEta[b+4].towers[16].tower_et();



Ey = sinphi[phi]*j;
Eys[b] = Ey;
Ex = cosphi[phi]*j;
Exs[b] = Ex;

}

// Step 3: Pack the outputs

OutputLinkPack(Exs, Eys, link_out);

}
