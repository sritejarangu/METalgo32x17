#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "algo_top.h"
#include "cordic.h"


using namespace std;

int main(){

	ap_uint<576> test_in[N_INPUT_LINKS];
	ap_uint<576> test_out[N_OUTPUT_LINKS];

test_in[0]  = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[1]  = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[2]  = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[3]  = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[4]  = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[5]  = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[6]  = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[7]  = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[8]  = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[9]  = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[10] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[11] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[12] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[13] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[14] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[15] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[16] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[17] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[18] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[19] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[20] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[21] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[22] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[23] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[24] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[25] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[26] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[27] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[28] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[29] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
test_in[30] = "0x00000000001000000210D43500C0000001C0401000000000000000000030E0380230000000000000012094250000000000400000024118460231244900000000000000000130E83A";
test_in[31] = "0x000002FE01C0300C0000000000816C5B0000000002306C1B0190F03C000000000000000000000000000000000030000002310C430200000000000000000000000191084201C14050";
// Run the algorithm

algo_top(test_in, test_out);





cout << hex << "link_out[0]: " << fixed << test_out[0] << endl;
cout << hex << "link_out[1]: " << fixed << test_out[1] << endl;

return 0;

}
