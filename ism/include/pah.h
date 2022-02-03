/**
 *  @file ism/include/pah.h
 *
 *  @brief Contains the template vectors [lambda, F(lambda)] of PAH emission
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __PAH_H__
#define __PAH_H__

// external includes
#include <vector>

namespace sed {

  namespace pah {

    static const std::vector< double > lpah =
      { 1.0000000000e+00, 1.9065787731e+00, 3.6350426182e+00, 6.9304950952e+00, 1.3213534836e+01, 
	2.5192645036e+01, 4.8031762264e+01, 9.1576338369e+01, 1.7459750286e+02, 3.3288389279e+02, 
	6.3466936390e+02, 1.2100471372e+03, 2.3070501862e+03, 4.3985729136e+03, 8.3862257491e+03, 
	1.5989000000e+04, 1.5990000000e+04, 1.6330000000e+04, 1.7060000000e+04, 1.7580000000e+04, 
	1.8250000000e+04, 1.9050000000e+04, 2.0400000000e+04, 2.4340000000e+04, 2.6660000000e+04, 
	3.0390000000e+04, 3.2360000000e+04, 3.2760000000e+04, 3.2790000000e+04, 3.2920000000e+04, 
	3.3080000000e+04, 3.3310000000e+04, 3.3310000000e+04, 3.3630000000e+04, 3.3770000000e+04, 
	3.9090000000e+04, 4.3000000000e+04, 4.5850000000e+04, 5.3310000000e+04, 5.5770000000e+04, 
	5.7750000000e+04, 5.8830000000e+04, 5.9930000000e+04, 6.0300000000e+04, 6.1800000000e+04, 
	6.2180000000e+04, 6.2480000000e+04, 6.4000000000e+04, 6.4740000000e+04, 6.5320000000e+04, 
	6.5680000000e+04, 6.6410000000e+04, 6.8670000000e+04, 7.0380000000e+04, 7.1650000000e+04, 
	7.3890000000e+04, 7.4750000000e+04, 7.5840000000e+04, 7.6100000000e+04, 7.7780000000e+04, 
	7.9560000000e+04, 8.0320000000e+04, 8.1320000000e+04, 8.4790000000e+04, 8.5720000000e+04, 
	8.6720000000e+04, 8.8950000000e+04, 8.9990000000e+04, 9.0920000000e+04, 9.2240000000e+04, 
	9.3440000000e+04, 9.4670000000e+04, 9.6760000000e+04, 9.9930000000e+04, 1.0180000000e+05, 
	1.0280000000e+05, 1.0360000000e+05, 1.0360000000e+05, 1.0450000000e+05, 1.0460000000e+05, 
	1.0580000000e+05, 1.0630000000e+05, 1.0720000000e+05, 1.0770000000e+05, 1.0780000000e+05, 
	1.0960000000e+05, 1.1040000000e+05, 1.1090000000e+05, 1.1180000000e+05, 1.1230000000e+05, 
	1.1260000000e+05, 1.1330000000e+05, 1.1490000000e+05, 1.1520000000e+05, 1.1640000000e+05, 
	1.1730000000e+05, 1.1870000000e+05, 1.2090000000e+05, 1.2380000000e+05, 1.2570000000e+05, 
	1.2680000000e+05, 1.2850000000e+05, 1.2910000000e+05, 1.3020000000e+05, 1.3090000000e+05, 
	1.3260000000e+05, 1.3430000000e+05, 1.3610000000e+05, 1.3630000000e+05, 1.3880000000e+05, 
	1.4030000000e+05, 1.4100000000e+05, 1.4120000000e+05, 1.4250000000e+05, 1.4310000000e+05, 
	1.4580000000e+05, 1.4760000000e+05, 1.4940000000e+05, 1.5090000000e+05, 1.5320000000e+05, 
	1.5410000000e+05, 1.5460000000e+05, 1.5630000000e+05, 1.5750000000e+05, 1.5730000000e+05, 
	1.5790000000e+05, 1.5950000000e+05, 1.5990000000e+05, 1.6160000000e+05, 1.6270000000e+05, 
	1.6340000000e+05, 1.6330000000e+05, 1.6450000000e+05, 1.6430000000e+05, 1.6580000000e+05, 
	1.6770000000e+05, 1.6900000000e+05, 1.7100000000e+05, 1.7310000000e+05, 1.7750000000e+05, 
	1.8750000000e+05, 2.0480000000e+05, 2.3810000000e+05, 2.6440000000e+05, 2.7840000000e+05, 
	3.0080000000e+05, 3.2360000000e+05, 3.5760000000e+05, 3.8850000000e+05, 4.1610000000e+05, 
	4.5200000000e+05, 4.8900000000e+05, 5.1790000000e+05, 5.4900000000e+05, 5.7470000000e+05, 
	6.1040000000e+05, 6.3130000000e+05, 7.0240000000e+05, 7.3840000000e+05, 7.7200000000e+05, 
	8.0770000000e+05, 8.4790000000e+05, 8.4790100000e+05, 1.5841226525e+06, 2.9595962007e+06, 
	5.5293759340e+06, 1.0330462721e+07, 1.9300272091e+07, 3.6058452837e+07, 6.7367548751e+07, 
	1.2586193438e+08, 2.3514625096e+08, 4.3932074947e+08, 8.2077736781e+08, 1.5334479155e+09, 
	2.8649212342e+09, 5.3524958984e+09, 1.0000000000e+10 };
    
    static const std::vector< double > fpah =
      { 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 6.9967034404e-08, 7.9230695174e-08, 9.7027192789e-08, 1.1191633162e-07, 
	1.3270688844e-07, 1.5845044618e-07, 2.0365424474e-07, 3.2426006578e-07, 3.7922194946e-07, 
	4.4349983895e-07, 6.5447560985e-07, 8.2204097648e-07, 1.2586166321e-06, 1.9719389552e-06, 
	2.6115206632e-06, 1.8918794789e-06, 1.2879335796e-06, 9.3517624491e-07, 6.7281159365e-07, 
	4.5383027728e-07, 4.4349983895e-07, 4.4760349244e-07, 5.6871340689e-07, 7.3095975521e-07, 
	1.0277641541e-06, 1.4187094660e-06, 2.0506590853e-06, 2.5579583363e-06, 6.2646017356e-06, 
	8.5682763746e-06, 1.1558286797e-05, 9.0969009598e-06, 7.0289987562e-06, 5.7133846871e-06, 
	5.2347201183e-06, 4.8740891280e-06, 4.6440134242e-06, 4.9078746998e-06, 5.5067259845e-06, 
	1.0420617776e-05, 1.3179334065e-05, 1.6553638487e-05, 1.8702233136e-05, 1.8488150443e-05, 
	1.5236789456e-05, 1.1992050990e-05, 9.7699753159e-06, 8.9308638051e-06, 9.3302539923e-06, 
	8.3155976237e-06, 5.6220339836e-06, 4.6977885653e-06, 3.8273084696e-06, 3.3181305299e-06, 
	2.9709371639e-06, 2.9369291274e-06, 2.7032951576e-06, 2.3653394026e-06, 2.1721682751e-06, 
	2.0553863427e-06, 2.1671724267e-06, 2.4316074564e-06, 2.5345068333e-06, 2.8306975334e-06, 
	2.5345068333e-06, 2.1276170721e-06, 2.2433315614e-06, 2.6417606473e-06, 2.8766926682e-06, 
	3.2726040723e-06, 3.9800952583e-06, 4.8294032103e-06, 6.7436258456e-06, 7.9779896116e-06, 
	8.7275725256e-06, 9.3087950036e-06, 8.5682763746e-06, 6.9967034404e-06, 5.5448967855e-06, 
	4.7521863931e-06, 4.6014368533e-06, 4.9418944619e-06, 5.5194202990e-06, 5.9277989633e-06, 
	6.1644376810e-06, 5.7529879778e-06, 5.2347201183e-06, 4.5592506266e-06, 3.8627221188e-06, 
	3.4189556438e-06, 2.9100031994e-06, 2.4428312465e-06, 2.2797827669e-06, 2.3328857293e-06, 
	2.3927287420e-06, 2.2850382090e-06, 2.1522537167e-06, 1.9810410118e-06, 1.8025755103e-06, 
	1.7656046603e-06, 1.8234483510e-06, 1.7615438822e-06, 1.6900263131e-06, 1.6251500534e-06, 
	1.7174870316e-06, 1.8616304180e-06, 1.6139625992e-06, 1.7656046603e-06, 1.9856077830e-06, 
	2.1821945483e-06, 2.0553863427e-06, 1.8918794789e-06, 1.7534503233e-06, 1.8488150443e-06, 
	2.0225229878e-06, 2.2433315614e-06, 2.4484625541e-06, 2.6846858047e-06, 2.7918589072e-06, 
	2.5697653355e-06, 2.1973207983e-06, 1.8530770032e-06, 1.6745320587e-06, 1.5735968098e-06, 
	1.5166782768e-06, 1.5027732868e-06, 1.4787456036e-06, 1.4089431387e-06, 1.3548570217e-06, 
	1.2470775654e-06, 1.1399701157e-06, 9.8603757845e-07, 8.5485698858e-07, 7.4454931167e-07, 
	6.2790431340e-07, 5.1510227319e-07, 4.5383027728e-07, 3.9526964753e-07, 3.4665180546e-07, 
	3.0053389949e-07, 2.6846858047e-07, 1.9810410118e-07, 1.7135369210e-07, 1.4753445831e-07, 
	1.2968611086e-07, 1.1038078198e-07, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 
	0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00 };

  } // endnamespace pah
  
} // endnamespace sed

#endif //__PAH_H__
