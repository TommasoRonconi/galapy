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
      { 1.0000000000e-01, 1.2251826030e-01, 1.5010724108e-01, 1.8390878036e-01, 2.2532183824e-01, 
	2.7606039630e-01, 3.3822439493e-01, 4.1438664459e-01, 5.0769930788e-01, 6.2202435959e-01, 
	7.6209342403e-01, 9.3370360501e-01, 1.1439574133e+00, 1.4015567213e+00, 1.7171629121e+00, 
	2.1038381265e+00, 2.5775858722e+00, 3.1580133685e+00, 3.8691430392e+00, 4.7404067403e+00, 
	5.8078638696e+00, 7.1156937738e+00, 8.7180242202e+00, 1.0681171607e+01, 1.3086385633e+01, 
	1.6033212015e+01, 1.9643612431e+01, 2.4067012211e+01, 2.9486484668e+01, 3.6126328040e+01, 
	4.4261348627e+01, 5.4228234324e+01, 6.6439489287e+01, 8.1400506430e+01, 9.9730484356e+01, 
	1.2218805443e+02, 1.4970267858e+02, 1.8341311742e+02, 2.2471456064e+02, 2.7531637034e+02, 
	3.3731282727e+02, 4.1326980776e+02, 5.0633097882e+02, 6.2034790663e+02, 7.6003946303e+02, 
	9.3118712773e+02, 1.1408742691e+03, 1.3977793067e+03, 1.7125348895e+03, 2.0981679537e+03, 
	2.5706388751e+03, 3.1495020284e+03, 3.8587150934e+03, 4.7276306026e+03, 5.7922107678e+03, 
	7.0965158659e+03, 8.6945277810e+03, 1.0652384179e+04, 1.3051115777e+04, 1.5990000000e+04, 
	1.6330000000e+04, 1.7060000000e+04, 1.7580000000e+04, 1.8250000000e+04, 1.9050000000e+04, 
	2.0400000000e+04, 2.4340000000e+04, 2.6660000000e+04, 3.0390000000e+04, 3.2360000000e+04, 
	3.2760000000e+04, 3.2790000000e+04, 3.2920000000e+04, 3.3080000000e+04, 3.3310000000e+04, 
	3.3630000000e+04, 3.3770000000e+04, 3.9090000000e+04, 4.3000000000e+04, 4.5850000000e+04, 
	5.3310000000e+04, 5.5770000000e+04, 5.7750000000e+04, 5.8830000000e+04, 5.9930000000e+04, 
	6.0300000000e+04, 6.1800000000e+04, 6.2180000000e+04, 6.2480000000e+04, 6.4000000000e+04, 
	6.4740000000e+04, 6.5320000000e+04, 6.5680000000e+04, 6.6410000000e+04, 6.8670000000e+04, 
	7.0380000000e+04, 7.1650000000e+04, 7.3890000000e+04, 7.4750000000e+04, 7.5840000000e+04, 
	7.6100000000e+04, 7.7780000000e+04, 7.9560000000e+04, 8.0320000000e+04, 8.1320000000e+04, 
	8.4790000000e+04, 8.5720000000e+04, 8.6720000000e+04, 8.8950000000e+04, 8.9990000000e+04, 
	9.0920000000e+04, 9.2240000000e+04, 9.3440000000e+04, 9.4670000000e+04, 9.6760000000e+04, 
	9.9930000000e+04, 1.0180000000e+05, 1.0280000000e+05, 1.0360000000e+05, 1.0450000000e+05, 
	1.0460000000e+05, 1.0580000000e+05, 1.0630000000e+05, 1.0720000000e+05, 1.0770000000e+05, 
	1.0780000000e+05, 1.0960000000e+05, 1.1040000000e+05, 1.1090000000e+05, 1.1180000000e+05, 
	1.1230000000e+05, 1.1260000000e+05, 1.1330000000e+05, 1.1490000000e+05, 1.1520000000e+05, 
	1.1640000000e+05, 1.1730000000e+05, 1.1870000000e+05, 1.2090000000e+05, 1.2380000000e+05, 
	1.2570000000e+05, 1.2680000000e+05, 1.2850000000e+05, 1.2910000000e+05, 1.3020000000e+05, 
	1.3090000000e+05, 1.3260000000e+05, 1.3430000000e+05, 1.3610000000e+05, 1.3630000000e+05, 
	1.3880000000e+05, 1.4030000000e+05, 1.4100000000e+05, 1.4120000000e+05, 1.4250000000e+05, 
	1.4310000000e+05, 1.4580000000e+05, 1.4760000000e+05, 1.4940000000e+05, 1.5090000000e+05, 
	1.5320000000e+05, 1.5410000000e+05, 1.5460000000e+05, 1.5630000000e+05, 1.5730000000e+05, 
	1.5750000000e+05, 1.5790000000e+05, 1.5950000000e+05, 1.5990000000e+05, 1.6160000000e+05, 
	1.6270000000e+05, 1.6330000000e+05, 1.6340000000e+05, 1.6430000000e+05, 1.6450000000e+05, 
	1.6580000000e+05, 1.6770000000e+05, 1.6900000000e+05, 1.7100000000e+05, 1.7310000000e+05, 
	1.7750000000e+05, 1.8750000000e+05, 2.0480000000e+05, 2.3810000000e+05, 2.6440000000e+05, 
	2.7840000000e+05, 3.0080000000e+05, 3.2360000000e+05, 3.5760000000e+05, 3.8850000000e+05, 
	4.1610000000e+05, 4.5200000000e+05, 4.8900000000e+05, 5.1790000000e+05, 5.4900000000e+05, 
	5.7470000000e+05, 6.1040000000e+05, 6.3130000000e+05, 7.0240000000e+05, 7.3840000000e+05, 
	7.7200000000e+05, 8.0770000000e+05, 8.4790000000e+05, 1.0383401336e+06, 1.2715535240e+06, 
	1.5571471353e+06, 1.9068856758e+06, 2.3351762323e+06, 2.8596617537e+06, 3.5019478327e+06, 
	4.2884927238e+06, 5.2516972611e+06, 6.4312395751e+06, 7.8757095879e+06, 9.6446106211e+06, 
	1.1810810568e+07, 1.4463543606e+07, 1.7712086095e+07, 2.1690258099e+07, 2.6561935950e+07, 
	3.2527802952e+07, 3.9833616302e+07, 4.8780330785e+07, 5.9736496266e+07, 7.3153439690e+07, 
	8.9583856988e+07, 1.0970458077e+08, 1.3434446167e+08, 1.6451851195e+08, 2.0146971774e+08, 
	2.4672024252e+08, 3.0213413089e+08, 3.6999409580e+08, 4.5309555238e+08, 5.5486177189e+08, 
	6.7948489958e+08, 8.3209864536e+08, 1.0189897613e+09, 1.2478570172e+09, 1.5281283429e+09, 
	1.8713492012e+09, 2.2916581903e+09, 2.8063694675e+09, 3.4366859863e+09, 4.2085729285e+09, 
	5.1538273108e+09, 6.3113878269e+09, 7.7289388836e+09, 9.4648749061e+09, 1.1590705831e+10, 
	1.4194002879e+10, 1.7382005951e+10, 2.1286041257e+10, 2.6066931150e+10, 3.1921619027e+10, 
	3.9091282185e+10, 4.7871266855e+10, 5.8623254654e+10, 7.1790161658e+10, 8.7914383826e+10, 
	1.0766014040e+11, 1.3184083567e+11, 1.6145256625e+11, 1.9771515416e+11, 2.4212239602e+11, 
	2.9650359834e+11, 3.6309893373e+11, 4.4465172232e+11, 5.4452143974e+11, 6.6682210695e+11, 
	8.1659176273e+11, 1.0000000000e+12 };
    
    static const std::vector< double > fpah =
      { 1.1866212103e-38, 3.9404064150e-38, 1.3084885539e-37, 4.3450906210e-37, 1.4428718118e-36, 
	4.7913363538e-36, 1.5910563827e-35, 5.2834120293e-35, 1.7544596769e-34, 5.8260244345e-34, 
	1.9346446748e-33, 6.4243637489e-33, 2.1333348762e-32, 7.0841531892e-32, 2.3524307865e-31, 
	7.8117037530e-31, 2.5940280953e-30, 8.6139745845e-30, 2.8604377215e-29, 9.4986395400e-29, 
	3.1542079184e-28, 1.0474160589e-27, 3.4781486475e-27, 1.1549868756e-26, 3.8353584568e-26, 
	1.2736053372e-25, 4.2292541184e-25, 1.4044060492e-24, 4.6636033111e-24, 1.5486401426e-23, 
	5.1425606582e-23, 1.7076872409e-22, 5.6707074679e-22, 1.8830686566e-21, 6.2530955538e-21, 
	2.0764619425e-20, 6.8952955563e-20, 2.2897169380e-19, 7.6034502271e-19, 2.5248734634e-18, 
	8.3843331854e-18, 2.7841808306e-17, 9.2454137089e-17, 3.0701193583e-16, 1.0194928178e-15, 
	3.3854240970e-15, 1.1241958860e-14, 3.7331109898e-14, 1.2396520780e-13, 4.1165057206e-13, 
	1.3669657518e-12, 4.5392755250e-12, 1.5073546842e-11, 5.0054642676e-11, 1.6621617192e-10, 
	5.5195311226e-10, 1.8328676122e-09, 6.0863932263e-09, 2.0211051938e-08, 6.7114727108e-08, 
	7.6000741356e-08, 9.3071739021e-08, 1.0735390058e-07, 1.2729690031e-07, 1.5199098470e-07, 
	1.9535198507e-07, 3.1104113547e-07, 3.6376241852e-07, 4.2541992693e-07, 6.2779496556e-07, 
	7.8852928780e-07, 1.2073073057e-06, 1.8915500132e-06, 2.5050582484e-06, 1.8147542772e-06, 
	8.9705243347e-07, 6.4538345648e-07, 4.3532923001e-07, 4.2541992693e-07, 4.2935628905e-07, 
	5.4552898278e-07, 7.0116112417e-07, 9.8586586266e-07, 1.3608737238e-06, 1.9670610034e-06, 
	2.4536794671e-06, 6.0092161902e-06, 8.2189782026e-06, 1.1087096528e-05, 8.7260526424e-06, 
	6.7424514613e-06, 5.4804703015e-06, 5.0213191858e-06, 4.6753898391e-06, 4.4546935040e-06, 
	4.7077980932e-06, 5.2822363397e-06, 9.9958062296e-06, 1.2642059461e-05, 1.5878805486e-05, 
	1.7939809568e-05, 1.7734454266e-05, 1.4615639709e-05, 1.1503177697e-05, 9.3716881495e-06, 
	8.5667842324e-06, 8.9498927013e-06, 7.9766002662e-06, 5.3928436415e-06, 4.5062764237e-06, 
	3.6712827073e-06, 3.1828621424e-06, 2.8498226159e-06, 2.8172009662e-06, 2.5930914229e-06, 
	2.2689129227e-06, 2.0836166955e-06, 1.9715955473e-06, 2.0788245101e-06, 2.4311839986e-06, 
	2.7153000566e-06, 2.4311839986e-06, 2.0408816867e-06, 2.1518789075e-06, 2.5340654558e-06, 
	2.7594201333e-06, 3.1391916366e-06, 3.8178409218e-06, 4.6325256074e-06, 6.4687121900e-06, 
	7.6527553327e-06, 8.3717804158e-06, 8.9293085193e-06, 8.2189782026e-06, 6.7114727108e-06, 
	5.3188510528e-06, 4.5584566454e-06, 4.4138526285e-06, 4.7404309906e-06, 5.2944131521e-06, 
	5.6861436699e-06, 5.9131354681e-06, 5.5184591068e-06, 5.0213191858e-06, 4.3733861843e-06, 
	3.7052526680e-06, 3.2795769748e-06, 2.7913727125e-06, 2.3432456996e-06, 2.1868441270e-06, 
	2.2377822703e-06, 2.2951856961e-06, 2.1918853233e-06, 2.0645139828e-06, 1.9002810113e-06, 
	1.7290909140e-06, 1.6936272341e-06, 1.7491128431e-06, 1.6897319995e-06, 1.6211299474e-06, 
	1.5588984623e-06, 1.6474711901e-06, 1.7857383630e-06, 1.5481670808e-06, 1.9046616115e-06, 
	1.6936272341e-06, 2.0932342332e-06, 1.9715955473e-06, 1.8147542772e-06, 1.6819683862e-06, 
	1.7734454266e-06, 2.1518789075e-06, 1.9400719146e-06, 2.5752407071e-06, 2.3486474388e-06, 
	2.6780447431e-06, 2.4650051370e-06, 2.1077438397e-06, 1.7775336405e-06, 1.6062673388e-06, 
	1.5094468612e-06, 1.4548486945e-06, 1.4415105615e-06, 1.4184624015e-06, 1.3515055350e-06, 
	1.2996243168e-06, 1.1962386458e-06, 1.0934975861e-06, 9.4584033123e-07, 8.2000750773e-07, 
	7.1419668272e-07, 6.0230688642e-07, 4.9410338444e-07, 4.3532923001e-07, 3.7915590898e-07, 
	3.3252004352e-07, 2.8828220065e-07, 2.5752407071e-07, 1.9002810113e-07, 1.6436821115e-07, 
	1.4152000286e-07, 1.2439926909e-07, 1.0588094985e-07, 5.4052511957e-08, 2.7593953899e-08, 
	1.4086788277e-08, 7.1913436068e-09, 3.6712004081e-09, 1.8741577615e-09, 9.5676261840e-10, 
	4.8842991063e-10, 2.4934479358e-10, 1.2729119314e-10, 6.4982499203e-11, 3.3173742021e-11, 
	1.6935285241e-11, 8.6455090295e-12, 4.4135557989e-12, 2.2531322012e-12, 1.1502300973e-12, 
	5.8719558312e-13, 2.9976493716e-13, 1.5303081313e-13, 7.8122644993e-14, 3.9881822071e-14, 
	2.0359778293e-14, 1.0393722017e-14, 5.3060232688e-15, 2.7087392642e-15, 1.3828187382e-15, 
	7.0593271487e-16, 3.6038056483e-16, 1.8397525540e-16, 9.3919866665e-17, 4.7946346564e-17, 
	2.4476739911e-17, 1.2495442085e-17, 6.3789570613e-18, 3.2564748741e-18, 1.6624392521e-18, 
	8.4867974531e-19, 4.3325331088e-19, 2.2117698982e-19, 1.1291145295e-19, 5.7641602851e-20, 
	2.9426194531e-20, 1.5022152087e-20, 7.6688493678e-21, 3.9149683937e-21, 1.9986019791e-21, 
	1.0202917288e-21, 5.2086169374e-22, 2.6590130679e-22, 1.3574333802e-22, 6.9297342078e-23, 
	3.5376481007e-23, 1.8059789465e-23, 9.2195714848e-24, 4.7066162387e-24, 2.4027403502e-24, 
	1.2266054630e-24, 6.2618541437e-25, 3.1966935172e-25, 1.6319207073e-25, 8.3309994545e-26, 
	4.2529978082e-26, 2.1711669116e-26, 1.1083865947e-26, 5.6583436163e-27, 2.8885997567e-27, 
	1.4746380073e-27, 7.5280670068e-28 };

  } // endnamespace pah
  
} // endnamespace sed

#endif //__PAH_H__
