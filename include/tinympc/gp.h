#include <Eigen/Eigen/Dense>

static Eigen::Matrix<double, 20, 10> X_x = (Eigen::Matrix<double, 20, 10>() <<
2.525210194289680132e-02,1.054893732070922852e+00,7.599338889121998175e-02,-4.687120765447609638e-02,2.354632616043089710e-01,-2.448737621307369577e-02,-1.004404053092002036e-01,-1.969468891620635986e-01,6.178320646286010742e-01,0.000000000000000000e+00,
-2.251791395246980146e-02,1.012061014771460793e-01,7.532731443643571334e-02,9.790093451738350605e-02,2.922848165035247803e-01,-9.477119892835599718e-03,-4.315099418163299005e-01,1.583226770162582120e-01,4.862392842769623358e-01,0.000000000000000000e+00,
2.560034692287445068e-01,1.078305721282958762e+00,4.932625219225880014e-02,6.076542660593980266e-02,-1.317088901996611994e-01,-4.607554897665969845e-02,3.046501576900482178e-01,-3.206796348094940186e-01,2.223450690507887961e-01,0.000000000000000000e+00,
3.747352063655852716e-01,8.806998133659361683e-01,5.496795848011969826e-02,-5.157397687435150146e-02,4.715345203876495361e-01,3.041370958089820165e-02,3.360834717750540296e-02,5.770052596926680127e-02,-2.818708866834640156e-02,0.000000000000000000e+00,
4.789977073669433039e-01,-7.929583638906470555e-02,4.043506830930709839e-02,1.232387721538543007e-01,8.120989799499510331e-02,-3.323925286531440038e-02,1.967552304267879831e-02,2.858565375208850165e-02,7.376379519700998477e-02,0.000000000000000000e+00,
3.806606233119964044e-01,1.225205421447753906e+00,5.473101139068600046e-02,1.291765421628952026e-01,-1.247422173619270047e-01,8.052053744903998715e-04,2.737577259540549535e-02,-3.018080294132232111e-01,-1.788078695535658680e-01,0.000000000000000000e+00,
3.162885606288909357e-01,9.017521739006041370e-01,5.037712678313249759e-02,7.723655086010599961e-03,-2.762356996536254883e-01,4.173275455832479996e-02,-1.674311310052870871e-01,-1.257086545228958130e-01,3.070139139890670082e-02,0.000000000000000000e+00,
2.128943353891372126e-01,5.723236501216880101e-02,4.719782620668409867e-02,3.781954525038599968e-03,3.468865752220153809e-01,-2.926944755017750130e-02,-1.984628289937972745e-01,1.197622343897818825e-01,5.140603184700012207e-01,0.000000000000000000e+00,
3.907436132431030273e-01,2.073203027248380140e-02,3.993448987603179928e-02,5.934532359242430249e-02,3.246031701564788818e-01,-5.526817589998239688e-02,-1.625695079565047940e-01,2.081824243068695068e-01,1.690512001514433982e-01,0.000000000000000000e+00,
5.518692135810852051e-01,4.971239715814589760e-02,4.469018056988709880e-02,4.423968959599700393e-03,3.341951221227640323e-02,1.575475744903079983e-02,-9.778903424739826544e-02,-1.435705721378326138e-01,-5.091053992509839837e-02,0.000000000000000000e+00,
5.338028669357299805e-01,5.943208932876586914e-01,2.987354621291159890e-02,2.215995825827120000e-02,-4.015657603740692139e-01,-3.248623013496389905e-02,-3.787758350372313343e-01,1.071962565183638832e-01,6.115486025810241699e-01,0.000000000000000000e+00,
7.083839923143379902e-02,1.100129127502441406e+00,6.661289185285559911e-02,8.671657741069788150e-02,-2.073245644569396695e-01,1.386250276118509983e-02,-1.365879178047179898e-01,-3.602059185504912775e-01,1.353314667940138938e-01,0.000000000000000000e+00,
4.171759784221648615e-01,-2.360611595213410116e-02,5.108126997947689402e-02,-2.848456986248490072e-02,-1.630377471446990689e-01,-1.098233368247739962e-02,6.993857026100150365e-02,1.675122976303100031e-01,-1.868413537740707120e-01,0.000000000000000000e+00,
3.952841460704802357e-01,4.524818658828735352e-01,4.338358715176580255e-02,-5.536650493741029910e-02,4.450365900993347168e-01,4.155770316719999830e-02,-1.506113409996031882e-01,4.829676449298849622e-02,3.141592070460309805e-02,0.000000000000000000e+00,
1.024060472846031050e-01,7.160696387290950427e-02,5.759761855006209630e-02,1.357473619282240085e-02,-2.869196236133574884e-01,-5.912186577916139774e-02,-9.234052267856000373e-04,2.298715859651564719e-01,-1.289369612932204923e-01,0.000000000000000000e+00,
2.469863444566726129e-01,1.218516007065770061e-02,4.269332438707350297e-02,1.566695608198639955e-02,-2.948193550109862726e-01,2.209597080945959954e-02,-1.925233453512190662e-01,1.691911667585372925e-01,-1.966827921569340096e-02,0.000000000000000000e+00,
1.391671299934386929e-01,7.392155528068542480e-01,6.487067788839340210e-02,1.865550689399239973e-02,-3.465391993522644043e-01,-2.002699859440320099e-02,-2.807072103023529053e-01,5.695895478129379963e-02,6.828495860099789705e-02,0.000000000000000000e+00,
1.754518300294876099e-01,1.019764542579650879e+00,5.397808924317359924e-02,-3.628390654921529596e-02,1.929001361131668091e-01,2.194768376648420158e-02,1.528149247169494074e-01,-3.952201306819915771e-01,-4.232682660222049364e-02,0.000000000000000000e+00,
-9.639078751206300866e-03,2.609251299872899880e-03,1.544445604085921964e-01,-7.139932364225380634e-02,-2.292802557349200160e-02,-7.408968806266783558e-01,-3.006224036216735840e-01,-1.078206226229666831e-01,8.403607010841369629e-01,0.000000000000000000e+00,
1.215276196599005959e-01,-1.442076116800307950e-01,5.242563784122460102e-02,6.448885798454280505e-02,-3.619535639882080075e-02,-1.314104534685600015e-03,1.158854439854621055e-01,4.081059694290160578e-01,2.476120591163634699e-01,0.000000000000000000e+00
).finished();

static Eigen::Matrix<double, 20, 3> Y_x = (Eigen::Matrix<double, 20, 3>() <<
5.582386255264280145e-02,-1.627518534660339078e-01,7.853822708129881702e-01,
-1.732333600521087091e-01,-2.508087456226339856e-02,7.248708009719847523e-01,
-4.035722017288208008e-01,-6.902384757995599918e-02,1.357974648475646751e+00,
-3.661954775452609667e-02,-2.388997375965118131e-01,1.335805773735046387e+00,
-1.132963746786116860e-01,-4.606480896472929520e-02,1.404734253883361816e+00,
-1.983494758605956754e-01,1.011072248220443032e-01,1.344862103462219238e+00,
-1.018607467412948053e-01,6.587640941143030338e-02,1.147483706474303977e+00,
-5.034810304641720163e-02,9.477403014898300171e-02,9.635381102561951794e-01,
-7.783040404319760408e-02,3.048740327358245850e-01,1.522274374961853027e+00,
-3.678511083126059789e-02,1.731926202774039872e-02,1.390193939208984375e+00,
1.889273971319198053e-01,-1.620316505432119886e-02,1.089550256729125977e+00,
2.475151419639580117e-02,-8.954030275344840306e-02,1.167394161224365234e+00,
-6.969685852527610082e-02,6.159329414367670230e-02,1.458480119705199973e+00,
2.690180391073219990e-02,-1.360644102096557062e-01,1.197959542274475098e+00,
2.912267856299869534e-02,-7.172435522079400119e-03,1.527319192886352317e+00,
3.466607630252829808e-02,4.879355430602999588e-03,1.139307022094726340e+00,
9.488518536090848055e-02,1.437393575906752707e-01,1.171738386154174805e+00,
-1.659252941608428678e-01,1.216307580471037986e-01,1.270860433578491211e+00,
6.648725867271423340e-01,-4.005653411149969617e-02,2.736298799514770508e+00,
-1.109584420919418057e-01,-1.509388089179992121e-01,1.081993579864501953e+00
).finished();

static Eigen::Matrix<double, 20, 10> X_y = (Eigen::Matrix<double, 20, 10>() <<
2.525210194289680132e-02,1.054893732070922852e+00,7.599338889121998175e-02,-4.687120765447609638e-02,2.354632616043089710e-01,-2.448737621307369577e-02,-1.004404053092002036e-01,-1.969468891620635986e-01,6.178320646286010742e-01,0.000000000000000000e+00,
-2.251791395246980146e-02,1.012061014771460793e-01,7.532731443643571334e-02,9.790093451738350605e-02,2.922848165035247803e-01,-9.477119892835599718e-03,-4.315099418163299005e-01,1.583226770162582120e-01,4.862392842769623358e-01,0.000000000000000000e+00,
2.560034692287445068e-01,1.078305721282958762e+00,4.932625219225880014e-02,6.076542660593980266e-02,-1.317088901996611994e-01,-4.607554897665969845e-02,3.046501576900482178e-01,-3.206796348094940186e-01,2.223450690507887961e-01,0.000000000000000000e+00,
3.747352063655852716e-01,8.806998133659361683e-01,5.496795848011969826e-02,-5.157397687435150146e-02,4.715345203876495361e-01,3.041370958089820165e-02,3.360834717750540296e-02,5.770052596926680127e-02,-2.818708866834640156e-02,0.000000000000000000e+00,
4.789977073669433039e-01,-7.929583638906470555e-02,4.043506830930709839e-02,1.232387721538543007e-01,8.120989799499510331e-02,-3.323925286531440038e-02,1.967552304267879831e-02,2.858565375208850165e-02,7.376379519700998477e-02,0.000000000000000000e+00,
3.806606233119964044e-01,1.225205421447753906e+00,5.473101139068600046e-02,1.291765421628952026e-01,-1.247422173619270047e-01,8.052053744903998715e-04,2.737577259540549535e-02,-3.018080294132232111e-01,-1.788078695535658680e-01,0.000000000000000000e+00,
3.162885606288909357e-01,9.017521739006041370e-01,5.037712678313249759e-02,7.723655086010599961e-03,-2.762356996536254883e-01,4.173275455832479996e-02,-1.674311310052870871e-01,-1.257086545228958130e-01,3.070139139890670082e-02,0.000000000000000000e+00,
2.128943353891372126e-01,5.723236501216880101e-02,4.719782620668409867e-02,3.781954525038599968e-03,3.468865752220153809e-01,-2.926944755017750130e-02,-1.984628289937972745e-01,1.197622343897818825e-01,5.140603184700012207e-01,0.000000000000000000e+00,
3.907436132431030273e-01,2.073203027248380140e-02,3.993448987603179928e-02,5.934532359242430249e-02,3.246031701564788818e-01,-5.526817589998239688e-02,-1.625695079565047940e-01,2.081824243068695068e-01,1.690512001514433982e-01,0.000000000000000000e+00,
5.518692135810852051e-01,4.971239715814589760e-02,4.469018056988709880e-02,4.423968959599700393e-03,3.341951221227640323e-02,1.575475744903079983e-02,-9.778903424739826544e-02,-1.435705721378326138e-01,-5.091053992509839837e-02,0.000000000000000000e+00,
5.338028669357299805e-01,5.943208932876586914e-01,2.987354621291159890e-02,2.215995825827120000e-02,-4.015657603740692139e-01,-3.248623013496389905e-02,-3.787758350372313343e-01,1.071962565183638832e-01,6.115486025810241699e-01,0.000000000000000000e+00,
7.083839923143379902e-02,1.100129127502441406e+00,6.661289185285559911e-02,8.671657741069788150e-02,-2.073245644569396695e-01,1.386250276118509983e-02,-1.365879178047179898e-01,-3.602059185504912775e-01,1.353314667940138938e-01,0.000000000000000000e+00,
4.171759784221648615e-01,-2.360611595213410116e-02,5.108126997947689402e-02,-2.848456986248490072e-02,-1.630377471446990689e-01,-1.098233368247739962e-02,6.993857026100150365e-02,1.675122976303100031e-01,-1.868413537740707120e-01,0.000000000000000000e+00,
3.952841460704802357e-01,4.524818658828735352e-01,4.338358715176580255e-02,-5.536650493741029910e-02,4.450365900993347168e-01,4.155770316719999830e-02,-1.506113409996031882e-01,4.829676449298849622e-02,3.141592070460309805e-02,0.000000000000000000e+00,
1.024060472846031050e-01,7.160696387290950427e-02,5.759761855006209630e-02,1.357473619282240085e-02,-2.869196236133574884e-01,-5.912186577916139774e-02,-9.234052267856000373e-04,2.298715859651564719e-01,-1.289369612932204923e-01,0.000000000000000000e+00,
2.469863444566726129e-01,1.218516007065770061e-02,4.269332438707350297e-02,1.566695608198639955e-02,-2.948193550109862726e-01,2.209597080945959954e-02,-1.925233453512190662e-01,1.691911667585372925e-01,-1.966827921569340096e-02,0.000000000000000000e+00,
1.391671299934386929e-01,7.392155528068542480e-01,6.487067788839340210e-02,1.865550689399239973e-02,-3.465391993522644043e-01,-2.002699859440320099e-02,-2.807072103023529053e-01,5.695895478129379963e-02,6.828495860099789705e-02,0.000000000000000000e+00,
1.754518300294876099e-01,1.019764542579650879e+00,5.397808924317359924e-02,-3.628390654921529596e-02,1.929001361131668091e-01,2.194768376648420158e-02,1.528149247169494074e-01,-3.952201306819915771e-01,-4.232682660222049364e-02,0.000000000000000000e+00,
-9.639078751206300866e-03,2.609251299872899880e-03,1.544445604085921964e-01,-7.139932364225380634e-02,-2.292802557349200160e-02,-7.408968806266783558e-01,-3.006224036216735840e-01,-1.078206226229666831e-01,8.403607010841369629e-01,0.000000000000000000e+00,
1.215276196599005959e-01,-1.442076116800307950e-01,5.242563784122460102e-02,6.448885798454280505e-02,-3.619535639882080075e-02,-1.314104534685600015e-03,1.158854439854621055e-01,4.081059694290160578e-01,2.476120591163634699e-01,0.000000000000000000e+00
).finished();

static Eigen::Matrix<double, 20, 3> Y_y = (Eigen::Matrix<double, 20, 3>() <<
5.582386255264280145e-02,-1.627518534660339078e-01,7.853822708129881702e-01,
-1.732333600521087091e-01,-2.508087456226339856e-02,7.248708009719847523e-01,
-4.035722017288208008e-01,-6.902384757995599918e-02,1.357974648475646751e+00,
-3.661954775452609667e-02,-2.388997375965118131e-01,1.335805773735046387e+00,
-1.132963746786116860e-01,-4.606480896472929520e-02,1.404734253883361816e+00,
-1.983494758605956754e-01,1.011072248220443032e-01,1.344862103462219238e+00,
-1.018607467412948053e-01,6.587640941143030338e-02,1.147483706474303977e+00,
-5.034810304641720163e-02,9.477403014898300171e-02,9.635381102561951794e-01,
-7.783040404319760408e-02,3.048740327358245850e-01,1.522274374961853027e+00,
-3.678511083126059789e-02,1.731926202774039872e-02,1.390193939208984375e+00,
1.889273971319198053e-01,-1.620316505432119886e-02,1.089550256729125977e+00,
2.475151419639580117e-02,-8.954030275344840306e-02,1.167394161224365234e+00,
-6.969685852527610082e-02,6.159329414367670230e-02,1.458480119705199973e+00,
2.690180391073219990e-02,-1.360644102096557062e-01,1.197959542274475098e+00,
2.912267856299869534e-02,-7.172435522079400119e-03,1.527319192886352317e+00,
3.466607630252829808e-02,4.879355430602999588e-03,1.139307022094726340e+00,
9.488518536090848055e-02,1.437393575906752707e-01,1.171738386154174805e+00,
-1.659252941608428678e-01,1.216307580471037986e-01,1.270860433578491211e+00,
6.648725867271423340e-01,-4.005653411149969617e-02,2.736298799514770508e+00,
-1.109584420919418057e-01,-1.509388089179992121e-01,1.081993579864501953e+00
).finished();

