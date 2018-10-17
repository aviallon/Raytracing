#include "includes.h"

inline int clamp(double x, double min, double max){
	if(x <= min)
		x = min;
	if(x >= max)
		x = max;
	return x;
}

//inline int clamp(double x, double min, double max){
//	return floor((x>max) ? max : ((x<min) ? min : x));
//}

void pVec(Vec v){
	std::cout << (std::string)v << std::endl;
}

static float sin_vals[256] = { 0.0, 0.024541228522912288, 0.049067674327418015, 0.07356456359966743, 0.0980171403295606, 0.1224106751992162, 0.14673047445536175, 0.17096188876030122, 0.19509032201612825, 0.2191012401568698, 0.24298017990326387, 0.26671275747489837, 0.29028467725446233, 0.3136817403988915, 0.33688985339222005, 0.3598950365349881, 0.3826834323650898, 0.40524131400498986, 0.4275550934302821, 0.44961132965460654, 0.47139673682599764, 0.49289819222978404, 0.5141027441932217, 0.5349976198870972, 0.5555702330196022, 0.5758081914178453, 0.5956993044924334, 0.6152315905806268, 0.6343932841636455, 0.6531728429537768, 0.6715589548470183, 0.6895405447370668, 0.7071067811865475, 0.7242470829514669, 0.7409511253549591, 0.7572088465064845, 0.773010453362737, 0.7883464276266062, 0.8032075314806448, 0.8175848131515837, 0.8314696123025452, 0.844853565249707, 0.8577286100002721, 0.8700869911087113, 0.8819212643483549, 0.8932243011955153, 0.9039892931234433, 0.9142097557035307, 0.9238795325112867, 0.9329927988347388, 0.9415440651830208, 0.9495281805930367, 0.9569403357322089, 0.9637760657954398, 0.970031253194544, 0.9757021300385286, 0.9807852804032304, 0.9852776423889412, 0.989176509964781, 0.99247953459871, 0.9951847266721968, 0.9972904566786902, 0.9987954562051724, 0.9996988186962042, 1.0, 0.9996988186962042, 0.9987954562051724, 0.9972904566786902, 0.9951847266721969, 0.99247953459871, 0.989176509964781, 0.9852776423889412, 0.9807852804032304, 0.9757021300385286, 0.970031253194544, 0.9637760657954398, 0.9569403357322089, 0.9495281805930367, 0.9415440651830208, 0.9329927988347388, 0.9238795325112867, 0.9142097557035307, 0.9039892931234434, 0.8932243011955152, 0.881921264348355, 0.8700869911087115, 0.8577286100002721, 0.8448535652497072, 0.8314696123025455, 0.8175848131515837, 0.8032075314806449, 0.7883464276266063, 0.7730104533627371, 0.7572088465064847, 0.740951125354959, 0.7242470829514669, 0.7071067811865476, 0.689540544737067, 0.6715589548470186, 0.6531728429537766, 0.6343932841636455, 0.6152315905806269, 0.5956993044924335, 0.5758081914178454, 0.5555702330196022, 0.5349976198870972, 0.5141027441932218, 0.49289819222978415, 0.47139673682599786, 0.4496113296546069, 0.42755509343028203, 0.4052413140049899, 0.3826834323650899, 0.35989503653498833, 0.33688985339222033, 0.3136817403988914, 0.2902846772544624, 0.2667127574748985, 0.24298017990326407, 0.21910124015687005, 0.1950903220161286, 0.17096188876030122, 0.1467304744553618, 0.12241067519921635, 0.09801714032956083, 0.07356456359966773, 0.049067674327417966, 0.024541228522912326, 1.2246467991473532e-16, -0.02454122852291208, -0.049067674327417724, -0.0735645635996675, -0.09801714032956059, -0.1224106751992161, -0.14673047445536158, -0.17096188876030097, -0.19509032201612836, -0.2191012401568698, -0.24298017990326382, -0.26671275747489825, -0.2902846772544621, -0.3136817403988912, -0.3368898533922201, -0.3598950365349881, -0.38268343236508967, -0.4052413140049897, -0.4275550934302818, -0.44961132965460665, -0.47139673682599764, -0.4928981922297839, -0.5141027441932216, -0.5349976198870969, -0.555570233019602, -0.5758081914178453, -0.5956993044924332, -0.6152315905806267, -0.6343932841636453, -0.6531728429537765, -0.6715589548470184, -0.6895405447370668, -0.7071067811865475, -0.7242470829514668, -0.7409511253549589, -0.7572088465064842, -0.7730104533627367, -0.7883464276266059, -0.803207531480645, -0.8175848131515838, -0.8314696123025452, -0.844853565249707, -0.857728610000272, -0.8700869911087113, -0.8819212643483549, -0.8932243011955152, -0.9039892931234431, -0.9142097557035305, -0.9238795325112865, -0.932992798834739, -0.9415440651830208, -0.9495281805930367, -0.9569403357322088, -0.9637760657954398, -0.970031253194544, -0.9757021300385285, -0.9807852804032303, -0.9852776423889411, -0.9891765099647809, -0.9924795345987101, -0.9951847266721969, -0.9972904566786902, -0.9987954562051724, -0.9996988186962042, -1.0, -0.9996988186962042, -0.9987954562051724, -0.9972904566786902, -0.9951847266721969, -0.9924795345987101, -0.9891765099647809, -0.9852776423889412, -0.9807852804032304, -0.9757021300385286, -0.970031253194544, -0.96377606579544, -0.9569403357322089, -0.9495281805930368, -0.9415440651830209, -0.9329927988347391, -0.9238795325112866, -0.9142097557035306, -0.9039892931234433, -0.8932243011955153, -0.881921264348355, -0.8700869911087115, -0.8577286100002722, -0.8448535652497072, -0.8314696123025455, -0.817584813151584, -0.8032075314806453, -0.7883464276266061, -0.7730104533627369, -0.7572088465064846, -0.7409511253549591, -0.724247082951467, -0.7071067811865477, -0.6895405447370672, -0.6715589548470187, -0.6531728429537771, -0.6343932841636459, -0.6152315905806274, -0.5956993044924332, -0.5758081914178452, -0.5555702330196022, -0.5349976198870973, -0.5141027441932219, -0.49289819222978426, -0.4713967368259979, -0.449611329654607, -0.42755509343028253, -0.4052413140049904, -0.3826834323650904, -0.359895036534988, -0.33688985339222, -0.3136817403988915, -0.2902846772544625, -0.2667127574748986, -0.24298017990326418, -0.21910124015687016, -0.19509032201612872, -0.17096188876030177, -0.1467304744553624, -0.12241067519921603, -0.0980171403295605, -0.07356456359966741, -0.04906767432741809, -0.024541228522912448 };

float st_sin(float x){
	int s_x = int(x * 128 / PI) % 256;
	
	return sin_vals[s_x];
}

float st_cos(float x){
	return st_sin(x + PI/2);
}

static float tan_vals[512] = {
1.2246467991473532e-16, 0.012272462379566355, 0.024548622108925482, 0.036832180994845636, 0.049126849769467205, 0.061436352581593676, 0.0737644315224496, 0.08611485119762818, 0.09849140335716448, 0.11089791159591321, 0.12333823613673883, 0.13581627870938784, 0.14833598753834748, 0.16090136245348918, 0.17351646013785574, 0.18618539952758367, 0.19891236737965837, 0.21170162402398365, 0.2245575093171296, 0.23748444881607042, 0.2504869601913056, 0.26356965989991815, 0.27673727014041444, 0.2899946261126062, 0.3033466836073424, 0.31679852695260374, 0.33035537734433384, 0.3440226015924267, 0.3578057213145244, 0.3717104226127437, 0.38574256627112147, 0.39990819851453735, 0.4142135623730952, 0.42866510969949956, 0.4432695138908644, 0.4580336833706724, 0.47296477589131986, 0.48807021372286274, 0.5033576997992947, 0.5188352348999761, 0.5345111359507919, 0.5503940555372643, 0.5664930027303442, 0.5828173653349762, 0.5993769336819238, 0.6161819260948661, 0.6332430161775691, 0.6505713620801532, 0.6681786379192988, 0.6860770675448633, 0.7042794608650446, 0.7227992529642062, 0.7416505462720356, 0.7608481560702515, 0.7804076596539438, 0.8003454494993202, 0.8206787908286604, 0.8414258840072547, 0.8626059322567398, 0.8842392152253504, 0.9063471690191476, 0.9289524733703679, 0.9520791467009256, 0.9757526499323768, 1.0000000000000002, 1.0248498941502275, 1.0503328462398598, 1.0764813364152659, 1.1033299757334756, 1.1309156874988269, 1.1592779073334354, 1.188458804282967, 1.2185035255879768, 1.2494604681335792, 1.281381580036555, 1.3143226963510801, 1.3483439134867203, 1.3835100076528737, 1.4198909034940923, 1.4575622000871051, 1.49660576266549, 1.5371103898618828, 1.5791725679602098, 1.6228973256934551, 1.6683992055835075, 1.7158033707956644, 1.7652468700941917, 1.8168800878924023, 1.8708684117893895, 1.9273941566300639, 1.9866587923433645, 2.048885533030752, 2.114322357548642, 2.1832455478841535, 2.2559638519291596, 2.332823403101351, 2.4142135623730954, 2.500573890994257, 2.5924025177380723, 2.6902662372796144, 2.794812772490478, 2.906785761665536, 3.0270432043177746, 3.1565803339407887, 3.2965582089383227, 3.448339762033026, 3.6135356813074297, 3.794063400088302, 3.9922237837700845, 4.2108020335028, 4.453202224414413, 4.723629327882303, 5.02733949212585, 5.37099043500373, 5.763142005118809, 6.214987771089044, 6.741452405414994, 7.362887641324249, 8.10778580367691, 9.017302360424734, 10.153170387608856, 11.612398861435265, 13.556669242352443, 16.27700795993542, 20.35546762498719, 27.150170665699672, 40.73548387208354, 81.48324020654685, -1.633123935319537e+16, -81.48324020654604, -40.73548387208334, -27.15017066569958, -20.355467624987142, -16.27700795993539, -13.55666924235242, -11.612398861435247, -10.153170387608842, -9.017302360424724, -8.107785803676903, -7.362887641324242, -6.7414524054149885, -6.2149877710890395, -5.763142005118804, -5.370990435003726, -5.027339492125846, -4.723629327882301, -4.453202224414411, -4.210802033502797, -3.9922237837700827, -3.7940634000883002, -3.613535681307428, -3.4483397620330245, -3.296558208938321, -3.1565803339407874, -3.0270432043177733, -2.9067857616655353, -2.7948127724904768, -2.690266237279613, -2.5924025177380714, -2.5005738909942563, -2.414213562373095, -2.33282340310135, -2.2559638519291587, -2.1832455478841513, -2.11432235754864, -2.0488855330307496, -1.986658792343365, -1.9273941566300643, -1.8708684117893888, -1.8168800878924019, -1.7652468700941912, -1.715803370795664, -1.668399205583507, -1.6228973256934547, -1.5791725679602087, -1.5371103898618816, -1.496605762665489, -1.4575622000871054, -1.4198909034940923, -1.3835100076528737, -1.34834391348672, -1.31432269635108, -1.2813815800365544, -1.2494604681335786, -1.2185035255879764, -1.1884588042829665, -1.1592779073334345, -1.130915687498827, -1.1033299757334758, -1.0764813364152659, -1.0503328462398598, -1.0248498941502273, -0.9999999999999999, -0.9757526499323765, -0.9520791467009252, -0.9289524733703675, -0.906347169019147, -0.8842392152253498, -0.8626059322567399, -0.8414258840072547, -0.8206787908286602, -0.8003454494993202, -0.7804076596539435, -0.7608481560702512, -0.7416505462720354, -0.7227992529642059, -0.7042794608650442, -0.6860770675448629, -0.6681786379192989, -0.6505713620801533, -0.6332430161775691, -0.616181926094866, -0.5993769336819237, -0.5828173653349761, -0.566493002730344, -0.550394055537264, -0.5345111359507916, -0.5188352348999757, -0.5033576997992942, -0.48807021372286286, -0.4729647758913199, -0.4580336833706724, -0.44326951389086433, -0.4286651096994995, -0.41421356237309503, -0.3999081985145372, -0.3857425662711212, -0.37171042261274345, -0.3578057213145241, -0.3440226015924263, -0.33035537734433396, -0.3167985269526038, -0.3033466836073424, -0.2899946261126061, -0.27673727014041427, -0.263569659899918, -0.25048696019130545, -0.23748444881607017, -0.2245575093171293, -0.21170162402398335, -0.198912367379658, -0.1861853995275837, -0.17351646013785577, -0.16090136245348918, -0.14833598753834742, -0.1358162787093877, -0.12333823613673868, -0.11089791159591302, -0.09849140335716425, -0.0861148511976279, -0.07376443152244928, -0.06143635258159376, -0.049126849769467254, -0.036832180994845636, -0.024548622108925444, -0.012272462379566276, 0.0, 0.012272462379566276, 0.024548622108925444, 0.036832180994845636, 0.049126849769467254, 0.06143635258159376, 0.07376443152244928, 0.0861148511976279, 0.09849140335716425, 0.11089791159591302, 0.12333823613673868, 0.1358162787093877, 0.14833598753834742, 0.16090136245348918, 0.17351646013785577, 0.1861853995275837, 0.198912367379658, 0.21170162402398335, 0.2245575093171293, 0.23748444881607017, 0.25048696019130545, 0.263569659899918, 0.27673727014041427, 0.2899946261126061, 0.3033466836073424, 0.3167985269526038, 0.33035537734433396, 0.3440226015924263, 0.3578057213145241, 0.37171042261274345, 0.3857425662711212, 0.3999081985145372, 0.41421356237309503, 0.4286651096994995, 0.44326951389086433, 0.4580336833706724, 0.4729647758913199, 0.48807021372286286, 0.5033576997992942, 0.5188352348999757, 0.5345111359507916, 0.550394055537264, 0.566493002730344, 0.5828173653349761, 0.5993769336819237, 0.616181926094866, 0.6332430161775691, 0.6505713620801533, 0.6681786379192989, 0.6860770675448629, 0.7042794608650442, 0.7227992529642059, 0.7416505462720354, 0.7608481560702512, 0.7804076596539435, 0.8003454494993202, 0.8206787908286602, 0.8414258840072547, 0.8626059322567399, 0.8842392152253498, 0.906347169019147, 0.9289524733703675, 0.9520791467009252, 0.9757526499323765, 0.9999999999999999, 1.0248498941502273, 1.0503328462398598, 1.0764813364152659, 1.1033299757334758, 1.130915687498827, 1.1592779073334345, 1.1884588042829665, 1.2185035255879764, 1.2494604681335786, 1.2813815800365544, 1.31432269635108, 1.34834391348672, 1.3835100076528737, 1.4198909034940923, 1.4575622000871054, 1.496605762665489, 1.5371103898618816, 1.5791725679602087, 1.6228973256934547, 1.668399205583507, 1.715803370795664, 1.7652468700941912, 1.8168800878924019, 1.8708684117893888, 1.9273941566300643, 1.986658792343365, 2.0488855330307496, 2.11432235754864, 2.1832455478841513, 2.2559638519291587, 2.33282340310135, 2.414213562373095, 2.5005738909942563, 2.5924025177380714, 2.690266237279613, 2.7948127724904768, 2.9067857616655353, 3.0270432043177733, 3.1565803339407874, 3.296558208938321, 3.4483397620330245, 3.613535681307428, 3.7940634000883002, 3.9922237837700827, 4.210802033502797, 4.453202224414411, 4.723629327882301, 5.027339492125846, 5.370990435003726, 5.763142005118804, 6.2149877710890395, 6.7414524054149885, 7.362887641324242, 8.107785803676903, 9.017302360424724, 10.153170387608842, 11.612398861435247, 13.55666924235242, 16.27700795993539, 20.355467624987142, 27.15017066569958, 40.73548387208334, 81.48324020654604, 1.633123935319537e+16, -81.48324020654685, -40.73548387208354, -27.150170665699672, -20.35546762498719, -16.27700795993542, -13.556669242352443, -11.612398861435265, -10.153170387608856, -9.017302360424734, -8.10778580367691, -7.362887641324249, -6.741452405414994, -6.214987771089044, -5.763142005118809, -5.37099043500373, -5.02733949212585, -4.723629327882303, -4.453202224414413, -4.2108020335028, -3.9922237837700845, -3.794063400088302, -3.6135356813074297, -3.448339762033026, -3.2965582089383227, -3.1565803339407887, -3.0270432043177746, -2.906785761665536, -2.794812772490478, -2.6902662372796144, -2.5924025177380723, -2.500573890994257, -2.4142135623730954, -2.332823403101351, -2.2559638519291596, -2.1832455478841535, -2.114322357548642, -2.048885533030752, -1.9866587923433645, -1.9273941566300639, -1.8708684117893895, -1.8168800878924023, -1.7652468700941917, -1.7158033707956644, -1.6683992055835075, -1.6228973256934551, -1.5791725679602098, -1.5371103898618828, -1.49660576266549, -1.4575622000871051, -1.4198909034940923, -1.3835100076528737, -1.3483439134867203, -1.3143226963510801, -1.281381580036555, -1.2494604681335792, -1.2185035255879768, -1.188458804282967, -1.1592779073334354, -1.1309156874988269, -1.1033299757334756, -1.0764813364152659, -1.0503328462398598, -1.0248498941502275, -1.0000000000000002, -0.9757526499323768, -0.9520791467009256, -0.9289524733703679, -0.9063471690191476, -0.8842392152253504, -0.8626059322567398, -0.8414258840072547, -0.8206787908286604, -0.8003454494993202, -0.7804076596539438, -0.7608481560702515, -0.7416505462720356, -0.7227992529642062, -0.7042794608650446, -0.6860770675448633, -0.6681786379192988, -0.6505713620801532, -0.6332430161775691, -0.6161819260948661, -0.5993769336819238, -0.5828173653349762, -0.5664930027303442, -0.5503940555372643, -0.5345111359507919, -0.5188352348999761, -0.5033576997992947, -0.48807021372286274, -0.47296477589131986, -0.4580336833706724, -0.4432695138908644, -0.42866510969949956, -0.4142135623730952, -0.39990819851453735, -0.38574256627112147, -0.3717104226127437, -0.3578057213145244, -0.3440226015924267, -0.33035537734433384, -0.31679852695260374, -0.3033466836073424, -0.2899946261126062, -0.27673727014041444, -0.26356965989991815, -0.2504869601913056, -0.23748444881607042, -0.2245575093171296, -0.21170162402398365, -0.19891236737965837, -0.18618539952758367, -0.17351646013785574, -0.16090136245348918, -0.14833598753834748, -0.13581627870938784, -0.12333823613673883, -0.11089791159591321, -0.09849140335716448, -0.08611485119762818, -0.0737644315224496, -0.061436352581593676, -0.049126849769467205, -0.036832180994845636, -0.024548622108925482, -0.012272462379566355 };

float st_tan(float x){
	
	int s_x = int(256 + x * 512 / (2*PI)) % 512;
	return tan_vals[s_x];
}