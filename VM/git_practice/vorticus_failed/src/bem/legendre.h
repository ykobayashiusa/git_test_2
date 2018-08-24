#ifndef BEM2_LEGENDRE_H
#define BEM2_LEGENDRE_H

#include "types.h"

namespace bem
{

constexpr real legendre_1_points[1] = 
{
	0.5
};

constexpr real legendre_1_weights[1] = 
{
	1
};


constexpr real legendre_2_points[2] = 
{
	0.211324865405187134,
	0.788675134594812866
};

constexpr real legendre_2_weights[2] = 
{
	0.5,
	0.5
};


constexpr real legendre_3_points[3] = 
{
	0.112701665379258298,
	0.5,
	0.887298334620741702
};

constexpr real legendre_3_weights[3] = 
{
	0.27777777777777779,
	0.44444444444444442,
	0.27777777777777779
};


constexpr real legendre_4_points[4] = 
{
	0.0694318442029737137,
	0.330009478207571871,
	0.669990521792428129,
	0.930568155797026231
};

constexpr real legendre_4_weights[4] = 
{
	0.173927422568726925,
	0.326072577431273047,
	0.326072577431273047,
	0.173927422568726925
};


constexpr real legendre_5_points[5] = 
{
	0.0469100770306680181,
	0.230765344947158446,
	0.5,
	0.769234655052841498,
	0.953089922969331926
};

constexpr real legendre_5_weights[5] = 
{
	0.118463442528094542,
	0.239314335249683235,
	0.284444444444444444,
	0.239314335249683235,
	0.118463442528094542
};


constexpr real legendre_6_points[6] = 
{
	0.033765242898423975,
	0.169395306766867759,
	0.380690406958401562,
	0.619309593041598494,
	0.830604693233132241,
	0.966234757101576025
};

constexpr real legendre_6_weights[6] = 
{
	0.0856622461895851783,
	0.180380786524069303,
	0.233956967286345519,
	0.233956967286345519,
	0.180380786524069303,
	0.0856622461895851783
};


constexpr real legendre_7_points[7] = 
{
	0.0254460438286207569,
	0.12923440720030277,
	0.297077424311301408,
	0.5,
	0.702922575688698537,
	0.87076559279969723,
	0.974553956171379188
};

constexpr real legendre_7_weights[7] = 
{
	0.0647424830844348514,
	0.139852695744638322,
	0.190915025252559462,
	0.208979591836734702,
	0.190915025252559462,
	0.139852695744638322,
	0.0647424830844348514
};


constexpr real legendre_8_points[8] = 
{
	0.0198550717512318564,
	0.101666761293186636,
	0.237233795041835505,
	0.40828267875217511,
	0.591717321247824946,
	0.762766204958164495,
	0.898333238706813364,
	0.980144928248768199
};

constexpr real legendre_8_weights[8] = 
{
	0.0506142681451881293,
	0.111190517226687241,
	0.156853322938943635,
	0.181341891689180995,
	0.181341891689180995,
	0.156853322938943635,
	0.111190517226687241,
	0.0506142681451881293
};


constexpr real legendre_9_points[9] = 
{
	0.0159198802461869571,
	0.0819844463366821152,
	0.193314283649704821,
	0.337873288298095542,
	0.5,
	0.662126711701904513,
	0.806685716350295179,
	0.918015553663317885,
	0.984080119753813043
};

constexpr real legendre_9_weights[9] = 
{
	0.0406371941807872061,
	0.090324080347428698,
	0.130305348201467719,
	0.156173538520001431,
	0.165119677500629891,
	0.156173538520001431,
	0.130305348201467719,
	0.090324080347428698,
	0.0406371941807872061
};


constexpr real legendre_10_points[10] = 
{
	0.0130467357414141283,
	0.0674683166555077318,
	0.160295215850487782,
	0.283302302935376393,
	0.42556283050918442,
	0.57443716949081558,
	0.716697697064623607,
	0.839704784149512218,
	0.932531683344492324,
	0.986953264258585872
};

constexpr real legendre_10_weights[10] = 
{
	0.033335672154344069,
	0.0747256745752902934,
	0.109543181257991021,
	0.134633359654998175,
	0.147762112357376435,
	0.147762112357376435,
	0.134633359654998175,
	0.109543181257991021,
	0.0747256745752902934,
	0.033335672154344069
};


constexpr real legendre_11_points[11] = 
{
	0.0108856709269715135,
	0.0564687001159523416,
	0.134923997212975322,
	0.240451935396594096,
	0.365228422023827548,
	0.5,
	0.634771577976172452,
	0.759548064603405848,
	0.865076002787024678,
	0.943531299884047714,
	0.989114329073028431
};

constexpr real legendre_11_weights[11] = 
{
	0.0278342835580868316,
	0.062790184732452306,
	0.0931451054638671311,
	0.116596882295995241,
	0.131402272255123326,
	0.136462543388950308,
	0.131402272255123326,
	0.116596882295995241,
	0.0931451054638671311,
	0.062790184732452306,
	0.0278342835580868316
};


constexpr real legendre_12_points[12] = 
{
	0.00921968287664037822,
	0.0479413718147625456,
	0.115048662902847654,
	0.206341022856691259,
	0.316084250500909936,
	0.437383295744265543,
	0.562616704255734401,
	0.683915749499090064,
	0.793658977143308686,
	0.884951337097152346,
	0.952058628185237454,
	0.990780317123359566
};

constexpr real legendre_12_weights[12] = 
{
	0.0235876681932559139,
	0.0534696629976592133,
	0.0800391642716731105,
	0.101583713361532962,
	0.116746268269177403,
	0.124573522906701387,
	0.124573522906701387,
	0.116746268269177403,
	0.101583713361532962,
	0.0800391642716731105,
	0.0534696629976592133,
	0.0235876681932559139
};


constexpr real legendre_13_points[13] = 
{
	0.00790847264070593248,
	0.0412008003885110385,
	0.0992109546333450609,
	0.178825330279829886,
	0.275753624481776538,
	0.384770842022432613,
	0.5,
	0.615229157977567387,
	0.724246375518223462,
	0.821174669720170058,
	0.900789045366654939,
	0.958799199611488961,
	0.992091527359294068
};

constexpr real legendre_13_weights[13] = 
{
	0.0202420023826579386,
	0.0460607499188642258,
	0.0694367551098936248,
	0.089072990380972869,
	0.103908023768444255,
	0.113141590131448616,
	0.116275776615436949,
	0.113141590131448616,
	0.103908023768444255,
	0.089072990380972869,
	0.0694367551098936248,
	0.0460607499188642258,
	0.0202420023826579386
};


constexpr real legendre_14_points[14] = 
{
	0.00685809565159384293,
	0.035782558168213241,
	0.0863993424651174902,
	0.156353547594157261,
	0.24237568182092295,
	0.340443815536055128,
	0.445972525646328166,
	0.554027474353671834,
	0.659556184463944817,
	0.75762431817907705,
	0.843646452405842684,
	0.91360065753488251,
	0.964217441831786815,
	0.993141904348406213
};

constexpr real legendre_14_weights[14] = 
{
	0.0175597301658759301,
	0.040079043579880104,
	0.0607592853439515926,
	0.0786015835790967732,
	0.092769198738968911,
	0.102599231860647802,
	0.107631926731578897,
	0.107631926731578897,
	0.102599231860647802,
	0.092769198738968911,
	0.0786015835790967732,
	0.0607592853439515926,
	0.040079043579880104,
	0.0175597301658759301
};


constexpr real legendre_15_points[15] = 
{
	0.0060037409897573113,
	0.0313633037996470243,
	0.0758967082947863969,
	0.137791134319914965,
	0.214513913695730585,
	0.302924326461218307,
	0.399402953001282757,
	0.5,
	0.600597046998717299,
	0.697075673538781748,
	0.785486086304269415,
	0.862208865680085035,
	0.924103291705213659,
	0.968636696200352976,
	0.993996259010242689
};

constexpr real legendre_15_weights[15] = 
{
	0.0153766209980586346,
	0.0351830237440540622,
	0.0535796102335859697,
	0.069785338963077162,
	0.0831346029084969601,
	0.0930805000077811057,
	0.0992157426635557893,
	0.101289120962780643,
	0.0992157426635557893,
	0.0930805000077811057,
	0.0831346029084969601,
	0.069785338963077162,
	0.0535796102335859697,
	0.0351830237440540622,
	0.0153766209980586346
};


constexpr real legendre_16_points[16] = 
{
	0.00529953250417503074,
	0.0277124884633836999,
	0.0671843988060841224,
	0.122297795822498501,
	0.191061877798678115,
	0.270991611171386315,
	0.359198224610370542,
	0.452493745081181287,
	0.547506254918818769,
	0.640801775389629458,
	0.729008388828613629,
	0.808938122201321885,
	0.877702204177501555,
	0.932815601193915933,
	0.9722875115366163,
	0.994700467495824969
};

constexpr real legendre_16_weights[16] = 
{
	0.0135762297058770482,
	0.0311267619693239468,
	0.0475792558412463928,
	0.0623144856277669384,
	0.074797994408288368,
	0.0845782596975012679,
	0.0913017075224617919,
	0.0947253052275342511,
	0.0947253052275342511,
	0.0913017075224617919,
	0.0845782596975012679,
	0.074797994408288368,
	0.0623144856277669384,
	0.0475792558412463928,
	0.0311267619693239468,
	0.0135762297058770482
};


constexpr real legendre_17_points[17] = 
{
	0.00471226234279131795,
	0.0246622391156161025,
	0.0598804231365070438,
	0.109242998051599316,
	0.171164420391654637,
	0.243654731456761531,
	0.32438411827306185,
	0.410757909252076059,
	0.5,
	0.589242090747923886,
	0.675615881726938206,
	0.756345268543238469,
	0.828835579608345308,
	0.890757001948400684,
	0.940119576863492901,
	0.975337760884383842,
	0.995287737657208682
};

constexpr real legendre_17_weights[17] = 
{
	0.0120741514342739657,
	0.0277297646869936014,
	0.0425180741585895888,
	0.055941923596701984,
	0.0675681842342627376,
	0.0770228805384051418,
	0.0840020510782250179,
	0.0882813526834963225,
	0.0897232351781032667,
	0.0882813526834963225,
	0.0840020510782250179,
	0.0770228805384051418,
	0.0675681842342627376,
	0.055941923596701984,
	0.0425180741585895888,
	0.0277297646869936014,
	0.0120741514342739657
};


constexpr real legendre_18_points[18] = 
{
	0.00421741578953455098,
	0.0220880252143011435,
	0.0536987667512221489,
	0.0981475205137384288,
	0.154156478469823388,
	0.22011458446302623,
	0.294124419268578685,
	0.374056887154247231,
	0.457612493479132354,
	0.542387506520867646,
	0.625943112845752769,
	0.705875580731421315,
	0.77988541553697377,
	0.845843521530176612,
	0.901852479486261571,
	0.946301233248777907,
	0.977911974785698801,
	0.995782584210465505
};

constexpr real legendre_18_weights[18] = 
{
	0.0108080067632416559,
	0.0248572744474848985,
	0.0382128651274445258,
	0.0504710220531435841,
	0.0612776033557392297,
	0.0703214573353253269,
	0.077342337563132621,
	0.0821382418729163649,
	0.0845711914815718002,
	0.0845711914815718002,
	0.0821382418729163649,
	0.077342337563132621,
	0.0703214573353253269,
	0.0612776033557392297,
	0.0504710220531435841,
	0.0382128651274445258,
	0.0248572744474848985,
	0.0108080067632416559
};


constexpr real legendre_19_points[19] = 
{
	0.003796578078207824,
	0.0198959239325849913,
	0.0484220481925910495,
	0.0886426717314285906,
	0.139516911332385307,
	0.199727347669159505,
	0.267714629312019503,
	0.341717950018185057,
	0.419820677179887303,
	0.5,
	0.580179322820112642,
	0.658282049981814943,
	0.732285370687980497,
	0.80027265233084055,
	0.860483088667614693,
	0.911357328268571409,
	0.951577951807409006,
	0.980104076067415009,
	0.99620342192179212
};

constexpr real legendre_19_weights[19] = 
{
	0.00973089411486323906,
	0.0224071133828497998,
	0.0345222713688206131,
	0.0457450108112249995,
	0.0557833227736669948,
	0.0643769812696681071,
	0.0713033510868033016,
	0.0763830210329298348,
	0.07948442169697717,
	0.0805272249243918492,
	0.07948442169697717,
	0.0763830210329298348,
	0.0713033510868033016,
	0.0643769812696681071,
	0.0557833227736669948,
	0.0457450108112249995,
	0.0345222713688206131,
	0.0224071133828497998,
	0.00973089411486323906
};


constexpr real legendre_20_points[20] = 
{
	0.00343570040745255767,
	0.0180140363610430954,
	0.0438827858743370269,
	0.0804415140888906088,
	0.126834046769924602,
	0.181973159636742488,
	0.244566499024586437,
	0.313146955642290226,
	0.386107074429177466,
	0.461736739433251331,
	0.538263260566748669,
	0.613892925570822534,
	0.686853044357709774,
	0.755433500975413619,
	0.818026840363257568,
	0.873165953230075398,
	0.919558485911109447,
	0.956117214125662973,
	0.98198596363895696,
	0.996564299592547442
};

constexpr real legendre_20_weights[20] = 
{
	0.00880700356957605894,
	0.0203007149001934693,
	0.0313360241670545339,
	0.0416383707883523774,
	0.0509650599086202208,
	0.059097265980759206,
	0.0658443192245883185,
	0.0710480546591910206,
	0.0745864932363018707,
	0.0763766935653629186,
	0.0763766935653629186,
	0.0745864932363018707,
	0.0710480546591910206,
	0.0658443192245883185,
	0.059097265980759206,
	0.0509650599086202208,
	0.0416383707883523774,
	0.0313360241670545339,
	0.0203007149001934693,
	0.00880700356957605894
};

constexpr const real *const legendre_points[21] =
{
	nullptr,
	legendre_1_points,
	legendre_2_points,
	legendre_3_points,
	legendre_4_points,
	legendre_5_points,
	legendre_6_points,
	legendre_7_points,
	legendre_8_points,
	legendre_9_points,
	legendre_10_points,
	legendre_11_points,
	legendre_12_points,
	legendre_13_points,
	legendre_14_points,
	legendre_15_points,
	legendre_16_points,
	legendre_17_points,
	legendre_18_points,
	legendre_19_points,
	legendre_20_points
};

constexpr const real *const legendre_weights[21] =
{
	nullptr,
	legendre_1_weights,
	legendre_2_weights,
	legendre_3_weights,
	legendre_4_weights,
	legendre_5_weights,
	legendre_6_weights,
	legendre_7_weights,
	legendre_8_weights,
	legendre_9_weights,
	legendre_10_weights,
	legendre_11_weights,
	legendre_12_weights,
	legendre_13_weights,
	legendre_14_weights,
	legendre_15_weights,
	legendre_16_weights,
	legendre_17_weights,
	legendre_18_weights,
	legendre_19_weights,
	legendre_20_weights
};

}

#endif

