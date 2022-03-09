#!/usr/bin/env python

"""
Model to simulate CO hydrogenation on Rh211 using kmos

Standard settings:
    p_COgas = 6.66
    p_H2gas = 13.33
    T = 523
"""

import numpy as np
from kmos.cli import main as cli_main
from kmos.types import Action, Bystander, Condition, Project, Site, Species
from sys import path
path.append('../..')
from tools.bep_processes import BEPProcessHolder

model_name = 'Rh211_model_without_lateral_interactions'
COMPILE = False

# Project
pt = Project()
pt.set_meta(
    author='Hector Prats & Michael Seibt & Martin Deimel',
    model_name=model_name,
    email='',
    model_dimension=2
)

# Species
pt.add_species(
    Species(name='empty', color='#ffffff'),
    Species(name='CO', color='#00ff00', representation="Atoms('CO',[[0,0,0],[0,0,1.2]])"),
    Species(name='OH'),
    Species(name='O'),
    Species(name='C'),
    Species(name='CH'),
    Species(name='CH2'),
    Species(name='CH3'),
    Species(name='CHO'),
    Species(name='CHOH'),
    Species(name='CHCO'),
    Species(name='CH2CO'),
    Species(name='CH3CO'),
    Species(name='CH3CHO'),
    Species(name='CH3CHOH'),
    Species(name='CH3CH2OH'),
)
pt.species_list.default_species = 'empty'

# Lattice
layer = pt.add_layer(name='Rh211')

layer.add_site(
    Site(name='s', pos='0.121 0.5 0.6', default_species='empty'),
    Site(name='t', pos='0.5 0.0 0.545', default_species='empty'),
    Site(name='f', pos='0.774 0.5 0.53', default_species='empty')
)

pt.lattice.representation = """[
Atoms(symbols='Rh3',
      pbc=np.array([False, False, False], dtype=bool),
      cell=np.array(
          [[  6.582, 0.000,  0.000  ],
           [  0.000, 2.686,  0.000  ],
           [  0.000, 0.000, 20.000  ]]
      ),
      positions=np.array(
          [[  0.000, 0.000, 10.776  ],
           [  2.195, 1.344, 10.000  ],
           [  4.389, 0.000,  9.224  ]]
      )
)
]"""

pt.lattice.cell = np.array(
    [[6.582, 0.000, 0.000],
     [0.000, 2.686, 0.000],
     [0.000, 0.000, 20.00]]
)

# Coordinates
s = pt.lattice.generate_coord('s.(0,0,0).Rh211')
t = pt.lattice.generate_coord('t.(0,0,0).Rh211')
f = pt.lattice.generate_coord('f.(0,0,0).Rh211')

s_N = pt.lattice.generate_coord('s.(0,1,0).Rh211')
s_S = pt.lattice.generate_coord('s.(0,-1,0).Rh211')
s_E = pt.lattice.generate_coord('s.(1,0,0).Rh211')

t_N = pt.lattice.generate_coord('t.(0,1,0).Rh211')
t_S = pt.lattice.generate_coord('t.(0,-1,0).Rh211')
t_W = pt.lattice.generate_coord('t.(-1,0,0).Rh211')
t_NW = pt.lattice.generate_coord('t.(-1,1,0).Rh211')

f_N = pt.lattice.generate_coord('f.(0,1,0).Rh211')
f_S = pt.lattice.generate_coord('f.(0,-1,0).Rh211')
f_W = pt.lattice.generate_coord('f.(-1,0,0).Rh211')

s_bystanders = [
]

t_bystanders = [
]

f_bystanders = [
]

# Parameters
pt.add_parameter(name='T', value=523., adjustable=True, min=450, max=700)
pt.add_parameter(name='p_COgas', value=6.66)
pt.add_parameter(name='p_H2gas', value=13.33)
pt.add_parameter(name='p_CH4gas', value=1.e-20)
pt.add_parameter(name='p_H2Ogas', value=1.e-20)
pt.add_parameter(name='p_CH3CHOgas', value=1.e-20)
pt.add_parameter(name='p_CH3CH2OHgas', value=1.e-20)

# Area of unit cell
pt.add_parameter(name='A', value='6.582*2.686*angstrom**2')

# Gas phase formation energies
pt.add_parameter(name='E_COgas', value='0.00')
pt.add_parameter(name='E_H2Ogas', value='0.00')
pt.add_parameter(name='E_CH4gas', value='-0.51')
pt.add_parameter(name='E_H2gas', value='0.67')
pt.add_parameter(name='E_CH3CHOgas', value='-0.58')
pt.add_parameter(name='E_CH3CH2OHgas', value='-0.65')

# Analytic MF solution to H coverage
pt.add_parameter(name='H_s_cov', value='(1/(1+sqrt(exp(-beta*(GibbsGas_H2gas-2*GibbsAds_H_s)*eV))))')
pt.add_parameter(name='H_t_cov', value='(1/(1+sqrt(exp(-beta*(GibbsGas_H2gas-2*GibbsAds_H_t)*eV))))')

# BEP offset for diffusion
pt.add_parameter(name='E_CO_diff', value='0.11')
pt.add_parameter(name='E_O_diff', value='0.82')
pt.add_parameter(name='E_OH_diff', value='0.44')
pt.add_parameter(name='E_C_diff', value='1.32')
pt.add_parameter(name='E_CH_t_t_diff', value='0.55')
pt.add_parameter(name='E_CH_diff', value='0.37')
pt.add_parameter(name='E_CH2_diff', value='0.61')
pt.add_parameter(name='E_CH3_diff', value='0.51')

# BEP slopes
pt.add_parameter(name='alpha', value='0.5')  # Default BEP slope

# Sensitivity Analysis Parameters
pt.add_parameter(name='sens_C_OH_f', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_C_H_f', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_CO_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_HCO_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_H_f', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH2_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH3_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_O_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_OH_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CHCO_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH2CO_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH3CO_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH3CHO_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH3CHOH_H_s', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_OH_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH2_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH3_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_OH_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH_CO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CHCO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH2CO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_H_CH3CO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH3CHO_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_CH3CHOH_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_O_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)

# Formation energies
pt.add_parameter(name='E_O_t', value='-0.19')
pt.add_parameter(name='E_O_H_t', value='1.13')
pt.add_parameter(name='E_CH_H_f', value='0.36')
pt.add_parameter(name='E_CH_t', value=-0.56)
pt.add_parameter(name='E_CH_CO_t', value=-0.86)
pt.add_parameter(name='E_CH_H_t', value=0.15)
pt.add_parameter(name='E_CH_OH_t', value=0.45)
pt.add_parameter(name='E_CH2_t', value=-0.07)
pt.add_parameter(name='E_CH2CO_t', value=-1.13)
pt.add_parameter(name='E_CH3_t', value=-0.23)
pt.add_parameter(name='E_CH3CH2OH_t', value=-1.13)
pt.add_parameter(name='E_CH3CHO_t', value=-0.99)
pt.add_parameter(name='E_CH3CHO_H_t', value=-0.03)
pt.add_parameter(name='E_CH3CHOH_t', value=-0.86)
pt.add_parameter(name='E_CH3CHOH_H_t', value=-0.18)
pt.add_parameter(name='E_CH3CO_t', value=-1.67)
pt.add_parameter(name='E_CHCO_t', value=-1.31)
pt.add_parameter(name='E_CHO_t', value=-0.7)
pt.add_parameter(name='E_CHOH_t', value=-0.22)
pt.add_parameter(name='E_CO_t', value=-1.7)
pt.add_parameter(name='E_H_t', value=0.02)
pt.add_parameter(name='E_H_CH2_t', value=0.44)
pt.add_parameter(name='E_H_CH2CO_t', value=-0.54)
pt.add_parameter(name='E_H_CH3_t', value=0.44)
pt.add_parameter(name='E_H_CH3CO_t', value=-0.77)
pt.add_parameter(name='E_H_CHCO_t', value=-0.58)
pt.add_parameter(name='E_H_CO_t', value=-0.29)
pt.add_parameter(name='E_H_OH_t', value=0.81)
pt.add_parameter(name='E_HCO_H_t', value=0.61)
pt.add_parameter(name='E_OH_t', value=-0.14)
pt.add_parameter(name='E_C_f', value=-0.48)
pt.add_parameter(name='E_C_H_f', value=0.54)
pt.add_parameter(name='E_C_OH_f', value=0.07)
pt.add_parameter(name='E_CH_f', value=-0.47)
pt.add_parameter(name='E_H_s', value=-0.02)
pt.add_parameter(name='E_CH_CO_s', value=-0.34)
pt.add_parameter(name='E_CH_H_s', value=0.21)
pt.add_parameter(name='E_CH2_s', value=-0.35)
pt.add_parameter(name='E_CH2_H_s', value=0.05)
pt.add_parameter(name='E_CH2CO_s', value=-1.51)
pt.add_parameter(name='E_CH3_s', value=-0.42)
pt.add_parameter(name='E_CH3_H_s', value=0.11)
pt.add_parameter(name='E_CH3CHO_s', value=-1.54)
pt.add_parameter(name='E_CH3CHO_H_s', value=-0.06)
pt.add_parameter(name='E_CH3CHOH_s', value=-1.07)
pt.add_parameter(name='E_CH3CHOH_H_s', value=-0.26)
pt.add_parameter(name='E_CH3CO_s', value=-1.94)
pt.add_parameter(name='E_CHCO_s', value=-1.51)
pt.add_parameter(name='E_CO_s', value=-1.85)
pt.add_parameter(name='E_H_CH2CO_s', value=-0.72)
pt.add_parameter(name='E_H_CH3CO_s', value=-1.05)
pt.add_parameter(name='E_H_CHCO_s', value=-0.68)
pt.add_parameter(name='E_H_OH_s', value=0.56)
pt.add_parameter(name='E_O_s', value=-0.11)
pt.add_parameter(name='E_O_H_s', value=1.06)
pt.add_parameter(name='E_OH_s', value=-0.73)

# Frequencies
pt.add_parameter(name='f_O_t', value='[359.5, 393.3, 507.0]')
pt.add_parameter(name='f_O_H_t', value='[130.2, 179.5, 293.7, 497.8, 975.7]')
pt.add_parameter(name='f_CH_H_f', value='[247.4, 286.7, 460.4, 579.8, 794.5, 958.4, 1325.6, 2936.3]')
pt.add_parameter(name='f_CH_t', value='[413.0, 437.0, 487.0, 710.0, 735.0, 3045.0]')
pt.add_parameter(name='f_CH_CO_t', value='[56.0, 56.0, 150.0, 347.0, 373.0, 460.0, 580.0, 767.0, 961.0, 1607.0, 3001.0]')
pt.add_parameter(name='f_CH_H_t', value='[270.0, 369.0, 382.0, 566.0, 681.0, 955.0, 1998.0, 3025.0]')
pt.add_parameter(name='f_CH_OH_t', value='[194.0, 314.0, 411.0, 461.0, 643.0, 671.0, 723.0, 755.0, 849.0, 3162.0, 3676.0]')
pt.add_parameter(name='f_CH2_t', value='[56.0, 306.0, 382.0, 468.0, 663.0, 790.0, 1356.0, 2737.0, 3004.0]')
pt.add_parameter(name='f_CH2CO_t', value='[56.0, 142.0, 169.0, 291.0, 350.0, 435.0, 565.0, 621.0, 871.0, 930.0, 1035.0, 1295.0, 1674.0, 2998.0, 3071.0]')
pt.add_parameter(name='f_CH3_t', value='[56.0, 114.0, 168.0, 622.0, 686.0, 703.0, 1382.0, 1417.0, 1576.0, 3026.0, 3093.0, 3099.0]')
pt.add_parameter(name='f_CH3CH2OH_t', value='[56.0, 56.0, 56.0, 56.0, 56.0, 255.0, 324.0, 421.0, 479.0, 732.0, 867.0, 1009.0, 1087.0, 1136.0, 1254.0, 1264.0, 1394.0, 1419.0, 1432.0, 1462.0, 1491.0, 2887.0, 2966.0, 2984.0, 3029.0, 3146.0, 3726.0]')
pt.add_parameter(name='f_CH3CHO_t', value='[56.0, 61.0, 93.0, 134.0, 227.0, 361.0, 430.0, 516.0, 737.0, 928.0, 1037.0, 1062.0, 1128.0, 1287.0, 1321.0, 1408.0, 1443.0, 2933.0, 2974.0, 3000.0, 3063.0]')
pt.add_parameter(name='f_CH3CHO_H_t', value='[56.0, 56.0, 123.0, 183.0, 208.0, 232.0, 385.0, 511.0, 578.0, 810.0, 900.0, 999.0, 1054.0, 1084.0, 1129.0, 1316.0, 1335.0, 1424.0, 1453.0, 2906.0, 2945.0, 2984.0, 3054.0]')
pt.add_parameter(name='f_CH3CHOH_t', value='[56.0, 71.0, 162.0, 228.0, 268.0, 318.0, 446.0, 541.0, 664.0, 773.0, 890.0, 1022.0, 1089.0, 1141.0, 1236.0, 1347.0, 1375.0, 1462.0, 1489.0, 2927.0, 2973.0, 3024.0, 3060.0, 3629.0]')
pt.add_parameter(name='f_CH3CHOH_H_t', value='[56.0, 56.0, 108.0, 137.0, 183.0, 249.0, 400.0, 460.0, 516.0, 598.0, 820.0, 924.0, 1014.0, 1096.0, 1167.0, 1240.0, 1354.0, 1385.0, 1436.0, 1485.0, 1652.0, 2967.0, 3012.0, 3039.0, 3102.0, 3525.0]')
pt.add_parameter(name='f_CH3CO_t', value='[56.0, 123.0, 156.0, 228.0, 286.0, 402.0, 469.0, 676.0, 945.0, 1012.0, 1091.0, 1407.0, 1435.0, 1454.0, 1575.0, 3053.0, 3076.0, 3126.0]')
pt.add_parameter(name='f_CHCO_t', value='[56.0, 195.0, 261.0, 342.0, 401.0, 522.0, 572.0, 755.0, 933.0, 1090.0, 1712.0, 3004.0]')
pt.add_parameter(name='f_CHO_t', value='[56.0, 72.0, 206.0, 336.0, 452.0, 660.0, 1238.0, 1262.0, 2832.0]')
pt.add_parameter(name='f_CHOH_t', value='[56.0, 144.0, 202.0, 314.0, 367.0, 569.0, 761.0, 1044.0, 1203.0, 1426.0, 2983.0, 3650.0]')
pt.add_parameter(name='f_CO_t', value='[60.0, 231.0, 256.0, 303.0, 470.0, 1747.0]')
pt.add_parameter(name='f_H_t', value='[463.0, 716.0, 982.0]')
pt.add_parameter(name='f_H_CH2_t', value='[56.0, 234.0, 339.0, 527.0, 615.0, 828.0, 969.0, 1378.0, 1913.0, 2972.0, 3028.0]')
pt.add_parameter(name='f_H_CH2CO_t', value='[56.0, 123.0, 221.0, 245.0, 344.0, 383.0, 484.0, 602.0, 863.0, 994.0, 1030.0, 1088.0, 1332.0, 1433.0, 1684.0, 3026.0, 3087.0]')
pt.add_parameter(name='f_H_CH3_t', value='[56.0, 56.0, 56.0, 143.0, 435.0, 740.0, 826.0, 1241.0, 1398.0, 1421.0, 1731.0, 2994.0, 3046.0, 3093.0]')
pt.add_parameter(name='f_H_CH3CO_t', value='[56.0, 56.0, 56.0, 106.0, 310.0, 367.0, 484.0, 656.0, 750.0, 963.0, 996.0, 1071.0, 1311.0, 1390.0, 1424.0, 1434.0, 1597.0, 2970.0, 2997.0, 3104.0]')
pt.add_parameter(name='f_H_CHCO_t', value='[104.0, 170.0, 250.0, 276.0, 322.0, 484.0, 565.0, 641.0, 863.0, 947.0, 1089.0, 1670.0, 1876.0, 3058.0]')
pt.add_parameter(name='f_H_CO_t', value='[71.0, 197.0, 315.0, 376.0, 537.0, 1047.0, 1445.0, 2114.0]')
pt.add_parameter(name='f_H_OH_t', value='[56.0, 273.0, 312.0, 399.0, 695.0, 831.0, 1111.0, 3615.0]')
pt.add_parameter(name='f_HCO_H_t', value='[137.0, 231.0, 290.0, 327.0, 493.0, 594.0, 680.0, 911.0, 1105.0, 1166.0, 2929.0]')
pt.add_parameter(name='f_OH_t', value='[56.0, 341.0, 396.0, 670.0, 718.0, 3682.0]')
pt.add_parameter(name='f_C_f', value='[303.0, 508.0, 546.0]')
pt.add_parameter(name='f_C_H_f', value='[157.0, 237.0, 505.0, 512.0, 1065.0]')
pt.add_parameter(name='f_C_OH_f', value='[56.0, 56.0, 295.0, 464.0, 500.0, 507.0, 669.0, 3698.0]')
pt.add_parameter(name='f_CH_f', value='[420.0, 428.0, 539.0, 654.0, 677.0, 2979.0]')
pt.add_parameter(name='f_H_s', value='[292.0, 965.0, 1253.0]')
pt.add_parameter(name='f_CH_CO_s', value='[56.0, 56.0, 150.0, 347.0, 373.0, 460.0, 580.0, 767.0, 961.0, 1607.0, 3001.0]')
pt.add_parameter(name='f_CH_H_s', value='[56.0, 380.0, 438.0, 559.0, 657.0, 894.0, 1574.0, 3012.0]')
pt.add_parameter(name='f_CH2_s', value='[56.0, 336.0, 522.0, 545.0, 577.0, 762.0, 1299.0, 2978.0, 3081.0]')
pt.add_parameter(name='f_CH2_H_s', value='[62.1, 278.1, 322.8, 504.8, 588.7, 852.4, 953.2, 1338.6, 1815.2, 3008.8, 3082.0]')
pt.add_parameter(name='f_CH2CO_s', value='[56.0, 56.0, 87.0, 180.0, 390.0, 505.0, 584.0, 644.0, 969.0, 1011.0, 1089.0, 1476.0, 1686.0, 3053.0, 3125.0]')
pt.add_parameter(name='f_CH3_s', value='[56.0, 56.0, 85.0, 442.0, 594.0, 624.0, 1121.0, 1403.0, 1407.0, 2983.0, 3062.0, 3076.0]')
pt.add_parameter(name='f_CH3_H_s', value='[56.0, 156.0, 193.0, 442.0, 648.0, 657.0, 1073.0, 1219.0, 1320.0, 1376.0, 1416.0, 2642.0, 3026.0, 3096.0]')
pt.add_parameter(name='f_CH3CHO_s', value='[56.0, 56.0, 99.0, 183.0, 261.0, 319.0, 390.0, 502.0, 697.0, 928.0, 1030.0, 1109.0, 1210.0, 1375.0, 1404.0, 1451.0, 1519.0, 2972.0, 3014.0, 3062.0, 3113.0]')
pt.add_parameter(name='f_CH3CHO_H_s', value='[56.0, 56.0, 123.0, 183.0, 208.0, 232.0, 385.0, 511.0, 578.0, 810.0, 900.0, 999.0, 1054.0, 1084.0, 1129.0, 1316.0, 1335.0, 1424.0, 1453.0, 2906.0, 2945.0, 2984.0, 3054.0]')
pt.add_parameter(name='f_CH3CHOH_s', value='[56.0, 56.0, 65.0, 165.0, 209.0, 240.0, 374.0, 458.0, 587.0, 837.0, 904.0, 1061.0, 1072.0, 1140.0, 1243.0, 1363.0, 1417.0, 1449.0, 1532.0, 2883.0, 2956.0, 2987.0, 3098.0, 3678.0]')
pt.add_parameter(name='f_CH3CHOH_H_s', value='[56.0, 56.0, 108.0, 137.0, 183.0, 249.0, 400.0, 460.0, 516.0, 598.0, 820.0, 924.0, 1014.0, 1096.0, 1167.0, 1240.0, 1354.0, 1385.0, 1436.0, 1485.0, 1652.0, 2967.0, 3012.0, 3039.0, 3102.0, 3525.0]')
pt.add_parameter(name='f_CH3CO_s', value='[56.0, 56.0, 103.0, 161.0, 262.0, 378.0, 501.0, 593.0, 900.0, 1021.0, 1100.0, 1301.0, 1428.0, 1445.0, 1476.0, 3040.0, 3074.0, 3100.0]')
pt.add_parameter(name='f_CHCO_s', value='[56.0, 203.0, 238.0, 327.0, 377.0, 528.0, 565.0, 669.0, 1000.0, 1134.0, 2086.0, 3128.0]')
pt.add_parameter(name='f_CO_s', value='[56.0, 245.0, 305.0, 332.0, 383.0, 1842.0]')
pt.add_parameter(name='f_H_CH2CO_s', value='[56.0, 123.0, 221.0, 245.0, 344.0, 383.0, 484.0, 602.0, 863.0, 994.0, 1030.0, 1088.0, 1332.0, 1433.0, 1684.0, 3026.0, 3087.0]')
pt.add_parameter(name='f_H_CH3CO_s', value='[56.0, 56.0, 56.0, 106.0, 310.0, 367.0, 484.0, 656.0, 750.0, 963.0, 996.0, 1071.0, 1311.0, 1390.0, 1424.0, 1434.0, 1597.0, 2970.0, 2997.0, 3104.0]')
pt.add_parameter(name='f_H_CHCO_s', value='[104.0, 170.0, 250.0, 276.0, 322.0, 484.0, 565.0, 641.0, 863.0, 947.0, 1089.0, 1670.0, 1876.0, 3058.0]')
pt.add_parameter(name='f_H_OH_s', value='[56.0, 185.0, 356.0, 409.0, 464.0, 561.0, 1052.0, 3642.0]')
pt.add_parameter(name='f_O_s', value='[288.0, 332.0, 441.0]')
pt.add_parameter(name='f_O_H_s', value='[106.0, 218.0, 415.0, 429.0, 938.0]')
pt.add_parameter(name='f_OH_s', value='[116.0, 214.0, 383.0, 621.0, 683.0, 3719.0]')

process_holder = BEPProcessHolder()
process_holder.add_site_bystanders('s', s, s_bystanders)
process_holder.add_site_bystanders('t', t, t_bystanders)
process_holder.add_site_bystanders('f', f, f_bystanders)
process_holder.add_param_list(pt.parameter_list)

# Processes
for site in [t, s, ]:
    process_holder.add_process(
        name='CO_ads_{}'.format(site.name),
        rate_constant='p_COgas*bar*A/2/sqrt(2*pi*umass*m_CO/beta)',
        condition_list=[Condition(coord=site, species='empty')],
        action_list=[Action(coord=site, species='CO')],
    )
    process_holder.add_process(
        name='CO_des_{}'.format(site.name),
        rate_constant='p_COgas*bar*A/2/sqrt(2*pi*umass*m_CO/beta)*exp(-beta*(GibbsGas_COgas-GibbsAds_CO_{})*eV)'.format(site.name),
        condition_list=[Condition(coord=site, species='CO')],
        action_list=[Action(coord=site, species='empty')],
    )

for site, name in [(t_W, 't_W'), (t_NW, 't_NW')]:
    process_holder.add_process(
        name='H_CO_s_react_{}'.format(name),
        rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_C_OH_f+sens_C_OH_f,GibbsAds_C_f+GibbsAds_OH_t)-(GibbsAds_CO_s+GibbsAds_H_t),0)*eV)',
        condition_list=[Condition(coord=s, species='CO'), Condition(coord=f_W, species='empty'), Condition(coord=site, species='empty')],
        action_list=[Action(coord=s, species='empty'), Action(coord=f_W, species='C'), Action(coord=site, species='OH')],
    )
    process_holder.add_process(
        name='C_OH_s_react_{}'.format(name),
        rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_C_OH_f+sens_C_OH_f,GibbsAds_CO_s+GibbsAds_H_t)-(GibbsAds_C_f+GibbsAds_OH_t),0)*eV)',
        condition_list=[Condition(coord=s, species='empty'), Condition(coord=f_W, species='C'), Condition(coord=site, species='OH')],
        action_list=[Action(coord=s, species='CO'), Action(coord=f_W, species='empty'), Action(coord=site, species='empty')],
    )

for site, name in [(f, 'f'), (f_S, 'f_S')]:
    process_holder.add_process(
        name='H_CO_t_react_{}'.format(name),
        rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_C_OH_f+sens_C_OH_f,GibbsAds_C_f+GibbsAds_OH_t)-(GibbsAds_CO_t+GibbsAds_H_t),0)*eV)',
        condition_list=[Condition(coord=t, species='CO'), Condition(coord=site, species='empty')],
        action_list=[Action(coord=t, species='OH'), Action(coord=site, species='C')],
    )
    process_holder.add_process(
        name='C_OH_t_react_{}'.format(name),
        rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_C_OH_f+sens_C_OH_f,GibbsAds_CO_t+GibbsAds_H_t)-(GibbsAds_C_f+GibbsAds_OH_t),0)*eV)',
        condition_list=[Condition(coord=site, species='C'), Condition(coord=t, species='OH')],
        action_list=[Action(coord=site, species='empty'), Action(coord=t, species='CO')],
    )

process_holder.add_process(
    name='H_C_f_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_C_H_f+sens_C_H_f,GibbsAds_CH_f)-(GibbsAds_C_f+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=f, species='C')],
    action_list=[Action(coord=f, species='CH')],
)
process_holder.add_process(
    name='CH_f_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_C_H_f+sens_C_H_f,GibbsAds_C_f+GibbsAds_H_s)-(GibbsAds_CH_f),0)*eV)',
    condition_list=[Condition(coord=f, species='CH')],
    action_list=[Action(coord=f, species='C')],
)

process_holder.add_process(
    name='H_CH_f_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_s+sens_CH_H_s,GibbsAds_CH2_s)-(GibbsAds_CH_f+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=f, species='CH'), Condition(coord=s_E, species='empty')],
    action_list=[Action(coord=f, species='empty'), Action(coord=s_E, species='CH2')],
)
process_holder.add_process(
    name='CH2_s_react',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_s+sens_CH_H_s,GibbsAds_CH_f+GibbsAds_H_s)-(GibbsAds_CH2_s),0)*eV)',
    condition_list=[Condition(coord=f, species='empty'), Condition(coord=s_E, species='CH2')],
    action_list=[Action(coord=f, species='CH'), Action(coord=s_E, species='empty')],
)

for site, name in [(t, 't'), (t_N, 't_N')]:
    process_holder.add_process(
        name='H_CH_f_{}_react'.format(name),
        rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_f+sens_CH_H_f,GibbsAds_CH2_t)-(GibbsAds_CH_f+GibbsAds_H_s),0)*eV)',
        condition_list=[Condition(coord=f, species='CH'), Condition(coord=site, species='empty')],
        action_list=[Action(coord=f, species='empty'), Action(coord=site, species='CH2')],
    )
    process_holder.add_process(
        name='CH2_{}_f_dis'.format(name),
        rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_f+sens_CH_H_f,GibbsAds_CH_f+GibbsAds_H_s)-GibbsAds_CH2_t,0)*eV)',
        condition_list=[Condition(coord=f, species='empty'), Condition(coord=site, species='CH2')],
        action_list=[Action(coord=f, species='CH'), Action(coord=site, species='empty')],
    )

process_holder.add_process(
    name='H_CH2_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH2_H_s+sens_CH2_H_s,GibbsAds_CH3_s)-(GibbsAds_CH2_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH2')],
    action_list=[Action(coord=s, species='CH3')],
)
process_holder.add_process(
    name='CH3_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH2_H_s+sens_CH2_H_s,GibbsAds_CH2_s+GibbsAds_H_s)-(GibbsAds_CH3_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3')],
    action_list=[Action(coord=s, species='CH2')],
)

process_holder.add_process(
    name='H_CH3_s_react',
    rate_constant='H_s_cov/(2*beta*h)*exp(-beta*max(max(GibbsAds_CH3_H_s+sens_CH3_H_s,GibbsGas_CH4gas)-(GibbsAds_CH3_s+GibbsAds_H_s),0)*eV)',
    tof_count={'CH4_formation': 1},
    condition_list=[Condition(coord=s, species='CH3')],
    action_list=[Action(coord=s, species='empty')],
)
process_holder.add_process(
    name='CH4_s_ads',
    rate_constant='(1-H_s_cov)/(2*beta*h)*exp(-beta*max(max(GibbsAds_CH3_H_s+sens_CH3_H_s,GibbsAds_CH3_s+GibbsAds_H_s)-(GibbsGas_CH4gas),0)*eV)',
    tof_count={'CH4_formation': -1},
    condition_list=[Condition(coord=s, species='empty')],
    action_list=[Action(coord=s, species='CH3')],
)

process_holder.add_process(
    name='H_O_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_O_H_s+sens_O_H_s,GibbsAds_OH_s)-(GibbsAds_O_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='O')],
    action_list=[Action(coord=s, species='OH')],
)
process_holder.add_process(
    name='OH_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_O_H_s+sens_O_H_s,GibbsAds_O_s+GibbsAds_H_s)-(GibbsAds_OH_s),0)*eV)',
    condition_list=[Condition(coord=s, species='OH')],
    action_list=[Action(coord=s, species='O')],
)

process_holder.add_process(
    name='H_OH_s_react',
    rate_constant='H_s_cov/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_OH_s+sens_H_OH_s,GibbsGas_H2Ogas)-(GibbsAds_OH_s+GibbsAds_H_s),0)*eV)',
    tof_count={'H2O_formation': 1},
    condition_list=[Condition(coord=s, species='OH')],
    action_list=[Action(coord=s, species='empty')],
)
process_holder.add_process(
    name='H2O_s_ads',
    rate_constant='(1-H_s_cov)/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_OH_s+sens_H_OH_s,GibbsAds_OH_s+GibbsAds_H_s)-(GibbsGas_H2Ogas),0)*eV)',
    tof_count={'H2O_formation': -1},
    condition_list=[Condition(coord=s, species='empty')],
    action_list=[Action(coord=s, species='OH')],
)

process_holder.add_process(
    name='CH_CO_s_react',
    rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_s+sens_CH_CO_s,GibbsAds_CHCO_s)-(GibbsAds_CH_f+GibbsAds_CO_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CO'), Condition(coord=f_W, species='CH')],
    action_list=[Action(coord=s, species='CHCO'), Action(coord=f_W, species='empty')],
)
process_holder.add_process(
    name='CHCO_s_dis',
    rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_s+sens_CH_CO_s,GibbsAds_CH_f+GibbsAds_CO_s)-(GibbsAds_CHCO_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CHCO'), Condition(coord=f_W, species='empty')],
    action_list=[Action(coord=s, species='CO'), Action(coord=f_W, species='CH')],
)

process_holder.add_process(
    name='H_CHCO_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_s+sens_H_CHCO_s,GibbsAds_CH2CO_s)-(GibbsAds_CHCO_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CHCO')],
    action_list=[Action(coord=s, species='CH2CO')],
)
process_holder.add_process(
    name='CH2CO_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_s+sens_H_CHCO_s,GibbsAds_CHCO_s+GibbsAds_H_s)-(GibbsAds_CH2CO_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH2CO')],
    action_list=[Action(coord=s, species='CHCO')],
)

process_holder.add_process(
    name='H_CH2CO_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_s+sens_H_CH2CO_s,GibbsAds_CH3CO_s)-(GibbsAds_CH2CO_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH2CO')],
    action_list=[Action(coord=s, species='CH3CO')],
)
process_holder.add_process(
    name='CH3CO_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_s+sens_H_CH2CO_s,GibbsAds_CH2CO_s+GibbsAds_H_s)-(GibbsAds_CH3CO_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3CO')],
    action_list=[Action(coord=s, species='CH2CO')],
)

process_holder.add_process(
    name='H_CH3CO_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_s+sens_H_CH3CO_s,GibbsAds_CH3CHO_s)-(GibbsAds_CH3CO_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3CO')],
    action_list=[Action(coord=s, species='CH3CHO')],
)
process_holder.add_process(
    name='CH3CHO_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_s+sens_H_CH3CO_s,GibbsAds_CH3CO_s+GibbsAds_H_s)-(GibbsAds_CH3CHO_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3CHO')],
    action_list=[Action(coord=s, species='CH3CO')],
)

process_holder.add_process(
    name='H_CH3CHO_s_react',
    rate_constant='H_s_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_s+sens_CH3CHO_H_s,GibbsAds_CH3CHOH_s)-(GibbsAds_CH3CHO_s+GibbsAds_H_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3CHO')],
    action_list=[Action(coord=s, species='CH3CHOH')],
)
process_holder.add_process(
    name='CH3CHOH_s_dis',
    rate_constant='(1-H_s_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_s+sens_CH3CHO_H_s,GibbsAds_CH3CHO_s+GibbsAds_H_s)-(GibbsAds_CH3CHOH_s),0)*eV)',
    condition_list=[Condition(coord=s, species='CH3CHOH')],
    action_list=[Action(coord=s, species='CH3CHO')],
)

for site in [s, t, ]:
    process_holder.add_process(
        name='CH3CHO_des_{}'.format(site.name),
        rate_constant='p_CH3CHOgas*bar*A/2/sqrt(2*pi*umass*m_CH3CHO/beta)*exp(-beta*(GibbsGas_CH3CHOgas-GibbsAds_CH3CHO_{})*eV)'.format(site.name),
        tof_count={'CH3CHO_formation': 1},
        condition_list=[Condition(coord=site, species='CH3CHO')],
        action_list=[Action(coord=site, species='empty')],
    )
    process_holder.add_process(
        name='CH3CHO_ads_{}'.format(site.name),
        rate_constant='p_CH3CHOgas*bar*A/2/sqrt(2*pi*umass*m_CH3CHO/beta)',
        tof_count={'CH3CHO_formation': -1},
        condition_list=[Condition(coord=site, species='empty')],
        action_list=[Action(coord=site, species='CH3CHO')],
    )

process_holder.add_process(
    name='H_CH3CHOH_s_react',
    rate_constant='H_s_cov/(2*beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_s+sens_CH3CHOH_H_s,GibbsGas_CH3CH2OHgas)-(GibbsAds_CH3CHOH_s+GibbsAds_H_s),0)*eV)',
    tof_count={'CH3CH2OH_formation': 1},
    condition_list=[Condition(coord=s, species='CH3CHOH')],
    action_list=[Action(coord=s, species='empty')],
)

process_holder.add_process(
    name='CH3CH2OH_s_ads',
    rate_constant='(1-H_s_cov)/(2*beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_s+sens_CH3CHOH_H_s,GibbsAds_CH3CHOH_s+GibbsAds_H_s)-(GibbsGas_CH3CH2OHgas),0)*eV)',
    tof_count={'CH3CH2OH_formation': -1},
    condition_list=[Condition(coord=s, species='empty')],
    action_list=[Action(coord=s, species='CH3CHOH')],
)

process_holder.add_process(
    name='H_CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CO_t+sens_H_CO_t,GibbsAds_CHO_t)-(GibbsAds_CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CO')],
    action_list=[Action(coord=t, species='CHO')],
)
process_holder.add_process(
    name='CHO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CO_t+sens_H_CO_t,GibbsAds_CO_t+GibbsAds_H_t)-(GibbsAds_CHO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHO')],
    action_list=[Action(coord=t, species='CO')],
)

process_holder.add_process(
    name='H_CHO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_HCO_H_t+sens_HCO_H_t,GibbsAds_CHOH_t)-(GibbsAds_CHO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHO')],
    action_list=[Action(coord=t, species='CHOH')],
)
process_holder.add_process(
    name='CHOH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_HCO_H_t+sens_HCO_H_t,GibbsAds_CHO_t+GibbsAds_H_t)-(GibbsAds_CHOH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHOH')],
    action_list=[Action(coord=t, species='CHO')],
)

for (site1, name1), (site2, name2) in [((t, 't'), (t_N, 't_N')), ((t, 't'), (t_S, 't_S'))]:
    process_holder.add_process(
        name='CHOH_{}_{}_dis'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CH_t+GibbsAds_OH_t)-(GibbsAds_CHOH_t),0)*eV)',
        condition_list=[Condition(coord=site1, species='CHOH'), Condition(coord=site2, species='empty')],
        action_list=[Action(coord=site1, species='CH'), Action(coord=site2, species='OH')],
    )
    process_holder.add_process(
        name='CH_OH_{}_{}_react'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CHOH_t)-(GibbsAds_CH_t+GibbsAds_OH_t),0)*eV)',
        condition_list=[Condition(coord=site1, species='CH'), Condition(coord=site2, species='OH')],
        action_list=[Action(coord=site1, species='CHOH'), Action(coord=site2, species='empty')],
    )

process_holder.add_process(
    name='H_CH_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_t+sens_CH_H_t,GibbsAds_CH2_t)-(GibbsAds_CH_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH')],
    action_list=[Action(coord=t, species='CH2')],
)
process_holder.add_process(
    name='CH2_t_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_t+sens_CH_H_t,GibbsAds_CH_t+GibbsAds_H_t)-(GibbsAds_CH2_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2')],
    action_list=[Action(coord=t, species='CH')],
)

process_holder.add_process(
    name='H_CH2_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2_t+sens_H_CH2_t,GibbsAds_CH3_t)-(GibbsAds_CH2_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2')],
    action_list=[Action(coord=t, species='CH3')],
)
process_holder.add_process(
    name='CH3_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2_t+sens_H_CH2_t,GibbsAds_CH2_t+GibbsAds_H_t)-(GibbsAds_CH3_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3')],
    action_list=[Action(coord=t, species='CH2')],
)

process_holder.add_process(
    name='H_CH3_t_react',
    rate_constant='H_t_cov/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_CH3_t+sens_H_CH3_t,GibbsGas_CH4gas)-(GibbsAds_CH3_t+GibbsAds_H_t),0)*eV)',
    tof_count={'CH4_formation': 1},
    condition_list=[Condition(coord=t, species='CH3')],
    action_list=[Action(coord=t, species='empty')],
)

process_holder.add_process(
    name='CH4_t_ads',
    rate_constant='(1-H_t_cov)/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_CH3_t+sens_H_CH3_t,GibbsAds_CH3_t+GibbsAds_H_t)-(GibbsGas_CH4gas),0)*eV)',
    tof_count={'CH4_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CH3')],
)

process_holder.add_process(
    name='H_OH_t_react',
    rate_constant='H_t_cov/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_OH_t+sens_H_OH_t,GibbsGas_H2Ogas)-(GibbsAds_OH_t+GibbsAds_H_t),0)*eV)',
    tof_count={'H2O_formation': 1},
    condition_list=[Condition(coord=t, species='OH')],
    action_list=[Action(coord=t, species='empty')],
)
process_holder.add_process(
    name='H2O_t_ads',
    rate_constant='(1-H_t_cov)/(2*beta*h)*exp(-beta*max(max(GibbsAds_H_OH_t+sens_H_OH_t,GibbsAds_OH_t+GibbsAds_H_t)-(GibbsGas_H2Ogas),0)*eV)',
    tof_count={'H2O_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='OH')],
)

for (site1, name1), (site2, name2) in [((t, 't'), (t_N, 't_N')), ((t, 't'), (t_S, 't_S'))]:
    process_holder.add_process(
        name='CH_CO_{}_{}_react'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CHCO_t)-(GibbsAds_CH_t+GibbsAds_CO_t),0)*eV)',
        condition_list=[Condition(coord=site1, species='CH'), Condition(coord=site2, species='CO')],
        action_list=[Action(coord=site1, species='CHCO'), Action(coord=site2, species='empty')],
    )
    process_holder.add_process(
        name='CHCO_{}_{}_dis'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CH_t+GibbsAds_CO_t)-(GibbsAds_CHCO_t),0)*eV)',
        condition_list=[Condition(coord=site1, species='CHCO'), Condition(coord=site2, species='empty')],
        action_list=[Action(coord=site1, species='CH'), Action(coord=site2, species='CO')],
    )

process_holder.add_process(
    name='H_CHCO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_t+sens_H_CHCO_t,GibbsAds_CH2CO_t)-(GibbsAds_CHCO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHCO')],
    action_list=[Action(coord=t, species='CH2CO')],
)
process_holder.add_process(
    name='CH2CO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_t+sens_H_CHCO_t,GibbsAds_CHCO_t+GibbsAds_H_t)-(GibbsAds_CH2CO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2CO')],
    action_list=[Action(coord=t, species='CHCO')],
)

process_holder.add_process(
    name='H_CH2CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_t+sens_H_CH2CO_t,GibbsAds_CH3CO_t)-(GibbsAds_CH2CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2CO')],
    action_list=[Action(coord=t, species='CH3CO')],
)
process_holder.add_process(
    name='CH3CO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_t+sens_H_CH2CO_t,GibbsAds_CH2CO_t+GibbsAds_H_t)-(GibbsAds_CH3CO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CO')],
    action_list=[Action(coord=t, species='CH2CO')],
)

process_holder.add_process(
    name='H_CH3CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_t+sens_H_CH3CO_t,GibbsAds_CH3CHO_t)-(GibbsAds_CH3CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CO')],
    action_list=[Action(coord=t, species='CH3CHO')],
)
process_holder.add_process(
    name='CH3CHO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_t+sens_H_CH3CO_t,GibbsAds_CH3CO_t+GibbsAds_H_t)-(GibbsAds_CH3CHO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHO')],
    action_list=[Action(coord=t, species='CH3CO')],
)

process_holder.add_process(
    name='H_CH3CHO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_t+sens_CH3CHO_H_t,GibbsAds_CH3CHOH_t)-(GibbsAds_CH3CHO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHO')],
    action_list=[Action(coord=t, species='CH3CHOH')],
)
process_holder.add_process(
    name='CH3CHOH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_t+sens_CH3CHO_H_t,GibbsAds_CH3CHO_t+GibbsAds_H_t)-(GibbsAds_CH3CHOH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHOH')],
    action_list=[Action(coord=t, species='CH3CHO')],
)

process_holder.add_process(
    name='H_CH3CHOH_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_t+sens_CH3CHOH_H_t,GibbsAds_CH3CH2OH_t)-(GibbsAds_CH3CHOH_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHOH')],
    action_list=[Action(coord=t, species='CH3CH2OH')],
)
process_holder.add_process(
    name='CH3CH2OH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_t+sens_CH3CHOH_H_t,GibbsAds_CH3CHOH_t+GibbsAds_H_t)-(GibbsAds_CH3CH2OH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CH2OH')],
    action_list=[Action(coord=t, species='CH3CHOH')],
)

process_holder.add_process(
    name='CH3CH2OH_t_des',
    rate_constant='p_CH3CH2OHgas*bar*A/2/sqrt(2*pi*umass*m_CH3CH2OH/beta)*exp(-beta*(GibbsGas_CH3CH2OHgas-GibbsAds_CH3CH2OH_t)*eV)',
    tof_count={'CH3CH2OH_formation': 1},
    condition_list=[Condition(coord=t, species='CH3CH2OH')],
    action_list=[Action(coord=t, species='empty')],
)
process_holder.add_process(
    name='CH3CH2OH_t_ads',
    rate_constant='p_CH3CH2OHgas*bar*A/2/sqrt(2*pi*umass*m_CH3CH2OH/beta)',
    tof_count={'CH3CH2OH_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CH3CH2OH')],
)

process_holder.add_process(
    name='OH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_O_H_t+sens_O_H_t,GibbsAds_O_t+GibbsAds_H_t)-GibbsAds_OH_t,0)*eV)',
    condition_list=[Condition(coord=t, species='OH')],
    action_list=[Action(coord=t, species='O')],
)
process_holder.add_process(
    name='H_O_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_O_H_t+sens_O_H_t,GibbsAds_OH_t)-(GibbsAds_O_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='O')],
    action_list=[Action(coord=t, species='OH')],
)

for ads in ['CO', 'O', 'OH', 'CH2', 'CH3']:
    for (site1, name1), (site2, name2) in [
            ((s, 's'), (s_N, 's_N')), ((s_N, 's_N'), (s, 's')), ((t, 't'), (t_S, 't_S')), ((t_S, 't_S'), (t, 't')),
            ((s, 's'), (t, 't')), ((s, 's'), (t_N, 't_N')), ((t, 't'), (s, 's')), ((t_N, 't_N'), (s, 's'))
    ]:
        process_holder.add_process(
            name='{}_diff_{}_{}'.format(ads, name1, name2),
            rate_constant='1/(beta*h)*exp(-beta*max((alpha*(GibbsAds_{}_{}-GibbsAds_{}_{})+E_{}_diff),0)*eV)'.format(ads, site2.name, ads, site1.name, ads),
            condition_list=[Condition(coord=site1, species=ads), Condition(coord=site2, species='empty')],
            action_list=[Action(coord=site1, species='empty'), Action(coord=site2, species=ads)]
        )

for (site1, name1), (site2, name2) in [((f, 'f'), (f_S, 'f_S')), ((f_S, 'f_S'), (f, 'f'))]:
    process_holder.add_process(
        name='C_diff_{}_{}'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max((alpha*(GibbsAds_C_{}-GibbsAds_C_{})+E_C_diff),0)*eV)'.format(site2.name, site1.name),
        condition_list=[Condition(coord=site1, species='C'), Condition(coord=site2, species='empty')],
        action_list=[Action(coord=site1, species='empty'), Action(coord=site2, species='C')]
    )

for (site1, name1), (site2, name2) in [((t, 't'), (t_S, 't_S')), ((t_S, 't_S'), (t, 't'))]:
    process_holder.add_process(
        name='CH_diff_{}_{}'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max((alpha*(GibbsAds_CH_{}-GibbsAds_CH_{})+E_CH_t_t_diff),0)*eV)'.format(site2.name, site1.name),
        condition_list=[Condition(coord=site1, species='CH'), Condition(coord=site2, species='empty')],
        action_list=[Action(coord=site1, species='empty'), Action(coord=site2, species='CH')]
    )

for (site1, name1), (site2, name2) in [
        ((f, 'f'), (f_S, 'f_S')), ((f_S, 'f_S'), (f, 'f')), ((f, 'f'), (t, 't')),
        ((f, 'f'), (t_N, 't_N')), ((t, 't'), (f, 'f')), ((t_N, 't_N'), (f, 'f'))
]:
    process_holder.add_process(
        name='CH_diff_{}_{}'.format(name1, name2),
        rate_constant='1/(beta*h)*exp(-beta*max((alpha*(GibbsAds_CH_{}-GibbsAds_CH_{})+E_CH_diff),0)*eV)'.format(site2.name, site1.name),
        condition_list=[Condition(coord=site1, species='CH'), Condition(coord=site2, species='empty')],
        action_list=[Action(coord=site1, species='empty'), Action(coord=site2, species='CH')]
    )

pt = process_holder.add_project_processes(pt)

# Export
pt.export_xml_file('{}.xml'.format(model_name))

if COMPILE:
    cli_main('export {}.xml {} -b otf -t'.format(model_name, model_name))
