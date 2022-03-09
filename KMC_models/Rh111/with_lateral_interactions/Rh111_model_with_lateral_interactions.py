#!/usr/bin/env python

"""
Model to simulate CO hydrogenation on Rh111 using kmos

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

model_name = 'Rh111_model_with_lateral_interactions'
COMPILE = False

# Project
pt = Project()
pt.set_meta(
    author='Martin Deimel & Michael Seibt',
    model_name=model_name,
    email='',
    model_dimension=2
)

# Species
pt.add_species(
    Species(name='empty', color='#ffffff'),
    Species(name='CO', color='#00ff00', representation="Atoms('CO',[[0,0,0],[0,0,1.2]])"),
    Species(name='OH'),
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
layer = pt.add_layer(name='Rh111')

layer.add_site(
    Site(name='t', pos='1.0 0.67 0.53', default_species='empty'),
)

pt.lattice.representation = """[
Atoms(symbols='Rh1',
      pbc=np.array([False, False, False], dtype=bool),
      cell=np.array(
          [[  2.729, 0.000,  0.000  ],
           [  1.365, 2.364,  0.000  ],
           [  0.000, 0.000, 10.000  ]]
      ),
      positions=np.array(
          [[  0.000, 0.000,  5.000  ]]
      )
)
]"""

pt.lattice.cell = np.array(
    [[2.729, 0.000, 0.000],
     [1.365, 2.364, 0.000],
     [0.000, 0.000, 10.00]]
)

# Coordinates
t = pt.lattice.generate_coord('t.(0,0,0).Rh111')
t_E = pt.lattice.generate_coord('t.(1,0,0).Rh111')
t_W = pt.lattice.generate_coord('t.(-1,0,0).Rh111')
t_NE = pt.lattice.generate_coord('t.(0,1,0).Rh111')
t_NW = pt.lattice.generate_coord('t.(-1,1,0).Rh111')
t_SE = pt.lattice.generate_coord('t.(0,-1,0).Rh111')
t_SW = pt.lattice.generate_coord('t.(1,-1,0).Rh111')
t_bystanders = [
    Bystander(coord=t_NE, allowed_species=[], flag='nn_t'),
    Bystander(coord=t_NW, allowed_species=[], flag='nn_t'),
    Bystander(coord=t_E, allowed_species=[], flag='nn_t'),
    Bystander(coord=t_W, allowed_species=[], flag='nn_t'),
    Bystander(coord=t_SE, allowed_species=[], flag='nn_t'),
    Bystander(coord=t_SW, allowed_species=[], flag='nn_t'),
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
area = np.linalg.norm(np.cross([  2.729, 0.000,  0.000  ], [  1.365, 2.364,  0.000  ]))
pt.add_parameter(name='A', value='{}*angstrom**2'.format(area))

# Gas phase formation energies
pt.add_parameter(name='E_COgas', value='0.00')
pt.add_parameter(name='E_H2Ogas', value='0.00')
pt.add_parameter(name='E_CH4gas', value='-0.51')
pt.add_parameter(name='E_H2gas', value='0.67')
pt.add_parameter(name='E_CH3CHOgas', value='-0.58')
pt.add_parameter(name='E_CH3CH2OHgas', value='-0.65')

# Analytic MF solution to H coverage
pt.add_parameter(name='H_t_cov', value='(1/(1+sqrt(exp(-beta*(GibbsGas_H2gas-2*GibbsAds_H_t)*eV))))')

# BEP offset for diffusion
pt.add_parameter(name='E_CO_diff', value='0.11')
pt.add_parameter(name='E_OH_diff', value='0.44')
pt.add_parameter(name='E_CH_diff', value='0.37')
pt.add_parameter(name='E_CH2_diff', value='0.61')
pt.add_parameter(name='E_CH3_diff', value='0.51')

# BEP slopes
pt.add_parameter(name='alpha', value='0.5')  # Default BEP slope
pt.add_parameter(name='alpha_ads', value='0.0')
pt.add_parameter(name='alpha_CO_to_C_OH', value='0.792')
pt.add_parameter(name='alpha_CH2_to_CH_H', value='0.761')
pt.add_parameter(name='alpha_CH3_to_CH2_H', value='0.757')
pt.add_parameter(name='alpha_CH4_to_CH3_H', value='0.789')
pt.add_parameter(name='alpha_H2O_to_OH_H', value='0.650')
pt.add_parameter(name='alpha_CO_CH_to_CHCO', value='0.204')
pt.add_parameter(name='alpha_CO_H_to_CHO', value='0.562')
pt.add_parameter(name='alpha_avg_hydrogenation_C', value='0.264')
pt.add_parameter(name='alpha_avg_hydrogenation_O', value='0.456')

# Sensitivity Analysis Parameters
pt.add_parameter(name='sens_H_CO_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
pt.add_parameter(name='sens_HCO_H_t', value=0.00, adjustable=True, min=-0.4, max=0.4)
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

# Formation energies
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

# Frequencies
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

# Interaction parameters
pt.add_parameter(name='I_CO_t_CO_t', value=0.19)
pt.add_parameter(name='I_CH2CO_t_CO_t', value=0.15)
pt.add_parameter(name='I_CH3CHOH_t_CO_t', value=0.12)
pt.add_parameter(name='I_CH3_t_CO_t', value=0.10)
pt.add_parameter(name='I_CHOH_t_CO_t', value=0.09)
pt.add_parameter(name='I_CHCO_t_CO_t', value=0.09)
pt.add_parameter(name='I_CH_t_CO_t', value=0.07)
pt.add_parameter(name='I_OH_t_CO_t', value=0.06)
pt.add_parameter(name='I_CHO_t_CO_t', value=0.05)
pt.add_parameter(name='I_CH3CO_t_CO_t', value=0.01)
pt.add_parameter(name='I_CH3CHO_t_CO_t', value=-0.08)
pt.add_parameter(name='I_CH3CH2OH_t_CO_t', value=-0.09)

process_holder = BEPProcessHolder()
process_holder.add_site_bystanders('t', t, t_bystanders)
process_holder.add_param_list(pt.parameter_list)
process_holder.interaction_energy_pattern = r'I_([A-Z0-9]+)_t_{species}_t|I_{species}_t_([A-Z0-9]+)_t'
process_holder.self_interaction_energy_pattern = r'I_{c1_species}_t_{c2_species}_t'

# Processes
process_holder.add_process(
    name='CO_ads_t',
    rate_constant='p_COgas*bar*A/sqrt(2*pi*umass*m_CO/beta)',
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CO')],
    alpha='alpha_ads',
)
process_holder.add_process(
    name='CO_des_t',
    rate_constant='p_COgas*bar*A/sqrt(2*pi*umass*m_CO/beta)*exp(-beta*(GibbsGas_COgas-GibbsAds_CO_t)*eV)',
    condition_list=[Condition(coord=t, species='CO')],
    action_list=[Action(coord=t, species='empty')],
    alpha='rev_alpha_ads',
)

process_holder.add_process(
    name='CH3CHO_ads_t',
    rate_constant='p_CH3CHOgas*bar*A/sqrt(2*pi*umass*m_CH3CHO/beta)',
    tof_count={'CH3CHO_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CH3CHO')],
    alpha='alpha_ads',
)
process_holder.add_process(
    name='CH3CHO_des_t',
    rate_constant='p_CH3CHOgas*bar*A/sqrt(2*pi*umass*m_CH3CHO/beta)*exp(-beta*(GibbsGas_CH3CHOgas-GibbsAds_CH3CHO_t)*eV)',
    tof_count={'CH3CHO_formation': 1},
    condition_list=[Condition(coord=t, species='CH3CHO')],
    action_list=[Action(coord=t, species='empty')],
    alpha='rev_alpha_ads',
)

process_holder.add_process(
    name='CH3CH2OH_t_ads',
    rate_constant='p_CH3CH2OHgas*bar*A/sqrt(2*pi*umass*m_CH3CH2OH/beta)',
    tof_count={'CH3CH2OH_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CH3CH2OH')],
    alpha='alpha_ads',
)
process_holder.add_process(
    name='CH3CH2OH_t_des',
    rate_constant='p_CH3CH2OHgas*bar*A/sqrt(2*pi*umass*m_CH3CH2OH/beta)*exp(-beta*(GibbsGas_CH3CH2OHgas-GibbsAds_CH3CH2OH_t)*eV)',
    tof_count={'CH3CH2OH_formation': 1},
    condition_list=[Condition(coord=t, species='CH3CH2OH')],
    action_list=[Action(coord=t, species='empty')],
    alpha='rev_alpha_ads',
)
for ads in ['CO', 'CH', 'OH', 'CH2', 'CH3']:
    for site1, site2 in [('t', 't_NW'), ('t', 't_NE'), ('t', 't_E')]:
        process_holder.add_process(
            name='{}_diff_{}_{}'.format(ads, site1, site2),
            rate_constant='1/(beta*h)*exp(-beta*E_{}_diff*eV)'.format(ads),
            condition_list=[Condition(coord=locals()[site1], species=ads), Condition(coord=locals()[site2], species='empty')],
            action_list=[Action(coord=locals()[site1], species='empty'), Action(coord=locals()[site2], species=ads)]
        )
        process_holder.add_process(
            name='{}_diff_{}_{}'.format(ads, site2, site1),
            rate_constant='1/(beta*h)*exp(-beta*E_{}_diff*eV)'.format(ads),
            condition_list=[Condition(coord=locals()[site2], species=ads), Condition(coord=locals()[site1], species='empty')],
            action_list=[Action(coord=locals()[site2], species='empty'), Action(coord=locals()[site1], species=ads)]
        )
process_holder.add_process(
    name='H_CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CO_t+sens_H_CO_t,GibbsAds_CHO_t)-(GibbsAds_CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CO')],
    action_list=[Action(coord=t, species='CHO')],
    alpha='alpha_CO_H_to_CHO',
)
process_holder.add_process(
    name='CHO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CO_t+sens_H_CO_t,GibbsAds_CO_t+GibbsAds_H_t)-(GibbsAds_CHO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHO')],
    action_list=[Action(coord=t, species='CO')],
    alpha='rev_alpha_CO_H_to_CHO',
)

process_holder.add_process(
    name='H_CHO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_HCO_H_t+sens_HCO_H_t,GibbsAds_CHOH_t)-(GibbsAds_CHO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHO')],
    action_list=[Action(coord=t, species='CHOH')],
    alpha='alpha_avg_hydrogenation_O'
)
process_holder.add_process(
    name='CHOH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_HCO_H_t+sens_HCO_H_t,GibbsAds_CHO_t+GibbsAds_H_t)-(GibbsAds_CHOH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHOH')],
    action_list=[Action(coord=t, species='CHO')],
    alpha='rev_alpha_avg_hydrogenation_O'
)
for site1, site2 in [('t', 't_NW'), ('t', 't_NE'), ('t', 't_E')]:
    process_holder.add_process(
        name='CH_OH_{}_{}_react'.format(site1, site2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CHOH_t)-(GibbsAds_CH_t+GibbsAds_OH_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site1], species='CH'), Condition(coord=locals()[site2], species='OH')],
        action_list=[Action(coord=locals()[site1], species='CHOH'), Action(coord=locals()[site2], species='empty')],
        alpha='rev_alpha_CO_to_C_OH'
    )
    process_holder.add_process(
        name='CHOH_{}_{}_dis'.format(site1, site2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CH_t+GibbsAds_OH_t)-(GibbsAds_CHOH_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site1], species='CHOH'), Condition(coord=locals()[site2], species='empty')],
        action_list=[Action(coord=locals()[site1], species='CH'), Action(coord=locals()[site2], species='OH')],
        alpha='alpha_CO_to_C_OH'
    )

    process_holder.add_process(
        name='CH_OH_{}_{}_react'.format(site2, site1),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CHOH_t)-(GibbsAds_CH_t+GibbsAds_OH_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site2], species='CH'), Condition(coord=locals()[site1], species='OH')],
        action_list=[Action(coord=locals()[site2], species='CHOH'), Action(coord=locals()[site1], species='empty')],
        alpha='rev_alpha_CO_to_C_OH'
    )
    process_holder.add_process(
        name='CHOH_{}_{}_dis'.format(site2, site1),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_OH_t+sens_CH_OH_t,GibbsAds_CH_t+GibbsAds_OH_t)-(GibbsAds_CHOH_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site2], species='CHOH'), Condition(coord=locals()[site1], species='empty')],
        action_list=[Action(coord=locals()[site2], species='CH'), Action(coord=locals()[site1], species='OH')],
        alpha='alpha_CO_to_C_OH'
    )

    process_holder.add_process(
        name='CH_CO_{}_{}_react'.format(site1, site2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CHCO_t)-(GibbsAds_CH_t+GibbsAds_CO_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site1], species='CH'), Condition(coord=locals()[site2], species='CO')],
        action_list=[Action(coord=locals()[site1], species='CHCO'), Action(coord=locals()[site2], species='empty')],
        alpha='alpha_CO_CH_to_CHCO',
    )
    process_holder.add_process(
        name='CHCO_{}_{}_dis'.format(site1, site2),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CH_t+GibbsAds_CO_t)-(GibbsAds_CHCO_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site1], species='CHCO'), Condition(coord=locals()[site2], species='empty')],
        action_list=[Action(coord=locals()[site1], species='CH'), Action(coord=locals()[site2], species='CO')],
        alpha='rev_alpha_CO_CH_to_CHCO',
    )

    process_holder.add_process(
        name='CH_CO_{}_{}_react'.format(site2, site1),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CHCO_t)-(GibbsAds_CH_t+GibbsAds_CO_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site2], species='CH'), Condition(coord=locals()[site1], species='CO')],
        action_list=[Action(coord=locals()[site2], species='CHCO'), Action(coord=locals()[site1], species='empty')],
        alpha='alpha_CO_CH_to_CHCO',
    )
    process_holder.add_process(
        name='CHCO_{}_{}_dis'.format(site2, site1),
        rate_constant='1/(beta*h)*exp(-beta*max(max(GibbsAds_CH_CO_t+sens_CH_CO_t,GibbsAds_CH_t+GibbsAds_CO_t)-(GibbsAds_CHCO_t),0)*eV)',
        condition_list=[Condition(coord=locals()[site2], species='CHCO'), Condition(coord=locals()[site1], species='empty')],
        action_list=[Action(coord=locals()[site2], species='CH'), Action(coord=locals()[site1], species='CO')],
        alpha='rev_alpha_CO_CH_to_CHCO',
    )
process_holder.add_process(
    name='H_CH_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_t+sens_CH_H_t,GibbsAds_CH2_t)-(GibbsAds_CH_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH')],
    action_list=[Action(coord=t, species='CH2')],
    alpha='rev_alpha_CH2_to_CH_H',
)
process_holder.add_process(
    name='CH2_t_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH_H_t+sens_CH_H_t,GibbsAds_CH_t+GibbsAds_H_t)-(GibbsAds_CH2_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2')],
    action_list=[Action(coord=t, species='CH')],
    alpha='alpha_CH2_to_CH_H',
)

process_holder.add_process(
    name='H_CH2_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2_t+sens_H_CH2_t,GibbsAds_CH3_t)-(GibbsAds_CH2_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2')],
    action_list=[Action(coord=t, species='CH3')],
    alpha='rev_alpha_CH3_to_CH2_H',
)
process_holder.add_process(
    name='CH3_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2_t+sens_H_CH2_t,GibbsAds_CH2_t+GibbsAds_H_t)-(GibbsAds_CH3_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3')],
    action_list=[Action(coord=t, species='CH2')],
    alpha='alpha_CH3_to_CH2_H',
)

process_holder.add_process(
    name='H_CH3_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3_t+sens_H_CH3_t,GibbsGas_CH4gas)-(GibbsAds_CH3_t+GibbsAds_H_t),0)*eV)',
    tof_count={'CH4_formation': 1},
    condition_list=[Condition(coord=t, species='CH3')],
    action_list=[Action(coord=t, species='empty')],
    alpha='rev_alpha_CH4_to_CH3_H',
)
process_holder.add_process(
    name='CH4_t_ads',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3_t+sens_H_CH3_t,GibbsAds_CH3_t+GibbsAds_H_t)-(GibbsGas_CH4gas),0)*eV)',
    tof_count={'CH4_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='CH3')],
    alpha='alpha_CH4_to_CH3_H',
)

process_holder.add_process(
    name='H_OH_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_OH_t+sens_H_OH_t,GibbsGas_H2Ogas)-(GibbsAds_OH_t+GibbsAds_H_t),0)*eV)',
    tof_count={'H2O_formation': 1},
    condition_list=[Condition(coord=t, species='OH')],
    action_list=[Action(coord=t, species='empty')],
    alpha='rev_alpha_H2O_to_OH_H',
)
process_holder.add_process(
    name='H2O_t_ads',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_OH_t+sens_H_OH_t,GibbsAds_OH_t+GibbsAds_H_t)-(GibbsGas_H2Ogas),0)*eV)',
    tof_count={'H2O_formation': -1},
    condition_list=[Condition(coord=t, species='empty')],
    action_list=[Action(coord=t, species='OH')],
    alpha='alpha_H2O_to_OH_H',
)

process_holder.add_process(
    name='H_CHCO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_t+sens_H_CHCO_t,GibbsAds_CH2CO_t)-(GibbsAds_CHCO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CHCO')],
    action_list=[Action(coord=t, species='CH2CO')],
    alpha='alpha_avg_hydrogenation_C',
)
process_holder.add_process(
    name='CH2CO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CHCO_t+sens_H_CHCO_t,GibbsAds_CHCO_t+GibbsAds_H_t)-(GibbsAds_CH2CO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2CO')],
    action_list=[Action(coord=t, species='CHCO')],
    alpha='rev_alpha_avg_hydrogenation_C',
)

process_holder.add_process(
    name='H_CH2CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_t+sens_H_CH2CO_t,GibbsAds_CH3CO_t)-(GibbsAds_CH2CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH2CO')],
    action_list=[Action(coord=t, species='CH3CO')],
    alpha='alpha_avg_hydrogenation_C',
)
process_holder.add_process(
    name='CH3CO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH2CO_t+sens_H_CH2CO_t,GibbsAds_CH2CO_t+GibbsAds_H_t)-(GibbsAds_CH3CO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CO')],
    action_list=[Action(coord=t, species='CH2CO')],
    alpha='rev_alpha_avg_hydrogenation_C',
)

process_holder.add_process(
    name='H_CH3CO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_t+sens_H_CH3CO_t,GibbsAds_CH3CHO_t)-(GibbsAds_CH3CO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CO')],
    action_list=[Action(coord=t, species='CH3CHO')],
    alpha='alpha_avg_hydrogenation_C',
)
process_holder.add_process(
    name='CH3CHO_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_H_CH3CO_t+sens_H_CH3CO_t,GibbsAds_CH3CO_t+GibbsAds_H_t)-(GibbsAds_CH3CHO_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHO')],
    action_list=[Action(coord=t, species='CH3CO')],
    alpha='rev_alpha_avg_hydrogenation_C',
)

process_holder.add_process(
    name='H_CH3CHO_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_t+sens_CH3CHO_H_t,GibbsAds_CH3CHOH_t)-(GibbsAds_CH3CHO_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHO')],
    action_list=[Action(coord=t, species='CH3CHOH')],
    alpha='alpha_avg_hydrogenation_O',
)
process_holder.add_process(
    name='CH3CHOH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHO_H_t+sens_CH3CHO_H_t,GibbsAds_CH3CHO_t+GibbsAds_H_t)-(GibbsAds_CH3CHOH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHOH')],
    action_list=[Action(coord=t, species='CH3CHO')],
    alpha='rev_alpha_avg_hydrogenation_O',
)

process_holder.add_process(
    name='H_CH3CHOH_t_react',
    rate_constant='H_t_cov/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_t+sens_CH3CHOH_H_t,GibbsAds_CH3CH2OH_t)-(GibbsAds_CH3CHOH_t+GibbsAds_H_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CHOH')],
    action_list=[Action(coord=t, species='CH3CH2OH')],
    alpha='alpha_avg_hydrogenation_C',
)
process_holder.add_process(
    name='CH3CH2OH_t_dis',
    rate_constant='(1-H_t_cov)/(beta*h)*exp(-beta*max(max(GibbsAds_CH3CHOH_H_t+sens_CH3CHOH_H_t,GibbsAds_CH3CHOH_t+GibbsAds_H_t)-(GibbsAds_CH3CH2OH_t),0)*eV)',
    condition_list=[Condition(coord=t, species='CH3CH2OH')],
    action_list=[Action(coord=t, species='CH3CHOH')],
    alpha='rev_alpha_avg_hydrogenation_C',
)

pt = process_holder.add_project_processes(pt)

# Export
pt.export_xml_file('{}.xml'.format(model_name))

if COMPILE:
    cli_main('export {}.xml {} -b otf -t'.format(model_name, model_name))
