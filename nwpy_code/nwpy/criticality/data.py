###############################################################################
"""
    Last edited on June 20, 2019
    
    @author: matz
    
    comments: Data about material molar masses, abundances, cross sections
    
"""
###############################################################################
# Molar masses of rock, water, and heavy metal isotopes
mm = {}
mm['h1']   = 1.00782503224
mm['c12']  = 12.00000
mm['c13']  = 13.00335483521
mm['o16']  = 15.999
mm['na23'] = 22.989769282
mm['mg24'] = 23.985041697
mm['mg25'] = 24.98583696
mm['mg26'] = 25.98259297
mm['al27'] = 26.98153841
mm['si28'] = 27.9769265350
mm['si29'] = 28.9764946653
mm['si30'] = 29.973770137
mm['p31']  = 30.97376163
mm['s32']  = 31.97207100
mm['s33']  = 32.97145876
mm['s34']  = 33.96786690
mm['s36']  = 35.96708076
mm['k39']  = 38.96370668
mm['k40']  = 39.96399848
mm['k41']  = 40.96182576
mm['ca40'] = 39.96259098
mm['ca42'] = 41.95861801
mm['ca43'] = 42.9587666
mm['ca44'] = 43.9554818
mm['ca46'] = 45.9536926
mm['ca48'] = 47.952534
mm['ti46'] = 45.9526316
mm['ti47'] = 46.9517631
mm['ti48'] = 47.9479463
mm['ti49'] = 48.9478700
mm['ti50'] = 49.9447912
mm['mn55'] = 54.9380451
mm['fe54'] = 53.9396090
mm['fe56'] = 55.9349363
mm['fe57'] = 56.9353928
mm['fe58'] = 57.9332744
mm['th232'] = 232.0380553
mm['u232']  = 232.0332966
mm['u233']  = 233.0449876
mm['u234']  = 234.04054
mm['u235']  = 235.0441616
mm['u236']  = 236.045766
mm['u237']  = 237.048379
mm['u238']  = 238.050992
mm['u239']  = 239.0546137
mm['u240']  = 240.056218
mm['u241']  = 241.0608484
mm['np235'] = 235.0440633
mm['np236'] = 236.04657
mm['np237'] = 237.0481734
mm['np238'] = 238.0509464
mm['np239'] = 239.0529390
mm['pu236'] = 236.045766
mm['pu237'] = 237.048379
mm['pu238'] = 238.0489747
mm['pu239'] = 239.0525963
mm['pu240'] = 240.0542007
mm['pu241'] = 241.0487444
mm['pu242'] = 242.058418
mm['pu243'] = 243.0620397
mm['pu244'] = 244.0646527
mm['pu246'] = 246.0698787
mm['am239'] = 239.0530245
mm['am240'] = 240.055300
mm['am241'] = 241.0568291
mm['am242'] = 242.0595492
mm['am243'] = 243.0613811
mm['am244'] = 244.0642848
mm['cm242'] = 242.0588358
mm['cm243'] = 243.0613891
mm['cm244'] = 244.0627526
mm['cm245'] = 245.0654912
mm['cm246'] = 246.0672237
mm['cm247'] = 247.070354
mm['cm248'] = 248.072349
mm['cm249'] = 249.075953
mm['cm250'] = 250.078357
###############################################################################
# Abundances of rock and water isotopes
a = {}
a['h1']   = 1.0000
a['c12']  = 0.9890
a['c13']  = 0.0110
a['o16']  = 1.0000
a['na23'] = 1.0000
a['mg24'] = 0.7900
a['mg25'] = 0.1000
a['mg26'] = 0.1100
a['al27'] = 1.0000
a['si28'] = 0.9220
a['si29'] = 0.0470
a['si30'] = 0.0310
a['p31']  = 1.0000
a['s32']  = 0.9499
a['s33']  = 0.0075
a['s34']  = 0.0425
a['s36']  = 0.0001
a['k39']  = 0.93258
a['k40']  = 0.00012
a['k41']  = 0.06730
a['ca40'] = 0.96941
a['ca42'] = 0.00647
a['ca43'] = 0.00135
a['ca44'] = 0.02086
a['ca46'] = 0.00004
a['ca48'] = 0.00187
a['ti46'] = 0.0825
a['ti47'] = 0.0744
a['ti48'] = 0.7372
a['ti49'] = 0.0541
a['ti50'] = 0.0518
a['mn55'] = 1.0000
a['fe54'] = 0.0585
a['fe56'] = 0.9175
a['fe57'] = 0.0212
a['fe58'] = 0.0028
###############################################################################
# CROSS SECTIONS AT 0.625 EV
# 't': total
# 'f': fission
# 'c': capture
# 's': elastic scattering
# 'tr': transport
xs = {}
xs['h1']    = {'t':20.9194, 's':20.8525, 'c':0.066967}
xs['c12']   = {'t':4.75893, 's':4.75816, 'c':7.777E-4}
xs['c13']   = {'t':5.74228, 's':5.74198, 'c':3.024E-4}
xs['o16']   = {'t':3.80231, 's':3.80228, 'c':3.433E-5}
xs['na23']  = {'t':3.33411, 's':3.22775, 'c':0.106374}
xs['mg24']  = {'t':3.76551, 's':3.75538, 'c':0.010129}
xs['mg25']  = {'t':2.63774, 's':2.59939, 'c':0.038345}
xs['mg26']  = {'t':2.84289, 's':2.83517, 'c':0.007718}
xs['al27']  = {'t':1.47434, 's':1.42732, 'c':0.047019}
xs['si28']  = {'t':1.99295, 's':1.95887, 'c':0.034073}
xs['si29']  = {'t':2.61178, 's':2.58761, 'c':0.024167}
xs['si30']  = {'t':2.48189, 's':2.46031, 'c':0.021574}
xs['p31']   = {'t':4.13796, 's':4.10440, 'c':0.033565}
xs['s32']   = {'t':1.07775, 's':0.96437, 'c':0.106377}
xs['s33']   = {'t':3.08757, 's':2.84669, 'c':0.070499}
xs['s34']   = {'t':2.12393, 's':2.07889, 'c':0.045035}
xs['s36']   = {'t':2.22931, 's':2.18617, 'c':0.043140}
xs['k39']   = {'t':2.48572, 's':2.05713, 'c':0.427045}
xs['k40']   = {'t':9.63628, 's':2.75422, 'c':6.045609}
xs['k41']   = {'t':2.84506, 's':2.55287, 'c':0.292196}
xs['ca40']  = {'t':2.71607, 's':2.63301, 'c':0.082561}
xs['ca42']  = {'t':1.36100, 's':1.22343, 'c':0.137574}
xs['ca43']  = {'t':6.51081, 's':4.16360, 'c':2.347166}
xs['ca44']  = {'t':3.50324, 's':3.32431, 'c':0.178918}
xs['ca46']  = {'t':3.05735, 's':2.90821, 'c':0.149141}
xs['ca48']  = {'t':3.94169, 's':3.72154, 'c':0.220158}
xs['ti46']  = {'t':2.82538, 's':2.70660, 'c':0.118778}
xs['ti47']  = {'t':3.56178, 's':3.23423, 'c':0.327550}
xs['ti48']  = {'t':5.67425, 's':3.99970, 'c':1.674553}
xs['ti49']  = {'t':0.91462, 's':0.54039, 'c':0.374228}
xs['ti50']  = {'t':3.75915, 's':3.72299, 'c':0.036163}
xs['mn55']  = {'t':4.78648, 's':2.10648, 'c':2.680011}
xs['fe54']  = {'t':2.61807, 's':2.16421, 'c':0.453872}
xs['fe56']  = {'t':12.6153, 's':12.0908, 'c':0.524492}
xs['fe57']  = {'t':1.10372, 's':0.60410, 'c':0.499616}
xs['fe58']  = {'t':7.74657, 's':7.48191, 'c':0.264663}
xs['th232'] = {'t':13.9525, 's':12.8068, 'f':0.00000, 'c':1.14568}
xs['u233']  = {'t':146.189, 's':11.5198, 'f':123.543, 'c':11.1265}
xs['u234']  = {'t':27.5985, 's':14.1970, 'f':0.00890, 'c':13.3925}
xs['u235']  = {'t':87.4908, 's':13.1251, 'f':65.1892, 'c':9.17644}
xs['u236']  = {'t':9.83997, 's':8.65483, 'f':0.01150, 'c':1.17364}
xs['u238']  = {'t':9.71562, 's':9.12833, 'f':3.93E-6, 'c':0.58729}
xs['np237'] = {'t':81.1333, 's':15.3470, 'f':0.00337, 'c':65.7831}
xs['pu238'] = {'t':46.0438, 's':38.6740, 'f':0.35730, 'c':7.01246}
xs['pu239'] = {'t':145.750, 's':11.1780, 'f':98.7192, 'c':35.8531}
xs['pu240'] = {'t':316.423, 's':1.81402, 'f':0.07627, 'c':314.533}
xs['pu241'] = {'t':65.4717, 's':12.0214, 'f':38.1565, 'c':15.2938}
xs['pu242'] = {'t':14.5659, 's':7.96817, 'f':0.00252, 'c':6.59529}
xs['pu244'] = {'t':10.6086, 's':10.3093, 'f':0.00030, 'c':0.29890}
xs['am241'] = {'t':982.080, 's':20.3627, 'f':2.38779, 'c':959.328}
xs['am243'] = {'t':48.0394, 's':6.18891, 'f':0.03066, 'c':41.8197}
xs['cm244'] = {'t':15.4059, 's':12.0101, 'f':0.17649, 'c':3.21931}
xs['cm245'] = {'t':225.273, 's':8.89549, 'f':200.809, 'c':15.5694}
###############################################################################
# CROSS SECTIONS AT 0.0253 eV
#xs = {}
#xs['h1']    = {'t':30.4009, 's':30.0683, 'c':0.332587}
#xs['c12']   = {'t':4.95134, 's':4.94748, 'c':0.003860}
#xs['c13']   = {'t':5.95553, 's':5.95403, 'c':0.001500}
#xs['o16']   = {'t':3.91352, 's':3.91335, 'c':1.698E-4}
#xs['na23']  = {'t':3.92136, 's':3.39336, 'c':0.528000}
#xs['mg24']  = {'t':3.87750, 's':3.82721, 'c':0.050289}
#xs['mg25']  = {'t':2.83739, 's':2.64702, 'c':0.190374}
#xs['mg26']  = {'t':2.92348, 's':2.88517, 'c':0.038310}
#xs['al27']  = {'t':1.68512, 's':1.45165, 'c':0.233463}
#xs['si28']  = {'t':2.16033, 's':1.99119, 'c':0.169141}
#xs['si29']  = {'t':2.74900, 's':2.62904, 'c':0.119961}
#xs['si30']  = {'t':2.60572, 's':2.49863, 'c':0.107085}
#xs['p31']   = {'t':4.33602, 's':4.16666, 'c':0.169361}
#xs['s32']   = {'t':1.51397, 's':0.97865, 'c':0.528215}
#xs['s33']   = {'t':3.41079, 's':2.88788, 'c':0.350075}
#xs['s34']   = {'t':2.33183, 's':2.10821, 'c':0.223618}
#xs['s36']   = {'t':2.36581, 's':2.21533, 'c':0.150482}
#xs['k39']   = {'t':4.21678, 's':2.08505, 'c':2.127420}
#xs['k40']   = {'t':37.5724, 's':2.78540, 'c':30.00980}
#xs['k41']   = {'t':4.05725, 's':2.59612, 'c':1.461130}
#xs['ca40']  = {'t':3.07542, 's':2.66289, 'c':0.410045}
#xs['ca42']  = {'t':1.91960, 's':1.23651, 'c':0.683094}
#xs['ca43']  = {'t':15.8743, 's':4.20917, 'c':11.66490}
#xs['ca44']  = {'t':4.24641, 's':3.35778, 'c':0.888633}
#xs['ca46']  = {'t':3.70652, 's':2.96613, 'c':0.740393}
#xs['ca48']  = {'t':4.84849, 's':3.75555, 'c':1.092930}
#xs['ti46']  = {'t':3.32225, 's':2.73250, 'c':0.589748}
#xs['ti47']  = {'t':4.89113, 's':3.26475, 'c':1.626380}
#xs['ti48']  = {'t':12.3559, 's':4.03801, 'c':8.317910}
#xs['ti49']  = {'t':2.40885, 's':0.54603, 'c':1.862820}
#xs['ti50']  = {'t':3.93521, 's':3.75567, 'c':0.179537}
#xs['mn55']  = {'t':15.3947, 's':2.11629, 'c':13.27840}
#xs['fe54']  = {'t':4.43640, 's':2.18248, 'c':2.253920}
#xs['fe56']  = {'t':14.7949, 's':12.1899, 'c':2.605060}
#xs['fe57']  = {'t':3.09099, 's':0.60915, 'c':2.481850}
#xs['fe58']  = {'t':8.85607, 's':7.54125, 'c':1.314830}
#xs['th232'] = {'t':20.3823, 's':13.0442, 'f':0.00000, 'c':7.33806}
#xs['u233']  = {'t':588.579, 's':12.1781, 'f':534.072, 'c':42.3290}
#xs['u234']  = {'t':118.259, 's':17.2836, 'f':0.06710, 'c':100.908}
#xs['u235']  = {'t':700.185, 's':14.1088, 'f':586.691, 'c':99.3843}
#xs['u236']  = {'t':14.0222, 's':8.84111, 'f':0.04711, 'c':5.13396}
#xs['u238']  = {'t':11.9234, 's':9.23968, 'f':1.85E-5, 'c':2.68368}
#xs['np237'] = {'t':191.322, 's':15.8714, 'f':0.02037, 'c':175.430}
#xs['pu238'] = {'t':585.184, 's':154.561, 'f':17.7677, 'c':412.855}
#xs['pu239'] = {'t':1025.60, 's':8.07126, 'f':747.393, 'c':270.139}
#xs['pu240'] = {'t':291.319, 's':1.73406, 'f':0.05632, 'c':289.529}
#xs['pu241'] = {'t':1386.61, 's':11.2591, 'f':1012.30, 'c':363.047}
#xs['pu242'] = {'t':29.9977, 's':8.71509, 'f':0.01382, 'c':21.2688}
#xs['pu244'] = {'t':12.0985, 's':10.3866, 'f':0.00172, 'c':1.71021}
#xs['am241'] = {'t':699.219, 's':11.8181, 'f':3.12241, 'c':684.279}
#xs['am243'] = {'t':88.3874, 's':7.88178, 'f':0.08134, 'c':80.4243}
#xs['cm244'] = {'t':28.6466, 's':12.3887, 'f':1.02240, 'c':15.2355}
#xs['cm245'] = {'t':2411.60, 's':10.2500, 'f':2054.33, 'c':347.019}
###############################################################################
# NEUTRONS PER FISSION (nu-bar) AT 1e-5 eV
nu = {}
nu['th232'] = 2.104700
nu['u233'] =  2.485240
nu['u234'] =  2.364900
nu['u235'] =  2.429850
nu['u236'] =  2.371200
nu['u238'] =  2.443040
nu['np237'] = 2.635810
nu['pu238'] = 2.900917
nu['pu239'] = 2.876938
nu['pu240'] = 2.897150
nu['pu241'] = 2.945300
nu['pu242'] = 2.893200
nu['pu244'] = 2.936700
nu['am241'] = 3.080270
nu['am243'] = 3.272833
nu['cm244'] = 3.463651
nu['cm245'] = 3.596482
###############################################################################
# ATOMIC MASS
#a = {}
#a['h1']    = 1
#a['c12']   = 12
#a['c13']   = 13
#a['o16']   = 16
#a['na23']  = 23
#a['mg24']  = 24
#a['mg25']  = 25
#a['mg26']  = 26
#a['al27']  = 27
#a['si28']  = 28
#a['si29']  = 29
#a['si30']  = 30
#a['p31']   = 31
#a['s32']   = 32
#a['s33']   = 33
#a['s34']   = 34
#a['s36']   = 36
#a['k39']   = 39
#a['k40']   = 40
#a['k41']   = 41
#a['ca40']  = 40
#a['ca42']  = 42
#a['ca43']  = 43
#a['ca44']  = 44
#a['ca46']  = 46
#a['ca48']  = 48
#a['ti46']  = 46
#a['ti47']  = 47
#a['ti48']  = 48
#a['ti49']  = 49
#a['ti50']  = 50
#a['mn55']  = 55
#a['fe54']  = 54
#a['fe56']  = 56
#a['fe57']  = 57
#a['fe58']  = 58
#a['th232'] = 232
#a['u233']  = 233
#a['u234']  = 234
#a['u235']  = 235
#a['u236']  = 236
#a['u238']  = 238
#a['np237'] = 237
#a['pu238'] = 238
#a['pu239'] = 239
#a['pu240'] = 240
#a['pu241'] = 241
#a['pu242'] = 242
#a['pu244'] = 244
#a['am241'] = 241
#a['am243'] = 243
#a['cm244'] = 244
#a['cm245'] = 245
################################################################################
