%Muscle attachment points
%Gluteus Medius, 1
glut_med1_o = [-0.041; 0.03; 0.121]; %Pelvis, right, origin
glut_med1_i = [-0.022; -0.012; 0.056]; %Femur, right, insertion
%Gluteus Medius, 2
glut_med2_o = [-0.086; 0.044; 0.077]; %Pelvis, right, origin
glut_med2_i = [-0.026; -0.006; 0.053]; %Femur, right, insertion
%Gluteus Medius, 3
glut_med3_o = [-0.122; 0.01; 0.065]; %Pelvis, right, origin
glut_med3_i = [-0.031; -0.005; 0.052]; %Femur, right, insertion
%Gluteus Minimus, 1
glut_min1_o = [-0.047; -0.008; 0.106]; %Pelvis, right, origin
glut_min1_i = [-0.007; -0.01; 0.056]; %Femur, right, insertion
%Gluteus Minimus, 2
glut_min2_o = [-0.063; -0.006; 0.099]; %Pelvis, right, origin
glut_min2_i = [-0.01; -0.01; 0.056]; %Femur, right, insertion
%Gluteus Minimus, 3
glut_min3_o = [-0.083; -0.006; 0.086]; %Pelvis, right, origin
glut_min3_i = [-0.014; -0.008; 0.055]; %Femur, right, insertion
%Gluteus Maximus, 1
glut_max1_o = [-0.12; 0.061; 0.07]; %Pelvis, right, origin
glut_max1_wr1 = [-0.129; 0.001; 0.089]; %Pelvis, right, wrapping point 1
glut_max1_wr2 = [-0.046; -0.025; 0.039]; %Femur, right, wrapping point 2
glut_max1_i = [-0.028; -0.057; 0.047]; %Femur, right, insertion
%Gluteus Maximus, 2
glut_max2_o = [-0.135; 0.018; 0.056]; %Pelvis, right, origin
glut_max2_wr1 = [-0.138; -0.052; 0.091]; %Pelvis, right, wrapping point 1
glut_max2_wr2 = [-0.043; -0.053; 0.029]; %Femur, right, wrapping point 2
glut_max2_i = [-0.016; -0.102; 0.042]; %Femur, right, insertion
%Gluteus Maximus, 3
glut_max3_o = [-0.156; -0.031; 0.006]; %Pelvis, right, origin
glut_max3_wr1 = [-0.153; -0.105; 0.04]; %Pelvis, right, wrapping point 1
glut_max3_wr2 = [-0.03; -0.104; 0.014]; %Femur, right, wrapping point 2
glut_max3_i = [-0.006; -0.142; 0.041]; %Femur, right, insertion
%Semimembranosus
semimem_o = [-0.119; -0.097; 0.072]; %Femur, right, origin
% semimem_wr1 = [-0.035; -0.035; -0.019]; %via point, don't use
semimem_i = [-0.027; -0.048; -0.02]; %Tibia, right, insertion
%Semitendinosus
semiten_o = [-0.126; -0.11; 0.06]; %Femur, right, origin
% semiten_wr1 = [-0.042; -0.029; -0.023]; %via point, don't use
semiten_wr2 = [-0.033; -0.053; -0.023]; %Tibia, right, wrapping point 2
semiten_wr3 = [-0.011; -0.075; -0.024]; %Tibia, right, wrapping point 3
semiten_i = [0.003; -0.096; -0.019]; %Tibia, right, insertion
%Bicep Femoris, Long Head
bifemlh_o = [-0.126; -0.103; 0.069]; %Pelvis, right, origin
bifemlh_wr1 = [-0.03; -0.036; 0.029]; %Tibia, right, wrapping point 1
bifemlh_i = [-0.023; -0.056; 0.034]; %Tibia, right, insertion
%Bicep Femoris, Short Head
bifemsh_o = [0.005; -0.211; 0.023]; %Femur, right, origin
bifemsh_wr1 = [-0.03; -0.036; 0.029]; %Tibia, right, wrapping point 1
bifemsh_i = [-0.023; -0.056; 0.034]; %Tibia, right, insertion
%Sartorius
sar_o = [-0.015; -0.001; 0.124]; %Pelvis, right, origin
sar_wr1 = [-0.003; -0.357; -0.042]; %Femur, right, wrapping point 1
sar_wr2 = [-0.006; -0.042; -0.04]; %Tibia, right, wrapping point 2
sar_wr3 = [0.006; -0.059; -0.038]; %Tibia, right, wrapping point 3
sar_i = [0.024; -0.084; -0.025]; %Tibia, right, insertion
%Adductor Longus
add_long_o = [-0.032; -0.084; 0.017]; %Pelvis, right, origin
add_long_i = [0.005; -0.211; 0.023]; %Femur, right, insertion
%Adductor Brevis
add_brev_o = [-0.059; -0.092; 0.016]; %Pelvis, right, origin
add_brev_i = [0.001; -0.12; 0.029]; %Femur, right, insertion
%Adductor Magnus, 1
add_mag1_o = [-0.073; -0.117; 0.026]; %Pelvis, right, origin
add_mag1_i = [-0.004; -0.121; 0.034]; %Femur, right, insertion
%Adductor Magnus, 2
add_mag2_o = [-0.083; -0.119; 0.031]; %Pelvis, right, origin
add_mag2_i = [0.005; -0.228; 0.023]; %Femur, right, insertion
%Adductor Magnus, 3
add_mag3_o = [-0.111; -0.114; 0.049]; %Pelvis, right, origin
add_mag3_i = [0.007; -0.384; -0.027]; %Femur, right, insertion
%Tensor Fasciae Latae
tfl_o = [-0.031; 0.021; 0.124]; %Pelvis, right, origin
tfl_wr1 = [0.029; -0.1; 0.06]; %Femur, right, wrapping point 1
tfl_wr2 = [0.005; -0.405; 0.036]; %Femur, right, wrapping point 2
tfl_i = [0.006; -0.049; 0.03]; %Tibia, right, insertion
%Pectineus
pect_o = [-0.043; -0.077; 0.045]; %Pelvis, right, origin
pect_i = [-0.012; -0.082; 0.025]; %Femur, right, insertion
%Gracilis
grac_o = [-0.074; -0.119; 0.028]; %Pelvis, right, origin
% grac_wr1 = [-0.027; -0.032; -0.038]; %via point, don't use
grac_wr2 = [-0.019; -0.052; -0.036]; %Tibia, right, wrapping point 2
grac_i = [0.006; -0.084; -0.023]; %Tibia, right, insertion
%Iliacus
iliacus_o = [-0.067; 0.036; 0.085]; %Pelvis, right, origin
iliacus_wr1 = [-0.026; -0.055; 0.081]; %Pelvis, right, wrapping point 1
% iliacus_wr2 = [-0.029; -0.08; 0.082]; %via point, don't use
iliacus_wr3 = [0.002; -0.054; 0.006]; %Femur, right, wrapping point 3
iliacus_i = [-0.019; -0.062; 0.013]; %Femur, right, insertion
%Psoas
psoas_o = [-0.065; 0.089; 0.029]; %Pelvis, right, origin
psoas_wr1 = [-0.024; -0.057; 0.076]; %Pelvis, right, wrapping point 1
% psoas_wr2 = [-0.029; -0.08; 0.082]; %via point, don't use
psoas_wr3 = [0.002; -0.051; 0.004]; %Femur, right, wrapping point 3
psoas_i = [-0.019; -0.06; 0.01]; %Femur, right, insertion
%Quadriceps Femoris
quad_fem_o = [-0.114; -0.115; 0.052]; %Pelvis, right, origin
quad_fem_i = [-0.038; -0.036; 0.037]; %Femur, right, insertion
%Gemelli
gem_o = [-0.113; -0.082; 0.071]; %Pelvis, right, origin
gem_i = [-0.014; -0.003; 0.044]; %Femur, right, insertion
%Piriformis
peri_o = [-0.14; 0; 0.024]; %Pelvis, right, origin
peri_wr1 = [-0.119; -0.028; 0.066]; %Pelvis, right, wrapping point 1
peri_i = [-0.015; -0.004; 0.044]; %Femur, right, insertion
%Rectus Femoris
rect_fem_o = [-0.03; -0.031; 0.097]; %Pelvis, right, origin
% rect_fem_wr1 = [0.033; -0.403; 0.002]; %via point, don't use
rect_fem_i = [0.062; 0.021; 0.001]; %Tibia, right, insertion
%Vastus Medialis
vas_med_o = [0.014; -0.21; 0.019]; %Femur, right, origin
vas_med_wr1 = [0.036; -0.277; 0.001]; %Femur, right, wrapping point 1
% vas_med_wr2 = [0.037; -0.405; -0.012]; %via point, don't use
% vas_med_wr3 = [0.027; -0.426; -0.013]; %via point, don't use
vas_med_i = [0.056; 0.022; -0.015]; %Tibia, right, insertion
%Vastus Intermedius
vas_int_o = [0.029; -0.192; 0.031]; %Femur, right, origin
vas_int_wr1 = [0.034; -0.208; 0.028]; %Femur, right, wrapping point 1
%vas_int_wr2 = [0.034; -0.403; 0.006]; %via point, don't use
vas_int_i = [0.055; 0.025; 0.002]; %Tibia, right, insertion
%Vastus Lateralis
vas_lat_o = [0.005; -0.185; 0.035]; %Femur, right, origin
vas_lat_wr1 = [0.027; -0.259; 0.041]; %Femur, right, wrapping point 1
% vas_lat_wr2 = [0.036; -0.403; 0.002]; %via point, don't use
% vas_lat_wr3 = [0.025; -0.424; 0.018]; %via point, don't use
vas_lat_i = [0.06; 0.02; 0.016]; %Tibia, right, insertion
%Medial Gastrocnemius
med_gas_o = [-0.019; -0.393; -0.024]; %Femur, right, origin
% med_gas_wr1 = [-0.03; -0.402; -0.026]; %via point, don't use
med_gas_i = [0; 0.031; -0.005]; %Calcn, right, insertion
%Lateral Gastrocenemius
lat_gas_o = [-0.022; -0.395; 0.027]; %Femur, right, origin
% lat_gas_wr1 = [-0.03; -0.402; 0.027]; %via point, don't use
lat_gas_i = [0; 0.031; -0.005]; %Calcn, right, insertion
%Soleus
soleus_o = [-0.002; -0.153; 0.007]; %Tibia, right, origin
soleus_i = [0; 0.031; -0.005]; %Calcn, right, insertion
%Tibialis Posterior
tib_post_o = [-0.009; -0.135; 0.002]; %Tibia, right, origin
tib_post_wr1 = [-0.014; -0.405; -0.023]; %Tibia, right, wrapping point 1
tib_post_wr2 = [0.042; 0.033; -0.029]; %Calcn, right, wrapping point 2
tib_post_i = [0.077; 0.016; -0.028]; %Calcn, right, insertion
%Flexor Digitorum Longus
flex_dig_o = [-0.008; -0.205; 0.002]; %Tibia, right, origin
flex_dig_wr1 = [-0.015; -0.405; -0.02]; %Tibia, right, wrapping point 1
flex_dig_wr2 = [0.044; 0.032; -0.028]; %Calcn, right, wrapping point 2
flex_dig_wr3 = [0.071; 0.018; -0.026]; %Calcn, right, wrapping point 3
flex_dig_wr4 = [0.166; -0.008; 0.012]; %Calcn, right, wrapping point 4
flex_dig_wr5 = [-0.002; -0.008; 0.015]; %Toes, right, wrapping point 5
flex_dig_wr6 = [0.028; -0.007; 0.022]; %Toes, right, wrapping point 6
flex_dig_i = [0.044; -0.006; 0.024]; %Toes, right, insertion
%Flexor Hallucis Longus
flex_hal_o = [-0.008; -0.233; 0.024]; %Tibia, right, origin
flex_hal_wr1 = [-0.019; -0.408; -0.017]; %Tibia, right, wrapping point 1
flex_hal_wr2 = [0.037; 0.028; -0.024]; %Calcn, right, wrapping point 2
flex_hal_wr3 = [0.104; 0.007; -0.026]; %Calcn, right, wrapping point 3
flex_hal_wr4 = [0.173; -0.005; -0.027]; %Calcn, right, wrapping point 4
flex_hal_wr5 = [0.016; -0.006; -0.026]; %Toes, right, wrapping point 5
flex_hal_i = [0.056; -0.01; -0.018]; %Toes, right, insertion
%Tibialis Anterior
tib_ant_o = [0.018; -0.162; 0.011]; %Tibia, right, origin
tib_ant_wr1 = [0.033; -0.395; -0.018]; %Tibia, right, wrapping point 1
tib_ant_i = [0.117; 0.018; -0.03]; %Calcn, right, insertion
%Peroneus Brevis
per_brev_o = [-0.007; -0.265; 0.033]; %Tibia, right, origin
per_brev_wr1 = [-0.02; -0.418; 0.028]; %Tibia, right, wrapping point 1
per_brev_wr2 = [-0.014; -0.429; 0.029]; %Tibia, right, wrapping point 2
per_brev_wr3 = [0.047; 0.027; 0.023]; %Calcn, right, wrapping point 3
per_brev_i = [0.068; 0.022; 0.034]; %Calcn, right, insertion
%Peroneus Longus
per_long_o = [0; -0.157; 0.036]; %Tibia, right, origin
per_long_wr1 = [-0.021; -0.42; 0.029]; %Tibia, right, wrapping point 1
per_long_wr2 = [-0.016; -0.432; 0.029]; %Tibia, right, wrapping point 2
per_long_wr3 = [0.044; 0.023; 0.022]; %Calcn, right, wrapping point 3
per_long_wr4 = [0.068; 0.011; 0.028]; %Calcn, right, wrapping point 4
per_long_wr5 = [0.085; 0.007; 0.012]; %Calcn, right, wrapping point 5
per_long_i = [0.12; 0.008; -0.018]; %Calcn, right, insertion
%Peroneus Tertius
per_tert_o = [0.001; -0.28; 0.023]; %Tibia, right, origin
per_tert_wr1 = [0.023; -0.407; 0.016]; %Tibia, right, wrapping point 1
per_tert_i = [0.086; 0.023; 0.03]; %Calcn, right, insertion
%Extensor Digitorum Longus
ext_dig_o = [0.003; -0.138; 0.028]; %Tibia, right, origin
ext_dig_wr1 = [0.029; -0.401; 0.007]; %Tibia, right, wrapping point 1
ext_dig_wr2 = [0.092; 0.039; 0]; %Calcn, right, wrapping point 2
ext_dig_wr3 = [0.162; 0.006; 0.013]; %Calcn, right, wrapping point 3
ext_dig_wr4 = [0; 0.005; 0.015]; %Toes, right, wrapping point 4
ext_dig_i = [0.044; 0; 0.025]; %Toes, right, insertion
%Extensor Hallucis Longus
ext_hal_o = [0.001; -0.177; 0.023]; %Tibia, right, origin
ext_hal_wr1 = [0.033; -0.398; -0.008]; %Tibia, right, wrapping point 1
ext_hal_wr2 = [0.097; 0.039; -0.021]; %Calcn, right, wrapping point 2
ext_hal_wr3 = [0.129; 0.031; -0.026]; %Calcn, right, wrapping point 3
ext_hal_wr4 = [0.173; 0.014; -0.028]; %Calcn, right, wrapping point 4
ext_hal_wr5 = [0.03; 0.004; -0.024]; %Toes, right, wrapping point 5
ext_hal_i = [0.056; 0.003; -0.019]; %Toes, right, insertion
%Erector Spinae
ercspn_o = [-0.14; 0.044; 0.044]; %Pelvis, right, origin
ercspn_i = [-0.055; 0.11; 0.024]; %Torso, right, insertion
%Internal Oblique
intobl_o = [-0.04; 0.07; 0.116]; %Pelvis, right, origin
intobl_i = [0.07; 0.16; 0.015]; %Torso, right, insertion
%External Oblique extobl
extobl_o = [-0.03; -0.064; 0.01]; %Pelvis, right, origin
extobl_i = [0.065; 0.11; 0.11]; %Torso, right, insertion

save AttachPoints.mat