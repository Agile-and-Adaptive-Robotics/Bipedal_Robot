ForceStrain = [0	4.2805	32.333	60.386	88.438	116.49	144.54	172.6	200.65	228.7	256.75	284.81	312.86	340.91	368.96	397.02	425.07	453.12	481.17	509.23;
0.034481774	0	3.8893	30.554	57.083	83.459	109.66	135.65	161.41	186.89	212.05	236.83	261.17	284.98	308.16	330.6	352.17	372.72	392.1	410.13;
0.068963547	0	0	7.5228	33.26	58.803	84.128	109.21	134.02	158.52	182.67	206.43	229.74	252.55	274.78	296.37	317.23	337.26	356.37	374.45;
0.103445321	0	0	0	13.869	38.914	63.732	88.299	112.59	136.57	160.2	183.46	206.28	228.64	250.46	271.69	292.27	312.1	331.1	349.2;
0.137927094	0	0	0	0	22.126	46.609	70.839	94.792	118.44	141.75	164.69	187.22	209.3	230.87	251.87	272.25	291.92	310.81	328.85;
0.172408868	0	0	0	0	7.5136	31.752	55.738	79.449	102.86	125.94	148.65	170.95	192.81	214.18	234.99	255.18	274.68	293.4	311.28;
0.206890641	0	0	0	0	0	18.55	42.343	65.86	89.075	111.96	134.48	156.59	178.26	199.43	220.04	240.03	259.33	277.83	295.47;
0.241372415	0	0	0	0	0	6.6021	30.228	53.576	76.62	99.329	121.67	143.6	165.08	186.05	206.45	226.22	245.27	263.51	280.89;
0.275854188	0	0	0	0	0	0	19.101	42.293	65.174	87.714	109.88	131.62	152.9	173.66	193.83	213.34	232.12	250.06	267.09;
0.310335962	0	0	0	0	0	0	8.7502	31.789	54.509	76.879	98.859	120.41	141.47	162	181.91	201.15	219.61	237.2	253.83;
0.344817735	0	0	0	0	0	0	0	21.901	44.457	66.649	88.438	109.78	130.61	150.88	170.52	189.44	207.56	224.77	241.01;
0.379299509	0	0	0	0	0	0	0	12.503	34.887	56.892	78.476	99.59	120.18	140.17	159.5	178.08	195.82	212.62	228.43;
0.413786996	0	0	0	0	0	0	0	3.495	25.697	47.503	68.868	89.739	110.06	129.75	148.75	166.96	184.3	200.66	216;
0.448268769	0	0	0	0	0	0	0	0	16.806	38.399	59.528	80.138	100.16	119.54	138.18	156	172.91	188.81	203.67;
0.482750543	0	0	0	0	0	0	0	0	8.146	29.511	50.387	70.716	90.432	109.46	127.72	145.12	161.58	177.01	191.34;
0.517232316	0	0	0	0	0	0	0	0	0	20.782	41.387	61.416	80.799	99.46	117.32	134.29	150.28	165.22	179.05;
0.55171409	0	0	0	0	0	0	0	0	0	12.162	32.479	52.188	71.218	89.491	106.93	123.44	138.95	153.39	166.72;
0.586218718	0	0	0	0	0	0	0	0	0	3.6114	23.621	42.991	61.648	79.514	96.51	112.56	127.58	141.51	154.34;
0.620671923	0	0	0	0	0	0	0	0	0	0	14.777	33.788	52.054	69.496	86.036	101.6	116.13	129.56	141.91;
0.655182265	0	0	0	0	0	0	0	0	0	0	5.915	24.55	42.407	59.409	75.482	90.559	104.58	117.52	129.39;
0.68963547	0	0	0	0	0	0	0	0	0	0	0	15.25	32.682	49.23	64.826	79.41	92.935	105.37	116.76;
0.724145812	0	0	0	0	0	0	0	0	0	0	0	5.8662	22.859	38.943	54.055	68.143	81.172	93.127	104.04;
0.758599017	0	0	0	0	0	0	0	0	0	0	0	0	12.92	28.531	43.155	56.75	69.29	80.771	91.215;
0.793109359	0	0	0	0	0	0	0	0	0	0	0	0	2.8523	17.985	32.119	45.224	57.285	68.306	78.345;
0.827562564	0	0	0	0	0	0	0	0	0	0	0	0	0	7.294	20.942	33.564	45.155	55.732	65.329;
0.862072906	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9.6187	21.767	32.903	43.051	52.263;
0.896526111	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9.8342	20.53	30.266	39.099;
0.931036453	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	8.0383	17.381	25.887;
0.965489658	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	4.3997	12.527;
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

RelativeStrain = ForceStrain(:,1);
Force2 = ForceStrain(:,20);

save ForceStrainTable.mat
