%% MATPOWER Case Format : Version 2
function mpc = case67acdc_scopf
mpc.version = '2';

%% stochastic data
%column_names%  dst         pa      pb
mpc.sdata = [
                'Normal'    0.0     0.0;				
                'Beta'		2.0		5.0;
                'Beta'		5.0		2.0;
];

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
	1	3	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	2	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	3	1	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	4	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	5	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	6	1	191.0	76.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	7	1	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	8	1	287.0	73.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	9	1	186.0	74.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	10	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	11	1	271.0	55.00000000000001	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	12	1	171.0	87.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	13	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	14	1	199.0	60.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	15	1	112.99999999999999	52.5	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	16	1	38.0	7.000000000000001	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	17	1	275.0	106.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	18	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	19	1	165.0	46.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	20	1	178.0	82.5	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	21	1	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	22	1	30.0	7.000000000000001	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	23	1	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	24	1	32.0	7.000000000000001	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	25	2	0.0	0.0	0.0	0.0	1	1.0	0.0	380.0	1	1.1	0.9				
	26	1	395.0	89.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	27	1	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	28	1	665.0	99.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	29	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	30	1	266.0	100.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	31	1	844.9999999999999	119.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	32	1	332.0	137.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	33	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	34	1	540.0	158.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	35	1	459.99999999999994	97.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	36	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	37	1	451.0	190.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	38	1	150.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	39	1	629.0	87.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	40	1	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	41	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	42	1	859.0	180.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	43	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	44	1	474.0	92.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	45	1	668.0	109.00000000000001	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	46	1	614.0	95.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	47	1	81.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	48	1	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	49	1	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	50	2	0.0	0.0	0.0	0.0	2	1.0	0.0	380.0	1	1.1	0.9				
	51	1	430.0	123.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	52	1	309.0	102.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	53	1	100.0	30.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	54	1	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	55	1	303.0	110.00000000000001	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	56	2	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	57	1	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	58	1	324.0	157.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	59	2	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	60	1	114.99999999999999	42.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	61	1	187.0	75.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	62	1	319.0	95.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	63	2	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	64	2	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	65	1	315.0	97.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	66	2	0.0	0.0	0.0	0.0	3	1.0	0.0	380.0	1	1.1	0.9				
	67	2	0.0	0.0	0.0	0.0	4	1.0	0.0	380.0	1	1.1	0.9				
];

%% generator data
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf
mpc.gen = [
	1	700.0	23.0	1000.0	-500.0	1.0526	100.0	1	1000.0	400.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	2	1500.0	100.0	100.0	100.0	1.0526	100.0	1	1500.0	1500.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	4	523.0	140.0	350.0	-350.0	1.0526	100.0	1	560.0	220.00000000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	5	1200.0	100.0	100.0	100.0	1.0526	100.0	1	1200.0	1200.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	10	436.00000000000006	105.0	350.0	-350.0	1.0526	100.0	1	560.0	220.00000000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	13	541.0	117.0	300.0	-300.0	1.0526	100.0	1	630.0	250.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	18	681.0	-35.0	400.0	-400.0	1.0263	100.0	1	720.0	300.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	25	469.00000000000006	59.0	250.0	-250.0	1.0395	100.0	1	560.0	220.00000000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	29	500.0	101.0	350.0	-350.0	1.0263	100.0	1	630.0	250.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	33	496.0	306.0	500.0	-500.0	1.0263	100.0	1	850.0	350.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	36	512.0	249.00000000000003	400.0	-400.0	1.0263	100.0	1	720.0	300.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	41	350.0	238.0	450.0	-450.0	1.0395	100.0	1	850.0	350.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	43	574.0	223.0	500.0	-250.0	1.0263	100.0	1	720.0	220.00000000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	50	581.0	150.0	400.0	-400.0	1.0395	100.0	1	720.0	300.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	56	496.0	56.00000000000001	250.0	-250.0	1.0395	100.0	1	560.0	220.00000000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	59	430.99999999999994	206.0	350.0	-350.0	1.0447	100.0	1	720.0	300.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	63	488.0	152.0	250.0	-300.0	1.0447	100.0	1	520.0	250.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	64	300.0	61.0	250.0	-400.0	1.0474	100.0	1	560.0	250.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	66	537.0	143.0	300.0	-400.0	1.0447	100.0	1	630.0	250.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	67	800.0	0.0	0.0	0.0	1.0526	100.0	1	800.0	800.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
];

%% branch data
%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax
mpc.branch = [
	1	5	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	1	7	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	1	8	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	1	14	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	2	3	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	2	9	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	2	12	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	3	4	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	3	10	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	3	12	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	3	9	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	4	14	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	4	19	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	5	6	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	5	7	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	5	8	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	6	7	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	7	15	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	7	16	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	8	9	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	10	11	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	10	22	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	11	12	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	11	13	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	12	13	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	13	53	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	14	15	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	14	18	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	16	17	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	16	18	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	17	24	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	18	24	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	19	20	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	19	23	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	20	21	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	21	25	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	21	22	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	21	23	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	22	25	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	18	20	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	24	49	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	25	43	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	26	27	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	26	31	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	26	40	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	27	28	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	28	35	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	28	37	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	29	39	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	29	44	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	30	31	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	31	27	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	30	26	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	32	40	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	41	40	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	43	44	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	33	51	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	33	34	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	34	51	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	35	33	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	35	36	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	35	47	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	29	35	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	36	37	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	36	38	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	37	38	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	39	40	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	39	43	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	41	42	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	42	43	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	42	49	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	43	49	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	44	45	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	44	48	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	45	46	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	45	50	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	47	48	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	46	48	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	47	50	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	47	51	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	47	59	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	52	53	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	52	54	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	63	55	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	22	56	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	54	65	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	55	57	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	58	61	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	56	59	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	57	58	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	56	58	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	58	60	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	62	66	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	61	62	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	52	64	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	62	63	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	59	60	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	63	57	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	65	66	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	66	54	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	66	64	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
	30	32	0.00207756	0.01800554	0.61242207	900.0	900.0	900.0	0	0	1	-29.999999999999996	29.999999999999996								
];

%%-----  OPF Data  -----%%
%% cost data
%    1    startup    shutdown    n    x1    y1    ...    xn    yn
%    2    startup    shutdown    n    c(n-1)    ...    c0
mpc.gencost = [
	2	0.0	0.0	2	24.804734	0.0
	2	0.0	0.0	2	10.0	0.0
	2	0.0	0.0	2	24.804734	0.0
	2	0.0	0.0	2	34.844643	0.0
	2	0.0	0.0	2	10.0	0.0
	2	0.0	0.0	2	32.306483	0.0
	2	0.0	0.0	2	18.157477	0.0
	2	0.0	0.0	2	31.550181	0.0
	2	0.0	0.0	2	22.503168000000002	0.0
	2	0.0	0.0	2	27.434444	0.0
	2	0.0	0.0	2	34.844643	0.0
	2	0.0	0.0	2	14.707625	0.0
	2	0.0	0.0	2	24.804734	0.0
	2	0.0	0.0	2	34.844643	0.0
	2	0.0	0.0	2	24.652994	0.0
	2	0.0	0.0	2	32.306483	0.0
	2	0.0	0.0	2	18.157477	0.0
	2	0.0	0.0	2	31.550181	0.0
	2	0.0	0.0	2	22.503168000000002	0.0
	2	0.0	0.0	2	10.0	0.0
];

%column_names% μ dst_id λvmax σ λvmin 
mpc.bus_data = {
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	191.0	1	1.61	19.1	1.61
	0.0	0	1.61	0.0	1.61
	287.0	0	1.61	28.700000000000003	1.61
	186.0	0	1.61	18.6	1.61
	0.0	0	1.61	0.0	1.61
	271.0	0	1.61	27.1	1.61
	171.0	0	1.61	17.1	1.61
	0.0	0	1.61	0.0	1.61
	199.0	0	1.61	19.900000000000002	1.61
	112.99999999999999	0	1.61	11.299999999999999	1.61
	38.0	0	1.61	3.8000000000000003	1.61
	275.0	0	1.61	27.5	1.61
	0.0	0	1.61	0.0	1.61
	165.0	0	1.61	16.5	1.61
	178.0	2	1.61	17.8	1.61
	0.0	0	1.61	0.0	1.61
	30.0	0	1.61	3.0	1.61
	0.0	0	1.61	0.0	1.61
	32.0	0	1.61	3.2	1.61
	0.0	0	1.61	0.0	1.61
	395.0	0	1.61	39.5	1.61
	0.0	0	1.61	0.0	1.61
	665.0	0	1.61	66.5	1.61
	0.0	0	1.61	0.0	1.61
	266.0	0	1.61	26.6	1.61
	844.9999999999999	0	1.61	84.5	1.61
	332.0	0	1.61	33.2	1.61
	0.0	0	1.61	0.0	1.61
	540.0	0	1.61	54.0	1.61
	459.99999999999994	0	1.61	46.0	1.61
	0.0	0	1.61	0.0	1.61
	451.0	0	1.61	45.1	1.61
	150.0	0	1.61	15.0	1.61
	629.0	0	1.61	62.900000000000006	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	859.0	0	1.61	85.9	1.61
	0.0	0	1.61	0.0	1.61
	474.0	0	1.61	47.400000000000006	1.61
	668.0	0	1.61	66.8	1.61
	614.0	0	1.61	61.400000000000006	1.61
	81.0	0	1.61	8.1	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	430.0	0	1.61	43.0	1.61
	309.0	0	1.61	30.900000000000002	1.61
	100.0	0	1.61	10.0	1.61
	0.0	0	1.61	0.0	1.61
	303.0	0	1.61	30.3	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	324.0	0	1.61	32.4	1.61
	0.0	0	1.61	0.0	1.61
	114.99999999999999	0	1.61	11.5	1.61
	187.0	0	1.61	18.7	1.61
	319.0	0	1.61	31.900000000000002	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
	315.0	0	1.61	31.5	1.61
	0.0	0	1.61	0.0	1.61
	0.0	0	1.61	0.0	1.61
};

%column_names% qrated prated qref λqmax λpmax pref λpmin λqmin 
mpc.gen_data = {
	1000.0	1000.0	23.0	1.61	1.61	700.0	1.61	1.61
	100.0	1500.0	100.0	1.61	1.61	1500.0	1.61	1.61
	350.0	560.0	140.0	1.61	1.61	523.0	1.61	1.61
	100.0	1200.0	100.0	1.61	1.61	1200.0	1.61	1.61
	350.0	560.0	105.0	1.61	1.61	436.0	1.61	1.61
	300.0	630.0	117.0	1.61	1.61	541.0	1.61	1.61
	400.0	720.0	-35.0	1.61	1.61	681.0	1.61	1.61
	250.0	560.0	59.0	1.61	1.61	469.0	1.61	1.61
	350.0	630.0	101.0	1.61	1.61	500.0	1.61	1.61
	500.0	850.0	306.0	1.61	1.61	496.0	1.61	1.61
	400.0	720.0	249.0	1.61	1.61	512.0	1.61	1.61
	450.0	850.0	238.0	1.61	1.61	350.0	1.61	1.61
	500.0	720.0	223.0	1.61	1.61	574.0	1.61	1.61
	400.0	720.0	150.0	1.61	1.61	581.0	1.61	1.61
	250.0	560.0	56.0	1.61	1.61	496.0	1.61	1.61
	350.0	720.0	206.0	1.61	1.61	431.0	1.61	1.61
	250.0	520.0	152.0	1.61	1.61	488.0	1.61	1.61
	250.0	560.0	61.0	1.61	1.61	300.0	1.61	1.61
	300.0	630.0	143.0	1.61	1.61	537.0	1.61	1.61
	0.0	800.0	0.0	1.61	1.61	800.0	1.61	1.61
};

%column_names% tap_fr tap_fr_max shift_fr shift_to_min λcmax tap_to_max shift_to_max g_shunt shift_fr_max tap_to_min tap_to shift_to tappable shiftable tap_fr_min shift_fr_min 
mpc.branch_data = {
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
	1	1	0	0	1.61	1	0	0	0	1	1	0	0	0	1	0
};

mpc.dcpol = 2;

%column_names% dVdcset Vtar Pacmax filter reactor Vdcset Vmmax xtf Imax λvmax status Pdcset islcc LossA Qacmin rc rtf xc busdc_i busac_i tm type_dc Q_g LossB basekVac LossCrec droop Pacmin Qacmax type_ac Vmmin P_g transformer λvmin bf LossCinv 
mpc.convdc = {
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	1	7	1	2	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	-550	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	2	40	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	1000	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	3	3	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	-550	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	4	23	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	-600	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	5	48	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	1000	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	6	54	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	50	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	7	57	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	-550	0	1.61	0.01	0.0
	0	1	2000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	8	27	1	3	0	0.0	500	0.0	0.005	-2000	1000	1	0.95	1000	0	1.61	0.01	0.0
	0	1	1000	0	0	1.0	1.05	0.01	1.1	1.61	1	0.0	0	0.0	-1000	0.01	0.01	0.01	9	67	1	1	0	0.0	500	0.0	0.005	-1000	1000	1	0.95	-800	0	1.61	0.01	0.0
};

%column_names% qrated voll prated status pmax qref load_idx load_bus pref qmax qmin pmin 
mpc.sc_load = {
	0.0	5000.0	0.0	1	0.0	0.0	1	1	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	2	2	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	3	3	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	4	4	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	5	5	0.0	0.0	0.0	0.0
	76.0	5000.0	191.0	1	191.0	76.0	6	6	191.0	76.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	7	7	0.0	0.0	0.0	0.0
	73.0	5000.0	287.0	1	287.0	73.0	8	8	287.0	73.0	0.0	0.0
	74.0	5000.0	186.0	1	186.0	74.0	9	9	186.0	74.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	10	10	0.0	0.0	0.0	0.0
	55.0	5000.0	271.0	1	271.0	55.0	11	11	271.0	55.0	0.0	0.0
	87.0	5000.0	171.0	1	171.0	87.0	12	12	171.0	87.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	13	13	0.0	0.0	0.0	0.0
	60.0	5000.0	199.0	1	199.0	60.0	14	14	199.0	60.0	0.0	0.0
	52.5	5000.0	113.0	1	113.0	52.5	15	15	113.0	52.5	0.0	0.0
	7.0	5000.0	38.0	1	38.0	7.0	16	16	38.0	7.0	0.0	0.0
	106.0	5000.0	275.0	1	275.0	106.0	17	17	275.0	106.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	18	18	0.0	0.0	0.0	0.0
	46.0	5000.0	165.0	1	165.0	46.0	19	19	165.0	46.0	0.0	0.0
	82.5	5000.0	178.0	1	178.0	82.5	20	20	178.0	82.5	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	21	21	0.0	0.0	0.0	0.0
	7.0	5000.0	30.0	1	30.0	7.0	22	22	30.0	7.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	23	23	0.0	0.0	0.0	0.0
	7.0	5000.0	32.0	1	32.0	7.0	24	24	32.0	7.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	25	25	0.0	0.0	0.0	0.0
	89.0	5000.0	395.0	1	395.0	89.0	26	26	395.0	89.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	27	27	0.0	0.0	0.0	0.0
	99.0	5000.0	665.0	1	665.0	99.0	28	28	665.0	99.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	29	29	0.0	0.0	0.0	0.0
	100.0	5000.0	266.0	1	266.0	100.0	30	30	266.0	100.0	0.0	0.0
	119.0	5000.0	845.0	1	845.0	119.0	31	31	845.0	119.0	0.0	0.0
	137.0	5000.0	332.0	1	332.0	137.0	32	32	332.0	137.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	33	33	0.0	0.0	0.0	0.0
	158.0	5000.0	540.0	1	540.0	158.0	34	34	540.0	158.0	0.0	0.0
	97.0	5000.0	460.0	1	460.0	97.0	35	35	460.0	97.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	36	36	0.0	0.0	0.0	0.0
	190.0	5000.0	451.0	1	451.0	190.0	37	37	451.0	190.0	0.0	0.0
	0.0	5000.0	150.0	1	150.0	0.0	38	38	150.0	0.0	0.0	0.0
	87.0	5000.0	629.0	1	629.0	87.0	39	39	629.0	87.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	40	40	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	41	41	0.0	0.0	0.0	0.0
	180.0	5000.0	859.0	1	859.0	180.0	42	42	859.0	180.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	43	43	0.0	0.0	0.0	0.0
	92.0	5000.0	474.0	1	474.0	92.0	44	44	474.0	92.0	0.0	0.0
	109.0	5000.0	668.0	1	668.0	109.0	45	45	668.0	109.0	0.0	0.0
	95.0	5000.0	614.0	1	614.0	95.0	46	46	614.0	95.0	0.0	0.0
	0.0	5000.0	81.0	1	81.0	0.0	47	47	81.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	48	48	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	49	49	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	50	50	0.0	0.0	0.0	0.0
	123.0	5000.0	430.0	1	430.0	123.0	51	51	430.0	123.0	0.0	0.0
	102.0	5000.0	309.0	1	309.0	102.0	52	52	309.0	102.0	0.0	0.0
	30.0	5000.0	100.0	1	100.0	30.0	53	53	100.0	30.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	54	54	0.0	0.0	0.0	0.0
	110.0	5000.0	303.0	1	303.0	110.0	55	55	303.0	110.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	56	56	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	57	57	0.0	0.0	0.0	0.0
	157.0	5000.0	324.0	1	324.0	157.0	58	58	324.0	157.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	59	59	0.0	0.0	0.0	0.0
	42.0	5000.0	115.0	1	115.0	42.0	60	60	115.0	42.0	0.0	0.0
	75.0	5000.0	187.0	1	187.0	75.0	61	61	187.0	75.0	0.0	0.0
	95.0	5000.0	319.0	1	319.0	95.0	62	62	319.0	95.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	63	63	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	64	64	0.0	0.0	0.0	0.0
	97.0	5000.0	315.0	1	315.0	97.0	65	65	315.0	97.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	66	66	0.0	0.0	0.0	0.0
	0.0	5000.0	0.0	1	0.0	0.0	67	67	0.0	0.0	0.0	0.0
};

%column_names% dcconv_id2 dcbranch_id2 dcconv_id1 branch_id2 prob gen_id2 branch_id1 dcconv_id3 dcbranch_id1 dcbranch_id3 gen_id1 branch_id3 gen_id3 
mpc.contingencies = {
	0	0	0	0	0.98	0	0	0	0	0	0	0	0
	0	0	0	0	0.005	0	7	0	0	0	0	0	0
	0	0	0	0	0.005	0	13	0	0	0	0	0	0
	0	0	0	0	0.005	0	26	0	0	0	0	0	0
	0	0	0	0	0.005	0	41	0	0	0	0	0	0
	0	0	0	0	0.005	0	42	0	0	0	0	0	0
	0	0	0	0	0.005	0	45	0	0	0	0	0	0
	0	0	0	0	0.005	0	54	0	0	0	0	0	0
	0	0	0	0	0.005	0	81	0	0	0	0	0	0
	0	0	0	0	0.005	0	82	0	0	0	0	0	0
	0	0	0	0	0.005	0	85	0	0	0	0	0	0
	0	0	0	0	0.005	0	93	0	0	0	0	0	0
	0	0	0	0	0.005	0	0	0	1	0	0	0	0
	0	0	0	0	0.005	0	0	0	2	0	0	0	0
};

%column_names% col_3 col_4 col_1 col_2 
mpc.areas = {
	3	4	1	2
};

%column_names% CostCorr Converter_id CostPrev 
mpc.conv_cost = {
	10	1	1
	10	2	1
	10	3	1
	10	4	1
	10	5	1
	10	6	1
	10	7	1
	10	8	1
	10	9	1
};

%column_names% c r status λcmax rateB fbusdc rateA l rateC tbusdc 
mpc.branchdc = {
	0	0.0012	1	1.61	1575	1	1575	0	1575	2
	0	0.0012	1	1.61	1575	3	1575	0	1575	4
	0	0.0012	1	1.61	1575	4	1575	0	1575	5
	0	0.0012	1	1.61	1575	6	1575	0	1575	7
	0	0.0012	1	1.61	1575	1	1575	0	1575	3
	0	0.0012	1	1.61	1575	2	1575	0	1575	8
	0	0.0012	1	1.61	1575	8	1575	0	1575	5
	0	0.0012	1	1.61	1575	4	1575	0	1575	6
	0	0.0012	1	1.61	1575	2	1575	0	1575	4
	0	0.0012	1	1.61	1575	5	1575	0	1575	7
	0	0.0012	1	1.61	1575	3	1575	0	1575	9
};

%column_names% Vdc λvmax Cdc basekVdc busdc_i grid λvmin Vdcmax Vdcmin Pdc 
mpc.busdc = {
	1	1.61	0	500	1	1	1.61	1.05	0.95	0
	1	1.61	0	500	2	1	1.61	1.05	0.95	0
	1	1.61	0	500	3	1	1.61	1.05	0.95	0
	1	1.61	0	500	4	1	1.61	1.05	0.95	0
	1	1.61	0	500	5	1	1.61	1.05	0.95	0
	1	1.61	0	500	6	1	1.61	1.05	0.95	0
	1	1.61	0	500	7	1	1.61	1.05	0.95	0
	1	1.61	0	500	8	1	1.61	1.05	0.95	0
	1	1.61	0	500	9	1	1.61	1.05	0.95	0
};

