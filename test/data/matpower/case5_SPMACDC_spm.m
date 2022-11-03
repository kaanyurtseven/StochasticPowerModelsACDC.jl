%% MATPOWER Case Format : Version 2
function mpc = pglib_opf_case5_pjm
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
	1	3	0.0	0.0	0.0	0.0	1	1.06	0.0	345.0	1	1.1	0.9				
	2	2	20.0	10.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	3	1	45.0	15.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	4	1	40.0	5.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	5	1	60.0	10.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
];

%% generator data
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf
mpc.gen = [
	1	0.0	0.0	500.0	-500.0	1.06	100.0	1	250.0	10.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	2	40.0	0.0	300.0	-300.0	1.0	100.0	1	300.0	10.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
];

%% branch data
%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax
mpc.branch = [
	1	2	0.02	0.06	0.06	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	1	3	0.08	0.24	0.05	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	3	0.06	0.18	0.04	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	4	0.06	0.18	0.04	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	5	0.04	0.12	0.03	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	3	4	0.01	0.03	0.02	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	4	5	0.08	0.24	0.05	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
];

%%-----  OPF Data  -----%%
%% cost data
%    1    startup    shutdown    n    x1    y1    ...    xn    yn
%    2    startup    shutdown    n    c(n-1)    ...    c0
mpc.gencost = [
	2	0.0	0.0	2	1.0	0.0
	2	0.0	0.0	2	2.0	0.0
];

%column_names% μ dst_id λvmax σ λvmin 
mpc.bus_data = {
	0.0	0	1.03643	0.0	1.03643
	20.0	0	1.03643	2.0	1.03643
	45.0	0	1.03643	4.5	1.03643
	40.0	0	1.03643	4.0	1.03643
	60.0	0	1.03643	6.0	1.03643
};

%column_names% λqmax λpmax λpmin λqmin 
mpc.gen_data = {
	1.03643	1.03643	1.03643	1.03643
	1.03643	1.03643	1.03643	1.03643
};

%column_names% λcmax c_rating_a 
mpc.branch_data = {
	1.03643	100.0
	1.03643	100.0
	1.03643	100.0
	1.03643	100.0
	1.03643	100.0
	1.03643	100.0
	1.03643	100.0
};

%column_names% λvmax μ λvmin σ dst_id 
mpc.convdc_sdata = {
	1.03643	0.0	1.03643	0.0	0
	1.03643	0.0	1.03643	0.0	0
	1.03643	0.0	1.03643	0.0	0
};

%column_names% pa dst pb 
mpc.sdata = {
	2.0	'Beta'	2.0
	0.0	'Normal'	0.0
};

mpc.dcpol = 2;

%column_names% λvmax μ λvmin σ dst_id 
mpc.busdc_sdata = {
	1.03643	0.0	1.03643	0.0	0
	1.03643	0.0	1.03643	0.0	0
	1.03643	0.0	1.03643	0.0	0
};

%column_names% dVdcset μ Vtar Pacmax filter reactor Vdcset Vmmax dst_id xtf Imax λvmax status σ Pdcset islcc LossA Qacmin rc rtf xc busdc_i busac_i tm type_dc Q_g LossB basekVac LossCrec droop Pacmin Qacmax type_ac Vmmin P_g transformer λvmin bf LossCinv 
mpc.convdc = {
	0	0.0	1	100	1	1	1.0079	1.1	0	0.01	1.1	1.03643	1	0.0	-58.6274	0	1.103	-50	0.01	0.01	0.01	1	2	1	1	-40	0.887	345	2.885	0.005	-100	50	1	0.9	-60	1	1.03643	0.01	2.885
	0	20.0	1	100	1	1	1.0	1.1	0	0.01	1.1	1.03643	1	2.0	21.9013	0	1.103	-50	0.01	0.01	0.01	2	3	1	2	0	0.887	345	2.885	0.007	-100	50	1	0.9	0	1	1.03643	0.01	2.885
	0	45.0	1	100	1	1	0.9978	1.1	0	0.01	1.1	1.03643	1	4.5	36.1856	0	1.103	-50	0.01	0.01	0.01	3	5	1	1	5	0.887	345	2.885	0.005	-100	50	1	0.9	35	1	1.03643	0.01	2.885
};

%column_names% λcmax 
mpc.branchdc_sdata = {
	1.03643
	1.03643
	1.03643
};

%column_names% col_1 col_2 
mpc.areas = {
	1	4
};

%column_names% c r status λcmax rateB fbusdc rateA l rateC tbusdc 
mpc.branchdc = {
	0	0.052	1	1.03643	100	1	100	0	100	2
	0	0.052	1	1.03643	100	2	100	0	100	3
	0	0.073	1	1.03643	100	1	100	0	100	3
};

%column_names% μ Vdc dst_id λvmax Cdc σ basekVdc busdc_i grid λvmin Vdcmax Vdcmin Pdc 
mpc.busdc = {
	0.0	1	0	1.03643	0	0.0	345	1	1	1.03643	1.1	0.9	0
	20.0	1	0	1.03643	0	2.0	345	2	1	1.03643	1.1	0.9	0
	45.0	1	0	1.03643	0	4.5	345	3	1	1.03643	1.1	0.9	0
};

