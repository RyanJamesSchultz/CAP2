% Simple script to tabulate case data into a structure.

% Define the structure table.
T=struct('Name',[], 'Filename',[], 'Mc',[], 'Mk',[], 'Ts',NaT(0),'Te',NaT(0),'md',[]);
i=1;

%%% Field-scale cases in CAPfr.

% UK-PNR1z.
T(i).Name='PNR1z';
T(i).Filename='data/PNR1z.mat';
T(i).Mc=-0.8;
T(i).Mk=-0.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2018,10,16, 10,26,24); T(i).Te( 1)=datetime(2018,10,16, 14,06,09); T(i).md( 1)=mean([3374 3376]); % S1.
T(i).Ts( 2)=datetime(2018,10,17, 08,12,49); T(i).Te( 2)=datetime(2018,10,17, 16,35,51); T(i).md( 2)=mean([3356 3358]); % S2.
T(i).Ts( 3)=datetime(2018,10,18, 09,55,20); T(i).Te( 3)=datetime(2018,10,18, 14,21,58); T(i).md( 3)=mean([3338 3340]); % S3.
T(i).Ts( 4)=datetime(2018,10,19, 11,32,45); T(i).Te( 4)=datetime(2018,10,19, 14,39,35); T(i).md( 4)=mean([3180 3182]); % S12.
T(i).Ts( 5)=datetime(2018,10,20, 08,06,18); T(i).Te( 5)=datetime(2018,10,20, 10,40,12); T(i).md( 5)=mean([3180 3182]); % S12.
T(i).Ts( 6)=datetime(2018,10,22, 07,46,22); T(i).Te( 6)=datetime(2018,10,22, 11,53,45); T(i).md( 6)=mean([3162 3164]); % S13.
T(i).Ts( 7)=datetime(2018,10,23, 13,32,55); T(i).Te( 7)=datetime(2018,10,23, 14,58,33); T(i).md( 7)=mean([3144 3146]); % S14.
T(i).Ts( 8)=datetime(2018,10,24, 08,09,04); T(i).Te( 8)=datetime(2018,10,24, 12,30,13); T(i).md( 8)=mean([3074 3076]); % S18.
T(i).Ts( 9)=datetime(2018,10,25, 07,39,26); T(i).Te( 9)=datetime(2018,10,25, 11,44,21); T(i).md( 9)=mean([2997 3000]); % S22.
T(i).Ts(10)=datetime(2018,10,25, 14,29,39); T(i).Te(10)=datetime(2018,10,25, 15,42,08); T(i).md(10)=mean([2997 3000]); % S22.
T(i).Ts(11)=datetime(2018,10,26, 07,19,27); T(i).Te(11)=datetime(2018,10,26, 11,41,12); T(i).md(11)=mean([2857 2859]); % S30.
T(i).Ts(12)=datetime(2018,10,27, 08,04,04); T(i).Te(12)=datetime(2018,10,27, 10,53,30); T(i).md(12)=mean([2840 2842]); % S31.
T(i).Ts(13)=datetime(2018,10,29, 08,53,08); T(i).Te(13)=datetime(2018,10,29, 11,39,46); T(i).md(13)=mean([2822 2824]); % S32.
T(i).Ts(14)=datetime(2018,10,30, 08,26,53); T(i).Te(14)=datetime(2018,10,30, 10,07,28); T(i).md(14)=mean([2699 2701]); % S39.
T(i).Ts(15)=datetime(2018,10,31, 09,28,27); T(i).Te(15)=datetime(2018,10,31, 13,34,31); T(i).md(15)=mean([2681 2683]); % S40.
T(i).Ts(16)=datetime(2018,11,02, 08,13,32); T(i).Te(16)=datetime(2018,11,02, 16,43,19); T(i).md(16)=mean([2769 2771]); % S35.
%%%
T(i).Ts(17)=datetime(2018,12,08, 09,24,47); T(i).Te(17)=datetime(2018,12,08, 12,14,03); T(i).md(17)=mean([2734 2736]); % S37.
T(i).Ts(18)=datetime(2018,12,10, 08,34,18); T(i).Te(18)=datetime(2018,12,10, 10,04,58); T(i).md(18)=mean([2734 2736]); % S37.
T(i).Ts(19)=datetime(2018,12,11, 08,32,44); T(i).Te(19)=datetime(2018,12,11, 10,36,32); T(i).md(19)=mean([2716 2718]); % S38.
T(i).Ts(20)=datetime(2018,12,13, 08,30,10); T(i).Te(20)=datetime(2018,12,13, 11,33,42); T(i).md(20)=mean([2699 2701]); % S39.
T(i).Ts(21)=datetime(2018,12,14, 12,27,18); T(i).Te(21)=datetime(2018,12,14, 14,04,44); T(i).md(21)=mean([2681 2683]); % S40.
T(i).Ts(22)=datetime(2018,12,15, 08,59,21); T(i).Te(22)=datetime(2018,12,15, 13,00,13); T(i).md(22)=mean([2663 2665]); % S41.
T(i).Ts(23)=datetime(2018,12,17, 08,13,24); T(i).Te(23)=datetime(2018,12,17, 13,12,47); T(i).md(23)=mean([2663 2665]); % S41.
i=i+1;

% UK-PNR2.
T(i).Name='PNR2';
T(i).Filename='data/PNR2.mat';
T(i).Mc=-1.0;
T(i).Mk=-0.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2019,08,15, 09,39,13); T(i).Te(1)=datetime(2019,08,15, 10,50,48); T(i).md(1)=mean([3106 3111]); % S1.
T(i).Ts(2)=datetime(2019,08,16, 08,01,38); T(i).Te(2)=datetime(2019,08,16, 10,25,14); T(i).md(2)=mean([3092 3097]); % S2.
T(i).Ts(3)=datetime(2019,08,17, 07,59,42); T(i).Te(3)=datetime(2019,08,17, 11,19,49); T(i).md(3)=mean([3077 3082]); % S3.
T(i).Ts(4)=datetime(2019,08,19, 08,12,10); T(i).Te(4)=datetime(2019,08,19, 10,49,39); T(i).md(4)=mean([3063 3068]); % S4.
T(i).Ts(5)=datetime(2019,08,20, 07,57,17); T(i).Te(5)=datetime(2019,08,20, 10,35,01); T(i).md(5)=mean([3048 3053]); % S5.
T(i).Ts(6)=datetime(2019,08,21, 08,36,50); T(i).Te(6)=datetime(2019,08,21, 15,53,28); T(i).md(6)=mean([3033 3038]); % S6.
T(i).Ts(7)=datetime(2019,08,23, 10,21,35); T(i).Te(7)=datetime(2019,08,23, 14,44,56); T(i).md(7)=mean([3019 3024]); % S7.
i=i+1;

% Utah-FORGE 2022.
T(i).Name='FORGE22';
T(i).Filename='data/FORGE2022.mat';
T(i).Mc=-1.3;
T(i).Mk=-1.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2022,04,17, 02,41,36); T(i).Te(1)=datetime(2022,04,17, 05,16,48); T(i).md(1)=mean([10787 10955])*0.3048; % S1.
T(i).Ts(2)=datetime(2022,04,19, 12,50,01); T(i).Te(2)=datetime(2022,04,19, 15,42,50); T(i).md(2)=mean([10560 10580])*0.3048; % S2.
T(i).Ts(3)=datetime(2022,04,21, 13,33,02); T(i).Te(3)=datetime(2022,03,21, 16,19,24); T(i).md(3)=mean([10120 10140])*0.3048; % S3.
i=i+1;

% Utah-FORGE 2024.
T(i).Name='FORGE24';
T(i).Filename='/Users/rschultz/Desktop/CAPfr/codes/data/FORGE2024.mat';
T(i).Mc=+0.25;
T(i).Mk=+0.30;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2024,04,03, 16,15,06); T(i).Te( 1)=datetime(2024,04,03, 19,56,56); T(i).md( 1)=mean([10120 10140])*0.3048; % S3a(r).
T(i).Ts( 2)=datetime(2024,04,03, 23,56,05); T(i).Te( 2)=datetime(2024,04,04, 03,01,50); T(i).md( 2)=mean([10070 10076])*0.3048; % S4a.
T(i).Ts( 3)=datetime(2024,04,04, 06,47,14); T(i).Te( 3)=datetime(2024,04,04, 09,20,52); T(i).md( 3)=mean([10020 10026])*0.3048; % S5a.
T(i).Ts( 4)=datetime(2024,04,04, 12,31,33); T(i).Te( 4)=datetime(2024,04,04, 22,46,44); T(i).md( 4)=mean([ 9959  9976])*0.3048; % S6a.
T(i).Ts( 5)=datetime(2024,04,05, 06,33,34); T(i).Te( 5)=datetime(2024,04,05, 15,12,12); T(i).md( 5)=mean([ 9798  9901])*0.3048; % S7a.
T(i).Ts( 6)=datetime(2024,04,05, 19,45,07); T(i).Te( 6)=datetime(2024,04,06, 04,01,09); T(i).md( 6)=mean([ 9545  9723])*0.3048; % S8a.
T(i).Ts( 7)=datetime(2024,04,06, 07,00,42); T(i).Te( 7)=datetime(2024,04,06, 13,12,28); T(i).md( 7)=mean([ 9490  9493])*0.3048; % S9a.
T(i).Ts( 8)=datetime(2024,04,07, 02,57,28); T(i).Te( 8)=datetime(2024,04,07, 05,20,44); T(i).md( 8)=mean([ 9270  9323])*0.3048; % S10a.
T(i).Ts( 9)=datetime(2024,04,11, 10,01,52); T(i).Te( 9)=datetime(2024,04,11, 11,12,48); T(i).md( 9)=mean([ ])*0.3048; % S1b.
T(i).Ts(10)=datetime(2024,04,13, 17,10,25); T(i).Te(10)=datetime(2024,04,13, 18,38,43); T(i).md(10)=mean([ ])*0.3048; % S2b.
T(i).Ts(11)=datetime(2024,04,15, 20,24,58); T(i).Te(11)=datetime(2024,04,15, 21,53,12); T(i).md(11)=mean([ ])*0.3048; % S3b.
T(i).Ts(12)=datetime(2024,04,17, 12,02,55); T(i).Te(12)=datetime(2024,04,17, 13,13,02); T(i).md(12)=mean([ ])*0.3048; % S4b.
i=i+1;


%%% Underground laboratories in CAPugl.

% Äspö HRL.
T(i).Name='Aspo';
T(i).Filename='data/Aspo.mat';
T(i).Mc=-4.9;
T(i).Mk=-4.9;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2015,06,03, 09,40,33); T(i).Te(1)=datetime(2015,06,03, 17,00,00); T(i).md(1)=24.80; % HF1.
T(i).Ts(2)=datetime(2015,06,04, 07,15,07); T(i).Te(2)=datetime(2015,06,04, 09,30,00); T(i).md(2)=22.00; % HF2.
T(i).Ts(3)=datetime(2015,06,04, 12,07,02); T(i).Te(3)=datetime(2015,06,04, 14,30,00); T(i).md(3)=19.00; % HF3.
T(i).Ts(4)=datetime(2015,06,09, 10,03,05); T(i).Te(4)=datetime(2015,06,09, 14,00,00); T(i).md(4)=13.65; % HF4.
T(i).Ts(5)=datetime(2015,06,10, 10,35,10); T(i).Te(5)=datetime(2015,06,10, 14,00,00); T(i).md(5)=11.80; % HF5.
T(i).Ts(6)=datetime(2015,06,11, 09,29,03); T(i).Te(6)=datetime(2015,06,12, 12,00,00); T(i).md(6)=04.68; % HF6.
i=i+1;

% SURF EGS Collab Experiment #1.
T(i).Name='SURF1';
T(i).Filename='data/SURF1.mat';
T(i).Mc=-3.5;
T(i).Mk=-3.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2018,05,21, 00,00,00); T(i).Te(1)=datetime(2018,05,22, 00,00,00); T(i).md(1)=142*0.3048; % N142-0521.
T(i).Ts(2)=datetime(2018,05,22, 00,00,00); T(i).Te(2)=datetime(2018,05,26, 00,00,00); T(i).md(2)=164*0.3048; % N164-0522.
T(i).Ts(3)=datetime(2018,06,25, 00,00,00); T(i).Te(3)=datetime(2018,06,26, 00,00,00); T(i).md(3)=164*0.3048; % N164-0625i.
T(i).Ts(4)=datetime(2018,07,18, 00,00,00); T(i).Te(4)=datetime(2018,07,21, 00,00,00); T(i).md(4)=128*0.3048; % N128-0718.
T(i).Ts(5)=datetime(2018,12,07, 00,00,00); T(i).Te(5)=datetime(2018,12,08, 00,00,00); T(i).md(5)=142*0.3048; % N142-1207.
T(i).Ts(6)=datetime(2018,12,20, 16,47,51); T(i).Te(6)=datetime(2018,12,22, 00,00,00); T(i).md(6)=142*0.3048; % N142-1220b.
i=i+1;

% Grimsel HF2 stimulation.
T(i).Name='GTS-HF2';
T(i).Filename='data/GTS_HF2.mat';
T(i).Mc=-5.1;
T(i).Mk=-5.1;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,05,17, 08,40,31); T(i).Te(1)=datetime(2017,05,17, 08,50,54); T(i).md(1)=mean([35.8 36.8]); % c1.
T(i).Ts(2)=datetime(2017,05,17, 08,57,42); T(i).Te(2)=datetime(2017,05,17, 09,07,44); T(i).md(2)=mean([35.8 36.8]); % c2.
i=i+1;

% Grimsel HF3 stimulation.
T(i).Name='GTS-HF3';
T(i).Filename='data/GTS_HF3.mat';
T(i).Mc=-4.6;
T(i).Mk=-4.6;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,05,16, 14,00,00); T(i).Te(1)=datetime(2017,05,16, 14,11,19); T(i).md(1)=mean([19.8 20.8]); % c1.
T(i).Ts(2)=datetime(2017,05,16, 14,21,22); T(i).Te(2)=datetime(2017,05,16, 14,31,38); T(i).md(2)=mean([19.8 20.8]); % c2.
i=i+1;

% Grimsel HF5 stimulation.
T(i).Name='GTS-HF5';
T(i).Filename='data/GTS_HF5.mat';
T(i).Mc=-5.0;
T(i).Mk=-5.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,05,17, 11,30,00); T(i).Te(1)=datetime(2017,05,17, 12,45,00); T(i).md(1)=mean([14 15]); % c1.
i=i+1;

% Grimsel HF6 stimulation.
T(i).Name='GTS-HF6';
T(i).Filename='data/GTS_HF6.mat';
T(i).Mc=-5.0;
T(i).Mk=-5.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,05,18, 07,30,00); T(i).Te(1)=datetime(2017,05,18, 12,30,00); T(i).md(1)=mean([38.4 39.4]); % c1.
i=i+1;

% Grimsel HF8 stimulation.
T(i).Name='GTS-HF8';
T(i).Filename='data/GTS_HF8.mat';
T(i).Mc=-4.7;
T(i).Mk=-4.7;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,05,18, 13,00,00); T(i).Te(1)=datetime(2017,05,18, 15,30,00); T(i).md(1)=mean([15.2 16.2]); % c1.
i=i+1;

% Grimsel HS1 stimulation.
T(i).Name='GTS-HS1';
T(i).Filename='data/GTS_HS1.mat';
T(i).Mc=-4.1;
T(i).Mk=-4.1;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,15, 10,30,00); T(i).Te(1)=datetime(2017,02,15, 15,30,00); T(i).md(1)=mean([39.75 40.75]); % c1.
i=i+1;

% Grimsel HS2 stimulation.
T(i).Name='GTS-HS2';
T(i).Filename='data/GTS_HS2.mat';
T(i).Mc=-4.3;
T(i).Mk=-4.3;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,08, 08,30,00); T(i).Te(1)=datetime(2017,02,08, 13,00,00); T(i).md(1)=mean([38 40]); % c1.
i=i+1;

% Grimsel HS3 stimulation.
T(i).Name='GTS-HS3';
T(i).Filename='data/GTS_HS3.mat';
T(i).Mc=-5.2;
T(i).Mk=-5.2;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,13, 11,13,35); T(i).Te(1)=datetime(2024,04,03, 19,56,56); T(i).md(1)=mean([34.3 35.3]); % c1.
i=i+1;

% Grimsel HS4 stimulation.
T(i).Name='GTS-HS4';
T(i).Filename='data/GTS_HS4.mat';
T(i).Mc=-5.0;
T(i).Mk=-5.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,09, 10,24,21); T(i).Te(1)=datetime(2024,02,09, 11,02,32); T(i).md(1)=mean([27.2 28.2]); % c1.
T(i).Ts(2)=datetime(2017,02,09, 11,36,30); T(i).Te(2)=datetime(2024,02,09, 12,14,33); T(i).md(2)=mean([27.2 28.2]); % c2.
T(i).Ts(3)=datetime(2017,02,09, 13,09,53); T(i).Te(3)=datetime(2024,02,09, 13,52,13); T(i).md(3)=mean([27.2 28.2]); % c3.
T(i).Ts(4)=datetime(2017,02,09, 15,08,54); T(i).Te(4)=datetime(2024,02,09, 15,52,06); T(i).md(4)=mean([27.2 28.2]); % c4.
i=i+1;

% Grimsel HS5 stimulation.
T(i).Name='GTS-HS5';
T(i).Filename='data/GTS_HS5.mat';
T(i).Mc=-5.0;
T(i).Mk=-5.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,10, 09,07,39); T(i).Te(1)=datetime(2017,02,10, 09,46,29); T(i).md(1)=mean([31.2 32.2]); % c1.
T(i).Ts(2)=datetime(2017,02,10, 10,21,05); T(i).Te(2)=datetime(2017,02,10, 10,53,37); T(i).md(2)=mean([31.2 32.2]); % c2.
T(i).Ts(3)=datetime(2017,02,10, 11,30,07); T(i).Te(3)=datetime(2017,02,10, 12,01,09); T(i).md(3)=mean([31.2 32.2]); % c3.
T(i).Ts(4)=datetime(2017,02,10, 13,17,48); T(i).Te(4)=datetime(2017,02,10, 14,26,39); T(i).md(4)=mean([31.2 32.2]); % c4.
i=i+1;

% Grimsel HS8 stimulation.
T(i).Name='GTS-HS8';
T(i).Filename='data/GTS_HS8.mat';
T(i).Mc=-4.7;
T(i).Mk=-4.7;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2017,02,14, 10,30,00); T(i).Te(1)=datetime(2017,02,14, 16,00,00); T(i).md(1)=mean([22.0 23.0]); % c1.
i=i+1;


%%% Field scale cases NLE.

% Basel.
T(i).Name='Basel';
T(i).Filename='data/BASEL.mat';
T(i).Mc=1.5;
T(i).Mk=1.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Cooper Basin, Habanero 1.
T(i).Name='CB1';
T(i).Filename='data/CB_HAB1.mat';
T(i).Mc=0.0;
T(i).Mk=1.6;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2003,11,09, 00,41,09); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Cooper Basin, Habanero 4.
T(i).Name='CB4';
T(i).Filename='data/CB_HAB4.mat';
T(i).Mc=2.0;
T(i).Mk=2.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Paralana
T(i).Name='Paralana';
T(i).Filename='data/Paralana.mat';
T(i).Mc=0.5;
T(i).Mk=1.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Pohang.
T(i).Name='Pohang';
T(i).Filename='data/Pohang.mat';
T(i).Mc=0.5;
T(i).Mk=0.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2015,09,20, 23,57,19); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Duvernay ESB-10.
T(i).Name='ESB10';
T(i).Filename='data/ESB10.mat';
T(i).Mc=0.8;
T(i).Mk=0.8;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Peace River.
T(i).Name='PeaceRiver';
T(i).Filename='data/PR.mat';
T(i).Mc=2.0;
T(i).Mk=2.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Paradox Valley.
T(i).Name='ParadoxValley';
T(i).Filename='data/ParaV.mat';
T(i).Mc=1.5;
T(i).Mk=2.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Guy-Greenbrier.
T(i).Name='GGB';
T(i).Filename='data/GGB.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Prague.
T(i).Name='Prague';
T(i).Filename='data/PRAGUE_cat.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Pawnee.
T(i).Name='Pawnee';
T(i).Filename='data/PAWNEE_cat.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;

% Youngstown.
T(i).Name='Youngstown';
T(i).Filename='data/Youngstown.mat';
T(i).Mc=1.4;
T(i).Mk=1.4;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(1900,01,01, 00,00,00); T(i).Te(1)=datetime(2900,01,01, 00,00,00); T(i).md(1)=mean([ ]);
i=i+1;
