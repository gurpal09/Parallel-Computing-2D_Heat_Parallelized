
% Serial Code Timing Results
%Averages found in Excel
% # of TimeSteps: 1000
% Size:	 Time:
serial_data = [
                10	0.0005
                100	0.032966667
                500	1.305233333
                1000 4.650533333
                2000	11.28463333
                3000	23.50126667
                4000	39.72456667
                5000	63.14916667
                6000	88.79743333
                7000	118.9147
                8000	155.0277
                9000	195.8286333
                10000	240.9138
                12500	376.0246333
                15000	539.5445
                17500	733.3160333
                20000	983.6465];

            
% MPI Weak Scaling Averages			
% Size:	Time:	
mpi_weak_data = [10	0.00224	
                100	0.032393333	
                500	2.180736667	
                1000	4.951223333	
                2000	10.43328667	
                3000	15.62732	
                4000	27.48549667	
                5000	41.15737667	
                6000	55.17076333	
                7000	71.71282333	
                8000	90.3176	
                9000	112.77675	
                10000	136.6736133	
                12500	206.0556567	
                15000	290.6245633	
                17500	390.0485667	
                20000	505.5194433];	
            
            
% MPI Strong Scaling Averages            
% #Procs:	Time:
mpi_strong_data = [1	490.3145333
                   2	250.13038
                   3	171.0467433
                   4	136.5394133
                   8	77.68160333];

%Comparing it to serial data
%matrix size: 10,000 & 1000 Timesteps
serial_s = zeros(5,2);
serial_s(1,1) = 1;
serial_s(2,1) = 2;
serial_s(3,1) = 3;
serial_s(4,1) = 4;
serial_s(5,1) = 8;
serial_s(:,2) = serial_data(13,2);


figure
title('Serial Code Results (1000 Timesteps)')
xlabel('Matrix Size NX')
ylabel('Time (s)')
hold on
plot(serial_data(:,1),serial_data(:,2),'b-*')

figure
title('MPI Weak Scaling Study')
xlabel('Matrix Size NX')
ylabel('Time (s)')
hold on
plot(mpi_weak_data(:,1),mpi_weak_data(:,2),'b-*')
plot(serial_data(:,1),serial_data(:,2),'r-*')
legend('MPI Code', 'Serial Code')

figure
title('MPI Strong Scaling Study')
xlabel('Number of Processors')
ylabel('Time (s)')
hold on
plot(mpi_strong_data(:,1),mpi_strong_data(:,2),'b-*')
plot(serial_s(:,1),serial_s(:,2),'r-*')
legend('MPI Code', 'Serial Code')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Speed Up Analysis

% Size:	Speed Up
speed_data = [
10	0.223214286
100	1.017699115
500	0.598528632
1000	0.939269554
2000	1.081599087
3000	1.503857774
4000	1.445291935
5000	1.534334104
6000	1.609501627
7000	1.658206921
8000	1.716472758
9000	1.736427352
10000	1.762694306
12500	1.824869258
15000	1.856499994
17500	1.880063397
20000	1.945813387];

figure
fname = {'10','100','1000','2000','5000','10000','15000','20000'}; 
x = [1:8]; 
y = [speed_data(1,2),speed_data(2,2),speed_data(4,2),speed_data(5,2),speed_data(8,2),speed_data(13,2),speed_data(15,2),speed_data(17,2)]; 
title('Speed-Up Analysis: Increasing Problem Size')
ylabel('Speed-Up (T_{serial}/T_{parallel})')
xlabel('Matrix Size NX')
hold on
bar(x,y)
set(gca, 'XTick', 1:length(fname),'XTickLabel',fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Speed Up Increasing # of Procs
% # of Procs	Speed-Up
proc_data = [
1	0.491345419
2	0.963152897
3	1.408467623
4	1.764426799
8	3.101297986];

figure
pname = {'1','2','3','4','8'}; 
x = [1:5]; 
y1 = proc_data(:,2); 
title('Speed-Up Analysis: Increasing # of Procs')
ylabel('Speed-Up (T_{serial}/T_{parallel})')
xlabel('# of Processors')
hold on
bar(x,y1)
set(gca, 'XTick', 1:length(pname),'XTickLabel',pname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Efficiency
% # of Procs
E = [1	0.491345419
2	0.481576448
3	0.469489208
4	0.4411067
8	0.387662248];

figure
title('Efficiency')
xlabel('Number of Processors')
ylabel('Efficiency (S/p)')
hold on
plot(E(:,1), E(:,2),'b-*')

