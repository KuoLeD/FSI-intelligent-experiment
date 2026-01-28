function output4 = Step4_DealCreatCoe(stepn,filenametxt,filenamecsv,filenamestorematlab)
%% 4 Process the experimental data in casename.txt and build the hydrodynamic-coefficient training database (DPQ)
% Note: In experimental data processing, the raw data should be evaluated first, and the stable segment data should be selected.
% After processing the data and obtaining the coefficient results, a preliminary assessment of the rationality of the coefficient results is required.
% Import the hydrodynamic coefficient results into the initial_data file.
output4=stepn;
%% Read txt data
input0=importdata(filenametxt);
tinterval=input0.data(23);
err_coe=input0.data(24);
err_coeType=input0.data(25);
ytrain_direction=input0.data(26);
ytrain_type=input0.data(27);
Coe_Zhengzhi=input0.data(31);
% Read csv data
data = readtable(filenamecsv);
case_No = table2array(data (:,1));
case_name = table2array(data (:,2));
case_record = table2array(data (:,3));
csv_input0 = table2array(data (:,4:13));
case_direction = table2array(data (:,14));
case_count = table2array(data (:,15));
csv_coe = table2array(data (:,16:20));
load('information.mat');
if sum(csv_input0(:,1))~=0
    CF_CL=1;
else
    CF_CL=0;
end

if sum(csv_input0(:,3))~=0
    IL_CL=1;
else
    IL_CL=0;
end
%% Process experimental data and obtain hydrodynamic coefficients
Train_data1=zeros(n/2,10); % forward
Train_data0=zeros(n/2,10); % backward
err_coe_ZF=zeros(n/2,5);
% if n0 ~= 1
%     train1 = load('Train_data1.mat');
%     Train_data1(1:(n0-1)/2,:) = train1.Train_data1;
%     train0 = load('Train_data0.mat');
%     Train_data0(1:(n0-1)/2,:) = train0.Train_data0;
% end

for im=1:length(m)
    i=m(im);
    filelist = dir([filenamestorematlab case_name{i} '.txt']);
    filedate = [filelist.datenum];
    [~,newestIdx] = max(filedate);
    if ~isempty(newestIdx)
        newestfilename = filelist(newestIdx).name;
    else
        fprintf('Missing file: %s\n', case_name{i});
    end
    newestfilename = filelist(newestIdx).name;
    data = readmatrix(newestfilename); % There is a file naming issue to be resolved here
    if mod(i,2)==1
        i1=(i+1)/2;
        Us = csv_input0(i,10);
        As = csv_input0(i,5);
        fs = csv_input0(i,6);
        Coe1 = DataDeal(data,Us,As,fs,L,D,M,tinterval,CF_CL,IL_CL); 
        Train_data1(i1,1:6) = csv_input0(i,5:10);
        Train_data1(i1,7:end) = Coe1(1:4); % output 4 hydrodynamic coefficients
        csv_coe(i,:) = Coe1;
    else
        i0=i/2;
        Us = csv_input0(i,10);
        As = csv_input0(i,5);
        fs = csv_input0(i,6);
        Coe0 = DataDeal(data,Us,As,fs,L,D,M,tinterval,CF_CL,IL_CL); 
        Train_data0(i0,1:6) = csv_input0(i,5:10);
        Train_data0(i0,7:end) = Coe0(1:4);
        csv_coe(i,:) = Coe0;
        % Evaluate whether the forward and backward datasets are close; use the forward direction as reference and compute error ratio
        err_coe_ZF(i0,1:5) = abs((Coe1-Coe0)./Coe1);
        switch err_coeType
            case 1 % based on CF Ce
                err_coe_cal=err_coe_ZF(i0,1);
            case 2 % based on CF Ca
                err_coe_cal=err_coe_ZF(i0,2);
            case 3 % based on IL Ce
                err_coe_cal=err_coe_ZF(i0,3);
            case 4 % based on IL Ca
                err_coe_cal=err_coe_ZF(i0,4);
            case 5 % combined CF and IL
                err_coe_cal=sum(err_coe_ZF(i0,1:4))/4;
            case 6 % combined CF
                err_coe_cal=sum(err_coe_ZF(i0,1:2))/2;
            case 7 % combined IL
                err_coe_cal=sum(err_coe_ZF(i0,3:4))/2;
            case 8 % based on Cd
                err_coe_cal=err_coe_ZF(i0,5);
        end
        if err_coe_cal <= err_coe % error controlled within 5%
            err_coe_ZF(i0,5) = 1;
            case_record(i-1:i) = [1;1];
        else
            case_record(i-1:i) = [0;0]; % mark as unqualified and redo
            % Check whether it has already been redone for the second time; process the data that has been redone twice
            if length(filelist) > 1
                Coe0_repeat = zeros(length(filelist),5);
                Coe1_repeat = zeros(length(filelist),5);
                filelist0 =filelist;
                filelist1 = dir([filenamestorematlab case_name{i-1} '.txt']);
                
                for ifile=1:length(filelist1) % number of forward cases (in this loop it must be forward)
                    filename1 = filelist1(ifile).name;
                    datafile1 = readmatrix(filename1);
                    Coe1_repeat(ifile,:) = DataDeal(datafile1,Us,As,fs,L,D,M,tinterval,CF_CL,IL_CL);
                end
                for ifile=1:length(filelist0) % number of backward cases (in this loop it must be backward)
                    filename0 = filelist0(ifile).name;
                    datafile0 = readmatrix(filename0);
                    Coe0_repeat(ifile,:) = DataDeal(datafile0,Us,As,fs,L,D,M,tinterval,CF_CL,IL_CL);
                end
                % Set up a loop to compute coefficient error under all combinations
                ii=1;
                err_coe_ZF1=zeros(length(filelist1)*length(filelist0),5);
                for ifile1=1:length(filelist1)
                    for ifile0=1:length(filelist0)
                        err_coe_ZF1(ii,1:5)=abs((Coe1_repeat(ifile1,:)-Coe0_repeat(ifile0,:))./Coe1_repeat(ifile1,:));
                        ii=ii+1;
                    end
                end
                % Choose the coefficients under the minimum discrepancy according to the selected reference coefficient
                switch err_coeType
                    case 1 % based on CF Ce
                        [~,err_min_idx]=min(err_coe_ZF1(:,1));
                    case 2 % based on CF Ca
                        [~,err_min_idx]=min(err_coe_ZF1(:,2));
                    case 3 % based on IL Ce
                        [~,err_min_idx]=min(err_coe_ZF1(:,3));
                    case 4 % based on IL Ca
                        [~,err_min_idx]=min(err_coe_ZF1(:,4));
                    case 5 % combined CF and IL
                        [~,err_min_idx]=min(sum(err_coe_ZF1(:,1:4),2)/4);
                    case 6 % combined CF
                        [~,err_min_idx]=min(sum(err_coe_ZF1(:,1:2),2)/2);
                    case 7 % combined IL
                        [~,err_min_idx]=min(sum(err_coe_ZF1(:,3:4),2)/2);
                    case 8 % based on Cd
                        [~,err_min_idx]=min(err_coe_ZF1(:,5));
                end
                % Process data based on the number of backward cases (using the embedded direction)
                coe1_n = ceil(err_min_idx/ifile); % selected index for forward coefficients
                coe0_n = mod(err_min_idx,ifile); % selected index for backward coefficients
                if coe0_n == 0
                    coe0_n = ifile;
                end
                % Extract positive values and then do nearest selection and maximum selection; reliability needs to be ensured here
                switch Coe_Zhengzhi
                    case 1 % based on CF Ce
                        Coe11=0;
                        Coe00=0;
                        for ifile1=1:length(filelist1)
                            if Coe1_repeat(ifile1,1)>Coe11
                                Coe11 = Coe1_repeat(ifile1,1);
                                coe1_n = ifile1;
                            end
                        end
                        for ifile0=1:length(filelist0)
                            if Coe0_repeat(ifile0,1)>Coe00
                                Coe00 = Coe0_repeat(ifile0,1);
                                coe0_n = ifile0;
                            end
                        end
                end
                Train_data1(i1,7:end) = Coe1_repeat(coe1_n,1:4);
                csv_coe(i-1,:) = Coe1_repeat(coe1_n,1:5); % keep forward coefficient results
                Train_data0(i0,7:end) = Coe0_repeat(coe0_n,1:4);
                csv_coe(i,:) = Coe0_repeat(coe0_n,1:5);   % keep backward coefficient results
                case_record(i-1:i) = [1;1]; % for those already repeated twice, mark as accepted
            end
        end        
    end
end
clear m
m=find(case_record==0);
save('information.mat','S','D','L','M','n0','n','m','N0');
%% Write the analyzed hydrodynamic coefficients into the csv file
switch ytrain_direction
    case 2 % use the mean of forward and backward as training data
        y_train1 = (csv_coe(1:2:end,:)+csv_coe(2:2:end,:))./2;
    case 1 % use forward data as reference
        y_train1 = csv_coe(1:2:end,:);
    case 0 % use backward data as reference
        y_train1 = csv_coe(2:2:end,:);
    case 3 % use the positive and larger ce between forward/backward data
        for iyt=1:2:length(csv_coe(:,1))
            if csv_coe(iyt,1) > csv_coe(iyt+1,1)
                y_train1((iyt+1)/2,:) = csv_coe(iyt,:); % forward data is larger, take forward
            else
                y_train1((iyt+1)/2,:) = csv_coe(iyt+1,:); % backward data is larger, take backward
            end
        end
end
switch ytrain_type
    case 1 % CF
        y_train = y_train1(:,[1,2]);
    case 2 % IL
        y_train = y_train1(:,[3,4]);
    case 3 % CF+IL
        y_train = y_train1;
    case 4 % Cd
        y_train = y_train1(:,5);
end

save('err_coe_ZF.mat','err_coe_ZF');
save('Train_data0.mat','Train_data0');
save('Train_data1.mat','Train_data1');
save('csv_coe.mat','csv_coe');
save('y_train.mat','y_train');

data = table(case_No, case_name, case_record, csv_input0(:,1),csv_input0(:,2),...
       csv_input0(:,3),csv_input0(:,4),csv_input0(:,5),csv_input0(:,6),...
       csv_input0(:,7),csv_input0(:,8),csv_input0(:,9),csv_input0(:,10),...
       case_direction,case_count,...
       csv_coe(:,1),csv_coe(:,2),csv_coe(:,3),csv_coe(:,4),csv_coe(:,5),...
       'VariableNames', {'Number', 'Name', 'Finished','A1','f1','A2','f2',...
       'A1non','f1non','A2non','f2non','Theta','Speed','Direction','Count',...
       'Cv1','Ca1','Cv2','Ca2','Cdm'});
writetable(data, filenamecsv);

return
