% this script fits the equation y=a*x+b where x is time, y is the ln of
% (H/H+L)

%open a tab delimited text file containing:
% header row, 16 columns:
% protein systematic name, protein common name, description [ln
% H/(L+H), st error, weight (# of unique peptides)] x 5
% Note that missing values are imported as NaNs
inputfile = '../../Python_analysis/AIR06-1/AIR06-1_xlink_exchangeM/evidence.txt_matlab2.txt';
[sys_name, ave1, ster1, weight1, ave2, ster2, weight2, ave3, ster3, weight3, ave4, ster4, weight4, ave5, ster5, weight5] = ...
    textread(inputfile, ['%s %f %f %f %f %f %f %f %f %f ' ...
'%f %f %f %f %f %f'], 'delimiter', '\t', 'headerlines',1, 'emptyvalue', NaN);
%open a text file to write results in
fileID1 = fopen('../../Python_analysis/AIR06-1/AIR06-1_xlink_exchangeM/data_fitted2.txt','w');
fileID2 = fopen('../../Python_analysis/AIR06-1/AIR06-1_xlink_exchangeM/data_notfitted2.txt','w');
% since the firs our measurement is unreliable becuase label recycling and
% GFP-nup pull I start the timecourse at 0 (t=1 hour), so I don't have to
% worry about what happened in the first hour
time = [1 2 3 4 5]';
fprintf(fileID1,'Protein\tp1\t95 percent CI\tp2\t95 percent CI\tR square\tTotal unique peptides\n');


for i=1:numel(ave1)
    % make column vectors with average, error and weight values for each
    % protein
    x = [ave1(i) ave2(i) ave3(i) ave4(i) ave5(i)]';
    %w = [weight1(i) weight2(i) weight3(i) weight4(i) weight5(i)]';
    e = [ster1(i) ster2(i) ster3(1) ster4(i) ster5(i)]';
    % the following line finds indices of the vector that has
    % non-NaN values and applies those indices to all the vectors
    % used in the calculation
    j = find(~isnan(x));
    x = x(j);
    %w = w(j);
    e = e(j);
    total_weight = sum(w);
    timemod = time(j);
    if numel(x)>3
        % only fit the curve if there are at least 4 time points present
        % fit a single expenential model, taking into account weights (the
        % quality of the data)
        [FO, G] = fit(timemod, x, 'poly1', 'Lower', [-1,-1], 'Upper', [0,1]);
        b=coeffvalues(FO);
        %[FO, G] = fit(timemod, x, 'exp1', 'Lower', [0,-1], 'Upper', [1,1]);
        %b=coeffvalues(FO);
        % get degrees of freedom adjusted r square value and 95% confidence
        % interval bounds for the fitted exponent coeffiient
        G.adjrsquare;
        a = confint(FO);
        % print the protein name and fitted parameters and goodness of fit
        % value to a text file and plot both raw data with error bars and
        % fitted curves
        fprintf(fileID1,'%s\t%4f\t%4f\t%4f\t%4f\t%4f\t%4f\n', sys_name{i}, b(1), (b(1)-a(1,1)), b(2), (b(2)-a(1,2)), G.adjrsquare, total_weight);
        %plot(FO, 'k-', timemod, x, 'k.')
        %hold on
        %errorbar(timemod, x, e,'kx')
        %options=fitoptions(FO);
    else
        fprintf(fileID2,'%s\n',sys_name{i});
        continue
    end
end
%hold off
fclose(fileID1);
fclose(fileID2);