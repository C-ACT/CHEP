%% DNC - The 'sym4' wavelet resembles the QRS complex, which makes it a good choice 
% for QRS detection. 
% To illustrate this more clearly, extract a QRS complex and plot the result 
% with a dilated and translated 'sym4' wavelet for comparison.

    % qrsEx = ecgsig(4560:4810);
    % [mpdict,~,~,longs] = wmpdictionary(numel(qrsEx),'lstcpt',{{'sym4',3}});
    % figure
    % plot(qrsEx)
    % hold on
    % plot(2*circshift(mpdict(:,11),[-2 0]),'r')
    % axis tight
    % legend('QRS Complex','Sym4 Wavelet')
    % title('Comparison of Sym4 Wavelet and QRS Complex')