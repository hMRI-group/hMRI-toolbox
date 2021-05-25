classdef hmri_test_utils
    
    % Class with general utility methods for unit testing the hMRI toolbox.
    
    methods(Static)
        
        function w_TEs=createDecaySignal(w_TE0,TEs,R2s)
            
            dims=size(w_TE0);
            
            % Account for 1D case
            if (length(dims)==2)&&(dims(2)==1), dims=dims(1); end
            
            TEs=reshape(TEs,[ones(1,length(dims)),length(TEs)]);
            
            w_TEs=w_TE0.*exp(-R2s.*TEs);
            
        end
        
    end
end