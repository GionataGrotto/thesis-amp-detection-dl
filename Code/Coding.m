classdef Coding
   
    methods(Static)
        
        
        
        function seq = hydrocoding(sequence)
           alphabet = 'ARNDCEQGHILKMFPSTWYV'; 
           values = [0.62,-2.53,-0.78,-0.90,0.29,-0.74,-0.85,0.48,-0.4,1.38,1.06...
               -1.50,0.64,1.19,0.12,-0.18,-0.05,0.81,0.26,1.08];
           for i = 1:length(sequence)
              if contains(alphabet,sequence(i))
                  seq(i) = values(find(alphabet == sequence(i)));
              end
           end
        end
        
        function seq = PIcoding(sequence)
           alphabet = 'ARNDCEQGHILKMFPSTWYV'; 
           values = [6.11,10.76,5.43,2.98,5.15,3.08,5.65,6.06,7.64,6.04,6.04...
               9.47,5.71,5.76,6.30,5.70,5.60,5.88,5.63,6.02];
           for i = 1:length(sequence)
              if contains(alphabet,sequence(i))
                  seq(i) = values(find(alphabet == sequence(i)));
              end
           end
        end
       
        function seq = pk1coding(sequence)
           alphabet = 'ARNDCEQGHOILKMFPSTWYV'; 
           values = [2.34,2.17,2.02,1.88,1.96,2.19,2.17,2.34,1.82,1.82...
               2.36,2.36,2.18,2.28,1.83,1.99,2.21,2.09,2.83,2.20,2.32];
           for i = 1:length(sequence)
              if contains(alphabet,sequence(i))
                  seq(i) = values(find(alphabet == sequence(i)));
              end
           end
        end
        
        function seq = pk2coding(sequence)
           alphabet = 'ARNDCEQGHOILKMFPSTWYV'; 
           values = [9.69,9.04,8.80,9.60,10.28,9.67,9.13,9.60,9.17,9.65,9.60...
               9.60,8.95,9.21,9.13,10.60,9.15,9.10,9.39,9.11,9.62];
           for i = 1:length(sequence)
              if contains(alphabet,sequence(i))
                  seq(i) = values(find(alphabet == sequence(i)));
              end
           end
        end
        
        function seq = weightcoding(sequence)
           alphabet = 'ARNDCEQGHOILKMFPSTWYV'; 
           values = [89.10,174.20,132.12,133.11,121.16,147.13,146.15,75.07,155.16...
               131.13,131.18,131.18,146.19,149.21,165.19,115.13,105.09...
               119.12,204.23,181.19,117.15];
           for i = 1:length(sequence)
              if contains(alphabet,sequence(i))
                  seq(i) = values(find(alphabet == sequence(i)));
              end
           end
        end
        
    end
end