classdef Extraction
   
    methods(Static)
        
        % translate amino acid sequence to numerical vector. The lenght of
        % the vector must be specified. If the sequence is shorter than the
        % set length then the vector will be filled with zeros
        function result_vector = convert_vector(sequence,padd)
            % I used the matlab function for the conversion. It will use
            % the number from 1 to 20 with the dataset that I used
            result_vector = aa2int(sequence);
            % padd vector
            while length(result_vector) < padd
               result_vector(end + 1) = 0; 
            end            
        end
        
        % function for the training of the basic dataset without 'augmentation'
        function [padd, l] = sequencePadding(dat1,dat2)
            index = 1;
            for i = 1:2:size(dat1,1)*2
               padd(i,1) = {Extraction.convert_vector(dat1(index).Sequence,200)};
               padd(i+1,1) = {Extraction.convert_vector(dat2(index).Sequence,200)};
               l(i) = 1;
               l(i+1) = 2;
            end
        end
        
        % translate all the sequence in numeric vector and put them into a
        % cell array
        function trad = aavec2int(dataset)
           for i = 1:length(dataset)
              trad(1,i) = {aa2int(dataset(i).Sequence)}; 
           end
        end
        
        
        
        % compute a frequency matrix from all the sequence of the input
        % dataset
        function freq = freq_vector(dataset)
           % find the longest sequence in the dataset analyzed
           maxLengthSequence = 0;
           for i = 1:size(dataset,1)
               if length(dataset{i,1}(1,:)) > maxLengthSequence
                   maxLengthSequence = length(dataset{i,1}(1,:));
               end
           end
           % create a matrix that fits all the element
           freq = zeros(20,maxLengthSequence);
           % for each sequence in the dataset
           for i = 1:size(dataset,1)
               % save the specific sequence
               sequence = dataset{i,1}(1,:);
               for j = 1:length(sequence)
                   % save the numeric representation of the amino acid
                   intAminoAcid = sequence(j);
                   % check for eventual errors
                   if intAminoAcid > 0
                       freq(intAminoAcid,j) = freq(intAminoAcid,j) + 1;
                   end
               end
           end
           % divide for the number of amino acid in the i-th position
           % (number of sequence)
           for i = 1:size(freq,2)
               freq(:,i) = freq(:,i)./sum(freq(:,i));
           end
           % I used methods that need the translated matrix
           freq = freq';
        end
        
        
        
        % Rappresentation functions
        
        
        % calculate the AAC matrix of a sequence
        % use the frequency matrix and the contact matrix which contain
        % value of the contact between amino acids (AAindex-contact 
        % Entry:TANS760101)
        function AAC = aac(freq,contact)
            % create a matrix of the right size
            AAC = zeros(size(freq,1),20);
            % for each row of the freq matrix
            for i = 1:size(freq,1)
               for j = 1:20
                  tmp = 0;
                  for k = 1:20
                     % frequence of the k amino acid at the position i
                     % multiplied to the contact value between k and j
                     % amino acid
                     tmp = tmp + freq(i,k) * contact(k,j);
                  end
                  AAC(i,j) = tmp;
               end
            end
        end
        
        % calcutate the PSSM of a vector of amino acid throught its
        % frequence matrix
        % use the dayhoff matrix calculated with the standard dayhoff
        % function
        function PSSM = pssm(frequence,d)
           % create a vector that fits the elements
           PSSM = zeros(size(frequence,1),20);
           %
           for i = 1:size(frequence,1)
              for j = 1:20
                  elem = 0;
                  for k = 1:j
                      % frequency of amino acid k at the position i
                      % multiplied to dayhoff matrix value between element
                      % amino acid k and j
                      elem = elem + frequence(i,k) * d(k,j);
                  end
                  % save
                  PSSM(i,j) = elem;
              end
           end
        end
        
        % use the frequency matrix and a matrix of physicochemical property
        % to return a matrix where each column in encoded with one of the
        % 57 property
        function PP = physicochem(freq,property)
           % matrix to return
           PP = zeros(size(freq,1),57);
           for i = 1:size(freq,1)
              for j = 1:57
                  tmp = 0;
                  for k = 1:20
                      tmp =  tmp + freq(i,k) * property(k,j);
                  end
                  PP(i,j) = tmp;
              end
           end
        end
        
        % multiply each vector with every column of the representation
        % matrix and return a matrix 
        function final = multiply(vector,rapp_matrix)
            final = zeros(length(vector),size(rapp_matrix,2));
            for i = 1:size(rapp_matrix,2)
                for j = 1:length(vector)
                    final(j,i) = rapp_matrix(j,i) * double(vector(j));
                end
            end
        end

        
        % Extraction function
        
        %Average block algorith
        
        % divide the matrix in n blocks and return a cell array with the
        % block in each cell
        % example if n = 3 it will return a cell array with 3 cell that
        % contains the 3 matrix blocks
        function sep_mat = divide_matrix(matrix_to_divide,n)
           % calculate the rows for block. Cast it to uint8
           rows_per_block = uint8(size(matrix_to_divide,1) / n);
           % prevent the creation of the last matrix to be empty. So
           % control that with the current number of rows per block don't
           % fill in n-1 matrix. If it does then subtract 1 for the
           % variable
           if (rows_per_block * (n - 1)) >= size(matrix_to_divide,1)
               rows_per_block = rows_per_block - 1;
           end
           % variable where it start to catch the matrix row
           count = 1;
           % index of return vector
           i = 1;
           % while the number of rows that I need for the division is
           % present
           for ind = 1:n-1
              % get the rows that i need from the matrix
              sep_mat(1,ind) = {matrix_to_divide(count:count+rows_per_block-1,:)};
              
              % update next block first row
              count = count + rows_per_block;
              i = ind;
           end
           % add the remaining rows of the vector that are less than
           % rows_per_block
           sep_mat(1,i+1) = {matrix_to_divide(count:end,:)};
        end
        
        % input is a cell array with blocks of matrix in each cell
        % calculate the means of each block column and then concatenate the
        % result of each block
        function mean_ret = get_means_blocks(blocks)
           % for the number of blocks in the array
           for i = 1:size(blocks,2)
              % save the specific block
              mat = blocks{1,i};
              % get the number of rows of the block
              l = size(mat,1);
              % for each column of the vector compute the mean
              for j = 1:size(mat,2)
                 v(j) = sum(mat(:,j)) / l;
              end
              % save it in a cell array
              mean(1,i) = {v};
           end
           % concatenate each vector
           mean_ret = [];
           for i = 1:length(mean)
               mean_ret = cat(2,mean_ret,mean{1,i}(1,:));
           end
           
        end
        
        % use the functions divide_matrix and get_means_blocks to return
        % with a single function the array needed
        % in the project I choose to concatenate the element from 2 to 6
        % blocks.
        function AV = AvBlock(matrix)
           AV = [];
           for i = 2:6
               tmp = Extraction.divide_matrix(matrix,i);
               tmp2 = Extraction.get_means_blocks(tmp);
               AV = cat(2,AV,tmp2);
           end
        end
        
        % Other way to calculate AvBlock. It was described in a paper but
        % it doesn't work with all the numbers because there are case where
        % the block doesn't have all the same number of rows
        function av = averageBlock(matrix,N)
            % number of rows for block
            M = uint8(size(matrix,1)/N);
            % for each block
            for i = 1:N
               % for each columns
               for j = 1:size(matrix,2)
                   % create a temporary variable
                   tmp = 0;
                   % sum all the element of the j column of the N block 
                   for z = 1:M
                       tmp = tmp + matrix(z+(i-1)*M,j);
                   end
                   % divide in for the number of rows and save it
                   av(j+size(matrix,2)*(i-1)) = tmp / double(M);
               end
            end
        end
        
        % Pseudo algorithm
        
        % create a normalizated matrix for the pseudo algorithm
        function pse = pseudo_matrix(matrix)
            % get the number of rows
            l = size(matrix,1);
            for i = 1:l
               % compute meand and standard deviation of each row
               meann(i) = mean(matrix(i,:));
               stdn(i) = std(matrix(i,:));
            end
            for i = 1:l
               for j = 1:size(matrix,2)
                   % compute each element
                   pse(i,j) = (matrix(i,j) - meann(i)) / stdn(i);
               end
            end
        end
        
        % return the vector of features of the pseudo algorithm. It calls
        % the pseudo_matrix function
        function ps = pseudo_alg(matrix)
           pse = Extraction.pseudo_matrix(matrix);
           % save number of rows and column of the matrix
           L = size(pse,1);
           S = size(pse,2);
           % the first S element 
           for k = 1:S
              % save the mean of each column
              ps(k) = mean(pse(:,k));
           end
           % distance between residue
           for j = 1:S
              for lag = 1:15
                 tmp = 0;
                 % compute pseudo algorith formula
                 for i = 1:L-lag
                     tmp = tmp + floor(pse(i,j)-pse(i+lag,j))^2;             
                 end
                 % divide the result to L - lag. But in case the sequence
                 % is shorter than 16 we have the case where L - lag is
                 % equal to 0. So in that case the function will divide it
                 % by 1
                 if L - lag == 0
                     divide = 1;
                 else
                     divide = L - lag;
                 end
                 ps(S+j+S*(lag-1)) = tmp / double(divide);
              end
           end
           ps(isnan(ps)) = 0;
        end
        
        
        
        % Extraction methods
        
        % Discrete wavelet algorithm
        
        % compute the discrete wavelet transform to a matrix, return a
        % vector of features
        function res = discrete_wavelet(matrix)
            % create empty array
            res = [];
            % for each column of the matrix
            for i = 1:size(matrix,2)
               % save the column
               tmp = matrix(:,1);
               % compute for iteration
               for j = 1:4
                   % calculate the Meyer wavelet decomposing the result in
                   % low frequency and high frequency. Save the low freq as
                   % tmp because we are going to reuse it
                   [tmp,h] = dwt(tmp,'dmey');
                   % save max, min, mean and standard deviation from low
                   % freq first, high freq then
                   res(end+1) = max(tmp);
                   res(end+1) = min(tmp);
                   res(end+1) = mean(tmp);
                   res(end+1) = std(tmp);
                   res(end+1) = max(h);
                   res(end+1) = min(h);
                   res(end+1) = mean(h);
                   res(end+1) = std(h);
                   % calculate the discrete cosine transform of the low 
                   % freq save the first 5 element of the result
                   dct_ = dct(tmp);
                   res = cat(2,res,dct_(1:5)');
               end
            end
        end
        
        % PseAAC with GM(1,1) algorithm
        
        % this function take as input a numerical vector and return a
        % monotone increasing vector
        function ag = AGO(seq)
           % vector of the same dimension of the input
           ag = zeros(1,length(seq));
           % for every single element of the vector
           for i = 1:length(seq)
               % the i-th element is equals to the sum of the first i
               % element of the input vector
               el = 0;
               for j = 1:i
                  el = el + seq(j); 
               end
               ag(i) = el;
            end
        end
        
        % compute the frequence of amino acid in the sequence
        % return a sequence
        function freq = aminoacidfreq(sequence)
            % size equals to the number of amino acids
            freq = zeros(1,20);
            for i = 1:length(sequence)
                % count the nu
                if sequence(i) ~= 0
                    freq(sequence(1,i)) = freq(sequence(1,i)) + 1;
                end
            end
            freq = freq ./ length(sequence);
        end
        
        % called in find_a_vector function. Use the increasing vector
        % calculated in AGO
        function Bm = b_matrix(sequence)
           % return a matrix with 2 column
           Bm = zeros(length(sequence)-1,2);
           % all the elements in the second column are equal to 1
           Bm(:,2) = 1;
           
           for i = 1:length(sequence)-1
               Bm(i,1) = -0.5 * (sequence(i) + sequence(i+1));
           end
        end
        
        % return the 2 coefficent of the grey model using the function
        % described in the paper
        function a = find_a_vector(sequence)
           % compute the increasing vector
           greymode = Extraction.AGO(sequence);
           % create the B matrix
           B = Extraction.b_matrix(greymode);
           % remove the first element from the increasing vector
           Y = sequence(2:end)';
           % find the vector
           a = inv(B' * B) * B';
           a = a * Y;
        end
        
        % compute the GM vector with the encoding of 5 physicochemical
        % property
        function seq = pseaac(sequence,p)
            % translate it to numerical vector
            s = int2aa(sequence);
            % compute the frequency
            seq = Extraction.aminoacidfreq(sequence);
            % encode in the specific psysicochemical property and add to
            % the final vector


            % Removed because redundand. Already present in the p matrix
            %{
            x = Coding.hydrocoding(s);
            a = Extraction.find_a_vector(x);
            seq = cat(2,seq,a');

            x = Coding.weightcoding(s);
            a = Extraction.find_a_vector(x);
            seq = cat(2,seq,a');
            %}

            x = Coding.PIcoding(s);
            a = Extraction.find_a_vector(x);
            seq = cat(2,seq,a');
            x = Coding.pk1coding(s);
            a = Extraction.find_a_vector(x);
            seq = cat(2,seq,a');
            x = Coding.pk2coding(s);
            a = Extraction.find_a_vector(x);
            seq = cat(2,seq,a');
            
            % encode the sequence with all the protscale property
            for i = 1:57
               x = Extraction.code_to_property(sequence,p,i);
               a = Extraction.find_a_vector(x);
               seq = cat(2,seq,a');
            end
        end
        
        % encode the sequence with physicochemical property
        function trad = code_to_property(vector,phy,n)
            trad = zeros(1,length(vector));
            for i = 1:length(vector)
               if vector(i) == 0
                   trad(i) = 0;
               else
                   trad(i) = phy(vector(i),n);
               end
            end
        end
        
        % q-statistic function applied to each pair of classifier. It
        % doesn't count the same classifier pair
        % use as input all the predictions of the classifiers and the true
        % lables of the test set
        function stat = qstatistic(savedPred, truelabels)
           % cast to categorical the labels for the comparison
           truelabels = categorical(truelabels);
           % matrix to save if the prediction is true or false
           mat = zeros(size(savedPred,1),size(savedPred,2));
           % for each predictors
           for i = 1:10
              % for each element predicted 
              for j = 1:size(savedPred,1)
                 tmp = 0;
                 % set to 1 if the prediction is true
                 if savedPred(j,i) == truelabels(j)
                    tmp = 1;
                 end
                 mat(j,i) = tmp;
              end
           end
           % create a matrix to save the q_statistic value for each pair
           stat = ones(10,10);
           % for each classifier
           for i = 1:10
              % triangle matrix
              for j = 1:i-1
                 % save the array predictions
                 c1 = mat(:,i);
                 c2 = mat(:,j);
                 % number of element predicted correctly
                 N11 = sum((c1==1)&(c2==1));
                 % number of element predicted wrongly
                 N00 = sum((c1==0)&(c2==0));
                 % number of elem where c1 = 0 and c2 = 1
                 N01 = sum((c1==0)&(c2==1));
                 % number of elem where c1 = 1 and c2 = 0
                 N10 = sum((c1==1)&(c2==0));
                 
                 % Other formula (correlation)
                 %tmp = double((N11*N00-N01*N10)/sqrt((N11+N10)*(N01+N00)*(N11+N01)*(N10+N00)));

                 % Q statistic formula (association)
                 tmp = (N11*N00-N01*N10)/(N11*N00+N01*N10);
                 % save the triangle matrix
                 stat(i,j) = tmp;
                 stat(j,i) = tmp;
              end
           end
        end
    end
end