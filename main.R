
library(stringr)
library(glmnet) # for logistic regression and Lasso and Elastic-Net Regularized Generalized Linear Models

# function to orthogonally encode the sequence
# We add X (unknown aminoacid) as a new class 
orthogonally_encoding<-function(sequence){
  
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
                   "R", "S", "T", "V", "W", "Y");
  
  # initialize the matrix
  encoded_sequence <- matrix(0, nrow = nchar(sequence), ncol = length(amino_acids)+1)# +1 because of the "_"
  
  for (i in 1:nchar(sequence)) {
    aa <- substring(sequence, i, i)
    idx <- match(aa, amino_acids)
    if (!is.na(idx) && idx !="_") {
      encoded_sequence[i, idx] <- 1
    }
  }
  
  return(encoded_sequence)
}


create_window<-function(protein, index, window_size){
  
  # protein - the protein for which the window will be determined
  # idx - the index of the aminoacid on which the window will be centered
  # window_size - the size of the window
  # returns the encoded aminoacids chain
  extension_char<-"_"
  
  
  protein_length <- nchar(protein);
  
  # the number of elements on the left and the right side from the center
  range<-floor(window_size/2);
  
  # we have to calculate the range of aminoacids and copy them
  extended_left="";
  extended_right="";
  start<-index-range;
  if(start<1){
     extended_left=toString(rep(extension_char, abs(start)+1, collapse = NULL));
     extended_left=gsub(',','',extended_left);
     extended_left=gsub(' ','',extended_left);
  }
  
  stop<-index+range;
  if(stop>protein_length){
    extended_right =toString(rep(extension_char, stop-protein_length+1), collapse = NULL);
    extended_right=gsub(',','',extended_right);
    extended_right=gsub(' ','',extended_right);
  }
  
  if(start<1){
    index = index+abs(start)+1;
    stop=stop+abs(start)+1;
  }
  
  temp_protein <-paste0(extended_left,protein,extended_right);
  window<-substring(temp_protein,start, stop);
  return(window);
  
}


sample_data<-function(protein, second_lv_struct, window_size){
  
  samples_df<-data.frame();
  current<-NULL;
  
  for(i in 1:nchar(protein)){
    
    window<-create_window(protein,i,window_size);
    encoded<-((orthogonally_encoding(window)));
    struct_str<-substring(second_lv_struct,i,i);
   
    current<-data.frame(
      struct_str<-struct_str,
      encoded<-t(as.vector(encoded))
    ) 
    
    samples_df<-rbind(current,samples_df);
  }

  return(samples_df);
}





# function to read the data from file
read_protein_data<-function(filepath){
  
  lines <- readLines(filepath); # read all lines
  protein_data <- data.frame(); # initialize variables
  current <- NULL;
  
  for( i in 1:length(lines)){
    
    line<-lines[i]; # read line
    
    if(i%%4 == 1){
      id<-substring(line, 2); # add id without '<<'
    } else if(i%%4 == 2){
      sequence<-line # read sequence
    } else if(i%%4 == 3){
      sec_lvl_struct<-line;  #read secondary structure
    } else if(i%%4 == 0){
    
      current<-data.frame(
        id,
        sequence,
        sec_lvl_struct
      )
      protein_data = rbind(protein_data,current);
      
      #current$sampled_data<-sample_data(current$sequence,current$structure,20)
      #protein_data<-append(current, protein_data); # add protein to the list
    }
    
    
    
    
  }
  return(protein_data);
}

create_samples<-function(window_size, filepath){
  # read data into the workspace
  protein_data = read_protein_data(filepath);
  samples <- data.frame(); # initialize variables
  
  # create samples 
  for(i in 1:length(protein_data)){
    current = protein_data[i,];
    protein = current$sequence;
    second_lv_struct = current$sec_lvl_struct;
    
    sample =data.frame(sample_data(protein, second_lv_struct, window_size));
    samples=rbind(sample,samples);
  }
  return(samples);
}

create_binary_classifier<-function(test_filepath, train_filepath, window_size){
  
  test_samples = create_samples(window_size, "test.txt");
  train_samples = create_samples(window_size,"train.txt");
  
  train_samples_CnotC = create_binary_samples(train_samples, "C");
  test_samples_CnotC = create_binary_samples(train_samples, "C");
  
  train_samples_HnotH = create_binary_samples(train_samples, "H");
  test_samples_HnotH = create_binary_samples(train_samples, "H");
  
  train_samplesEnotE = create_binary_samples(train_samples, "E");
  test_samplesEnotE = create_binary_samples(train_samples, "E");
  
  
  
  
}

create_binary_samples_OVR<-function(data, className){
  
  binary_data<-data.frame();
  for(i in 1:nrow(data)){
    row = data[i,];
    general_class = row$struct_str;
    if(general_class==className){
      binary_class=1;
    }
    else{
      binary_class=-1;
    }
    
    current = data.frame(binary_class,row);
    binary_data=rbind(current, binary_data);
  }  
  
  return(binary_data);
}

create_binary_samples_OVO<-function(data, className1, className2){
  
  binary_data<-data.frame();
  for(i in 1:nrow(data)){
    row = data[i,];
    general_class = row$struct_str;
    if(general_class==className1){
      binary_class=1;
    }
    else if(general_class==className2){
      binary_class=-1;
    }
    else{
      binary_class=0;
    }
    if(binary_class!=0){
    current = data.frame(binary_class,row);
    binary_data=rbind(current, binary_data);
    }
  }  
  
  return(binary_data);
}


create_binary_model<-function(dataset){
  df=dataset[,3:ncol(dataset)];
  model<-glm(as.factor(dataset$binary_class)~., family = binomial(), df)
}





