# Author: Julia Merta 
# Year: 2023


library(stringr)
library(glmnet) # for logistic regression and Lasso and Elastic-Net Regularized Generalized Linear Models
library(parallel) 
library(doParallel)
# function to orthogonally encode the sequence
# We add X (unknown aminoacid) as a new class 
# @input - sequence - original sequence
# @output - encoded_sequence - orthogonally encoded sequence
orthogonally_encoding<-function(sequence){
  
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
                   "R", "S", "T", "V", "W", "Y");
  
  # initialize the matrix
  encoded_sequence <- matrix(0, nrow = nchar(sequence), ncol = length(amino_acids)+1)# +1 because of the "_"
  
  # for each aminoacid in the sequence
  for (i in 1:nchar(sequence)) {
    aa <- substring(sequence, i, i)
    idx <- match(aa, amino_acids)
    if (!is.na(idx) && idx !="_") {
      encoded_sequence[i, idx] <- 1
    }
  }
  
  return(encoded_sequence)
}

# Function to create the window for orthogonal encoding
# @input - protein - input protein structure
# @input - index - the index of the aminoacid in the middle 
# @input - window_size - the size of the window
# @output - window - the window of aminoacids (with the supplementary characters if needed)
create_window<-function(protein, index, window_size){
  
  # character for the supplementary characters in the ends of the sequence
  extension_char<-"_"
  
  
  protein_length <- nchar(protein);
  
  # the number of elements on the left and the right side from the center
  range<-floor(window_size/2);
  
  # we have to calculate the range of aminoacids and copy them
  extended_left="";
  extended_right="";
  
  # calculate if we need to add additional characters on the left side of the window
  start<-index-range;
  if(start<1){
    extended_left=toString(rep(extension_char, abs(start)+1, collapse = NULL));
    extended_left=gsub(',','',extended_left);
    extended_left=gsub(' ','',extended_left);
  }
  
  # calculate if we need to add additional characters on the right side of the window
  stop<-index+range;
  if(stop>protein_length){
    extended_right =toString(rep(extension_char, stop-protein_length+1), collapse = NULL);
    extended_right=gsub(',','',extended_right);
    extended_right=gsub(' ','',extended_right);
  }
  
  # if we added characters we need to 
  # calculate the index of the central element
  if(start<1){
    index = index+abs(start)+1;
    stop=stop+abs(start)+1;
  }
  
  # create new protein with additional characters if needed
  temp_protein <-paste0(extended_left,protein,extended_right);
  
  # create the window
  window<-substring(temp_protein,start, stop);
  return(window);
  
}

# The function to create the samples for the protein.
# for each secondary level structure element the input vector 
# of the size (window_size * 21) is created
# @input - protein - the protein sequence
# @input - second_lv_struct - the secondary level structure
# @input - window_size - the size of the window
# @output - samples_df - the dataframe with the sampled data
sample_protein<-function(protein, second_lv_struct, window_size){
  
  samples_df<-data.frame();
  current<-NULL;
  
  for(i in 1:nchar(protein)){
    
    window<-create_window(protein,i,window_size);
    encoded<-((orthogonally_encoding(window)));
    struct_str<-substring(second_lv_struct,i,i);
    HnotH<-ifelse(struct_str=="H", 1,-1);
    EnotE<-ifelse(struct_str=="E", 1,-1);
    CnotC<-ifelse(struct_str=="C", 1,-1);
    HE<-ifelse(struct_str=="H", 1, ifelse(struct_str=="E",-1, 0));
    EC<-ifelse(struct_str=="E", 1,ifelse(struct_str=="C",-1, 0));
    CH<-ifelse(struct_str=="C", 1,ifelse(struct_str=="H",-1, 0));
    
    current<-data.frame(
      struct_str<-struct_str,
      HnotH<-HnotH,
      EnotE<-EnotE,
      CnotC<-CnotC,
      HE<-HE,
      EC<-EC,
      CH<-CH,
      encoded<-t(as.vector(encoded))
    ) 
    
    samples_df<-rbind(current,samples_df);
  }
  
  return(samples_df);
}

# function to read the data from file
# @input - filepath - the filepath of the file with the dataset
# @warning - data should be in the following format
# > 1ais - the index of the protein
# NLAFALSELDRITAQ - the sequence of the aminoacids
# CHHHHHHHHHHHEEE - the secondary level structure
# (empty line)
read_protein_data<-function(filepath){
  
  lines <- readLines(filepath); # read all lines
  protein_data <- data.frame(); # initialize variables
  current <- NULL;
  
  for( i in 1:length(lines)){
    
    line<-lines[i]; # read line
    
    if(i%%4 == 1){
      id<-substring(line, 2); # add id without '<<'
    } else if(i%%4 == 2){
      sequence<-line          # read the sequence
    } else if(i%%4 == 3){
      sec_lvl_struct<-line;   #read the secondary level structure
    } else if(i%%4 == 0){
      current<-data.frame(
        id,
        sequence,
        sec_lvl_struct
      )
      protein_data = rbind(protein_data,current);
    }
  }
  return(protein_data);
}


# The function to read the file and create the samples
# @input - window_size - the size of the window
# @input - filepath - the file-path to the file with data
# @output - samples - sampled data
create_samples<-function(window_size, filepath){
  # read data into the workspace
  protein_data = read_protein_data(filepath);
  samples <- data.frame(); # initialize variables
  
  
  cl <- makeCluster(12);  # Tworzy klaster z 12 rdzeniami procesora
  registerDoParallel(cl);
  clusterExport(cl=cl, c('sample_protein', 'create_window','orthogonally_encoding'));
  
  samples <-foreach(i = 1:nrow(protein_data), .combine = rbind) %dopar% {
    current = protein_data[i, ]
    protein = current$sequence;
    second_lv_struct = current$sec_lvl_struct;
    sample =data.frame(sample_protein(protein, second_lv_struct, window_size));
    sample;
  }
  
  stopCluster(cl);
  return(samples);
}

# The function to create the binary model (logisitic regression)
# @input - dataset - the sampled dataset
# @input - maxit - the max number of interations
# @input - epsilon - the epsilon parameter
# @output - model - the trained logisitic regression model
create_binary_model<-function(dataset,target_col, maxit=80, epsilon=1e-6){
  control <- glm.control(maxit = maxit, epsilon = epsilon)
  df=dataset[,8:ncol(dataset)];
  model<-glm(as.factor(target_col)~., family = binomial(), data = df);
  return(model);
}


# The function to save a column from the dataframe to a FASTA format
# @input - data - the dataframe to save 
# @input - column - column id
# @input - id - the sequence id 
# @input - filepath - filepath to the result file
save_as_fasta<-function(column, id, filepath){
  writeLines(paste0(">", id), filepath);
  aminos<-as.character(column);
  aminos_str<-paste(aminos, collapse="");
  write(aminos_str, filepath, append = TRUE);
}


