
library(stringr)
library(glmnet) # for logistic regression and Lasso and Elastic-Net Regularized Generalized Linear Models

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
    
    current<-data.frame(
      struct_str<-struct_str,
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
  
  # create samples 
  for(i in 1:length(protein_data)){
    current = protein_data[i,];
    protein = current$sequence;
    second_lv_struct = current$sec_lvl_struct;
    
    sample =data.frame(sample_protein(protein, second_lv_struct, window_size));
    samples=rbind(sample,samples);
  }
  return(samples);
}


# The function to create the samples for one vs rest binary classifier
# @input - data - sampled protein data
# @input - classname - the name of the class (for instance "H")
# @return - binary_data - the encoded data
# @warning - the target values are 1 if class is the same as the classname
#            and (-1) if it is different
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

# The function to create the samples for one vs one binary classifier
# @input - data - sampled protein data
# @input - className1 - the name of the first class (for instance "H")
# @input - className2 - the name of the second class (for instance "E")
# @return - binary_data - the encoded data
# @warning - the target values are 1 if class is the same as the className1
#            and (-1) if it is the same as the className2, other classes are
#            not added to the data
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

# The function to create the binary model (logisitic regression)
# @input - dataset - the sampled dataset
# @input - maxit - the max number of interations
# @input - epsilon - the epsilon parameter
# @output - model - the trained logisitic regression model
create_binary_model<-function(dataset, maxit=80, epsilon=1e-6){
  control <- glm.control(maxit = maxit, epsilon = epsilon)
  df=dataset[,3:ncol(dataset)];
  model<-glm(as.factor(dataset$binary_class)~., family = binomial(), df);
  return(model);
}

# Function to save the pred_class column to the file
# @input - filepath - the filepath for the result file
# @input - dataframe - the dataframe which column will be saved
save_dataframe<-function(filepath, dataframe){
  
  # Zapisywanie kolumny "Imię" do pliku tekstowego za pomocą write.table()
  write.table(dataframe$pred_class, filepath, row.names = FALSE, col.names = FALSE)
}


# The function to load data from file and save it as a one string
# @input - filepath - the filepath with the file to read
save_as_fasta<-function(filepath){
  
  lines <- readLines(filepath);
  sequence <- "";
  
  for( i in 1:length(lines)){
    line<-lines[i]; # read line
    sequence<-paste0(sequence, line);
  }
  sequence<-gsub("\"","", sequence);
  path<-gsub(".txt",".fasta",filepath);
  writeLines(sequence, paste0("corrected",path));
}

# The function to load filepaths from the file
# and save these files as one string
# @input - filepath - the file with filepaths
save_all_as_fasta<-function(filepath){
  filepaths<-readLines(filepath)
  
  for( i in 1:length(filepaths)){
    line<-filepaths[i]; # read line
    res<-save_as_fasta(line);
  }
  
}