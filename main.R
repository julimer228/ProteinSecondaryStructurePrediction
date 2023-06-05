
library(stringr)



# function to orthogonally encode the sequence
# We add X (unknown aminoacid) as a new class 
orthogonally_encoding<-function(sequence){
  
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
                   "R", "S", "T", "V", "W", "Y","X");
  
  # initialize the matrix
  encoded_sequence <- matrix(0, nrow = nchar(sequence), ncol = length(amino_acids))
  
  for (i in 1:nchar(sequence)) {
    aa <- substring(sequence, i, i)
    idx <- match(aa, amino_acids)
    if (!is.na(idx)) {
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
  
  protein_length <- length(protein);
  
  # the number of elements on the left and the right side from the center
  range<-floor(window_size/2);
  
  # we need to create the auxiliary protein
  # we need to complement the protein on the right and left with the first and last amino    # acid 
  first<-substring(protein, 1,1);
  last<-substring(protein, protein_length, protein_length);
  temp_protein <-paste0(first,protein,last);
  
  # new index in the temp_protein
  index = index + range;
  
  # we have to calculate the range of aminoacids and copy them
  start<-index-range;
  stop<-index+range;
  
  window<-substring(temp_protein,start, stop);
  return(window);
  
}


sample_data<-function(protein, second_lv_struct, window_size){
  
  sample_list<-list();
  sample<-NULL; 
  
  for(i in 1:nchar(protein)){
    
    vector<-create_window(protein,i,window_size);
    encoded_vector<-orthogonally_encoding(vector);
    struct_str<-substring(second_lv_struct,i,i);
    
    sample$vector<-encoded_vector;
    sample$struct<-struct_str;
    
    sample_list<-append(sample,sample_list);
  }
  
  return(sample_list);
}





# function to read the data from file
read_protein_data<-function(filepath){
  lines <- readLines(filepath); # read all lines
  protein_data <- list(); # initialize variables
  current <- NULL;
  
  for( i in 1:length(lines)){
    
    line<-lines[i]; # read line
    
    if(i%%4 == 1){
      current$id<-substring(line, 2); # add id without '<<'
    } else if(i%%4 == 2){
      current$sequence<-line # read sequence
    } else if(i%%4 == 3){
      current$structure<-line;  #read secondary structure    
    } else if(i%%4 == 0){
      current$sampled_data<-sample_data(current$sequence,current$structure,20)
      protein_data<-append(current, protein_data); # add protein to the list
    }
    
  }
  return(protein_data);
}

protein_data = read_protein_data("test.txt");


