


### The following functions up until decision_tree are internal functions ###

Node<-function(){# data that is smallar than the split value is always put to the left, and the reset to the right

  nd = list(
    left_child=0,
    right_child=0,
    split_val = NULL,
    split_var = NULL,
    left_prediction=-1,
    right_prediction = -1

  )

}


Tree<-function(){
  tr=list(
    nodeList=NULL,
    gene_id=c(),
    marker_types=c()
  )
}







Data_obj<-function(){
  obj = list(
    data=NULL,
    labels=NULL,
    belongs_to=0,
    entropy=-1,
    is_left=NULL # is this data object a left child or right child
  )
}


maxEntropy<-function(data_list){
  index<-1
  max_entropy <- 0
  for(i in 1:length(data_list)){
    if(max_entropy<data_list[[i]]$entropy){
      max_entropy<-data_list[[i]]$entropy
      index<-i
    }

  }
  return(index)
}


get.split_index<-function(data_list){
  index<-1
  entropy <- c()
  for(i in 1:length(data_list)){
    entropy<-c(entropy,data_list[[i]]$entropy)
  }
  next_split <- which(entropy!=0)
  entropy<-entropy[next_split]
  for(i in 2:length(entropy)){
    if(length(entropy)<2){
      return(next_split[1])

    }
    entropy[i]<-entropy[i-1]+entropy[i]
  }

  random_val <- runif(1, min = 0, max = max(entropy))
  retIndex <- which(entropy > random_val)[1]
  return(next_split[retIndex])


}


getEntropy<-function(labels){#generalizable to multiclass
  entropy<-0
  data_length<-length(labels)
  n.unique<-unique(labels)
  for(i in n.unique){
    current_val <- i
    ratio <- length(which(labels == current_val))/data_length
    entropy<- entropy - ratio*(log2(ratio))
  }

  return(entropy)
}

row.predict<-function(node_list, data_row){
  current_index<-1

  while(TRUE){
    gene<-as.character(node_list[[current_index]]$split_var)
    val<-as.numeric(data_row[gene])
    if(val<=node_list[[current_index]]$split_val){
      #go left
      if(node_list[[current_index]]$left_child !=0){
        current_index<-node_list[[current_index]]$left_child
      }else {
        return(node_list[[current_index]]$left_prediction)
      }
    }else{# go right
      if(node_list[[current_index]]$right_child !=0){
        current_index<-node_list[[current_index]]$right_child
      }else {
        return(node_list[[current_index]]$right_prediction)
      }
    }
  }
}

tree.predict<-function(node_list, df){
  retPrediction<-c()
  for(i in 1:NROW(df)){
    retPrediction<-c(retPrediction,row.predict(node_list = node_list, df[i,]))
  }
  return(retPrediction)
}


get.Error<-function(node_list, df,labels){# only works for binary classes
  predictions<-tree.predict(node_list, df)
  FP<-length(which(predictions==0 & labels==1))
  FN<-length(which(predictions==1 & labels == 0))
  TP <- length(which(predictions == 0 & labels == 0))
  TN <- length(which(predictions == 1 & labels == 1))
  percision <- TP/(TP + FP)
  recall <- TP/(TP+FN)

  #class_0_error <-length(which(predictions==1 & labels==0)) / length(which(labels == 0))
  #class_1_error <-length(which(predictions== 0& labels==1)) / length(which(labels == 1))

  score<-length(which(predictions != labels))/NROW(df)
  score<-1/score



  return(c(percision,recall, score))
}


get.geneid<-function(node_list){
  names<-c()
  for(i in 1:length(node_list)){
    names<-c(names,as.character(node_list[[i]]$split_var))
  }
  return(names)
}

get.prediction<-function(node_list, index){
  # as long as the returned value is -1 means this left and right
  left_prediction<-NULL
  right_prediction<-NULL
  if(node_list[[index]]$left_child == 0  ){ # no left child
    left_prediction<-node_list[[index]]$left_prediction
  }else{
    left_prediction<-get.prediction(node_list, node_list[[index]]$left_child)
  }

  if(node_list[[index]]$right_child == 0  ){ # no left child
    right_prediction<-node_list[[index]]$right_prediction
  }else{
    right_prediction<-get.prediction(node_list, node_list[[index]]$right_child)
  }

  if(left_prediction == right_prediction){
    return(left_prediction)
  }else{
    # -1 means left and right have different predicitons
    return(-1)
  }

}

has.zero<-function(node_list, index){

  if(node_list[[index]]$left_child == 0 && node_list[[index]]$left_prediction == 0){
    return(TRUE)
  }else if(node_list[[index]]$right_child == 0 && node_list[[index]]$right_prediction == 0){
    return(TRUE)
  }

  if(node_list[[index]]$left_child!=0 && has.zero(node_list, node_list[[index]]$left_child)){ # check if left child has 0
    return(TRUE)
  }else if(node_list[[index]]$right_child!=0 && has.zero(node_list, node_list[[index]]$right_child)){# check if right child has 0
    return(TRUE)
  }
  # else, no child is found in any branch. Thus return false
  return(FALSE)
}

get.type<-function(node_list, index){
  if(get.prediction(node_list, index) != (-1)){
    return("X")
  }
  # to be reviewed
  zero_on_left<-NULL
  zero_on_right<-NULL
  if(node_list[[index]]$left_child==0 && node_list[[index]]$left_prediction == 0){
    zero_on_left<-TRUE
  }else if(node_list[[index]]$left_child==0 && node_list[[index]]$left_prediction == 1){
    zero_on_left<-FALSE
  }else{
    zero_on_left<- has.zero(node_list, node_list[[index]]$left_child)
  }

  if(node_list[[index]]$right_child ==0 && node_list[[index]]$right_prediction ==0){
    zero_on_right<-TRUE
  }else if(node_list[[index]]$right_child ==0 && node_list[[index]]$right_prediction ==1){
    zero_on_right<-FALSE
  }else{
    zero_on_right<-has.zero(node_list, node_list[[index]]$right_child)
  }

  if(zero_on_left && zero_on_right){
    return("?")
  }else if(zero_on_left){
    return("-")
  }else if(zero_on_right){
    return("+")
  }
  #
  return("Something went wrong")

}

get.markers_type<-function(node_list, gene_id){
  node_index<-1
  while(node_list[[node_index]]$split_var != gene_id){
    node_index<-node_index+1
  }
  return(get.type(node_list, node_index))

}



get.markers_types<-function(node_list, genes){
  types<-c()
  for(i in 1:length(genes)){
    types<-c(types,get.markers_type(node_list,genes[i]))
  }
  return(types)
}





# creates a decision tree object given the training data and labels
# training data is a dataFrame and training_labels is a vector representing the
# classes. This implementation only accepts binary classes, and labels can only take
# the values of 0 and 1, where 0 represents cells within the cluster of interest
# and 1 represents cells outside
decision_tree<-function(training_data, training_labels, n_gene = 2){

  tree <- list()
  data_list<-list()
  data_list[[1]]<-Data_obj()
  data_list[[1]]$data<-training_data
  data_list[[1]]$labels<-training_labels
  data_index<-1 # keeps tract of the largest error
  current_index<-1 # keeps of tract of howmany data are in the list

  # stores the error change and split value each time a node is added
  error_change<-c()
  split_values<-c()

  breakpoint <-1
  j<-1
  while( j <= n_gene){
    # finds the split that results in lowest errors rate




    training_data<-data_list[[data_index]]$data
    training_labels<-data_list[[data_index]]$labels

    df.rf <- randomForest( y=as.factor(training_labels), training_data, ntree = 1,maxnode=2)
    tree_matrix<-getTree(df.rf, labelVar = TRUE)


    if(length(unique(tree_matrix$prediction))<=2){
      next
    }


    currentNode<-Node()
    currentNode$split_val<-tree_matrix[1,4]
    currentNode$split_var<-tree_matrix[1,3]
    predicted_labels<-predict(df.rf,training_data, type = "class")



    left<-Data_obj()
    left$data<-training_data[which(training_data[,currentNode$split_var]<=currentNode$split_val),]
    left$labels<-training_labels[which(training_data[,currentNode$split_var]<=currentNode$split_val)] # might not be data frame
    left$belongs_to<-j
    left$entropy<-getEntropy(left$labels)
    left$is_left<-TRUE
    left$prediction<-round(mean(as.numeric(left$labels)))# needs to be chaged if classification is multi-class

    right<-Data_obj()
    right$data<-training_data[which(training_data[,currentNode$split_var]>currentNode$split_val),]
    right$labels<-training_labels[which(training_data[,currentNode$split_var]>currentNode$split_val)]
    right$belongs_to<-j
    right$entropy<-getEntropy(right$labels)
    right$is_left<-FALSE
    right$prediction<-round(mean(as.numeric(right$labels)))# need to be changed if classificiation is multi-class

    if(is.na(right$prediction)){
      right$prediction<-left$prediction
    }else if(is.na(left$prediction)){
      left$prediction<-right$prediction
    }

    if(is.na(right$prediction)||is.na(left$prediction)||right$prediction == left$prediction){
      if(breakpoint<10){
        breakpoint<-breakpoint+1
        next
      }else{
        breakpoint<-1
      }

    }


    currentNode$left_prediction<-left$prediction
    currentNode$right_prediction<-right$prediction





    if(j>1){
      parent<-data_list[[data_index]]$belongs_to
      if(data_list[[data_index]]$is_left){
        tree[[parent]]$left_child<-j
      }else{
        tree[[parent]]$right_child<-j
      }

    }

    data_list[[data_index]]$entropy<-0


    tree[[j]]<-currentNode
    current_index<-current_index+1
    data_list[[current_index]]<-left
    current_index<-current_index+1
    data_list[[current_index]]<-right

    error_change<-c(error_change, get.Error(tree,data_list[[1]]$data,data_list[[1]]$labels))

    data_index<-get.split_index(data_list) # decides which place to split

    # stores the split value of the current node
    split_values<-c(split_values, currentNode$split_val)

    j<-j+1

  }

  gene_id<-get.geneid(tree)
  markers_tyeps<-get.markers_types( tree, gene_id)
  #markers_tyeps<-NULL
  return(list(Tree = tree,
              genes = gene_id,
              types = markers_tyeps,
              error_change = error_change,
              split_values = split_values))

}



### The above functions are all internal functions###

##################################
# END OF DECISION TREE ALGORITHM
##################################












#
# This function is used to perform a filtering on the gene list
# It is used inside select_markers3 to prevent stacking one good gene
# marker with manying bad ones marked by X
#
# @Returns the filtered dataframe
#
# @Param df is the gene list dataframe with the set of gene markers found
# by select_markers3
# @Param sga_hits the data Frame to be printed
#
# @Author Simon Hu
genelist_filter<-function(df, n.genes){
  i <-1
  while(i <NROW(df)){
    currentRow<-df[i,]
    next.pattern<-df[(i+1),1:n.genes]
    current.pattern<-df[i,1:n.genes]
    next.genes<-df[(i+1),(n.genes+1):(2*n.genes)]
    current.genes<-df[i,(n.genes+1):(2*n.genes)]

    next.non_X.gene<-next.genes[which(next.pattern!='X')]
    current.non_X.gene<-current.genes[which(current.pattern!='X')]
    if(length(next.non_X.gene) == length(current.non_X.gene)&&all(current.non_X.gene == next.non_X.gene)){
      df<-df[-(i+1),]
      next
    }
    i<-i+1
  }

  return(df)
}





#
#' Title
#' @description This function finds the gene markers for a given cluster.
#' @param df A dataframe represents the single cell data matrix, with the genes as the features and cells as the rows. NAs are not allowed
#' @param cluster_index A one column dataframe that contains the cluster labels assigned to each cell.
#' @param cluster_i The cluster which we wish to extract the markers for
#' @param n_genes The number of genes used in each decision tree
#' @param n_trees The number of trees one wishes to construct for each cluster, default is the n_genes times the number of features in the data matrix
#' @param reduced_form Indicates whether or not the output should be in its reduced form. The full form contains the recall, percision and 1/error rate at each iteration.
#' Whereas the reduced form only contains the final recall, percision, and 1/error rate
#'
#' @return The dataframe with the list of markers
#' @importFrom randomForest randomForest
#' @importFrom randomForest getTree
#'
#' @export
#'
#' @examples ### Suppose you want to extract the markers for cluster 0
#' df <- t(pbmc_small@scale.data)
#' labels <- data.frame(as.numeric(pbmc_small@meta.data$res.1))
#' #Make sure the labels are not factors as they cause problems during the mapping to binary labels. Use either character or integers
#' markers<-selectMarkersRF(df, cluster_index = labels, cluster_i = 0)


selectMarkersRF<-function(df,cluster_index,cluster_i, n_genes = 2, n_trees = NCOL(df)*2, reduced_form = TRUE){



  # convert cluster assignments to a binary class
  # for example if there are 4 clusters and we are trying to identify cluster 2 then the
  # cluster id of cluster 2 would be 0 and clusters 1, 3, and 4 would be assigned 1
  colnames(cluster_index)[1]<-"cluster_assignment"
  cluster_index$cluster_assignment[cluster_index$cluster_assignment != cluster_i] <- -1
  cluster_index$cluster_assignment[cluster_index$cluster_assignment == cluster_i] <- 0
  cluster_index$cluster_assignment[cluster_index$cluster_assignment == -1] <- 1

  cluster_index$cluster_assignment<-as.integer(cluster_index$cluster_assignment)

  ###### Compute sampling size#########
  sample_size<-cluster_index %>% group_by(cluster_assignment) %>% count()
  sample_size<-min(sample_size[,2])

  ###### CREATE A BALANCED DATASET by SUBSAMPLING ########
  cell_in_cluster_0<-which(cluster_index == 0)
  cell_in_cluster_1<-which(cluster_index == 1)
  random_numbers0<-sample(1:length(cell_in_cluster_0),sample_size,replace = TRUE)
  random_numbers1<-sample(1:length(cell_in_cluster_1),sample_size, replace = TRUE)

  cell_in_cluster_0<-cell_in_cluster_0[random_numbers0]
  cell_in_cluster_1<-cell_in_cluster_1[random_numbers1]


  ####### NOW do a train test split, into training and testing set########
  train_size <- floor(sample_size*0.80)

  cluster_0_train_row_location <- cell_in_cluster_0[1:train_size]

  cluster_1_train_row_location <- cell_in_cluster_1[1:train_size]

  train_set_row_location<-c(cluster_0_train_row_location,cluster_1_train_row_location)

  training_data <- df[train_set_row_location,]

  training_labels <- data.frame(c(rep(0,length(cluster_0_train_row_location)),rep(1,length(cluster_1_train_row_location))))
  testing_data<-df[-train_set_row_location,]
  testing_labels<-data.frame(cluster_index[-train_set_row_location,])

  #######################################

  l<-list()
  errors<-c()
  marker.types<-c()


  df.reduced<-training_data
  i<-1
  while(i <= n_trees){

    ###################################
    # Create decision tree
    tree<- decision_tree(training_data, training_labels[,1], n_gene = n_genes)
    training_error<-get.Error(tree$Tree, training_data, training_labels[,1])
    testing_error<- get.Error(tree$Tree, testing_data, testing_labels[,1])



    print("Tree created")

    #remove trees that split at the same gene twice
    current_name<-get.geneid(tree$Tree)
    if(length(unique(current_name)) < length(current_name)){
      next
    }

    # Report the lowest confidence score from the training and testing data
    # in other words use the largest error rate between the training and testing error
    # for the following computations
    if(training_error[3]<testing_error[3]){
      errors<-rbind(errors,c(i,training_error))
    }else{
      errors<-rbind(errors,c(i,testing_error))
    }

    l[[i]]<-tree
    i<-i+1
  }

  errors<-errors[order(-errors[,NCOL(errors)]),]
  l<-l[errors[,1]]

  ordered_marker.types<-c()
  gene_id <-c()
  filtered_errors<-c()
  error_change<-c()

  split_vals<-c()

  print("start extracting genes")
  i<-1
  while(i<=NROW(errors)){
    currentId<-get.geneid(l[[i]]$Tree)
    ########################################################3
    ordered_marker.types<-rbind(ordered_marker.types, get.markers_types(l[[i]]$Tree, currentId))
    gene_id<-rbind(gene_id,currentId)
    filtered_errors<-rbind(filtered_errors, errors[i,])

    error_change<-rbind(error_change,l[[i]]$error_change)

    split_vals <-rbind(split_vals,l[[i]]$split_values)

    i<-i+1

  }




  print("finished extracting genes")

  colnames(split_vals) <- rep(c("split_val"), n_genes)
  colnames(gene_id) <- rep(c("gene"), n_genes)
  colnames(ordered_marker.types) <-rep(c("marker_type"), n_genes)
  colnames(error_change)<-rep(c("percision_on_itr","recall_on_iter","invers_error_itr"),n_genes)

  errors<-errors[,-1]
  filtered_errors<-filtered_errors[,-1]
  # append markers
  retFrame<-data.frame(ordered_marker.types,gene_id, split_vals,filtered_errors)
  if(reduced_form == FALSE){
      retFrame<-data.frame(ordered_marker.types,gene_id, split_vals,error_change,filtered_errors)
  }


  print("retFrame constructed")

  retFrame<-retFrame %>% group_by(.dots = colnames(retFrame)[1:(NCOL(gene_id)*2)]) %>% summarise_all(mean)


  print("grouping worked properly")
  #retFrame<-aggregate(retFrame$errors, by = list(retFrame$gene_id),FUN = mean)

  retFrame<-retFrame[order(-retFrame[,NCOL(retFrame)]),]

  print("ordering worked properly")

  colnames(retFrame)[NCOL(retFrame)] <- "inverse_error"

  colnames(retFrame)[NCOL(retFrame)-1] <- "recall"

  colnames(retFrame)[NCOL(retFrame)-2] <-"percision"

  ################################################################## Break

  retFrame<-data.frame(retFrame)
  retFrame<-genelist_filter(retFrame, n_genes)
  return(retFrame)


}

# input is a dataframe outputted from select_best_markers3
# returns an evaluation score on the that clustering solution
# assumes the last column is the cluster size
# To be reviewed later

clustering_evaluation_index<-function(data_frame){
  score<-0
  total_sample_size<-sum(unique(data_frame$Cluster))
  cluster_index<-unique(data_frame$cluster_index)
  for(i in cluster_index){
    cluster_size<-data_frame[,NCOL(data_frame)][which(data_frame$cluster_index == i)[1]]
    size_ratio<-cluster_size/total_sample_size
    score<-score+min(mean(data_frame$conf_index[which(data_frame$cluster_index == i)][1:3]),10) * size_ratio
  }
  return(score)
}



#
# Finds the markers for all of the clusters given the cluster labels

#
#' Title
#'
#' @param seurat The Seurat Object
#' @param labels A dataframe that contains labels for the clustering solution. The dataframe should only have 1 column
#' @param output.graphs Indicates whether or not you wish to output the feature plots of the markers found to a pdf file
#' @param output.table Indicates whether or not you wish to output the table of markers to a pdf file
#' @param n_genes The number of genes used in each decision tree
#' @param topn_markers An integer indicating how many markers you wish to return for each cluster
#' @param graph.name The name of the feature plot if one wishes to output them to a pdf file
#' @param table.name The name of the list of markers if one wishes to output it into a pdf file
#' @param n.trees The number of trees one wishes to construct for each cluster, if zero, then the value will be determined automatically as the number of most variable genes indicated in the seurat object
#'
#' @return A dataframe that contains the list of markers and their associated attributes
#' @export
#'
#' @examples
#' data <- t(pbmc_small@scale.data)
#' labels <- data.frame(as.numeric(pbmc_small@meta.data$res.1))
#' marker_list<-getAllMarkers(data, labels = labels)
#'
getAllMarkers<-function(df , labels, specify_clusters = NULL ,output.graphs = FALSE,n_genes = 2, topn_markers = 10, graph.name = "Feature plot for cluster ", table.name = "Table of computed confidence values for resolution: ", n_trees=NCOL(df)*n_genes, method = "RF"){
  Clusters<-c()
  if(is.null(specify_clusters)){
    Clusters <-unique(labels)[,1]
  }else{
    Clusters <-specify_clusters
  }


  Table.merge <- c()
  for(i in Clusters){
    retVal <- c()
    if (method == "RF"){
      retVal<-selectMarkersRF(df, labels,i, n_genes = n_genes, n_trees = n_trees)
    }

    retVal<-retVal[1:topn_markers,]
    retVal<-cbind(retVal,rep(i,NROW(retVal)), rep(length(which(labels[,1] == i)),NROW(retVal)))
    colnames(retVal)[(NCOL(retVal)-1):NCOL(retVal)]<-c("cluster_index", "cluster_size")
    Table.merge<-rbind(Table.merge,retVal)

    if(output.graphs){
      plot.name<-paste(graph.name , i, ".pdf")
      pdf(plot.name)
      for(j in 1:NROW(retVal)){
        num_of_genes_used = (NCOL(retVal)-4)/2
        name_index <- c((num_of_genes_used+1):(2*num_of_genes_used))
        genenames <- as.character(unlist(retVal[j,name_index]))

        FeaturePlot(object = seurat, features.plot = genenames, cols.use = c("lightgrey", "blue"))
      }
      dev.off()

    }


  }




  return(Table.merge)
}






increase_sample_size<-function(df,labels, cluster_index, multiplication_factor){

    df.sampled<-df[which(labels==cluster_index),]
    labels.sampled<-data.frame(labels[which(labels==cluster_index),])
    colnames(labels.sampled)<-colnames(labels)

    for(i in 1:multiplication_factor){
        df<-rbind(df,df.sampled)
        labels<-rbind(labels,labels.sampled)
    }
    return(list(df = df,
                labels = labels))

}



