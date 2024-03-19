
create_graph <- function(all_networks, de.genes = NULL) {
  nodes <-
    data.table(
      id = c(all_networks$target, all_networks$source),
      tool = c(all_networks$tool, all_networks$tool),
      keys = c(all_networks$key, all_networks$key)
    )
  
  # Necessary because key is reserved keywort in df.
  colnames(nodes)<-c('id', 'tool', 'key')
  nodes <- unique(nodes)
  de.genes$gene<-toupper(de.genes$gene)
  if (!is.null(de.genes)) {
    nodes <-
      merge(nodes,
            de.genes,
            by.x = "id",
            by.y = "gene",
            all.x = TRUE
      )
  }
  edges <- all_networks
  colnames(edges) <-
    c(
      "target",
      "source",
      "weight",
      "interaction",
      "effect",
      "tool",
      "key"
    )
  return(list(nodes, edges))
}

read_metacell_files<-function(selections){
  metacell.data<-list()
  
  for (k in unique(selections$key)){
    all_files <-
      list.files(
        path = file.path(opt$results.path, k),
        full.names = TRUE,
        recursive = F
      )
    files <-
      c(grep("out(.)*.tsv", all_files, value = TRUE, perl = T))
    data<- lapply(files, function(x) fread(x))
    condition_names<-as.character(sapply(files, function(x) str_replace(str_replace(basename(x), 'out_', ''), '.tsv', '')))
    
    for (f in files){
      print(f)
      data<-fread(f)
      condition_name<-as.character( str_replace(str_replace(basename(f), 'out_', ''), '.tsv', ''))
      dd<-data[, 2:ncol(data)]
      rownames(dd)<-data$Gene
      meta<-data.table(k = k, condition = rep( condition_name,ncol(dd)))
      print(group.var)
      meta[[group.var]]<-selections[condition==condition_name]$variables
      colnames(meta)<-c('key', 'condition', group.var)
      rownames(meta)<-colnames(dd)
      seu.obj<-CreateSeuratObject(counts = dd, project = "Cancer", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 1,
                                  names.delim = "_", meta.data = meta)
      metacell.data[[paste0(k, condition_name)]]<-seu.obj
    }
  }
  metacell.data<-merge(metacell.data[[1]], metacell.data[[2:length(metacell.data)]])
  metacell.data<-JoinLayers(metacell.data)
  colnames(metacell.data@meta.data)<-c("orig.ident" , "nCount_RNA", "nFeature_RNA" , "key","condition", group.var)
  return (metacell.data)
}

read_all_networks<-function(selections, res.path){
  an <- list()
  for (k in unique(selections$key)) {
    all_files <-
      list.files(
        path = file.path(res.path, k),
        full.names = TRUE,
        recursive = TRUE
      )
    network_files <-
      c(grep("aggregated_filtered_network", all_files, value = TRUE))
    all_networks <- read_files(network_files)
    an[[k]] <- all_networks
  }
  all_networks <- rbindlist(an)
  colnames(all_networks) <-
    c(
      "target",
      "source",
      "weight",
      "interaction",
      "effect",
      "tool",
      "key"
    )
  return(all_networks)
}

read_files <- function(file_list) {
  
  all_data <- list()
  
  for (file in file_list) {
    tryCatch(
      {
        data <- fread(file)
        if (!is.null(data) && nrow(data) > 0) {
          data$tool <- basename(dirname(file))
          data$key <- basename(dirname(dirname(file)))
          all_data[[basename(dirname(file))]] <- data
        } else {
          print("Data empty")
        }
      },
      error = function(e) {
        cat("Error reading", file, ": ", conditionMessage(e), "\n")
      }
    )
  }
  all_data <- rbindlist(all_data)
  return(all_data)
}

make_color_map<-function(all_networks, selections){
  my.vars<- merge(unique(all_networks[, .(interaction, key)]), selections, by = c('key'), allow.cartesian=TRUE)
  
  n.colors<-length(unique(c(my.vars$interaction, my.vars$condition)))
  color.pal<-rcartocolor::carto_pal(n=max(3, 2*n.colors), 'Prism')
  color.pal<-color.pal[seq(1, 2*n.colors, 2)]
  
  
  color.df<-data.table(color.pal[1:n.colors], unique(c(my.vars$interaction, my.vars$condition)))
  colnames(color.df)<-c('color', 'interaction')
  color.df<-color.df %>% pivot_longer(-c('interaction'))
  combos<-unique(all_networks[, .(interaction, tool)])
  colnames(combos)<-c('interaction', 'tool')
  combos<-merge(combos, color.df, by = 'interaction', allow.cartesian=T, all.y=T)
  combos<-combos %>% group_by(interaction) %>% mutate(row_count = row_number()) %>% as.data.table()
  combos$value<-sapply(1:nrow(combos), function(i) lighten(combos$value[i], (combos$row_count[i]-1)*(1/max(combos$row_count)))) 
  combos<-as.data.table(combos)
  combos<-combos[, .(tool, interaction, value)]
  
  cl<-list()
  possible_styles<-c('solid', 'dotted', 'dashed')
  dashes<- ceiling(sqrt(length(unique(all_networks$key))-1))
  possible_styles<-c('solid', 'dotted', rep('dashed', dashes*dashes))
  
  for (i in 1:length(unique(all_networks$key))){
    c<-combos
    c$arrow_style<-possible_styles[i]
    c$key<-unique(all_networks$key)[i]
    cl[[i]]<-c
  }
  combos<-rbindlist(cl)
  print(combos)
  return(combos)
}


