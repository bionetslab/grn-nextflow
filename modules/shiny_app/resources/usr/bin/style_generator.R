# Function to generate the style JSON based on tools and colors
generate_style_json <- function(tools, colors, conditions, keys, arrow_styles, dash_styles) {
  edge_styles <- c()
  
  for (i in seq_along(tools)) {
    tool <- tools[i]
    color <- colors[i]
    condition <-conditions[i]
    key <- keys[i]
    arrow_style<-arrow_styles[i]
    
    if (length(unique(keys))>3){
      print(unique(keys))
      print('many arrow stke')
      arrow_style<-'solid'
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect < 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "tee",
        "line-style": "', arrow_style, '"
        }
    }')
      edge_styles <- c(edge_styles, edge_style)
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect >= 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "triangle",
        "line-style": "', arrow_style, '"

      }
    }')

      edge_styles <- c(edge_styles, edge_style)
    }
    else{
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect < 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "tee",
        "line-style": "', arrow_style, '"
        }
    }')
      edge_styles <- c(edge_styles, edge_style)
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect >= 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "triangle",
        "line-style": "', arrow_style, '"

      }
    }')
      e<-'{
        "selector": "edge[key]",
        "css": {
          "label": "data(key)",
          "width": 3
        }
      }'
      edge_styles <- c(edge_styles, edge_style)
    }
  }
  if (length(unique(keys))>3){
    
  e<-'{
        "selector": "edge[key]",
        "css": {
          "label": "data(key)",
          "width": 3
        }
      }'
  print(e)
  edge_styles<-c(edge_styles, e)
  
  }
  
  style_json <- paste0('[', 
                       '{
                         "selector": "node",
                         "css": {
                           "border-width": "2px",
                           "width": "data(size)",
                           "height": "data(size)",
                           "content": "data(id)"
                         }
                       },
                       {
                         "selector": "edge",
                         "css": {
                           "curve-style": "bezier",
                           "control-point-step-size": 40,
                           "target-arrow-shape": "triangle"
                         }
                       },',
                       paste(edge_styles, collapse = ","),
                       ',{
                         "selector": "node:selected",
                         "css": {
                           "overlay-opacity": 0.3,
                           "overlay-color": "gray",
                           "background-color": "blue"
                         }
                       },
                       {
                         "selector": "node[value <= 0]",
                         "css": {
                           "background-color": "mapData(value, -5, 0, green, white)"
                         }
                       },
                       {
                         "selector": "node[value > 0]",
                         "css": {
                           "background-color": "mapData(value, 0, 10, white, red)"
                         }
                       }
                     ]')
  
  return(style_json)
}

