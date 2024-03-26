[{
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
                       },{
      "selector": "edge[tool = 'boostdiff'][interaction = 'BRCA_M0'][key = 'NSCLC_BRCA'][effect < 0]",
      "css": {
        "line-color": "#38A6A5",
        "target-arrow-color": "#38A6A5",
        "target-arrow-shape": "tee",
        "line-style": "solid"
        }
    },{
      "selector": "edge[tool = 'boostdiff'][interaction = 'BRCA_M0'][key = 'NSCLC_BRCA'][effect >= 0]",
      "css": {
        "line-color": "#38A6A5",
        "target-arrow-color": "#38A6A5",
        "target-arrow-shape": "triangle",
        "line-style": "solid"

      }
    },{
      "selector": "edge[tool = 'boostdiff'][interaction = 'NSCLC_M0'][key = 'NSCLC_BRCA'][effect < 0]",
      "css": {
        "line-color": "#5F4690",
        "target-arrow-color": "#5F4690",
        "target-arrow-shape": "tee",
        "line-style": "solid"
        }
    },{
      "selector": "edge[tool = 'boostdiff'][interaction = 'NSCLC_M0'][key = 'NSCLC_BRCA'][effect >= 0]",
      "css": {
        "line-color": "#5F4690",
        "target-arrow-color": "#5F4690",
        "target-arrow-shape": "triangle",
        "line-style": "solid"

      }
    },{
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
                           "background-color": "mapData(value, -2, 0, green, white)"
                         }
                       },
                       {
                         "selector": "node[value > 0]",
                         "css": {
                           "background-color": "mapData(value, 0, 5, white, red)"
                         }
                       }
                     ]
