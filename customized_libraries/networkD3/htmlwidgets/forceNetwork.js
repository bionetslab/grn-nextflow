HTMLWidgets.widget({

  name: "forceNetwork",

  type: "output",

  initialize: function(el, width, height) {

    d3.select(el).append("svg")
        .attr("width", width)
        .attr("height", height);

    return d3.forceSimulation();
  },

  resize: function(el, width, height, force) {

    d3.select(el).select("svg")
        .attr("width", width)
        .attr("height", height);

    force.force("center", d3.forceCenter(width / 2, height / 2))
        .restart();
  },

  renderValue: function(el, x, force) {

  // Compute the node radius  using the javascript math expression specified
    function nodeSize(d) {
            if(options.nodesize){
                    return eval(options.radiusCalculation);

            }else{
                    return 6}

    }


    // alias options
    var options = x.options;

    // convert links and nodes data frames to d3 friendly format
    var links = HTMLWidgets.dataframeToD3(x.links);
    var nodes = HTMLWidgets.dataframeToD3(x.nodes);

    links.sort(function(a,b) {
      if (a.source > b.source) {return 1;}
      else if (a.source < b.source) {return -1;}
      else {
          if (a.target > b.target) {return 1;}
          if (a.target < b.target) {return -1;}
          else {return 0;}
      }
    });

    for (var i=0; i<links.length; i++) {
      if (i != 0 &&
          links[i].source == links[i-1].source &&
          links[i].target == links[i-1].target) {
              links[i].linknum = links[i-1].linknum + 1;
          }
      else {links[i].linknum = 1;};
    };

    console.log(links)
  

    // create linkedByIndex to quickly search for node neighbors
    // adapted from: http://stackoverflow.com/a/8780277/4389763
    var linkedByIndex = {};
    links.forEach(function(d) {
      linkedByIndex[d.source + "," + d.target] = 1;
      linkedByIndex[d.target + "," + d.source] = 1;
    });
    function neighboring(a, b) {
      return linkedByIndex[a.index + "," + b.index];
    }

    // get the width and height
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var color = eval(options.colourScale);

    // set this up even if zoom = F
    var zoom = d3.zoom();

    // create d3 force layout
    force
      .nodes(d3.values(nodes))
      .force("link", d3.forceLink(links).distance(options.linkDistance))
      .force("center", d3.forceCenter(width / 2, height / 2))
      .force("charge", d3.forceManyBody().strength(options.charge))
      .on("tick", tick);

    force.alpha(1).restart();

      var drag = d3.drag()
        .on("start", dragstart)
        .on("drag", dragged)
        .on("end", dragended)
      function dragstart(d) {
        if (!d3.event.active) force.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
      }
      function dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
      }
      function dragended(d) {
        if (!d3.event.active) force.alphaTarget(0);
        d.fx = null;
        d.fy = null;
      }

    // select the svg element and remove existing children
    var svg = d3.select(el).select("svg");
    svg.selectAll("*").remove();
    // add two g layers; the first will be zoom target if zoom = T
    //  fine to have two g layers even if zoom = F
    svg = svg
        .append("g").attr("class","zoom-layer")
        .append("g")

    // add zooming if requested
    if (options.zoom) {
      function redraw() {
        d3.select(el).select(".zoom-layer")
          .attr("transform", d3.event.transform);
      }
      zoom.on("zoom", redraw)

      d3.select(el).select("svg")
        .attr("pointer-events", "all")
        .call(zoom);

    } else {
      zoom.on("zoom", null);
    }


    var linkColoursArr = d3.nest().key(function(d) { return d.colour; }).entries(links);


    // <marker id="arrow-#009e73" viewBox="0 -5 10 10" refX="25" refY="-2" markerWidth="6" markerHeight="6" orient="auto" style="fill: rgb(0, 158, 115);"><path d="M0,-4L8,0L0,4"></path></marker>
    svg.append("defs").selectAll("marker")
      .data(linkColoursArr)
    .enter().append("marker")
      .attr("id", function(d) { return "arrow-" + d.key; })
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 15.75)
      .attr("refY", -0.75)
      .attr("markerWidth", 6)
      .attr("markerHeight", 6)
      .attr("orient", "auto")
      .style("fill", "none")
      .style("stroke", function(d) { return d.key; })
      .style("stroke-width", "1.5")
      .style("opacity", options.opacity)
    .append("path")
      .attr("d", "M0,-3L4,0L0,3");
      

    var path = svg.selectAll(".link")
      .data(links)
    .enter().append("path")
      .attr("class", function(d) { return "link " + d.value; })
      .attr("marker-end", function(d) { return "url(#arrow-" + d.colour + ")"; })
      .style("fill", "none")
      .style("stroke", function(d) { return d.colour ; })
      .style("stroke-width", "1")
      .style("opacity", options.opacity);
  

    // // draw links
    // var link = svg.selectAll(".link")
    //   .data(links)
    //   .enter().append("line")
    //   .attr("class", "link")
    //   .style("stroke", function(d) { return d.colour ; })
    //   //.style("stroke", options.linkColour)
    //   .style("opacity", options.opacity)
    //   .style("stroke-width", 1)
    //   .on("mouseover", function(d) {
    //       d3.select(this)
    //         .style("opacity", 1);
    //   })
    //   .on("mouseout", function(d) {
    //       d3.select(this)
    //         .style("opacity", options.opacity);
    //   });

    // if (options.arrows) {
    //   link.style("marker-end",  function(d) { return "url(#arrow-" + d.colour + ")"; });

    //   var linkColoursArr = d3.nest().key(function(d) { return d.colour; }).entries(links);

    //   svg.append("defs").selectAll("marker")
    //       .data(linkColoursArr)
    //       .enter().append("marker")
    //         .attr("id", function(d) { return "arrow-" + d.key; })
    //         .attr("viewBox", "0, -5, 10, 10")
    //         .attr("refX", 0)
    //         .attr("markerWidth", 4)
    //         .attr("markerHeight", 4)
    //         .attr("orient", "auto")
    //         .style("fill", "context-fill")
    //         .style("fill", function(d) { return d.key; })
    //         .style("opacity", options.opacity)
    //       .append("path")
    //         .attr("d", "M0,-5 L10,0 L0,5");
    // }

    // draw nodes
    var node = svg.selectAll(".node")
      .data(force.nodes())
      .enter().append("g")
      .attr("class", "node")
      .style("fill", function(d) { return color(d.group); })
      .style("opacity", options.opacity)
      .on("mouseover", mouseover)
      .on("mouseout", mouseout)
      .on("click", click)
      .call(drag);

    node.append("circle")
      .attr("r", function(d){return nodeSize(d);})
      .style("stroke", "#fff")
      .style("opacity", options.opacity)
      .style("stroke-width", "1.5px");

    node.append("svg:text")
      .attr("class", "nodetext")
      .attr("dx", 12)
      .attr("dy", ".35em")
      .text(function(d) { return d.name })
      .style("font", options.fontSize + "px " + options.fontFamily)
      .style("opacity", options.opacityNoHover)
      .style("pointer-events", "none");

    function tick() {
      node.attr("transform", function(d) {
        if(options.bounded){ // adds bounding box
            d.x = Math.max(nodeSize(d), Math.min(width - nodeSize(d), d.x));
            d.y = Math.max(nodeSize(d), Math.min(height - nodeSize(d), d.y));
        }

        return "translate(" + d.x + "," + d.y + ")"});

        path.attr("d", function(d) {
          var dx = d.target.x - d.source.x,
              dy = d.target.y - d.source.y,
              dr = 75/d.linknum;  //linknum is defined above
          return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;
        });

      // function idx(d, type) {
      //   var linkWidthFunc = eval("(" + options.linkWidth + ")");
      //   var a = d.target.x - d.source.x;
			//   var b = d.target.y - d.source.y;
			//   var c = Math.sqrt(Math.pow(a, 2) + Math.pow(b, 2));

      //   var m = - (a / b)

      //   var x1_withoutOffset = d.source.x + (nodeSize(d.source) * a) / c
      //   var y1_withoutOffset = d.source.y + (nodeSize(d.source) * b) / c

      //   var t = y1_withoutOffset - m * x1_withoutOffset

      //   var x_offset = nodeSize(d.source)
      //   var x_1 = x1_withoutOffset + x_offset

      //   var y_offset = y1_withoutOffset - (m * (x_1) + t)
      //   var y_1 = y1_withoutOffset - y_offset

      //   // if (type == "x1") return (d.source.x + ((nodeSize(d.source)) * a) / c);
      //   // if (type == "y1") return (d.source.y + ((nodeSize(d.source) * b) / c));
      //   // if (options.arrows) {
      //   //   if (type == "x2") return (d.target.x - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * a) / c));
      //   //   if (type == "y2") return (d.target.y - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * b) / c));
      //   // } else {
      //   //   if (type == "x2") return (d.target.x - ((nodeSize(d.target) * a) / c));
      //   //   if (type == "y2") return (d.target.y - ((nodeSize(d.target) * b) / c));
      //   // }  

      //   if (d.value == 1) {
      //     if (type == "x1") return (d.source.x);
      //     if (type == "y1") return (d.source.y);
      //     if (options.arrows) {
      //       if (type == "x2") return (d.target.x - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * a) / c));
      //       if (type == "y2") return (d.target.y - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * b) / c));
      //     } else {
      //       if (type == "x2") return (d.target.x - ((nodeSize(d.target) * a) / c));
      //       if (type == "y2") return (d.target.y - ((nodeSize(d.target) * b) / c));
      //     }  
      //   } else {
      //     if (type == "x1") return (x_1);
  		// 	  // if (type == "y1") return (d.source.y + ((nodeSize(d.source) * b) / c));
  		// 	  if (type == "y1") return (y_1);
  		// 	  if (options.arrows) {
  		// 	    if (type == "x2") return (d.target.x - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * a) / c));
  		// 	    if (type == "y2") return (d.target.y - ((((5 * linkWidthFunc(d)) + nodeSize(d.target)) * b) / c));
  		// 	  } else {
  		// 	    if (type == "x2") return (d.target.x - ((nodeSize(d.target) * a) / c));
  		// 	    if (type == "y2") return (d.target.y - ((nodeSize(d.target) * b) / c));
  		// 	  }
      //   }
		  // }

      // link
      //   .attr("x1", function(d) { return idx(d, "x1"); })
      //   .attr("y1", function(d) { return idx(d, "y1"); })
      //   .attr("x2", function(d) { return idx(d, "x2"); })
      //   .attr("y2", function(d) { return idx(d, "y2"); });
    }

    function mouseover(d) {
      // unfocus non-connected links and nodes
      //if (options.focusOnHover) {
        var unfocusDivisor = 4;

        link.transition().duration(200)
          .style("opacity", function(l) { return d != l.source && d != l.target ? +options.opacity / unfocusDivisor : +options.opacity });

        node.transition().duration(200)
          .style("opacity", function(o) { return d.index == o.index || neighboring(d, o) ? +options.opacity : +options.opacity / unfocusDivisor; });
      //}

      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", function(d){return nodeSize(d)+5;});
      d3.select(this).select("text").transition()
        .duration(750)
        .attr("x", 13)
        .style("stroke-width", ".5px")
        .style("font", options.clickTextSize + "px ")
        .style("opacity", 1);
    }

    function mouseout() {
      node.style("opacity", +options.opacity);
      link.style("opacity", +options.opacity);

      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", function(d){return nodeSize(d);});
      d3.select(this).select("text").transition()
        .duration(1250)
        .attr("x", 0)
        .style("font", options.fontSize + "px ")
        .style("opacity", options.opacityNoHover);
    }

    function click(d) {
      return eval(options.clickAction)
    }

    // add legend option
    if(options.legend){
        var legendRectSize = 18;
        var legendSpacing = 4;
        var legend = d3.select(el).select('svg').selectAll('.legend')
          .data(color.domain())
          .enter()
          .append('g')
          .attr('class', 'legend')
          .attr('transform', function(d, i) {
            var height = legendRectSize + legendSpacing;
            var offset =  height * color.domain().length / 2;
            var horz = legendRectSize;
            var vert = i * height+4;
            return 'translate(' + horz + ',' + vert + ')';
          });

        legend.append('rect')
          .attr('width', legendRectSize)
          .attr('height', legendRectSize)
          .style('fill', color)
          .style('stroke', color);

        legend.append('text')
          .attr('x', legendRectSize + legendSpacing)
          .attr('y', legendRectSize - legendSpacing)
          .text(function(d) { return d; });
    }

    // make font-family consistent across all elements
    d3.select(el).selectAll('text').style('font-family', options.fontFamily);
  },
});