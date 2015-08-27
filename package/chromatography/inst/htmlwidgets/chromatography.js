HTMLWidgets.widget({

    name: 'chromatography',

    type: 'output',

    initialize: function(el, w, h) {
	  
      
        var margin  = {top: 10, right: 10, bottom: 100, left: 40},
            margin2 = {top: 430, right: 10, bottom: 20, left: 40},
            width   = w - margin.left - margin.right,
            height  = h - margin.top - margin.bottom,
            height2 = h - margin2.top - margin2.bottom;
    
        var widthScale   = d3.scale.linear()
    	    			             .range([0,width]);
            width2Scale  = d3.scale.linear()
    				                 .range([0,width]),	
            heightScale  = d3.scale.linear()
        				             .range([height,0]),
    	      height2Scale = d3.scale.linear()
    				                 .range([height2,0]);   
                             
        var line = d3.svg.line()
                		 .x(function(d,i){return widthScale(i)})
    		    		     .y(function(d){return heightScale(d)});
      
        var svg = d3.select(el).append("svg")
                               .attr("width", width + margin.left + margin.right)
                               .attr("height", height + margin.top + margin.bottom);

        var focus = svg.append("g")
                       .attr("class", "focus")
            		       .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
        
        svg.append("defs").append("clipPath")
                  		    .attr("id", "clip")
        			            .append("rect")
        			            .attr("width", width)
        			            .attr("height", height);	
                          
	      var context = svg.append("g")
				    	           .attr("class", "context")
					               .attr("transform", "translate(" + margin2.left + "," + margin2.top + ")");

    	  var brush = d3.svg.brush().on("brushend", brushed);
	
	      function brushed() {
		        widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
		        //console.log(brush.extent())
		        focus.selectAll("g").selectAll("path").attr("d", line);
		    }

        return {
            svg: svg,
            line: line,
		        context: context,
		        brush: brush,
            focus: focus,
            widthScale: widthScale,
            heightScale: heightScale,
            width: w,
            height: h,
		        height2: height2
        }
    },

    resize: function(el, width, height, instance) {
        if (instance.lastValue) {
            this.renderValue(el, instance.lastValue, instance);
        }
/*    
     d3.select(el).selectAll("svg")
          .attr("width", width)
          .attr("height", height);

     instance.size([width, height]).resume();
*/
    },

    renderValue: function(el, x, instance) {
      
        instance.lastValue = x;
        console.log(x);
        
	      var data = x["data"];
        var domain_y = x["helperdat"]["max_y"];
        var domain_x = x["helperdat"]["max_x"];
        var intrex = HTMLWidgets.dataframeToD3(x["helperdat"]["helper_intrex"])
        
        console.log(data)
        var svg = instance.svg;
        var line = instance.line;
        var focus = instance.focus;
	      var context = instance.context;
	      var brush = instance.brush;
        var widthScale = instance.widthScale;
        var heightScale = instance.heightScale;
	      var height2 = instance.height2;
           
        
        
	      widthScale.domain([0,domain_x]);
	      width2Scale.domain([0,domain_x]);
	      heightScale.domain([0,domain_y]);
	      height2Scale.domain([0,domain_y]);
        
        //visualise introns/exons     
        //TO DO
        //R must generage a readable structure for the d3 data function 
        //map color to name so that introns have different color maybe add label
        context.selectAll("rect").data(intrex).enter()
                .append("rect")
                .attr("x",function(d){return widthScale(d["start"]);})
                .attr("y",0).attr("rx",5).attr("ry",5).attr("opacity",0.5)
                .attr("width",function(d){return widthScale(d["end"]-d["start"]);})
                .attr("height",55)
                .attr("fill","rgba(0, 255, 0, 0.6)");
                   
    
						
	      brush.x(width2Scale);
	      var group_a = focus.append("g");
	      var group_c = focus.append("g");
	      var group_g = focus.append("g");	
	      var group_t = focus.append("g");
	
	      var group_ac = context.append("g");
	      var group_cc = context.append("g");
	      var group_gc = context.append("g");	
	      var group_tc = context.append("g");
	      
  
	      var linec = d3.svg.line()
		                  .x(function(d,i){return widthScale(i)})
		                  .y(function(d){return height2Scale(d)});
			
	      group_a.selectAll("path")
			         .data([data["A"]])
			         .enter()
			         .append("path")
               .attr("class","path")
			         .attr("d",line)
			         .attr("fill","none")
			         .attr("stroke","#33CC33")
			         .attr("stroke-width",0.75);				
	      group_c.selectAll("path")
			         .data([data["C"]])
			         .enter()
			         .append("path")
			         .attr("d",line)
			         .attr("fill","none")
			         .attr("stroke","#0000FF")
			         .attr("stroke-width",0.75);			
 	      group_g.selectAll("path")
			         .data([data["G"]])
			         .enter()
			         .append("path")
			         .attr("d",line)
			         .attr("fill","none")
			         .attr("stroke","#000000")
			         .attr("stroke-width",0.75);				
	      group_t.selectAll("path")
			         .data([data["T"]])
			         .enter()
			         .append("path")
			         .attr("d",line)
			         .attr("fill","none")
			         .attr("stroke","#FF0000")
			         .attr("stroke-width",0.75);
			
	      group_ac.selectAll("path")
			          .data([data["A"]])
			          .enter()
			          .append("path")
			          .attr("d",linec)
			          .attr("fill","none")
			          .attr("stroke","#33CC33")
			          .attr("stroke-width",0.5);
	      group_cc.selectAll("path")
			          .data([data["C"]])
			          .enter()
			          .append("path")
			          .attr("d",linec)
			          .attr("fill","none")
			          .attr("stroke","#0000FF")
			          .attr("stroke-width",0.5);			
	      group_gc.selectAll("path")
			          .data([data["G"]])
			          .enter()
			          .append("path")
			          .attr("d",linec)
			          .attr("fill","none")
			          .attr("stroke","#000000")
			          .attr("stroke-width",0.5);				
	      group_tc.selectAll("path")
			          .data([data["T"]])
			          .enter()
			          .append("path")
		          	.attr("d",linec)
			          .attr("fill","none")
			          .attr("stroke","#FF0000")
			          .attr("stroke-width",0.5);
		
	      context.append("g")
	             .attr("class", "x brush")
	             .call(brush)
	             .selectAll("rect")
	             .attr("y", -6)
	             .attr("height", height2 + 7);

    }

});
