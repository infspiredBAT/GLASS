HTMLWidgets.widget({

    name: 'chromatography',
    type: 'output',

    initialize: function(el, w, h) {
        //console.log("I'm being initialized.")
        var instanceCounter = 0;
        var intrex = "";
        var max_x = 0;
        var max_y = 0;
        var margin  = {top: 10,  right: 10, bottom: 100, left: 40},
            margin2 = {top: 430, right: 10, bottom: 20,  left: 40},
            width   = w - margin.left - margin.right,
            height  = h - margin.top  - margin.bottom,
            height2 = h - margin2.top - margin2.bottom;
        var widthScale   = d3.scale.linear().range([0,width]);
            width2Scale  = d3.scale.linear().range([0,width]),  //remains constant, to be used with context
            heightScale  = d3.scale.linear().range([height,0]),
    	      height2Scale = d3.scale.linear().range([height2,0]);
        var line = d3.svg.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return heightScale(d)});
        //lines in the brush tool
        //var noise_area = d3.svg.line()
        //        .x(function(d){return widthScale(d[0]);})
        //        .y(function(d){return heightScale(d[1]);});
        var noise_area = d3.svg.area()
                .x(function(d){return widthScale(d[0]);})
                .y0(heightScale(0)+2)
                .y1(function(d){return heightScale(d[1]);});
        var linec = d3.svg.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return height2Scale(d)});
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
        
        var label_pos = {}        //map for pisitioning labels representing called bases 
        label_pos["reference"] = 10;
        label_pos["call"]      = 28;
        label_pos["sec_fwd"]   = 38;
        label_pos["call_rev"]  = 58;
        label_pos["sec_rev"]   = 68;
        label_pos["user_pri"]  = 58;  // + rev
        label_pos["user_sec"]  = 70;  // + rev
        function redraw()  {
            widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            var w = brush.extent()[1]-brush.extent()[0] ;
            focus.selectAll("g").selectAll(".path").attr("d", line);
            focus.selectAll("g").selectAll(".area").attr("d",noise_area);
            focus.selectAll(".peak_label").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".q").attr("x",function(d){return widthScale(d["trace_peak"])-9;});
            focus.selectAll(".line").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                    .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            //conditional visibility
            if(w<260){
                focus.selectAll(".peak_label").attr("visibility","visible");
                if(w==0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                focus.selectAll(".peak_label").attr("visibility","hidden");
                if(w<800){focus.selectAll(".q").attr("visibility","visible")}
                if(w<2000){focus.selectAll(".short").attr("visibility","visible");}
            }
        }
        function brushed() { redraw(); }
        function reHeight(domain_y){
            heightScale.domain([0,domain_y]);
            focus.selectAll("g").selectAll(".path").attr("d", line);
            focus.selectAll("g").selectAll(".area").attr("d",noise_area);
        }
        //setting brush programmatically
        function setBrush(start,end){
            context.call(brush.extent([start,end]));
            redraw();
        }
        function setPeakLabel(calls,label,offset){
            focus.append("g").selectAll("text.seq").data(calls).enter() //reference
  	            .append("text")
                .attr("class",function(d){
                    if(label.indexOf("sec") > -1){
                        return "peak_label short sec";
                    }else if(label.indexOf("user")> -1){
                        return "peak_label short user";
                    }else{return "peak_label short";}    
                })
                .text(function(d){
                    if(label.indexOf("sec") > -1){
                        return d[label].toLowerCase();
                    }else if(label.indexOf("user") > -1 && d[label]==="low qual"){
                        return "N";}
                    else {
                        return d[label];}  
                 })
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos[label]+offset)
                .attr("fill", "black")
                .attr("opacity", function(d){
                    if (d[label]==="low qual") {return 0.3;}
                    else                       {return 0.8;}})
                .attr("font-family", "sans-serif")
                .attr("font-size",function(){if(label.indexOf("user")>-1){return "12px";}else{return "11px";}})
                .attr("stroke",function(d) {
                    if      (d[label] === "A"){ return "#33CC33"; }
                    else if (d[label] === "C"){ return "#0000FF"; }
                    else if (d[label] === "G"){ return "#000000"; }
                    else if (d[label] === "T"){ return "#FF0000"; }
                    else if (d[label] === "-"){ return "white"; }
                    else    {                   return "orange"; }});  
        }
        function showVarInMap(choices){
            //genomic
            //console.log("changing choices");
            context.selectAll("lines.choices").data(choices).enter()
    			      .append("line")
                .attr("class","varInMinimap context")
      			    .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y1",1)
      			    .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y2",24)
      			    .attr("stroke-width",3)
      			    .attr("stroke",function(d) {
      			        if      (d["reference"] === "A"){ return "#33CC33"; }
      			        else if (d["reference"] === "C"){ return "#0000FF"; }
      			        else if (d["reference"] === "G"){ return "#000000"; }
      			        else if (d["reference"] === "T"){ return "#FF0000"; }
      			        else if (d["reference"] === "-"){ return "white"; }
      			        else    {                         return "yellow"; }});
  			    // user
  			    context.selectAll("lines.choices").data(choices).enter()
  				      .append("line")
                .attr("class","varInMinimap context")
  				      .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y1",26)
  				      .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y2",49)
  				      .attr("stroke-width",3)
  				      .attr("stroke",function(d) {
  				          if      (d["user_pri"] === "A"){ return "#33CC33"; }
  				          else if (d["user_pri"] === "C"){ return "#0000FF"; }
  				          else if (d["user_pri"] === "G"){ return "#000000"; }
  				          else if (d["user_pri"] === "T"){ return "#FF0000"; }
  				          else if (d["user_pri"] === "-"){ return "white"; }
  				          else    {                        return "yellow";  }});
/*
  			context.selectAll("text.choices.coord").data(choices).enter()
  				.append("text").attr("class","varInMinimap")
  				.attr("x",function(d){return width2Scale(d["trace_peak"])+4;})
  				.attr("y",30)
  				.attr("opacity",0.6)
  				.text(function(d){return d["id"];})
  				.attr("fill","black");
*/
            focus.append("g").selectAll("variance_indicator").data(choices).enter()  //variance indicator
  		         .append("line").attr("class","peak_label short line varind")
  		         .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y1",180)
  		         .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y2",280)
  		         .attr("stroke-width",8).attr("stroke","rgba(255,0,255,0.15)").attr("stroke-dasharray",2);
        }

        Shiny.addCustomMessageHandler("zoom",
            function(message) {
                setBrush(Number(message)-100,Number(message)+120);
            }
        );
        Shiny.addCustomMessageHandler("opac_f",
            function(message){
                //console.log(message);
                focus.selectAll(".line_f").attr("opacity",Number(message));
//                focus.selectAll(".fwd").attr("opacity",Number(message));
            }
        );
        Shiny.addCustomMessageHandler("opac_r",
            function(message){
                //console.log(message);
                focus.selectAll(".line_r").attr("opacity",Number(message));
//                focus.selectAll(".rev").attr("opacity",Number(message));
            }
        );
        function callShiny(message){
            //console.log(message);
            Shiny.onInputChange("posClick", {id: message});
        };
        function brushZoomIn(){
            return function(event) {
                event.preventDefault();
                var ext = brush.extent();
                full = ext[1] - ext[0];
                if(full >= 20){
                    ten = full/10;
                    setBrush(ext[0]+ten,ext[1]-ten);}
            };
        };
        function brushZoomOut(){
            return function(event) {
                event.preventDefault();
                var ext = brush.extent();
                full = ext[1] - ext[0];
                ten = full/10;
                from = (ext[0]-ten);
                if(from < 0){from=0};
                to   = (ext[1]+ten);
                //if(to > width){to = width}; must get the scales right to set boundaries
              
                setBrush(from,to);}
        };
        function brushMoveLeft(){
            return function(event) {
                event.preventDefault();
                var ext = brush.extent();
                full = ext[1] - ext[0];
                ten = full/10;
                from = (ext[0]-ten);
                if(from < 0){from=0};
                to   = (ext[1]-(ext[0]-from));
                setBrush(from,to);}
          
        };
        function brushMoveRight(){
            return function(event) {
                event.preventDefault();
                var ext = brush.extent();
                full = ext[1] - ext[0];
                ten = full/10;
                to = ext[1] + ten;
                //if(to > width){to=width};must set scales
                from  = (ext[0]+(to-ext[1]));
                setBrush(from,to);}
        };
        d3.select('body').call(d3.keybinding()
            .on('←', brushMoveLeft())
            .on('↑', brushZoomIn())
            .on('→', brushMoveRight())
            .on('↓', brushZoomOut())
            );
        transpose = function(a) {
        
          // Calculate the width and height of the Array
          var w = a.length ? a.length : 0,
            h = a[0] instanceof Array ? a[0].length : 0;
        
          // In case it is a zero matrix, no transpose routine needed.
          if(h === 0 || w === 0) { return []; }
        
          /**
           * @var {Number} i Counter
           * @var {Number} j Counter
           * @var {Array} t Transposed data is stored in this array.
           */
          var i, j, t = [];
        
          // Loop through every item in the outer array (height)
          for(i=0; i<h; i++) {
        
            // Insert a new row (array)
            t[i] = [];
        
            // Loop through every item per item in outer array (width)
            for(j=0; j<w; j++) {
        
              // Save transposed data.
              t[i][j] = a[j][i];
            }
          }
        
          return t;
        };
        //passing arguments
        //this enables to access vars and functions from the render function as instance.*
        return {
            svg:     svg,
            line:    line,
            noise_area: noise_area,
            transpose: transpose,
            linec:   linec,
            context: context,
	          brush:   brush,
            focus:   focus,
            widthScale:  widthScale,
            width2Scale: width2Scale,
            heightScale: heightScale,
            width:   w,
            height:  h,
	          height2: height2,
            instanceCounter: instanceCounter,
            reHeight: reHeight,
            setBrush: setBrush,
            setPeakLabel: setPeakLabel,
            showVarInMap: showVarInMap,
            callShiny: callShiny
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
    //function called everytime input parameters change
    renderValue: function(el, x, instance) {

        instance.lastValue = x;
        //the render function behaves differently when called repeatedly
        //the first run is actually still a part initialization step
        if(x.new_sample){
  		    //console.log(x)
  		    var intens = x["intens"];
            var intens_rev = "";
            var rev = 0;  //offset on labels in case we have alternative reference
            if(x["intens_rev"] !== null){
                var intens_rev = x["intens_rev"];
                rev = 30;
            }
        if(instance.instanceCounter>=1){ //cleanup after previous sample
            instance.focus.selectAll(".area").remove();
            instance.focus.selectAll(".line_f").remove();
            instance.focus.selectAll(".line_r").remove();
            instance.focus.selectAll(".peak_label").remove();
            instance.context.selectAll(".context").remove();
        }
        instance.instanceCounter = instance.instanceCounter+1;
  			var intens_guide_line = x["intens_guide_line"];
  			var calls       = HTMLWidgets.dataframeToD3(x["calls"]);
  			var choices     = HTMLWidgets.dataframeToD3(x["choices"]);
  			var domain_y    = x["helperdat"]["max_y"];
  			instance.max_y  = domain_y;
  			var domain_x    = x["helperdat"]["max_x"];
  			instance.max_x  = domain_x;
  			var intrex      = HTMLWidgets.dataframeToD3(x["helperdat"]["helper_intrex"])
  			instance.intrex = intrex;

  			var svg     = instance.svg;
  			var line    = instance.line;
        var noise_area = noise_area;
  			var linec   = instance.linec;
  			var focus   = instance.focus;
  			var context = instance.context;
  			var brush   = instance.brush;
  			var widthScale  = instance.widthScale;
  			var heightScale = instance.heightScale;
        var width2Scale = instance.width2Scale;
  			var height2 = instance.height2;

  			widthScale.domain([0,domain_x]);
  			width2Scale.domain([0,domain_x]);
  			heightScale.domain([0,domain_y]);
  			height2Scale.domain([0,domain_y]);

  			//visualise fullseq width
  			context
  				.append("line").attr("class","context")
  				.attr("x1",0)
  				.attr("y1",25)
  				.attr("x2",widthScale(domain_x))
  				.attr("y2",25)
  				.attr("stroke-width",3).attr("stroke","rgba(180,180,180,1.0)");
  			//visualise introns/exons
/*
            //lines
  			context.selectAll("lines.intrex").data(intrex).enter()
  				.append("line")
  				.attr("x1",function(d){return widthScale(d["start"]);})
  				.attr("y1",0)
  				.attr("x2",function(d){return widthScale(d["start"]);})
  				.attr("y2",50)
  				.attr("stroke-width",2).attr("stroke","rgba(20,20,20,0.6)").attr("stroke-dasharray",2);
//  				.on("mouseover", function(){d3.select(this).style("fill", "white");})
//  				.on("mouseout",  function(){d3.select(this).style("fill", "rgba(200,200,200,0.2)");});
*/
            //intron/exon boxes
  			context.selectAll("rect").data(intrex).enter()
  				.append("rect").attr("class","context")
  				.attr("x",function(d){return widthScale(d["start"]);})
  				.attr("y",0).attr("rx",3).attr("ry",3)
  				.attr("width",function(d){return widthScale(d["end"]-d["start"]);})
  				.attr("height",50)
  				.attr("fill",function(d) {
  				    if (/exon/.test(d["attr"])){ return "rgba(190,190,190,1.0)";
  				    } else {                     return "rgba(230,230,230,1.0)"; }
  				});
  			context.selectAll("text.intrex.name").data(intrex).enter()
  				.append("text").attr("class","context")
  				.attr("x",function(d){return widthScale(d["start"]);})
  				.attr("y",-2)
  				.attr("opacity",0.6)
  				.text(function(d){return d["attr"];})
  				.attr("fill","black");
  			context.selectAll("text.intrex.start").data(intrex).enter()
  				.append("text").attr("class","context")
  				.attr("x",function(d){return widthScale(d["start"]);})
  				.attr("y",60)
  				.attr("opacity",0.6)
  				.text(function(d){return d["id"];}) 
  				.attr("fill","black");

  			brush.x(width2Scale);

  			var group_a = focus.append("g");  
  			var group_c = focus.append("g");
  			var group_g = focus.append("g");
  			var group_t = focus.append("g");

        //forward strand
  			group_a.selectAll("path").data([intens["A"]]).enter()
  			    .append("path").attr("class","path line_f")
  				.attr("d",line)
  				.attr("fill","none")
  				.attr("stroke","#33CC33").attr("stroke-width",0.75);
  			group_c.selectAll("path").data([intens["C"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line)
  				.attr("fill","none")
  				.attr("stroke","#0000FF").attr("stroke-width",0.75);
  			group_g.selectAll("path").data([intens["G"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line)
  				.attr("fill","none")
  				.attr("stroke","#000000").attr("stroke-width",0.75);
  			group_t.selectAll("path").data([intens["T"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line)
  				.attr("fill","none")
  				.attr("stroke","#FF0000").attr("stroke-width",0.75);
        
        

        var a = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
        //var a = [x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]
        
        var group_noise_f = focus.append("g");
        group_noise_f.selectAll("path").data([a]).enter()
            .append("path").attr("class","area").attr("d",noise_area)
            .attr("fill","#000000").attr("stroke","none").attr("opacity",0.3);
        //reverse strand
        if(intens_rev != ""){
            var group_a_r = focus.append("g");
    	    	var group_c_r = focus.append("g");
  		    	var group_g_r = focus.append("g");
  		    	var group_t_r = focus.append("g");
            group_a_r.selectAll("path").data([intens_rev["A"]]).enter()
        			.append("path").attr("class","path line_r")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#33CC33").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,3,10,3");
      			group_c_r.selectAll("path").data([intens_rev["C"]]).enter()
      				.append("path").attr("class","path line_r")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#0000FF").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,3,10,3");
      			group_g_r.selectAll("path").data([intens_rev["G"]]).enter()
      				.append("path").attr("class","path line_r")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#000000").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,3,10,3");
      			group_t_r.selectAll("path").data([intens_rev["T"]]).enter()
      			    .append("path").attr("class","path line_r")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#FF0000").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,3,10,3");
            }

            //trace peak labels
            focus.append("g").selectAll("qualities").data(calls).enter()  //quality box
  		        .append("rect").attr("class","peak_label q")
  		        .attr("x",function(d){return (widthScale(d["trace_peak"])-9);})
  		        .attr("y",0).attr("rx",1).attr("ry",1)
  		        .attr("width",18)
  		        .attr("height",function(d){return d["quality"];})
  		        .attr("fill", "rgba(200,200,200,0.3)");
            focus.append("g").selectAll("text.qualities").data(calls).enter() //quality number
  		        .append("text").attr("class","peak_label")
  		        .text(function(d){return d["quality"];})
  		        .attr("text-anchor", "middle")
  		        .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
  		        .attr("y",-2)
  		        .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "10px");
              instance.setPeakLabel(calls,"reference",0); 
              instance.setPeakLabel(calls,"call",0);
              instance.setPeakLabel(calls,"sec_fwd",0);
              if(rev!==0){
                  instance.setPeakLabel(calls,"call_rev",0);
                  instance.setPeakLabel(calls,"sec_rev",0);
              }
              instance.setPeakLabel(calls,"user_pri",rev);
              instance.setPeakLabel(calls,"user_sec",rev);
            focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa
                .append("text").attr("class","peak_label short")
                .text(function(d){
                    if   (d["ord_in_cod"] == 1) {return d["AA_ref"].toUpperCase();}
                    else {                       return ".";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",88+rev)
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "11px")
                .attr("stroke","#000000");
            focus.append("g").selectAll("text.seq.codon").data(calls).enter() //codon stuff
                .append("text").attr("class","peak_label")
                .text(function(d){
                    if   (d["coding_seq"] > 0){return d["coding_seq"]+" : "+d["codon"]+"."+d["ord_in_cod"];}
                    else {                     return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(110+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
                .append("text").attr("class","peak_label")
                .text(function(d){return d["gen_coord"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(120+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.append("g").selectAll("text.exon_intron").data(calls).enter() //intrex
                .append("text").attr("class","peak_label")
                .text(function(d){return d["exon_intron"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(130+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.selectAll(".peak_label").attr("visibility","hidden")
/*
            //horizontal line on top of the chrom
                focus
      				.append("line")
      				.attr("x1",0)
      				.attr("y1",intens_guide_line)
      				.attr("x2",1200)
      				.attr("y2",intens_guide_line)
      				.attr("stroke-width",1).attr("stroke","rgba(0,0,0,0.6)").attr("stroke-dasharray",2);
*/
  			context.append("g")
//				.attr("class", "x brush")
      			.call(brush).attr("class","context")
  				.selectAll("rect")
  				.attr("y", -14)
  				.attr("height", 78) //height2 + 10)
  				.attr("rx",3)
  				.attr("ry",3)
  				.attr("fill","rgba(255,255,255,0.3)")
  				.attr("stroke-width",2).attr("stroke","red").attr("stroke-dasharray","2,6")
  				.attr("opacity",0.5);

            //In SVG, z-index is defined by the order the element appears in the document
            //http://stackoverflow.com/questions/17786618/how-to-use-z-index-in-svg-elements
            //the vars on the minimap are shown last sto that they stay on to
            instance.showVarInMap(choices);
	          //zooming in so that the first view is not ugly dense graph

            if (typeof choices[0] !== 'undefined') {
                from = choices[0]["trace_peak"]-100;
                to   = choices[0]["trace_peak"]+120;
                if(from < 0) {from = 0;to = 220}
                instance.setBrush(from,to);
            }else{
                instance.setBrush(200,1000);
            }

        }else{
            //console.log(x)
  			if(x["helperdat"]["max_y"]!= instance.max_y){
  				instance.reHeight(x["helperdat"]["max_y"]);
  			}else if(x.choices != instance.choices){
                var choices = HTMLWidgets.dataframeToD3(x["choices"]);
                var calls   = HTMLWidgets.dataframeToD3(x["calls"]);
                var rev = 0;  //offset on labels in case we have alternative reference
                if(x["intens_rev"] !== null){
                    var intens_rev = x["intens_rev"];
                    rev = 30;
                }
                instance.focus.selectAll(".peak_label short line").remove();
                instance.context.selectAll(".varInMinimap").remove();
                instance.showVarInMap(choices);
                instance.choices = x.choices;
                instance.focus.selectAll(".user").remove();
                instance.focus.selectAll(".sec").remove();
                instance.setPeakLabel(calls,"sec_fwd",0);
                if(rev!==0)instance.setPeakLabel(calls,"sec_rev",0);
                instance.setPeakLabel(calls,"user_pri",rev);
                instance.setPeakLabel(calls,"user_sec",rev);
  
                instance.focus.selectAll(".varind").remove();
                instance.focus.append("g").selectAll("variance_indicator").data(choices).enter() //variance indicator
                    .append("line").attr("class","peak_label short line varind")
                    .attr("x1",function(d){return instance.widthScale(d["trace_peak"]);})
                    .attr("y1",180)
                    .attr("x2",function(d){return instance.widthScale(d["trace_peak"]);})
                    .attr("y2",280)
                    .attr("stroke-width",8).attr("stroke","rgba(255,0,255,0.15)").attr("stroke-dasharray",2);

  			} else { console.log(x) }
        }
    }
});
