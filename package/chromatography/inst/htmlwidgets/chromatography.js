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
            half_height = height/1.6;
            height2 = h - margin2.top - margin2.bottom;
        var widthScale   = d3.scale.linear().range([0,width]),
            width2Scale  = d3.scale.linear().range([0,width]),  //remains constant, to be used with context
            heightScale  = d3.scale.linear().range([height,0]),
    	      height2Scale = d3.scale.linear().range([height2,0]),
            heightScale_fwd_split = d3.scale.linear().range([half_height,(2*half_height -  height)]),
            heightScale_rev_split = d3.scale.linear().range([height,half_height]);
            //heightScale_fwd = heightScale;
            //heightScale_rev = heightScale;
            heightScale_fwd = heightScale_fwd_split;
            heightScale_rev = heightScale_rev_split;
        var line_fwd = d3.svg.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return heightScale_fwd(d)});
        var line_rev = d3.svg.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return heightScale_rev(d)});
        //lines in the brush tool
        //var noise_area = d3.svg.line()
        //        .x(function(d){return widthScale(d[0]);})
        //        .y(function(d){return heightScale(d[1]);});
        var mult = 2;

        var noise_area_fwd = d3.svg.area()
                .x(function(d){return widthScale(d[0]);})
                .y0(function(d){return (heightScale_fwd(0)+2);})
                .y1(function(d){return heightScale_fwd(d[1]*mult);});
        var noise_area_rev = d3.svg.area()
                .x(function(d){return widthScale(d[0]);})
                .y0(function(d){return heightScale_rev(0)+2;})
                .y1(function(d){return heightScale_rev(d[1]*mult);});
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
         var join = false;

        var label_pos = {}        //map for pisitioning labels representing called base

        label_pos["reference"]      =  10;
        label_pos["aa"]             =  25;
        label_pos["aa_sample"]      =  40;
        label_pos["user_sample"]    =  55;
        label_pos["user_mut"]       =  70;
        label_pos["aa_mut"]         =  83;
        label_pos["codon"]          = 105;
        label_pos["gen_coord"]      = 115;
//        label_pos["intrex"]         = 120;
        label_pos["call"]           = 130;
        label_pos["mut_call_fwd"]   = 145;
        label_pos["call_rev"]       = 170;
        label_pos["mut_call_rev"]   = 185;
        function brushed() { redraw(); }
        //setting brush programmatically
        function setBrush(start,end){
            context.call(brush.extent([start,end]));
            redraw();
        }
        function reHeight(domain_y){
            heightScale.domain([0,domain_y-1]);
            heightScale_fwd_split.domain([0,domain_y-1]);
            heightScale_rev_split.domain([0,domain_y-1]);
            redraw();
        }
        function redraw()  {
            widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            var w = brush.extent()[1]-brush.extent()[0] ;
            focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
            focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);

            focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd);
            focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev);
            focus.selectAll(".scope").attr("x",function(d){return widthScale(d["trace_peak"])-12;});
            focus.selectAll(".peak_label").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".qual_fwd").attr("x",function(d){return widthScale(d["trace_peak"])-9;});
            focus.selectAll(".qual_rev").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".line").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                    .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            //conditional visibility
            if(w<260){
                focus.selectAll(".peak_label").attr("visibility","visible");
                if(w==0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                focus.selectAll(".peak_label").attr("visibility","hidden");
                if(w<800){
                  focus.selectAll(".qual_fwd").attr("visibility","visible");
                  focus.selectAll(".qual_rev").attr("visibility","visible");
                }
                if(w<2000){focus.selectAll(".short").attr("visibility","visible");}
            }
        }

        function setPeakLabel(calls,label,offset){
            //Bind data
            var text = focus.append("g").selectAll("text").data(calls)
            //Ender
            text = text.enter().append("text");
            //Update  - not working correctly still have to remove old data
            text.attr("class",function(d){
                    if(label.indexOf("user") > -1){
                        return "peak_label short user ".concat("id").concat(d["id"]);
                    }else if(label.indexOf("call") > -1){
                        return "peak_label short call";
                    }else if(label.indexOf("mut")> -1){
                        return "peak_label short mut";
                    }else{return "peak_label short";}
                })
                .text(function(d){
                    if(label.indexOf("mut") > -1){
                        return d[label].toLowerCase();
//                    }else if(label.indexOf("user") > -1 && d[label]==="low qual"){
//                        return "N";
                    } else {
                        return d[label];}
                 })
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"],d["trace_peak"]);})
                .attr("y",label_pos[label]+offset)
                .attr("fill", "black")
                .attr("opacity", 0.8)
                .attr("font-family", "sans-serif")
                .attr("font-size",function(){if(label.indexOf("user")>-1){return "12px";}else{return "11px";}})
                .attr("stroke",function(d) {
                    if      (d[label] === "A"){ return "#33CC33"; }
                    else if (d[label] === "C"){ return "#0000FF"; }
                    else if (d[label] === "G"){ return "#000000"; }
                    else if (d[label] === "T"){ return "#FF0000"; }
                    else if (d[label] === "-"){ return "black"; }
                    else    {                   return "orange"; }});
        }
        function showVarInMinimap(choices){
            //genomic
            //console.log("changing choices");
            context.selectAll("lines.choices").data(choices).enter()
    			      .append("line")
                .attr("class","minimap context")
      			    .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y1",4)
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
                .attr("class","minimap context")
  				      .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y1",26)
  				      .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y2",38)
  				      .attr("stroke-width",3)
  				      .attr("stroke",function(d) {
  				          if      (d["user_sample"] === "A"){ return "#33CC33"; }
  				          else if (d["user_sample"] === "C"){ return "#0000FF"; }
  				          else if (d["user_sample"] === "G"){ return "#000000"; }
  				          else if (d["user_sample"] === "T"){ return "#FF0000"; }
  				          else if (d["user_sample"] === "-"){ return "white"; }
  				          else    {                           return "yellow";  }});
              context.selectAll("lines.choices").data(choices).enter()
    			      .append("line")
                .attr("class","minimap context")
  				      .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y1",38)
  				      .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				      .attr("y2",46)
  				      .attr("stroke-width",3)
  				      .attr("stroke",function(d) {
  				          if      (d["user_mut"] === "A"){ return "#33CC33"; }
  				          else if (d["user_mut"] === "C"){ return "#0000FF"; }
  				          else if (d["user_mut"] === "G"){ return "#000000"; }
  				          else if (d["user_mut"] === "T"){ return "#FF0000"; }
  				          else if (d["user_mut"] === "-"){ return "white"; }
  				          else    {                        return "yellow";  }});
            focus.append("g").selectAll("variance_indicator").data(choices).enter()  //variance indicator
  		         .append("line").attr("class","peak_label short line var_noise_indic")
  		         .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y1",140)
  		         .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y2",400)
  		         .attr("stroke-width",20).attr("stroke","rgba(255,0,0,0.15)").attr("stroke-dasharray","2,8");
        }
        function showNoiseInMinimap(noisy_neighbors){
            context.selectAll("lines.noisy_neighbors").data(noisy_neighbors).enter()
    			      .append("line")
                .attr("class","minimap context")
      			    .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y1",-3)
      			    .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y2",3)
      			    .attr("stroke-width",3)
      			    .attr("stroke", "brown");
            focus.append("g").selectAll("noise_indicator").data(noisy_neighbors).enter()  //noise indicator
  		         .append("line").attr("class","peak_label short line var_noise_indic")
  		         .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y1",140)
  		         .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y2",400)
  		         .attr("stroke-width",10).attr("stroke","brown").attr("opacity",0.2).attr("stroke-dasharray","1,3");
        }
        Shiny.addCustomMessageHandler("goto",
            function(message) {
                setBrush(Number(message)-180,Number(message)+200);
                focus.selectAll(".scope").attr("opacity",0);
                focus.selectAll(".".concat("scope").concat(message))
                .transition().attr("opacity", 1);
            }
        );
        Shiny.addCustomMessageHandler("input_change",
            function(message){
                if(message<brush.extent()[0] || message>brush.extent()[1]){
                   setBrush(Number(message)-100,Number(message)+120);
                }
                focus.selectAll(".scope").attr("opacity",0);
                focus.selectAll(".".concat("scope").concat(message))
                       .transition().attr("opacity", 1);
            }
        );
        Shiny.addCustomMessageHandler("s2n_min",
            function(message) {
                mult = Number(message);
                //console.log(mult);
                redraw();
            }
        );
        Shiny.addCustomMessageHandler("show",
            function(message){
                //console.log(message);
                if(message==="TRUE"){
                    focus.selectAll(".call").attr("opacity",0.8);
                    redraw();
                }else if(message==="FALSE"){
                    focus.selectAll(".call").attr("opacity",0);
                    redraw();
                }
            }
        );
        function joinView(join){
            if(join==="FALSE"){
                heightScale_fwd = heightScale_fwd_split;
                heightScale_rev = heightScale_rev_split;
                redraw();
                join = false;
            }else if(join==="TRUE"){
                split_peak_offset = 100;
                heightScale_fwd = heightScale;
                heightScale_rev = heightScale;
                redraw();
                join = true;
            }
        }
        
        Shiny.addCustomMessageHandler("join",
            function(message){
                joinView(message);
                }
        );
        Shiny.addCustomMessageHandler("opac_f",
            function(message){
                //console.log(message);
                focus.selectAll(".line_f").attr("opacity",Number(message));
                focus.selectAll("g").selectAll(".area_fwd").attr("opacity",Number(message)/5);
                focus.selectAll(".qual_fwd").attr("opacity",Number(message)*0.8);
//                focus.selectAll(".fwd").attr("opacity",Number(message));
            }
        );
        Shiny.addCustomMessageHandler("opac_r",
            function(message){
                //console.log(message);
                focus.selectAll(".line_r").attr("opacity",Number(message));
                focus.selectAll("g").selectAll(".area_rev").attr("opacity",Number(message)/5);
                focus.selectAll(".qual_rev").attr("opacity",Number(message)*0.8);
//                focus.selectAll(".rev").attr("opacity",Number(message));
            }
        );
        function callShiny(id,trace_peak){
            //console.log(message);
            Shiny.onInputChange("pos_click", {id: id});
/*
            setBrush(Number(trace_peak)-100,Number(trace_peak)+120);
            focus.append("g").selectAll("position_indicator")  //position indicator
  		         .append("line").attr("class","peak_label short line var_noise_indic")
  		         .attr("x1",function(d){return widthScale(Number(trace_peak));})
  		         .attr("y1",0)
  		         .attr("x2",function(d){return widthScale(Number(trace_peak));})
  		         .attr("y2",180)
  		         .attr("stroke-width",8).attr("stroke","rgba(0,0,255,0.1)").attr("stroke-dasharray",2);
*/
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
            .on('a', brushMoveLeft())
            .on('w', brushZoomIn())
            .on('d', brushMoveRight())
            .on('s', brushZoomOut())
            );


        //passing arguments
        //this enables to access vars and functions from the render function as instance.*
        return {
            svg:     svg,
            line_fwd: line_fwd,
            line_rev: line_rev,
            noise_area_fwd: noise_area_fwd,
            noise_area_rev: noise_area_rev,
            label_pos:  label_pos,
            linec:   linec,
            context: context,
	        brush:   brush,
            join:    join,
            joinView:joinView,
            focus:   focus,
            redraw:  redraw,
            widthScale:  widthScale,
            width2Scale: width2Scale,
            heightScale: heightScale,
            heightScale_fwd_split: heightScale_fwd_split,
            heightScale_rev_split: heightScale_rev_split,
            width:   w,
            height:  h,
	        height2: height2,
            instanceCounter: instanceCounter,
            reHeight: reHeight,
            setBrush: setBrush,
            setPeakLabel: setPeakLabel,
            showVarInMinimap: showVarInMinimap,
            showNoiseInMinimap: showNoiseInMinimap,
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
                rev = 1;
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
  			var noisy_neighbors     = HTMLWidgets.dataframeToD3(x["noisy_neighbors"]);
  			var domain_y    = x["intrexdat"]["max_y"];
  			instance.max_y  = domain_y;
  			var domain_x    = x["intrexdat"]["max_x"];
  			instance.max_x  = domain_x;
  			var intrex      = HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"])
  			instance.intrex = intrex;

  			var svg     = instance.svg;
  			var line_fwd   = instance.line_fwd;
            var line_rev   = instance.line_rev;
            var noise_area_fwd = instance.noise_area_fwd;
            var noise_area_rev = instance.noise_area_rev;
  			var linec   = instance.linec;
  			var focus   = instance.focus;
  			var context = instance.context;
  			var brush   = instance.brush;
  			var widthScale  = instance.widthScale;
  			var heightScale = instance.heightScale;
            var width2Scale = instance.width2Scale;

  			widthScale.domain([0,domain_x]);
  			width2Scale.domain([0,domain_x]);
  			heightScale.domain([0,domain_y]);
  	        instance.heightScale_fwd_split.domain([0,domain_y]);
            instance.heightScale_rev_split.domain([0,domain_y]);

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
  				           if (d["splicevar"] != ''){   return "rgba(255,205,205,1.0)";
  				    } else if (/exon/.test(d["attr"])){ return "rgba(200,200,200,1.0)";
  				    } else {                            return "rgba(230,230,230,1.0)"; }
//  				           if (/exon/.test(d["attr"])){                         return "rgba(200,200,200,1.0)";
//  				    } else if (/intron9/.test(d["attr"]) & d["length"] == 133){ return "rgba(255,205,205,1.0)";
//  				    } else {                                                    return "rgba(230,230,230,1.0)"; }
  				});
  			context.selectAll("text.intrex.name").data(intrex).enter()
  				.append("text").attr("class","context")
  				.attr("x",function(d){return widthScale(d["start"]);})
  				.attr("y",-4)
  				.attr("opacity",0.8)
  				.attr("fill","black")
  				.text(function(d) {
  				    return d["attr"]+d["splicevar"];
//  				    if (/intron9/.test(d["attr"]) & d["length"] == 133){ return d["attr"]+" | "+"beta variant";
//  				    } else {                                             return d["attr"]; }
  				});
  			context.selectAll("text.intrex.start").data(intrex).enter()
  				.append("text").attr("class","context")
  				.attr("x",function(d){return widthScale(d["start"]);})
  				.attr("y",62)
  				.attr("opacity",0.8)
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
  				.attr("d",line_fwd)
  				.attr("fill","none")
  				.attr("stroke","#33CC33").attr("stroke-width",0.75);
  			group_c.selectAll("path").data([intens["C"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line_fwd)
  				.attr("fill","none")
  				.attr("stroke","#0000FF").attr("stroke-width",0.75);
  			group_g.selectAll("path").data([intens["G"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line_fwd)
  				.attr("fill","none")
  				.attr("stroke","#000000").attr("stroke-width",0.75);
  			group_t.selectAll("path").data([intens["T"]]).enter()
  				.append("path").attr("class","path line_f")
  				.attr("d",line_fwd)
  				.attr("fill","none")
  				.attr("stroke","#FF0000").attr("stroke-width",0.75);

        //noise indocator fwd

        var a_noise_fwd = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
        var group_noise_fwd = focus.append("g");
        var group_noise_rev = focus.append("g");
        group_noise_fwd.selectAll("path").data([a_noise_fwd]).enter()
            .append("path").attr("class","area area_fwd").attr("d",noise_area_fwd)
            .attr("fill","#000000").attr("stroke","none").attr("opacity",0.15);

        //reverse strand
        if(intens_rev != ""){
            //Noise indicator rev
            var a_noise_rev = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_rev"]]);
            group_noise_rev.selectAll("path").data([a_noise_rev]).enter()
            .append("path").attr("class","area area_rev").attr("d",noise_area_rev)
            .attr("fill","#440000").attr("stroke","none").attr("opacity",0.15);
            var group_a_r = focus.append("g");
    	    	var group_c_r = focus.append("g");
  		    	var group_g_r = focus.append("g");
  		    	var group_t_r = focus.append("g");

            group_a_r.selectAll("path").data([intens_rev["A"]]).enter()
        			.append("path").attr("class","path line_r")
      				.attr("d",line_rev)
      				.attr("fill","none")
      				.attr("stroke","#33CC33").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,1,10,1");
      			group_c_r.selectAll("path").data([intens_rev["C"]]).enter()
      				.append("path").attr("class","path line_r")
      				.attr("d",line_rev)
      				.attr("fill","none")
      				.attr("stroke","#0000FF").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,1,10,1");
      			group_g_r.selectAll("path").data([intens_rev["G"]]).enter()
      				.append("path").attr("class","path line_r")
      				.attr("d",line_rev)
      				.attr("fill","none")
      				.attr("stroke","#000000").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,1,10,1");
      			group_t_r.selectAll("path").data([intens_rev["T"]]).enter()
      			    .append("path").attr("class","path line_r")
      				.attr("d",line_rev)
      				.attr("fill","none")
      				.attr("stroke","#FF0000").attr("stroke-width",0.75)
                    .attr("stroke-dasharray","20,3,10,1,10,1");
            }
            
            //on single strand always show "join view"
            if(intens_rev != ""){
                instance.joinView(true);
            }else{
                instance.joinView(instance.join);
            }
            

            //trace peak labels
            focus.append("g").selectAll("scope").data(calls).enter()        //scope (position indicator)
                 .append("rect").attr("class",function(d){return "scope ".concat("scope").concat(d["trace_peak"]);})
                 .attr("x",function(d){return (widthScale(d["trace_peak"])-12)})
                 .attr("y",-14).attr("rx",2).attr("ry",2)
                 .attr("width",24).attr("height",instance.height-90)
                 .attr("fill","rgba(155, 155, 255, 0.12)").attr("opacity",0);
            if(rev==0){
                focus.append("g").selectAll("qualities").data(calls).enter()  //quality box
      		        .append("rect").attr("class","peak_label qual_fwd q")
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
            }else{
                focus.append("g").selectAll("qualities.fwd").data(calls).enter()  //quality box
        	        .append("rect").attr("class","peak_label qual_fwd q")
      		        .attr("x",function(d){return (widthScale(d["trace_peak"])-9);})
      		        .attr("y",0).attr("rx",1).attr("ry",1)
      		        .attr("width",9)
      		        .attr("height",function(d){return d["quality_fwd"];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
                focus.append("g").selectAll("qualities.rev").data(calls).enter()  //quality box
                  .append("rect").attr("class","peak_label qual_rev q")
      		        //.attr("x",function(d){return (widthScale(d["trace_peak"]) + 900);})
      		        .attr("y",0).attr("rx",1).attr("ry",1)
      		        .attr("width",9)
      		        .attr("height",function(d){return d["quality_rev"];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
                focus.append("g").selectAll("text.qualities").data(calls).enter() //quality number
      		        .append("text").attr("class","peak_label")
      		        .text(function(d){return d["quality"];})
      		        .attr("text-anchor", "middle")
      		        .attr("x",function(d){return widthScale(d["trace_peak"]);})
                    .on("click",function(d,i){instance.callShiny(d["id"]);})
      		        .attr("y",-2)
      		        .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "10px");
            }


            focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa
              .append("text").attr("class","peak_label short")
              .text(function(d){
//                    if   (d["ord_in_cod"] == 1) {return d["aa_ref"].toUpperCase()+""+d["codon"];}
                    if   (d["ord_in_cod"] == 1) {return d["aa_ref"]+""+d["codon"];}
                    else {                       return "";}})
              .attr("text-anchor", "right")
              .attr("x",function(d){return widthScale(d["trace_peak"]);})
              .on("click",function(d,i){instance.callShiny(d["id"]);})
              .attr("y",instance.label_pos["aa"])
              .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
              .attr("stroke","#000000");
              instance.setPeakLabel(calls,"reference",0);
              instance.setPeakLabel(calls,"call",0);
              instance.setPeakLabel(calls,"mut_call_fwd",0);
              if(rev!==0){
                  instance.setPeakLabel(calls,"call_rev",0);
                  instance.setPeakLabel(calls,"mut_call_rev",0);
              }
              //default
              focus.selectAll(".call").attr("opacity",0);
              instance.setPeakLabel(calls,"user_sample",rev);
              instance.setPeakLabel(calls,"user_mut",rev);
            focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa_sample
                .append("text").attr("class","peak_label short aa_sample")
                .text(function(d){
                    if   (d["sample_ord_in_cod"] == 1) {return d["aa_sample"];}
                    else {                       return "";}})
                .attr("text-anchor", "right")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",instance.label_pos["aa_sample"]+rev)
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa_mut
                .append("text").attr("class","peak_label short aa_mut")
                .text(function(d){
                    if   (d["mut_ord_in_cod"] == 1) {return d["aa_mut"];}
                    else {                       return "";}})
                .attr("text-anchor", "right")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",instance.label_pos["aa_mut"]+rev)
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            focus.append("g").selectAll("text.seq.codon").data(calls).enter() //codon stuff
                .append("text").attr("class","peak_label")
                .text(function(d){
//                    if   (d["coding_seq"] > 0){return d["coding_seq"]+" : "+d["codon"]+"."+d["ord_in_cod"];}
                    if   (d["coding_seq"] > 0){return d["coding_seq"];}
                    else {                     return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(instance.label_pos["codon"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
                .append("text").attr("class","peak_label")
                .text(function(d){return d["gen_coord"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(instance.label_pos["gen_coord"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
/*
            focus.append("g").selectAll("text.exon_intron").data(calls).enter() //intrex
                .append("text").attr("class","peak_label")
                .text(function(d){return d["exon_intron"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){instance.callShiny(d["id"]);})
                .attr("y",(instance.label_pos["intrex"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.selectAll(".peak_label").attr("visibility","hidden")
*/
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
  				.attr("y", -16)
  				.attr("height", 80) //height2 + 10)
  				.attr("rx",3)
  				.attr("ry",3)
  				.attr("fill","rgba(255,255,255,0.3)")
  				.attr("stroke-width",2).attr("stroke","red").attr("stroke-dasharray","3,6")
  				.attr("opacity",0.6);

            //In SVG, z-index is defined by the order the element appears in the document
            //http://stackoverflow.com/questions/17786618/how-to-use-z-index-in-svg-elements
            //the vars on the minimap are shown last sto that they stay on to
            instance.showVarInMinimap(choices);
            instance.showNoiseInMinimap(noisy_neighbors);
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
            console.log("render");
  			if(x["intrexdat"]["max_y"]!= instance.max_y){
  				instance.reHeight(x["intrexdat"]["max_y"]);
          instance.max_y = x["intrexdat"]["max_y"];
  			}else if(x.choices != instance.choices){
                var choices = HTMLWidgets.dataframeToD3(x["choices"]);
                var noisy_neighbors = HTMLWidgets.dataframeToD3(x["noisy_neighbors"]);
                var calls   = HTMLWidgets.dataframeToD3(x["calls"]);
                var rev = 0;  //offset on labels in case we have alternative reference
                if(x["intens_rev"] !== null){
                    var intens_rev = x["intens_rev"];
                    rev = 1;
                }
                instance.focus.selectAll(".peak_label short line").remove();
                instance.context.selectAll(".minimap").remove();
                instance.showVarInMinimap(choices);
                instance.choices = x.choices;
                instance.showNoiseInMinimap(noisy_neighbors);
                instance.noisy_neighbors = x.noisy_neighbors;
                instance.focus.selectAll(".user").remove();
                instance.setPeakLabel(calls,"user_sample",rev);
                instance.setPeakLabel(calls,"user_mut",rev);

                instance.focus.selectAll(".var_noise_indic").remove();
                instance.focus.append("g").selectAll("variance_indicator").data(choices).enter() //variance indicator
                    .append("line").attr("class","peak_label short line var_noise_indic")
                    .attr("x1",function(d){return instance.widthScale(d["trace_peak"]);})
                    .attr("y1",140)
                    .attr("x2",function(d){return instance.widthScale(d["trace_peak"]);})
                    .attr("y2",400)
                    .attr("stroke-width",20).attr("stroke","rgba(255,0,0,0.15)").attr("stroke-dasharray","2,8");
                instance.focus.append("g").selectAll("noise_indicator").data(noisy_neighbors).enter()  //noise indicator
      		         .append("line").attr("class","peak_label short line var_noise_indic")
      		         .attr("x1",function(d){return instance.widthScale(d["trace_peak"]);})
      		         .attr("y1",140)
      		         .attr("x2",function(d){return instance.widthScale(d["trace_peak"]);})
      		         .attr("y2",400)
      		         .attr("stroke-width",10).attr("stroke","brown").attr("opacity",0.2).attr("stroke-dasharray","1,3");

                instance.focus.selectAll("g").selectAll(".aa_sample").remove();
                instance.focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa_sample
                    .append("text").attr("class","peak_label short aa_sample")
                    .text(function(d){
                        if   (d["sample_ord_in_cod"] == 1) {return d["aa_sample"];}
                        else {                             return "";}})
                    .attr("text-anchor", "right")
                    .attr("x",function(d){return instance.widthScale(d["trace_peak"]);})
                    .on("click",function(d,i){instance.callShiny(d["id"]);})
                    .attr("y",instance.label_pos["aa_sample"]+rev)
                    .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                    .attr("stroke","#000000");
                instance.focus.selectAll("g").selectAll(".aa_mut").remove();
                instance.focus.append("g").selectAll("text.seq.aa").data(calls).enter() //aa_mut
                    .append("text").attr("class","peak_label short aa_mut")
                    .text(function(d){
                        if   (d["mut_ord_in_cod"] == 1) {return d["aa_mut"];}
                        else {                       return "";}})
                    .attr("text-anchor", "right")
                    .attr("x",function(d){return instance.widthScale(d["trace_peak"]);})
                    .on("click",function(d,i){instance.callShiny(d["id"]);})
                    .attr("y",instance.label_pos["aa_mut"]+rev)
                    .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                    .attr("stroke","#000000");
                var w = instance.brush.extent()[1]-instance.brush.extent()[0] ;
                var focus = instance.focus;
                //conditional visibility !duplicate code! remove when setpeak label is fixed and instance.focus.selectAll(".user").remove(); is no longer required
                if(w<260){
                    focus.selectAll(".peak_label").attr("visibility","visible");
                    if(w==0){focus.selectAll(".peak_label").attr("visibility","hidden");}
                }else{
                    focus.selectAll(".peak_label").attr("visibility","hidden");
                    if(w<800){
                      focus.selectAll(".qual_fwd").attr("visibility","visible");
                      focus.selectAll(".qual_rev").attr("visibility","visible");
                    }
                    if(w<2000){focus.selectAll(".short").attr("visibility","visible");}
                }

  			} else { console.log(x) }
        }
    }
});
