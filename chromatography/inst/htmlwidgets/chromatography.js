HTMLWidgets.widget({

    name: 'chromatography',
    type: 'output',

    initialize: function(el, w, h) {
        var instanceCounter = 0;
        var intrex = "";
        var max_x = 0;
        var max_y = 0;
        //var margin  = {top: 10,  right: 10, bottom: 100, left: 40},   //minimap on bottom
        //    margin2 = {top: 430, right: 10, bottom: 20,  left: 40},   //minimap on top
        var margin  = {top: 55, right: 10,bottom: 20, left:10},
            margin2 = {top: 20, right: 10,bottom: 420,  left:10},
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
        var call_opacity = 0;
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
        var svg = d3.select(el).append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom);
        var brush = d3.svg.brush().on("brushend", brushed);
        var brush_fw = d3.svg.brush().on("brushend",brushed_fw);
        var brush_rv = d3.svg.brush().on("brushend",brushed_rv);
        var brush_fw_g;
        var brush_rv_g;
        var brush_fw_extent;
        var brush_rv_extent;
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
        var join = "FALSE";

        var label_pos = {};        //map for pisitioning labels representing called base
        label_pos["reference"]      =  10 + margin.top - 10;
        label_pos["aa"]             =  25 + margin.top - 10;
        label_pos["aa_sample"]      =  40 + margin.top - 10;
        label_pos["user_sample"]    =  55 + margin.top - 10;
        label_pos["user_mut"]       =  70 + margin.top - 10;
        label_pos["aa_mut"]         =  83 + margin.top - 10;
        label_pos["codon"]          = 105 + margin.top - 10;
        label_pos["gen_coord"]      = 115 + margin.top - 10;
//        label_pos["intrex"]         = 120;
        label_pos["call"]           = 130 + margin.top - 10;
        label_pos["mut_call_fwd"]   = 145 + margin.top - 10;
        label_pos["call_rev"]       = 170 + margin.top - 10;
        label_pos["mut_call_rev"]   = 185 + margin.top - 10;
        label_pos["qual"]           =   0 + margin.top - 10;
        //variables grouped at one place so its easier to read their order and determine their visibility (who's on top)
        var scope_g    = focus.append("g");
        var quals_g    = focus.append("g");
        var quals_g_r  = focus.append("g");
        var quals_txt  = focus.append("g");
        var gLine_a    = focus.append("g");
      	var gLine_c    = focus.append("g");
  		var gLine_g    = focus.append("g");
  		var gLine_t    = focus.append("g");
        var gLine_a_r  = focus.append("g");
        var gLine_c_r  = focus.append("g");
  		var gLine_g_r  = focus.append("g");
  		var gLine_t_r  = focus.append("g");
        var gNoise_fwd = focus.append("g");
        var gNoise_rev = focus.append("g");
        var text_ref          = focus.append("g");
        var text_call         = focus.append("g");
        var text_call_rev     = focus.append("g");
        var text_user_sample  = focus.append("g");
        var text_user_mut     = focus.append("g");
        var text_mut_call_fwd = focus.append("g");
        var text_mut_call_rev = focus.append("g");
        var var_ind_focus     = focus.append("g");
        var full_line         = context.append("g");
        var exon_boxes        = context.append("g");
        var intrex_txt        = context.append("g");
        var intrex_num        = context.append("g");
        var varim_gen         = context.append("g");
        var varim_user_s      = context.append("g");
        var varim_user_m      = context.append("g");
        var noisyn_con        = context.append("g");  //noisy neighbours
        var noisyn_foc        = focus.append("g");
        var aa_ref            = focus.append("g");
        var aa_sample         = focus.append("g");
        var aa_mut            = focus.append("g");

        function resetHandlers_fw(brush_fw_g){
            var oldMousedown = brush_fw_g.on('mousedown.brush');
            brush_fw_g.on('mousedown.brush', function() {
                //console.log("md"+brush_fw.extent());
                brush_fw_g.on('mouseup.brush', function() {
                    clearHandlers();
                });
                //console.log(d3.event.target.className.baseVal);
                if(d3.event.target.className.baseVal == "hook"){
                    brush_fw_g.on('mousemove.brush', function() {
                        clearHandlers();
                        oldMousedown.call(this);
                        if(brush_fw.extent()[0]>10){
                            console.log("jump detected");
                            brush_fw.extent(brush_fw_extent);
                            brush_fw(brush_fw_g);
                            brush_fw.event(d3.select(".brush_fw").transition().delay(1));
                        }
                        //console.log(this);
                        //console.log("brush from hook" + brush_fw.extent());
                        //brushg.on('mousemove.brush').call(this);
                    });
                }else{
                    brush_fw_g.on('mousemove.brush', function() {
                        clearHandlers();
                        //oldMousedown.call(this);
                        //brushg.on('mousemove.brush').call(this);
                    });
                }
                function clearHandlers() {
                    brush_fw_g.on('mousemove.brush', null);
                    brush_fw_g.on('mouseup.brush', null);
                }
            })
        }
        function resetHandlers_rv(brush_rv_g){
            var oldMousedown = brush_rv_g.on('mousedown.brush');
            brush_rv_g.on('mousedown.brush', function() {
                console.log("md"+brush_rv.extent());
                brush_rv_g.on('mouseup.brush', function() {
                    clearHandlers();
                });
                //console.log(d3.event.target.className.baseVal);
                if(d3.event.target.className.baseVal == "hook"){
                    brush_rv_g.on('mousemove.brush', function() {
                        clearHandlers();
                        oldMousedown.call(this);
                        //console.log("brush from hook" + brush_rv.extent());
                        //brushg.on('mousemove.brush').call(this);
                    });
                }else{
                    brush_rv_g.on('mousemove.brush', function() {
                        clearHandlers();
                        //oldMousedown.call(this);
                        //brushg.on('mousemove.brush').call(this);
                    });
                }
                function clearHandlers() {
                    brush_rv_g.on('mousemove.brush', null);
                    brush_rv_g.on('mouseup.brush', null);
                }
            })
        }
        function finish_fwBrushInit(to,rev){
            console.log("brush fw finish");
            brush_fw.x(widthScale);
            if(brush_fw_g==undefined){
                 brush_fw_g = focus.append("g")
                    .attr("class","brush_fw")
                    .call(brush_fw);
            };
            brush_fw_g.selectAll(".resize").append("rect")
                 .attr("class","hook")
                 .attr("fill", "red")
                 .attr("width", function(d,i){return i ? 0 : 2;})
                 .attr("x", function(d, i) {
                     return i ? -10 : -3;
                 })
                 .attr("rx", 2);

            brush_fw_g.selectAll("rect")
                .attr("y",150)
                .attr("height",function (d){if(rev!=0){return 120;}else{return 280;}});

            brush_fw_g.selectAll(".extent")
                .attr("fill","red")
                .attr("opacity",0.09);
            console.log(widthScale);
            console.log("brush fw extent in init",to);
            brush_fw.extent([0,to]);
            brush_fw_extent = brush_fw.extent();
            //brush_fw(brush_fw_g);
            //brush_fw.event(d3.select(".brush_fw"));
            //resetHandlers_fw(brush_fw_g);
        }
        function finish_rvBrushInit(from,to){
            if(to == 0){
                brush_rv_g = undefined;
            }else{
                brush_rv.x(widthScale);
                if(brush_rv_g==undefined){
                     brush_rv_g = focus.append("g")
                        .attr("class","brush_rv")
                        .call(brush_fw);
                };
                brush_rv_g.selectAll(".resize").append("rect")
                     .attr("class","hook")
                     .attr("fill", "red")
                     .attr("width", function(d,i){return i ? 2 : 0;})
                     .attr("x", function(d, i) {
                         return i ? 0 : -3;
                     })
                     .attr("rx", 2);

                brush_rv_g.selectAll("rect")
                    .attr("y",310)
                    .attr("height",120);
                brush_rv_g.selectAll(".extent")
                    .attr("fill","red")
                    .attr("opacity",0.09);
                brush_rv.extent([from,to]);
                brush_rv_extent = brush_rv.extent();
                //brush_rv(brush_rv_g);
                //brush_rv.event(d3.select(".brush_rv"));
                resetHandlers_rv(brush_rv_g);
                //console.log("setting brush rv" + brush_rv_extent);
            }
        }

        function brushed() { redraw(); }
        function brushed_fw() {
            console.log("1brushing fw to: " + brush_fw_extent);
            //brush_fw.empty() ? width2Scale.domain() : brush_fw.extent();
            brush_fw_extent = brush_fw.extent();
            Shiny.onInputChange("brush_fw", {coord: brush_fw_extent[1]});
            console.log("2brushing fw to: " + brush_fw_extent);

        }
        function brushed_rv() {
            brush_rv.empty() ? width2Scale.domain() : brush_rv.extent();
            brush_rv_extent = brush_rv.extent();
            //console.log("setting brush rv" + brush_rv_extent);
        }

        //setting brush programmatically
        function setBrush(start,end){
            //context.call(brush.extent([start,end]));
            //console.log(start+en);
            brush.extent([start,end]);
            brush(d3.select(".brush").transition());
            brush.event(d3.select(".brush").transition().delay(100));
            redraw();
        }
        function reHeight(domain_y){
            heightScale.domain([0,domain_y-1]);
            heightScale_fwd_split.domain([0,domain_y-1]);
            heightScale_rev_split.domain([0,domain_y-1]);
            redraw();
        }
        function redraw()  {
            console.log("redraw");
            widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            //console.log("old"+old);
            //console.log("extent"+brush_fw.extent());
            //console.log("setting brush_fw_extent in redraw: "+ brush_fw_extent);
            //brush_fw.extent(brush_fw_extent);
            //brush_fw(brush_fw_g);
            //brush_fw.event(d3.select(".brush_fw").transition().delay(1));
            //resetHandlers_fw(brush_fw_g);
            //if(brush_rv_g!=undefined){ //this is not good
            //    brush_rv.extent(brush_rv_extent);
            //    brush_rv(brush_rv_g);
            //    brush_rv.event(d3.select(".brush_rv").transition().delay(1));
            //    resetHandlers_rv(brush_rv_g);
            //}
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
                if(w===0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                focus.selectAll(".peak_label").attr("visibility","hidden");
                if(w<800){
                  focus.selectAll(".qual_fwd").attr("visibility","visible");
                  focus.selectAll(".qual_rev").attr("visibility","visible");
                }
                if(w<2000){focus.selectAll(".short").attr("visibility","visible");}
            }
        }
        function joinView(j){
            j = typeof j !== 'undefined' ? j : join;
            focus.selectAll(".rev")
                 .attr("y",(j=="FALSE")*110 + label_pos["call_rev"]);
            focus.selectAll(".mut_rev")
                 .attr("y",(j=="FALSE")*110 + label_pos["mut_call_rev"]);
            if(j==="FALSE"){
                heightScale_fwd = heightScale_fwd_split;
                heightScale_rev = heightScale_rev_split;
                redraw();
            }else if(j==="TRUE"){
                split_peak_offset = 100;
                heightScale_fwd = heightScale;
                heightScale_rev = heightScale;
                redraw();
            }
        }
        function updateLine(data,base,rev){
            switch(base) {
                case "A": if(rev){var g = gLine_a_r;}else{var g = gLine_a;} var col = "#33CC33"; break;
                case "C": if(rev){var g = gLine_c_r;}else{var g = gLine_c;} var col = "#0000FF"; break;
                case "G": if(rev){var g = gLine_g_r;}else{var g = gLine_g;} var col = "#000000"; break;
                case "T": if(rev){var g = gLine_t_r;}else{var g = gLine_t;} var col = "#FF0000"; break;
            }
            //console.log("updating line + ",rev);
            if(rev){var c = "line_r";var l = line_rev;}
            else{var c = "line_f";var l = line_fwd;}

            var line = g.selectAll("path").data(data); //join
            line.enter().append("path");                //enter
            line.attr("class","path "+c)
                .attr("d",l)
                .attr("fill","none").attr("stroke",col)
                .attr("stroke-width",0.75);         // on reverse attr("stroke-dasharray","20,3,10,1,10,1");
            line.exit().remove();                      //exit
        }
        function setNoiseArea(fwd,rev){
            var gnf = gNoise_fwd.selectAll("path").data([fwd]);
            gnf.enter().append("path");
            gnf.attr("class","area area_fwd").attr("d",noise_area_fwd)
               .attr("fill","#000000").attr("stroke","none").attr("opacity",0.15);
            gnf.exit().remove();
            if(typeof rev !== 'undefined'){
                var gnr = gNoise_rev.selectAll("path").data([rev]);
                gnr.enter().append("path");
                gnr.attr("class","area area_rev").attr("d",noise_area_rev)
                   .attr("fill","#440000").attr("stroke","none").attr("opacity",0.15);
                gnr.exit().remove();
            }
        }
        function setPeakLabel(calls,label,opacity){
            opacity = typeof opacity !== 'undefined' ? opacity : 0.8;
            switch(label) {
                case "reference":    var t = text_ref;          c = "ref"; break;
                case "call":         var t = text_call;         c = "call"; break;
                case "call_rev":     var t = text_call_rev;     c = "call rev"; break
                case "user_sample":  var t = text_user_sample;  c = "user"; break;
                case "user_mut":     var t = text_user_mut;     c = "user"; break;
                case "mut_call_fwd": var t = text_mut_call_fwd; c = "call call_fwd"; break;
                case "mut_call_rev": var t = text_mut_call_rev; c = "call mut_rev mut"; break ;
            }
            var text = t.selectAll("text").data(calls);         //Join
            text.enter().append("text");                        //Enter
            text.attr("class",function(d){
                    if(label.indexOf("user") > -1) {return "peak_label short user ".concat("id").concat(d["id"]);}
                    return("peak_label short "+c);
                })
                .text(function(d){
                    if(label.indexOf("mut") > -1){
                        if(typeof(d[label])=='undefined'){ console.log(label);console.log(d);}
                        return d[label].toLowerCase();
                    } else {
                        return d[label];}
                 })
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",function(d){
                                if(label.indexOf("rev") > -1){
                                    return((join=="FALSE")*110 + label_pos[label] );
                                }else{
                                    return(label_pos[label] );
                                }})
                .attr("fill", "black")
                .attr("opacity", opacity)
                .attr("font-family", "sans-serif")
                .attr("font-size",function(){if(label.indexOf("user")>-1){return "12px";}else{return "11px";}})
                .attr("stroke",function(d) {
                    if      (d[label] === "A"){ return "#33CC33"; }
                    else if (d[label] === "C"){ return "#0000FF"; }
                    else if (d[label] === "G"){ return "#000000"; }
                    else if (d[label] === "T"){ return "#FF0000"; }
                    else if (d[label] === "-"){ return "black"; }
                    else    {                   return "orange"; }});

            text.exit().remove();
        }
        function showVarInMinimap(choices){
            //genomic
            //console.log("changing choices");
            var g = varim_gen.selectAll("line").data(choices);
            g.enter().append("line").attr("class","enter");
            g.attr("class","minimap context")
      			.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			.attr("y1",4)
      			.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      			.attr("y2",24)
      			.attr("stroke-width",3)
      			.attr("stroke",function(d) { return get_color(d["reference"])});
            g.exit().remove();
  			    // user
            var us = varim_user_s.selectAll("line").data(choices);
            us.enter().append("line").attr("class","enter");
            us.attr("class","minimap context")
  		        .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y1",26)
  				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y2",38)
  				.attr("stroke-width",3)
  				.attr("stroke",function(d) { return get_color(d["user_sample"])});
            us.exit().remove();
            var um = varim_user_m.selectAll("line").data(choices);
            //um.attr("class","update");
            um.enter().append("line").attr("class","enter");
            um.attr("class","minimap context")
  				.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y1",38)
  				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y2",46)
  				.attr("stroke-width",3)
  				.attr("stroke",function(d) { return get_color(d["user_mut"])});
            um.exit().remove();
            var v = var_ind_focus.selectAll("line").data(choices);
            v.enter().append("line").attr("class","enter");
            v.attr("class","peak_label short line var_noise_indic")
  		        .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		        .attr("y1",140)
  		        .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		        .attr("y2",400)
  		        .attr("stroke-width",20).attr("stroke","rgba(255,0,0,0.15)").attr("stroke-dasharray","2,8");
            v.exit().remove();

        }
        function showNoisyNbr(noisy_neighbors){
            var nnc = noisyn_con.selectAll("line").data(noisy_neighbors);
            nnc.enter().append("line");
            nnc.attr("class","minimap context")
      			    .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y1",-3)
      			    .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      			    .attr("y2",3)
      			    .attr("stroke-width",3)
      			    .attr("stroke", "brown");
            nnc.exit().remove();
        /*  var nnf = noisyn_foc.selectAll("line").data(noisy_neighbors);
            nnf.enter().append("line");
            nnf.attr("class","peak_label short line var_noise_indic")
  		         .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y1",140)
  		         .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		         .attr("y2",400)
  		         .attr("stroke-width",10).attr("stroke","brown").attr("opacity",0.2).attr("stroke-dasharray","1,3");
            nnf.exit().remove();
        */
        }
        function setCodingLabel(calls){
            var aa_r = aa_ref.selectAll("text").data(calls);
            aa_r.enter().append("text");
            aa_r.attr("class","peak_label short")
                .text(function(d){
//                    if   (d["ord_in_cod"] == 1) {return d["aa_ref"].toUpperCase()+""+d["codon"];}
                    if   (d["ord_in_cod"] == 1) {return d["aa_ref"]+""+d["codon"];}
                    else {                       return "";}})
                .attr("text-anchor", "right")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            aa_r.exit().remove();
            var aa_s = aa_sample.selectAll("text").data(calls);
            aa_s.enter().append("text"); //aa_sample
            aa_s.attr("class","peak_label short aa_sample")
                .text(function(d){
                    if   (d["sample_ord_in_cod"] == 1) {return d["aa_sample"];}
                    else {                       return "";}})
                .attr("text-anchor", "right")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa_sample"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            aa_s.exit().remove();
            var aa_m = aa_mut.selectAll("text").data(calls);
            aa_m.enter().append("text");
            aa_m.attr("class","peak_label short aa_mut")
                .text(function(d){
                    if   (d["mut_ord_in_cod"] == 1) {return d["aa_mut"];}
                    else {                       return "";}})
                .attr("text-anchor", "right")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa_mut"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            aa_m.exit().remove();
        }
        function setQualityLabels(calls,rev){
            q = rev ? "quality_fwd" : "quality";
            var qt = quals_txt.selectAll("text").data(calls);
            qt.enter().append("text");
            qt.attr("class","peak_label")
      		        .text(function(d){return d["quality"];})
      		        .attr("text-anchor", "middle")
      		        .attr("x",function(d){return widthScale(d["trace_peak"]);})
                    .on("click",function(d,i){callShiny(d["id"]);})
      		        .attr("y",label_pos["qual"] -2)
      		        .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "10px");
            qt.exit().remove();
            var qg = quals_g.selectAll("rect").data(calls);
            qg.enter().append("rect");
            qg.attr("class","peak_label qual_fwd q")
      		        .attr("x",function(d){return (widthScale(d["trace_peak"])-9);})
      		        .attr("y",label_pos["qual"]).attr("rx",1).attr("ry",1)
      		        .attr("width",9+(!rev)*9)
      		        .attr("height",function(d){return d[q];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
            qg.exit().remove();
            if(rev!=0){
                var qg_r = quals_g_r.selectAll("rect").data(calls);
                qg_r.enter().append("rect");
                qg_r.attr("class","peak_label qual_rev q")
      		        //.attr("x",function(d){return (widthScale(d["trace_peak"]) + 900);})
      		        .attr("y",label_pos["qual"]).attr("rx",1).attr("ry",1)
      		        .attr("width",9)
      		        .attr("height",function(d){return d["quality_rev"];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
                qg_r.exit().remove();
            }
        }
        function setIntrexBoxes(intrex){
            var eb = exon_boxes.selectAll("rect").data(intrex);
            eb.enter().append("rect");
  			eb.attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",0).attr("rx",3).attr("ry",3)
  				.attr("width",function(d){return width2Scale(d["end"]-d["start"]);})
  				.attr("height",50)
  				.attr("fill",function(d) {
  				           if (d["splicevar"] != ''){   return "rgba(255,205,205,1.0)";
  				    } else if (/exon/.test(d["attr"])){ return "rgba(200,200,200,1.0)";
  				    } else {                            return "rgba(230,230,230,1.0)"; }
  				});
            eb.exit().remove();
            var iet = intrex_txt.selectAll("text").data(intrex);
            iet.enter().append("text");
            iet.attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",-4)
  				.attr("opacity",0.8)
  				.attr("fill","black")
  				.text(function(d) {
  				    return d["attr"]+d["splicevar"];
  				});
            iet.exit().remove();
            var ien = intrex_num.selectAll("text").data(intrex);
            ien.enter().append("text");
  		    ien.attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",62)
  				.attr("opacity",0.8)
  				.text(function(d){return d["id"];})
  				.attr("fill","black");
            ien.exit().remove();
        }
        Shiny.addCustomMessageHandler("goto",
            function(message) {
                /*console.log("goto message");*/
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
        Shiny.addCustomMessageHandler("join",
            function(message){
                join = message;
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
            //console.log(id);
            Shiny.onInputChange("pos_click", {id: id});
        };
        function get_color(col){
            if      (col === "A"){ return "#33CC33"; }
            else if (col === "C"){ return "#0000FF"; }
  			else if (col === "G"){ return "#000000"; }
  			else if (col === "T"){ return "#FF0000"; }
  			else if (col === "-"){ return "white"; }
  	        return "yellow";
        }
/*        function brushZoomIn(){
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
*/
        //passing arguments
        //this enables to access vars and functions from the render function as instance.*
        return {
            label_pos:  label_pos,
            full_line: full_line,
            scope_g: scope_g,
            context: context,
	        brush:   brush,
            finish_fwBrushInit: finish_fwBrushInit,
            finish_rvBrushInit: finish_rvBrushInit,
	        join:    join,
            joinView:joinView,
            focus:   focus,
            setNoiseArea: setNoiseArea,
            updateLine: updateLine,
            widthScale:  widthScale,
            width2Scale: width2Scale,
            heightScale: heightScale,
            heightScale_fwd_split: heightScale_fwd_split,
            heightScale_rev_split: heightScale_rev_split,
            width:   w,
            height:  h,
	        height2: height2,
            call_opacity: call_opacity,
            instanceCounter: instanceCounter,
            reHeight: reHeight,
            setBrush: setBrush,
            setPeakLabel: setPeakLabel,
            showVarInMinimap: showVarInMinimap,
            showNoisyNbr: showNoisyNbr,
            setCodingLabel: setCodingLabel,
            setQualityLabels: setQualityLabels,
            setIntrexBoxes: setIntrexBoxes
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
        //the first run is actually still a part of the initialization step
        if(x.new_sample){
  		    var intens = x["intens"];
            var intens_rev = "";
            var rev = 0;  //offset on labels in case we have alternative reference
            if(x["intens_rev"] !== null){
                intens_rev = x["intens_rev"];
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
            var show_calls  = x["show_calls"];
            if(show_calls){
                instance.call_opacity = 0.8;
            }else{
                instance.call_opacity = 0;
            }
  			var domain_y    = x["intrexdat"]["max_y"];
  			instance.max_y  = domain_y;
  			var domain_x    = x["intrexdat"]["max_x"];
  			instance.max_x  = domain_x;
  			var intrex      = HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"]);
  			instance.intrex = intrex;

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
  			instance.full_line
  				.append("line").attr("class","context")
  				.attr("x1",0)
  				.attr("y1",25)
  				.attr("x2",widthScale(domain_x))
  				.attr("y2",25)
  				.attr("stroke-width",3).attr("stroke","rgba(180,180,180,1.0)");

            //intron/exon boxes
            instance.setIntrexBoxes(intrex);
  			brush.x(width2Scale);
            //instance.finish_fwBrushInit(x["brush_fw"],rev);
            //if(rev!=0){
            //    instance.finish_rvBrushInit(x["brush_rv"],domain_x);
            //}else{
            //    instance.finish_rvBrushInit(0,0);
            //}
            console.log("brush_rv: " + x["brush_rv"] );

            context.append("g")
//				.attr("class", "x brush")
                .attr("class","brush context")
                .call(brush)
      			.selectAll("rect")
  				.attr("y", -16)
                .attr("height", 80) //height2 + 10)
                .attr("rx",3)
  				.attr("ry",3)
  				.attr("fill","rgba(255,255,255,0.3)")
  				.attr("stroke-width",2).attr("stroke","red").attr("stroke-dasharray","3,6")
  				.attr("opacity",0.6);


            instance.updateLine([intens["A"]],"A",false);
            instance.updateLine([intens["C"]],"C",false);
            instance.updateLine([intens["G"]],"G",false);
            instance.updateLine([intens["T"]],"T",false);

            //noise indicator
            var a_noise_fwd = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
            //reverse strand
            if(intens_rev != ""){
                var a_noise_rev = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_rev"]]);
                instance.updateLine([intens_rev["A"]],"A",true);
                instance.updateLine([intens_rev["C"]],"C",true);
                instance.updateLine([intens_rev["G"]],"G",true);
                instance.updateLine([intens_rev["T"]],"T",true);
            }
            instance.setNoiseArea(a_noise_fwd,a_noise_rev);
            //on single strand always show "join view"
            if(intens_rev != ""){
                instance.joinView();
            }else{
                instance.joinView("TRUE");
            }

            instance.scope_g.selectAll("scope").data(calls).enter()        //scope (position indicator)
                 .append("rect").attr("class",function(d){return "scope ".concat("scope").concat(d["trace_peak"]);})
                 .attr("x",function(d){return (widthScale(d["trace_peak"])-12)})
                 .attr("y",20).attr("rx",2).attr("ry",2)
                 .attr("width",24).attr("height",instance.height-90)
                 .attr("fill","rgba(155, 155, 255, 0.12)").attr("opacity",0);
            focus.append("g").selectAll("text.seq.codon").data(calls).enter() //codon stuff
                .append("text").attr("class","peak_label")
                .text(function(d){
//                    if   (d["coding_seq"] > 0){return d["coding_seq"]+" : "+d["codon"]+"."+d["ord_in_cod"];}
                    if   (d["coding_seq"] > 0){return d["coding_seq"];}
                    else {                     return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .attr("y",(instance.label_pos["codon"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
                .append("text").attr("class","peak_label")
                .text(function(d){return d["gen_coord"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .attr("y",(instance.label_pos["gen_coord"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            instance.setPeakLabel(calls,"reference");
            instance.setPeakLabel(calls,"call");
            instance.setPeakLabel(calls,"mut_call_fwd");
            if(rev!==0){
                instance.setPeakLabel(calls,"call_rev");
                instance.setPeakLabel(calls,"mut_call_rev");
            }
            //console.log(x);
            //console.log(x["qual_present"]);
            if(x["qual_present"]){
                instance.setQualityLabels(calls,rev);
            }
            //default
            focus.selectAll(".call").attr("opacity",instance.call_opacity);
            instance.setPeakLabel(calls,"user_sample");
            instance.setPeakLabel(calls,"user_mut");
            instance.setCodingLabel(calls);
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

            instance.showVarInMinimap(choices);
            instance.showNoisyNbr(noisy_neighbors);
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
            //console.log("render");
            var calls   = HTMLWidgets.dataframeToD3(x["calls"]);
            var rev = 0;
            if(x["intens_rev"] !== null){
                rev = 1;
            }
            var a_noise_fwd = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
            var intens = x["intens"];
            instance.updateLine([intens["A"]],"A",false);
            instance.updateLine([intens["C"]],"C",false);
            instance.updateLine([intens["G"]],"G",false);
            instance.updateLine([intens["T"]],"T",false);

            if(rev){
                var a_noise_rev = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_rev"]]);
                var intens_rev = x["intens_rev"];
                instance.updateLine([intens_rev["A"]],"A",true);
                instance.updateLine([intens_rev["C"]],"C",true);
                instance.updateLine([intens_rev["G"]],"G",true);
                instance.updateLine([intens_rev["T"]],"T",true);
            }
            instance.setNoiseArea(a_noise_fwd,a_noise_rev);
            if(x["qual_present"]){
                instance.setQualityLabels(calls,rev);
            }
            if(instance.max_x != x["intrexdat"]["max_x"]){
                instance.max_x = x["intrexdat"]["max_x"];
                var domain_x    = x["intrexdat"]["max_x"];
                instance.width2Scale.domain([0,domain_x]);
                var intrex      = HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"]);
                instance.setIntrexBoxes(intrex);
            }
  			if(x["intrexdat"]["max_y"]!= instance.max_y){
  				instance.reHeight(x["intrexdat"]["max_y"]);
                instance.max_y = x["intrexdat"]["max_y"];
  			}else if(x.choices != instance.choices){
                var choices = HTMLWidgets.dataframeToD3(x["choices"]);
                var noisy_neighbors = HTMLWidgets.dataframeToD3(x["noisy_neighbors"]);
                var show_calls  = x["show_calls"];
                if(show_calls){ instance.call_opacity = 0.8; }
                else{ instance.call_opacity = 0; }
                instance.focus.selectAll(".call").attr("opacity",instance.call_opacity);
                instance.showVarInMinimap(choices);
                instance.choices = x.choices;
                instance.showNoisyNbr(noisy_neighbors);
                instance.setPeakLabel(calls,"user_sample");
                instance.setPeakLabel(calls,"user_mut");
                instance.setPeakLabel(calls,"reference")
                instance.setPeakLabel(calls,"mut_call_fwd",instance.call_opacity);
                if(rev){
                    instance.setPeakLabel(calls,"mut_call_rev",instance.call_opacity);
                }
                instance.setCodingLabel(calls);
  			} else { console.log(x) }
        }
    }
});
