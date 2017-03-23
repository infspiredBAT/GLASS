HTMLWidgets.widget({

    name: 'chromatography',
    type: 'output',

    initialize: function(el, w, h) {
        var instanceCounter = 0;
        var intrex = "";
        var max_x = 0;
        var max_y = 0;
        var margin  = {top: 65, right: 10,bottom: 20, left:10},
            margin2 = {top: 20, right: 10,bottom: 420,  left:10},
            width   = w - margin.left - margin.right,
            height  = h - margin.top  - margin.bottom,
            half_height = height/1.6;
            height2 = h - margin2.top - margin2.bottom;
        var widthScale   = d3.scaleLinear().range([0,width]),
            width2Scale  = d3.scaleLinear().range([0,width]),  //remains constant, to be used with context
            heightScale  = d3.scaleLinear().range([height,0]),
    	    height2Scale = d3.scaleLinear().range([height2,0]),
            heightScale_fwd_split = d3.scaleLinear().range([half_height,(2*half_height -  height)]),
            heightScale_rev_split = d3.scaleLinear().range([height,half_height]);
            heightScale_fwd = heightScale_fwd_split;
            heightScale_rev = heightScale_rev_split;

        var line_fwd = d3.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return heightScale_fwd(d)});
        var line_rev = d3.line()
            .x(function(d,i){return widthScale(i)})
            .y(function(d){return heightScale_rev(d)});
        var call_opacity = 0
        var mult = 2;

        var noise_area_fwd = d3.area()
                .x(function(d){return widthScale(d[0]);})
                .y0(function(d){return (heightScale_fwd(0)+2);})
                .y1(function(d){return heightScale_fwd(d[1]*mult);});
        var noise_area_rev = d3.area()
                .x(function(d){return widthScale(d[0]);})
                .y0(function(d){return heightScale_rev(0)+2;})
                .y1(function(d){return heightScale_rev(d[1]*mult);});
        var svg = d3.select(el).append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom);
        var brush = d3.brushX().extent([[0, -16], [width, 68]])
                                    .on("brush",brushed)
                                    .on("start",brush_start)
                                    .on("end",brush_end);
                                  //.on("brush",redrawLines)
                                  //.on("start",hideLabels)
                                  //.on("end",brushed);


        var zoom_height = height - 25;
        var zoom = d3.zoom()
                     .scaleExtent([1, Infinity])
                     .translateExtent([[0, 0], [width, zoom_height]])
                     .extent([[0, 0], [width, zoom_height]])
                     .on("start",zoom_start)
                     .on("zoom", zoomed)
                     .on("end",zoom_end);
        // var brush_fw = d3.svg.brush().on("brushend",brushed_fw);
        // var brush_rv = d3.svg.brush().on("brushend",brushed_rv);
        var filt_fwd_x = 0;
        var filt_rev_x = 0;
        var filt_single_rev = false;
        var filt_fwd_mini = svg.append("line")
                               .attr("x1", 10)
                               .attr("y1", 18)
                               .attr("x2", 250)
                               .attr("y2", 18)
                               .attr("clip-path","url(#clip)")
                               .attr("stroke-width", 2)
                               .style("stroke-dasharray", ("3, 3"))
                               .attr("stroke", "red")
                               .attr("opacity",1);
        var filt_rev_mini = svg.append("line")
                               //.attr("x1", 10)
                               .attr("y1", 72)
                               //.attr("x2", 250) // don't know before init
                               .attr("y2", 72)
                               .attr("clip-path","url(#clip)")
                               .attr("stroke-width", 2)
                               .style("stroke-dasharray", ("3, 3"))
                               .attr("stroke", "red")
                               .attr("opacity",1);
        // var brush_fw_g;
        // var brush_rv_g;

        // var brush_fw_extent;
        // var brush_rv_extent;
        // var frame = svg.append("rect").attr("x", 0 + margin.right)
        //                               .attr("y",0+margin.top+30)
        //                               .attr("rx",3)
        //                 			  .attr("ry",3)
        //                               .attr("width",width)
        //                               .attr("height",height-20)
        //                               .attr("stroke","rgb(70,130,180)") //steal blue
        //                               .attr("fill-opacity",0)
        //                               .attr("fill","white")
        //                               .attr("stroke-width",2);

        var context = svg.append("g")
            .attr("class", "context")
            .attr("transform", "translate(" + margin2.left + "," + margin2.top + ")");
        var focus = svg.append("g")
            .attr("class", "focus")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
        var clip = svg.append("defs").append("clipPath")
            .attr("id", "clip");
        var clip_rect = clip.append("rect")
            .attr("width", width)
            .attr("height", height);

        var defs = svg.append("defs");

        var filter = defs.append("filter")
                         .attr("id", "drop-shadow")
                         .attr("height", "130%");
        filter.append("feGaussianBlur")
              .attr("in", "SourceAlpha")
              .attr("stdDeviation", 5)
              .attr("result", "blur");
        filter.append("feOffset")
              .attr("in", "blur")
              .attr("dx", 5)
              .attr("dy", 5)
              .attr("result", "offsetBlur");
        var feMerge = filter.append("feMerge");
        feMerge.append("feMergeNode")
               .attr("in", "offsetBlur");
        feMerge.append("feMergeNode")
               .attr("in", "SourceGraphic");

        //filter for Grayscale to be used on filtered part of traces
        var filter_gs   = svg.append("filter")
                           .attr("id","monochrome");
        var colormatrix = filter_gs.append("feColorMatrix")
                            .attr("type","matrix")
                            .attr("values","2 0.5 0.5 0 0 0.5 2 0.5 0 0 0.5 0.5 2 0 0 0 0 0 1 0");


        var clip_fwd = svg.append("clipPath")
                        .attr("id", "rect_fwd_clip")
                        .append("rect")
                        .attr("x", 0)
                        .attr("y", 0)
                        .attr("width", 0)
                        .attr("height",height)
                        .attr("clip-path","url(#clip)");
        var clip_rev = svg.append("clipPath")
                        .attr("id", "rect_rev_clip")
                        .append("rect")
                        .attr("x", 0)
                        .attr("y", 0)
                        .attr("width", 0)
                        .attr("height",height);
        var join = "FALSE";
        var show_qual = "hidden";

        var label_pos = {};        //map for pisitioning labels representing called base
        label_pos["reference"]      =  10 + margin.top - 10;
        label_pos["aa"]             =  24 + margin.top - 10;
        label_pos["aa_sample"]      =  43 + margin.top - 10;
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
        var aa_ref_line       = focus.append("g");
        var aa_sample_line    = focus.append("g");
        var aa_mut_line       = focus.append("g");
        var gLine_a_gray    = focus.append("g");
      	var gLine_c_gray    = focus.append("g");
  		var gLine_g_gray    = focus.append("g");
  		var gLine_t_gray    = focus.append("g");
        var gLine_a_r_gray  = focus.append("g");
        var gLine_c_r_gray  = focus.append("g");
  		var gLine_g_r_gray  = focus.append("g");
  		var gLine_t_r_gray  = focus.append("g");
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
        var varim_ref         = context.append("g");
        var varim_user_s      = context.append("g");
        var varim_user_m      = context.append("g");
        // var noisyn_con        = context.append("g");  //noisy neighbours
        // var noisyn_foc        = focus.append("g");


        var aa_ref            = focus.append("g");
        var aa_sample         = focus.append("g");
        var aa_mut            = focus.append("g");
        var filt_line_fwd = svg.append("line")
                               .attr("class","filt_indic")
                               .attr("x1", 0)
                               .attr("y1", 0)
                               .attr("x2", 0)
                               .attr("y2", 0)
                               .attr("stroke-width", 2)
                               .attr("stroke", "red")
                               .attr("opacity",1);
        var filt_line_rev = svg.append("line")
                               .attr("x1", 0)
                               .attr("y1", 0)
                               .attr("x2", 0)
                               .attr("y2", 0)
                               .attr("stroke-width", 2)
                               .attr("stroke", "red")
                               .attr("opacity",1);
        focus.append("rect")
                 .attr("class", "zoom")
                 .attr("y",25 +113)
                 .attr("width", width)
                 .attr("height", zoom_height - 108)
                 .call(zoom).on("mousewheel.zoom", null)
             .on("DOMMouseScroll.zoom", null) // disables older versions of Firefox
             .on("wheel.zoom", null)

        // function finish_fwBrushInit(to,rev,height){
        //     //console.log("brush fw finish:" + to + " " + rev);
        //
        //     brush_fw.x(widthScale);
        //     if(brush_fw_g==undefined){
        //          brush_fw_g = focus.append("g")
        //             .attr("class","brush_fw excl")
        //             .call(brush_fw);
        //     };
        //     brush_fw_g.selectAll(".resize").append("rect")
        //          .attr("class","hook")
        //          .attr("fill", "red")
        //          .attr("width", function(d,i){return i ? 0 : 2;})
        //          .attr("x", function(d, i) {
        //              return i ? -10 : -3;
        //          })
        //          .attr("rx", 2);
        //
        //     brush_fw_g.selectAll("rect")
        //         .attr("y",138)
        //         //.attr("height",function (d){if(rev!=0){return 120;}else{return 280;}});
        //         .attr("height",height);
        //
        //     brush_fw_g.selectAll(".extent")
        //         .attr("fill","white")
        //         .attr("opacity",0.6);
        //     //console.log(widthScale);
        //     //console.log("brush fw extent in init",to);
        //     brush_fw.extent([0,to]);
        //     brush_fw_extent = brush_fw.extent();
        //     //brush_fw(brush_fw_g);
        //     //brush_fw.event(d3.select(".brush_fw"));
        //     //resetHandlers_fw(brush_fw_g);
        //     brush_fw_mini.attr("x2",width2Scale(to));
        // }

        // function finish_rvBrushInit(from,to,rev = null){
        //     //console.log("brush fw finish:" + from + " " + to);
        //     mod = 0;
        //     if(!rev){
        //         mod = 160;
        //     }
        //     if(to == 0){
        //         if(brush_rv_g!= undefined){
        //             brush_rv_g.selectAll("rect")
        //                 .attr("height",0);
        //     }
        //     }else{
        //         brush_rv.x(widthScale);
        //         if(brush_rv_g==undefined){
        //              brush_rv_g = focus.append("g")
        //                 .attr("class","brush_rv excl")
        //                 .call(brush_fw);
        //         };
        //         brush_rv_g.selectAll(".resize").append("rect")
        //              .attr("class","hook")
        //              .attr("fill", "red")
        //              .attr("width", function(d,i){return i ? 2 : 0;})
        //              .attr("x", function(d, i) {
        //                  return i ? 0 : -3;
        //              })
        //              .attr("rx", 2);
        //
        //         brush_rv_g.selectAll("rect")
        //             .attr("y",295 - mod)
        //             .attr("height",140 + mod);
        //         brush_rv_g.selectAll(".extent")
        //             .attr("fill","white")
        //             .attr("opacity",0.6);
        //         brush_rv.extent([from,to]);
        //         brush_rv_extent = brush_rv.extent();
                //  brush_rv_mini.attr("x1",width2Scale(brush_rv_extent[0]))
                //               .attr("x2",width2Scale.domain()[1]);
        //         //brush_rv(brush_rv_g);
        //         //brush_rv.event(d3.select(".brush_rv"));
        //         //console.log("setting brush rv" + brush_rv_extent);
        //     }
        // }

       function brush_start(){
           hideLabels();
       }
       function brushed(){
           var s = d3.event.selection || width2Scale().domain();
           var t = s.map(width2Scale.invert, width2Scale);
           widthScale.domain(t);
           svg.select(".zoom").call(zoom.transform, d3.zoomIdentity
              .scale(width / (s[1] - s[0]))
              .translate(-s[0], 0));
          redrawLines();

       }
       function brush_end() {
           var s = [0,0];
           if(d3.event == null){
               var s = d3.brushSelection(d3.selectAll(".brush"));
               if (s == null){console.log(d3.brushSelection(d3.selectAll(".brush")));console.log('return');return;}
           }
           else if (d3.event.sourceEvent && d3.event.sourceEvent.type === "zoom") {return;}
           else {s = d3.event.selection || width2Scale.range();}
           //console.log(s);
           var t = s.map(width2Scale.invert, width2Scale)
           widthScale.domain(t);
           updateFilt();
           redraw();
       }
       function zoom_start(){
           if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return;
           var t = d3.event.transform;
           hideLabels();
           console.log("zoomestart");
       }
       function zoomed() {
           if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return;
           var t = d3.event.transform;
           widthScale.domain(t.rescaleX(width2Scale).domain());
           //focus.select(".area").attr("d", area);
           var s = widthScale.range().map(t.invertX,t)
           context.select(".selection").attr("x",s[0])
                                       .attr("width",s[1]-s[0])
           //context.select(".brush").call(brush.move, widthScale.range().map(t.invertX, t));
           redrawLines();
           console.log("zoome");
       }
       function zoom_end() {
           if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return;
           if(d3.event.sourceEvent==null) return;
           var t = d3.event.transform;
           widthScale.domain(t.rescaleX(width2Scale).domain());
           //focus.select(".area").attr("d", area);
           context.select(".brush").call(brush.move, widthScale.range().map(t.invertX, t));
           console.log("zoomeend");
       }
        // function brushed_fw() {
        //
        //     if (!d3.event.sourceEvent) return; // only transition after input
        //     //resetHandlers_fw();
        //     var extent0 = brush_fw.extent();
        //     var extent1 = [0,extent0[1]];
        //     Shiny.onInputChange("brush_fw", {coord: extent1[1]});
        //     brush_fw_extent = extent1;
        //     brush_fw_mini.attr("x2",width2Scale(extent1[1]));
        //     //console.log(extent1[1]);
        //     //console.log(width2Scale(extent1[1]));
        //     d3.select(this).transition()
        //         .call(brush_fw.extent(extent1))
        //         .call(brush_fw.event);
        //
        // }

        // function brushed_rv() {
        //     if (!d3.event.sourceEvent) return; // only transition after input
        //     //resetHandlers_fw();
        //     var extent0 = brush_rv.extent();
        //     var extent1 = [extent0[0],width2Scale.domain()[1]];
        //     Shiny.onInputChange("brush_rv", {coord: extent1[0]});
        //     brush_rv_mini.attr("x1",width2Scale(extent1[0]))
        //                  .attr("x2",width2Scale.domain()[1]);
        //     brush_rv_extent = extent1;
        //     d3.select(this).transition()
        //         .call(brush_rv.extent(extent1))
        //         .call(brush_rv.event);
        // }

        // //setting brush programmatically
        function setBrush(start,end){
            brush.move(d3.selectAll(".brush"),[width2Scale(start),width2Scale(end)]);
        }

        function reHeight(domain_y){
            heightScale.domain([0,domain_y-1]);
            heightScale_fwd_split.domain([0,domain_y-1]);
            heightScale_rev_split.domain([0,domain_y-1]);
            redraw();
        }

        function rescaleWidth(width_in){
                var width_new   = width_in - margin.left - margin.right;
                width = width_new;



                var sel = d3.brushSelection(d3.select(".brush").node());
                //console.log("sel:"+ sel);
                //if a brush selection exists redraw it in new scale
                if(sel != null){
                    var sin = sel.map(width2Scale.invert, width2Scale);
                    width2Scale.range([0,width_new]);
                    widthScale.range([0,width_new]);
                    context.selectAll(".brush_context").remove();
                    brush = d3.brushX().extent([[0, -16], [width, 68]])
                                                 .on("brush",brushed)
                                                 .on("start",brush_start)
                                                 .on("end",brush_end);
                    context.append("g").attr("class","brush context brush_context")
                                    .call(brush)
                                    .selectAll(".selection")
                                    .attr("fill","none")
                                    .attr("stroke-width",3).attr("stroke","rgb(70,130,180)")//.attr("stroke-dasharray","3,6")
                                    .attr("y", -16)
                                    .attr("height", 82) //height2 + 10)
                                    .attr("rx",3)
                                    .attr("ry",3)
                                    .attr("opacity",1)
                                    .attr("shape-rendering","geometricPrecision")
                                    .style("filter", "url(#drop-shadow)");
                    //var nis = sin.map(width2Scale, width2Scale.invert);
                    setBrush(sin[0],sin[1]);

                }else{
                    width2Scale.range([0,width_new]);
                    widthScale.range([0,width_new]);
                }
                clip_rect.attr("width",width_new);
                context.selectAll(".fullSeqW").attr("x2",width_new);
                d3.selectAll(".zoom").attr("width",width_new);
        }

        function hideLabels(){
            //console.log("hiding labels");
            focus.selectAll(".peak_label").attr("visibility","hidden");
            focus.selectAll(".var_noise_indic").attr("visibility","visible");
        }

        function redrawLines(){
            //widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            //widthScale.domain(brush.extent());
            focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
            focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");
            focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
            focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
            focus.selectAll(".var_noise_indic").attr("stroke-width",widthScale(12)-widthScale(0));
            focus.selectAll(".var_noise_indic").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                    .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            updateFilt();
            // focus.selectAll(".ref").attr("x",function(d){return widthScale(d["trace_peak"]);});
            // focus.selectAll(".user").attr("x",function(d){return widthScale(d["trace_peak"]);});
            // focus.selectAll(".peak_label").attr("x",function(d){return widthScale(d["trace_peak"]);});
        }

        function redraw()  {

            var wsd = widthScale.domain();
            var wsr = widthScale.range();
            var bp = (wsd[1]-wsd[0])/(wsr[1]*11);  //base per
            var bpt = bp*1000;

            console.log('bpt',bpt);

            focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
            focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
            focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
            focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");
        //     focus.selectAll(".scope").attr("x",function(d){return widthScale(d["trace_peak"])-12;});
            focus.selectAll(".peak_label").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".aa_ref_line").attr("x1",function(d){if(d["ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                                                else if(d["ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                           .attr("x2",function(d){if(d["ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                                                else if(d["ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                           .attr("stroke-width",function(d){if(d["ord_in_cod"]==0){return 0}
                                                                else{return (widthScale(10)-widthScale(0))}});
            focus.selectAll(".aa_sample_line").attr("x1",function(d){if(d["sample_ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                                                else if(d["sample_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["sample_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                              .attr("x2",function(d){if(d["sample_ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                                                else if(d["sample_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["sample_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                              .attr("stroke-width",function(d){if(d["sample_ord_in_cod"]==0||  !d["sample_ord_in_cod"]){return 0}
                                                                else{return (widthScale(10)-widthScale(0))}});
            focus.selectAll(".aa_mut_line").attr("x1",function(d){if(d["mut_ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                                                else if(d["mut_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["mut_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                           .attr("x2",function(d){if(d["mut_ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                                                else if(d["mut_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                                                else if(d["mut_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                                                else {return widthScale(d["trace_peak"]);}})
                                           .attr("stroke-width",function(d){if(d["mut_ord_in_cod"]==0||!d["mut_ord_in_cod"]){return 0}
                                                                else{return (widthScale(10)-widthScale(0))}});
            focus.selectAll(".qual_fwd").attr("x",function(d){return widthScale(d["trace_peak"])-9;});
            focus.selectAll(".qual_rev").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".line").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                   .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".var_noise_indic").attr("stroke-width",widthScale(12)-widthScale(0));
        //     //conditional visibility
            //console.log(w);
            if(bpt<22){
                focus.selectAll(".peak_label").attr("visibility","visible");
                if(bp===0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                focus.selectAll(".peak_label").attr("visibility","hidden");
                focus.selectAll(".ref").attr("visibility","visible");
                focus.selectAll(".user").attr("visibility","visible");
                focus.selectAll(".var_noise_indic").attr("visibility","visible");

                if(bpt<70){
                  focus.selectAll(".qual_fwd").attr("visibility","visible");
                  focus.selectAll(".qual_rev").attr("visibility","visible");
                  focus.selectAll(".aa").attr("visibility","visible");
                  focus.selectAll(".coding_ten").attr("visibility","visible");
                  //focus.selectAll(".aa").attr("visibility","visible");
                }
                //if(bpt<25){focus.selectAll(".short").attr("visibility","visible");}
            }
            //console.log(show_qual); //setting quality bars on/off
            focus.selectAll(".q").attr("visibility",show_qual);

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
                focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
                focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
                focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
                focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");

            }else if(j==="TRUE"){
                split_peak_offset = 100;
                heightScale_fwd = heightScale;
                heightScale_rev = heightScale;
                focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
                focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
                focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
                focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");
            }

        }

        function updateLine(data,base,rev){
            switch(base) {
                case "A": if(rev){var g1 = gLine_a_r; var g2 = gLine_a_r_gray;}
                          else   {var g1 = gLine_a;   var g2 = gLine_a_gray;  } var col = "#00A100"; break;
                case "C": if(rev){var g1 = gLine_c_r; var g2 = gLine_c_r_gray;}
                          else   {var g1 = gLine_c;   var g2 = gLine_c_gray;  } var col = "#2985EA"; break;
                case "G": if(rev){var g1 = gLine_g_r; var g2 = gLine_g_r_gray;}
                          else{var g1 = gLine_g;      var g2 = gLine_g_gray;  } var col = "#6C6A6C"; break;
                case "T": if(rev){var g1 = gLine_t_r; var g2 = gLine_t_r_gray;}
                          else{var g1 = gLine_t;      var g2 = gLine_t_gray;  } var col = "#EA2929"; break;
            }
            if(rev){var c = "line_r";var l = line_rev; cl = "url(#rect_rev_clip)"}
            else{var c = "line_f";   var l = line_fwd; cl = "url(#rect_fwd_clip)"}

            var line = g1.selectAll("path").data(data); //UPDATE
            line.exit().remove();                      //EXIT
            line.enter().append("path")                //enter
                .merge(line).attr("class","path "+c)
                .attr("d",l)
                .attr("clip-path","url(#clip)")
                .attr("fill","none").attr("stroke",col)
                .attr("stroke-width",1.5);         // on reverse attr("stroke-dasharray","20,3,10,1,10,1");
            var line = g2.selectAll("path").data(data); //UPDATE
            line.exit().remove();                      //EXIT
            line.enter().append("path")                //enter
                .merge(line).attr("class","path "+c)
                .attr("d",l)
                .attr("clip-path", cl)
                .attr("filter","url(#monochrome)")
                .attr("fill","none").attr("stroke",col)
                .attr("stroke-width",2);         // on reverse attr("stroke-dasharray","20,3,10,1,10,1");
        }
        function updateFilt(fwd_x,rev_x,single_rev){
            if(fwd_x != undefined){filt_fwd_x = fwd_x;}
            if(rev_x != undefined){filt_rev_x = rev_x;}
            if(single_rev!=undefined){filt_single_rev = single_rev;}
            var fwd_width = 0;
            if(widthScale(filt_fwd_x)>0){fwd_width = widthScale(filt_fwd_x)};
            clip_fwd.attr("x",0).attr("width",fwd_width);
            filt_fwd_mini.attr("x2",width2Scale(filt_fwd_x));
            if(fwd_width!=0){
                //console.log("drawing line");
                filt_line_fwd.attr("x1",widthScale(filt_fwd_x))
                             .attr("x2",widthScale(filt_fwd_x))
                             .attr("y1",200)
                             .attr("y2",340);
            }else{
                filt_line_fwd.attr("x1",0)
                         .attr("x2",0)
                         .attr("y1",0)
                         .attr("y2",0);
            }
            // console.log("filt_rev_x ", filt_rev_x);
            // console.log("widthScale(filt_rev_x) ", widthScale(filt_rev_x));
            // console.log("width",width);
            // console.log("width - widthScale(filt_rev_x)",width - widthScale(filt_rev_x));
            if(filt_single_rev){
                if(filt_rev_x == 0  ){
                    filt_rev_mini.attr("x1",0).attr("x2",0);
                }else{
                    filt_rev_mini.attr("x1",width2Scale(filt_rev_x)).attr("x2",width);
                }
                //console.log("single_rev true");
                //console.log("filt_rev_x ",filt_rev_x);
                if(width - widthScale(filt_rev_x)<0){
                    clip_fwd.attr("x",0).attr("width",0);
                }else{
                    //console.log("filt_rev_x",filt_rev_x);
                    //console.log("width filt_rev_x",(width - widthScale(filt_rev_x)));
                    clip_fwd.attr("x",widthScale(filt_rev_x)).attr("width",(width - widthScale(filt_rev_x)));
                }
            }else{
                if(filt_rev_x == 0  ){
                    clip_rev.attr("x",0).attr("width",0);
                    filt_rev_mini.attr("x1",0).attr("x2",0);
                    filt_line_rev.attr("x1",0)
                                 .attr("x2",0)
                                 .attr("y1",0)
                                 .attr("y2",0);

                }else{
                    if(width - widthScale(filt_rev_x)<0){
                        clip_rev.attr("x",0).attr("width",0);
                    }else{
                        clip_rev.attr("x",widthScale(filt_rev_x)).attr("width",(width - widthScale(filt_rev_x)));
                    }
                    filt_rev_mini.attr("x1",width2Scale(filt_rev_x)).attr("x2",width);

                    filt_line_rev.attr("x1",(widthScale(filt_rev_x)-1))
                                     .attr("x2",(widthScale(filt_rev_x)-1))
                                     .attr("y1",350)
                                     .attr("y2",490);


                }
            }
        }
        function setNoiseArea(fwd,rev){
            var gnf = gNoise_fwd.selectAll("path").data([fwd]);
            gnf.enter().append("path")
                .merge(gnf).attr("class","area area_fwd").attr("d",noise_area_fwd)
               .attr("fill","#000000").attr("stroke","none").attr("opacity",0.15).attr("clip-path","url(#clip)");
            gnf.exit().remove();
            if(typeof rev !== 'undefined'){
                var gnr = gNoise_rev.selectAll("path").data([rev]);
                gnr.enter().append("path")
                   .merge(gnr).attr("class","area area_rev").attr("d",noise_area_rev)
                   .attr("fill","#440000").attr("stroke","none").attr("opacity",0.15).attr("clip-path","url(#clip)");
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
            text.exit().remove();                               //EXIT
            text.enter().append("text")                        //Enter
                .merge(text).attr("class",function(d){
                    if(label.indexOf("user") > -1) {return "peak_label short user ".concat("id").concat(d["id"]);}
                    return("peak_label short "+c);
                })
                .text(function(d){
                    if(label.indexOf("mut") > -1){
                        if(typeof(d[label])=='undefined'){ console.log(label);console.log(d);}
                        return d[label].toLowerCase();
                    } else {
                        if(label=="reference" & d[label]!="NA"){
                                if(d["exon_intron"].indexOf("exon")>-1){
                                    return d[label];
                                }
                        }
                        return d[label].toLowerCase();}
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
                    if      (d[label] === "A"){ return "#00A100"; }
                    else if (d[label] === "C"){ return "#2985EA"; }
                    else if (d[label] === "G"){ return "#6C6A6C"; }
                    else if (d[label] === "T"){ return "#EA2929"; }
                    else if (d[label] === "-"){ return "black"; }
                    else    {                   return "orange"; }});
        }

        function showVarInMinimap(choices){
            //reference
            //console.log("showVarInMinimap");
            var r = varim_ref.selectAll("line").data(choices);
            r.exit().remove();
            r.enter().append("line").attr("class","enter")
                .merge(r).attr("class","minimap context")
      			.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      			.attr("y1",4)
      			.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      			.attr("y2",24)
      			.attr("stroke-width",3)
      			.attr("stroke",function(d) { return get_color(d["reference"])});
  			    // user
            var us = varim_user_s.selectAll("line").data(choices);
            us.exit().remove();
            us.enter().append("line").attr("class","enter")
                .merge(us).attr("class","minimap context")
  		        .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y1",26)
  				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y2",38)
  				.attr("stroke-width",3)
  				.attr("stroke",function(d) { return get_color(d["user_sample"])});
            var um = varim_user_m.selectAll("line").data(choices);
            //um.attr("class","update");
            um.exit().remove();
            um.enter().append("line").attr("class","enter")
                .merge(um).attr("class","minimap context")
  				.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y1",38)
  				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
  				.attr("y2",46)
  				.attr("stroke-width",3)
  				.attr("stroke",function(d) { return get_color(d["user_mut"])});


            //console.log(choices);
            var v = var_ind_focus.selectAll("line").data(choices);
            v.exit().remove();
            //strandness in choices 0 = undetermined (should not occure); 1 = forward strand; 2 = reverse strand; 3 = both, 4 = undetermined indel, 5 = indel forward only, 6 = indel reverse only, 7 = indel both strands
            v.enter().append("line").attr("class","enter")
                .merge(v).attr("class","peak_label short line var_noise_indic")
  		        .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  		        .attr("y1",function(d){if(d["strand"]==2){return 280;}else{return 140;}})
  		        .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  		        .attr("y2",function(d){if(d["strand"]==1){return 270;}else{return 420;}})
  		        .attr("stroke-width",widthScale(12)-widthScale(0))
                .attr("stroke","rgba(255,0,0,0.12)").attr("clip-path","url(#clip)");

        }

        // function showNoisyNbr(noisy_neighbors){
        //     var nnc = noisyn_con.selectAll("line").data(noisy_neighbors);
        //     nnc.enter().append("line");
        //     nnc.attr("class","minimap context")
     //  			    .attr("x1",function(d){return width2Scale(d["trace_peak"]);})
     //  			    .attr("y1",-3)
     //  			    .attr("x2",function(d){return width2Scale(d["trace_peak"]);})
     //  			    .attr("y2",3)
     //  			    .attr("stroke-width",3)
     //  			    .attr("stroke", "brown");
        //     nnc.exit().remove();
        // /*  var nnf = noisyn_foc.selectAll("line").data(noisy_neighbors);
        //     nnf.enter().append("line");
        //     nnf.attr("class","peak_label short line var_noise_indic")
  // 		         .attr("x1",function(d){return widthScale(d["trace_peak"]);})
  // 		         .attr("y1",140)
  // 		         .attr("x2",function(d){return widthScale(d["trace_peak"]);})
  // 		         .attr("y2",400)
  // 		         .attr("stroke-width",10).attr("stroke","brown").attr("opacity",0.2).attr("stroke-dasharray","1,3");
        //     nnf.exit().remove();
        // */
        // }

        function setCodingLabel(calls){
            //console.log(calls);
            var aa_r = aa_ref.selectAll("text").data(calls);
            aa_r.exit().remove();
            aa_r.enter().append("text")
                .merge(aa_r).attr("class","peak_label short aa")
                .text(function(d){
//                    if   (d["ord_in_cod"] == 1) {return d["aa_ref"].toUpperCase()+""+d["codon"];}
                    //if   (d["ord_in_cod"] == 1) {return d["aa_ref"]+""+d["codon"];}
                    if   (d["ord_in_cod"] == 2) {return d["aa_ref"]+d["codon"];}
                    else {                       return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");

            var aa_s = aa_sample.selectAll("text").data(calls);
            aa_s.enter().append("text") //aa_sample
                .merge(aa_s).attr("class","peak_label short aa_sample aa")
                .text(function(d){
                    if   (d["sample_ord_in_cod"] == 2) {return d["aa_sample"];}
                    else {                       return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa_sample"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            aa_s.exit().remove();
            var aa_m = aa_mut.selectAll("text").data(calls);
            aa_m.enter().append("text")
                .merge(aa_m).attr("class","peak_label short aa_mut aa")
                .text(function(d){
                    if   (d["mut_ord_in_cod"] == 2) {return d["aa_mut"];}
                    else {                       return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){callShiny(d["id"]);})
                .attr("y",label_pos["aa_mut"])
                .attr("fill", "black").attr("opacity", 0.6).attr("font-family", "sans-serif").attr("font-size", "10px")
                .attr("stroke","#000000");
            aa_m.exit().remove();
            //console.log(calls);
            var aa_r_l = aa_ref_line.selectAll("line").data(calls);
            aa_r_l.exit().remove();
            //strandness in choices 0 = undetermined (should not occure); 1 = forward strand; 2 = reverse strand; 3 = both, 4 = undetermined indel, 5 = indel forward only, 6 = indel reverse only, 7 = indel both strands
            aa_r_l.enter().append("line").attr("class","enter")
                .merge(aa_r_l).attr("class","peak_label short aa_ref_line aa")
  		        .attr("y1",function(d){return label_pos["aa"]-25}) //-11
  		        .attr("y2",function(d){return label_pos["aa"]-10}) //-9
  		        .attr("x1",function(d){
                                if(d["ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                else if(d["ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                else if(d["ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                else {return widthScale(d["trace_peak"]);}
                })
                .attr("x2",function(d){
                                if(d["ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                else if(d["ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                else if(d["ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                else {return widthScale(d["trace_peak"]);}
                })
                .attr("stroke-width",function(d){
                                if(d["ord_in_cod"]==0){return 0}
                                else{return (widthScale(10)-widthScale(0))}
                })
                .attr("stroke","rgba(0,0,0,0.07)");

                var aa_s_l = aa_sample_line.selectAll("line").data(calls);
                aa_s_l.exit().remove();
                //strandness in choices 0 = undetermined (should not occure); 1 = forward strand; 2 = reverse strand; 3 = both, 4 = undetermined indel, 5 = indel forward only, 6 = indel reverse only, 7 = indel both strands
                aa_s_l.enter().append("line").attr("class","enter")
                    .merge(aa_s_l).attr("class","peak_label short aa_sample_line aa")
      		        .attr("y1",function(d){return label_pos["aa_sample"]+3}) //+2
      		        .attr("y2",function(d){return label_pos["aa_sample"]+15}) //+4
      		        .attr("x1",function(d){
                                    if(d["sample_ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                    else if(d["sample_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                    else if(d["sample_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                    else {return widthScale(d["trace_peak"]);}
                    })
                    .attr("x2",function(d){
                                    if(d["sample_ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                    else if(d["sample_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                    else if(d["sample_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                    else {return widthScale(d["trace_peak"]);}
                    })
                    .attr("stroke-width",function(d){
                                    if(d["sample_ord_in_cod"]==0|| !d["sample_ord_in_cod"]){return 0}
                                    else{return (widthScale(10)-widthScale(0))}
                    })
                    .attr("stroke","rgba(0,0,0,0.07)");
                var aa_m_l = aa_mut_line.selectAll("line").data(calls);
                aa_m_l.exit().remove();
                    //strandness in choices 0 = undetermined (should not occure); 1 = forward strand; 2 = reverse strand; 3 = both, 4 = undetermined indel, 5 = indel forward only, 6 = indel reverse only, 7 = indel both strands
                aa_m_l.enter().append("line").attr("class","enter")
                    .merge(aa_m_l).attr("class","peak_label short aa_mut_line aa")
          		    .attr("y1",function(d){return label_pos["aa_mut"]-23}) //-9
          		    .attr("y2",function(d){return label_pos["aa_mut"]-10}) //-11
          		    .attr("x1",function(d){
                                    if(d["mut_ord_in_cod"]==1){return (widthScale(d["trace_peak"] +2));}
                                    else if(d["mut_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                    else if(d["mut_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                    else {return widthScale(d["trace_peak"]);}
                    })
                    .attr("x2",function(d){
                                    if(d["mut_ord_in_cod"]==1){return (widthScale(d["trace_peak"]+2))}
                                    else if(d["mut_ord_in_cod"]==3){return (widthScale(d["trace_peak"]-2))}
                                    else if(d["mut_ord_in_cod"]==4){return (widthScale(d["trace_peak"]-4))}
                                    else {return widthScale(d["trace_peak"]);}
                    })
                    .attr("stroke-width",function(d){
                                    if(d["mut_ord_in_cod"]==0 || !d["mut_ord_in_cod"]){return(0);}
                                    else{return (widthScale(10)-widthScale(0))}
                    })
                    .attr("stroke","rgba(0,0,0,0.07)");

        }

        function setQualityLabels(calls,rev){
            q = rev ? "quality_fwd" : "quality";
            var qt = quals_txt.selectAll("text").data(calls);
            qt.enter().append("text")
                    .merge(qt).attr("class","peak_label q")
      		        .text(function(d){return d["quality"];})
      		        .attr("text-anchor", "middle")
      		        .attr("x",function(d){return widthScale(d["trace_peak"]);})
                    .on("click",function(d,i){callShiny(d["id"]);})
      		        .attr("y",label_pos["qual"] -2)
      		        .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "10px");
            qt.exit().remove();
            var qg = quals_g.selectAll("rect").data(calls);
            qg.enter().append("rect")
            .merge(qg).attr("class","peak_label qual_fwd q")
      		        .attr("x",function(d){return (widthScale(d["trace_peak"])-9);})
      		        .attr("y",label_pos["qual"]).attr("rx",1).attr("ry",1)
      		        .attr("width",9+(!rev)*9)
      		        .attr("height",function(d){return d[q];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
            qg.exit().remove();
            if(rev!=0){
                var qg_r = quals_g_r.selectAll("rect").data(calls);
                qg_r.enter().append("rect")
                .merge(qg_r).attr("class","peak_label qual_rev q")
      		        //.attr("x",function(d){return (widthScale(d["trace_peak"]) + 900);})
      		        .attr("y",label_pos["qual"]).attr("rx",1).attr("ry",1)
      		        .attr("width",9)
      		        .attr("height",function(d){return d["quality_rev"];})
      		        .attr("fill", "rgba(200,200,200,0.3)");
                qg_r.exit().remove();
            }
        }

        function setIntrexBoxes(intrex){
            //console.log("setIntrexBoxes");
            var eb = exon_boxes.selectAll("rect").data(intrex);
            eb.exit().remove();
            eb.enter().append("rect")
                .merge(eb).attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",0).attr("rx",3).attr("ry",3)
  				.attr("width",function(d){return width2Scale(d["end"]-d["start"]);})
  				.attr("height",50)
  				.attr("fill",function(d) {
  				           if (d["splicevar"] != ''){   return "rgba(255,205,205,1.0)";
  				    } else if (/exon/.test(d["attr"])){ return "rgba(200,200,200,1.0)";
  				    } else {                            return "rgba(230,230,230,1.0)"; }
  				});
            var iet = intrex_txt.selectAll("text").data(intrex);
            iet.exit().remove();
            iet.enter().append("text")
                .merge(iet).attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",-5)
  				.attr("opacity",0.8)
  				.attr("fill","black")
  				.text(function(d) {
  				    return d["attr"]+d["splicevar"];
  				});
            var ien = intrex_num.selectAll("text").data(intrex);
            ien.exit().remove();
            ien.enter().append("text")
  		        .merge(ien).attr("class","context")
  				.attr("x",function(d){return width2Scale(d["start"]);})
  				.attr("y",64)
  				.attr("opacity",0.8)
  				.text(function(d){return d["id"];})
  				.attr("fill","black");
        }

        Shiny.addCustomMessageHandler("goto",
            function(message) {
                /*console.log("goto message");*/
                setBrush(Number(message)-180,Number(message)+200);
                //focus.selectAll(".scope").attr("opacity",0);
                //focus.selectAll(".".concat("scope").concat(message))
                //.transition().attr("opacity", 1);
            }
        );

        // Shiny.addCustomMessageHandler("input_change",
        //     function(message){
        //         if(message<brush.extent()[0] || message>brush.extent()[1]){
        //            setBrush(Number(message)-100,Number(message)+120);
        //         }
        //         focus.selectAll(".scope").attr("opacity",0);
        //         focus.selectAll(".".concat("scope").concat(message))
        //                .transition().attr("opacity", 1);
        //     }
        // );

        Shiny.addCustomMessageHandler("s2n_min",
            function(message) {
                mult = Number(message);
                //console.log(mult);
                redraw();
            }
        );

        // Shiny.addCustomMessageHandler("show",
        //     function(message){
        //         //console.log(message);
        //         if(message==="TRUE"){
        //             focus.selectAll(".call").attr("opacity",0.8);
        //             redraw();
        //         }else if(message==="FALSE"){
        //             focus.selectAll(".call").attr("opacity",0);
        //             redraw();
        //         }
        //     }
        // );

        Shiny.addCustomMessageHandler("join",
            function(message){
                join = message;
                joinView(message);
                }
        );

//         Shiny.addCustomMessageHandler("opac_f",
//             function(message){
//                 //console.log(message);
//                 focus.selectAll(".line_f").attr("opacity",Number(message));
//                 focus.selectAll("g").selectAll(".area_fwd").attr("opacity",Number(message)/5);
//                 focus.selectAll(".qual_fwd").attr("opacity",Number(message)*0.8);
// //                focus.selectAll(".fwd").attr("opacity",Number(message));
//             }
//         );

//         Shiny.addCustomMessageHandler("opac_r",
//             function(message){
//                 //console.log(message);
//                 focus.selectAll(".line_r").attr("opacity",Number(message));
//                 focus.selectAll("g").selectAll(".area_rev").attr("opacity",Number(message)/5);
//                 focus.selectAll(".qual_rev").attr("opacity",Number(message)*0.8);
// //                focus.selectAll(".rev").attr("opacity",Number(message));
//             }
//         );

        function callShiny(id,trace_peak){
            //console.log(id);
            Shiny.onInputChange("pos_click", {id: id});
        };

        function get_color(col){
            if      (col === "A"){ return "#00A100"; }
            else if (col === "C"){ return "#2985EA"; }
  			else if (col === "G"){ return "#6C6A6C"; }
  			else if (col === "T"){ return "#EA2929"; }
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
            instanceCounter: instanceCounter,
            width:   w,
            height:  h,
	        height2: height2,
            call_opacity: call_opacity,
            label_pos:  label_pos,
            widthScale:  widthScale,
            width2Scale: width2Scale,
            heightScale: heightScale,
            heightScale_fwd_split: heightScale_fwd_split,
            heightScale_rev_split: heightScale_rev_split,
            context: context,
            focus:   focus,
            full_line: full_line,
            updateFilt: updateFilt,
            rescaleWidth: rescaleWidth,
            // scope_g: scope_g,
            //
	        brush:   brush,
            redraw:  redraw,
            // finish_fwBrushInit: finish_fwBrushInit,
            // finish_rvBrushInit: finish_rvBrushInit,
	        join:    join,
            show_qual: show_qual,
            joinView:joinView,
            //
            setNoiseArea: setNoiseArea,
            updateLine: updateLine,
            reHeight: reHeight,
            setBrush: setBrush,
            setPeakLabel: setPeakLabel,
            showVarInMinimap: showVarInMinimap,
            // showNoisyNbr: showNoisyNbr,
            setCodingLabel: setCodingLabel,
            setQualityLabels: setQualityLabels,
            setIntrexBoxes: setIntrexBoxes
        }
    },

    resize: function(el, width, height, instance) {


        d3.select(el).selectAll("svg")
          .attr("width", width);
        instance.rescaleWidth(width);

        if (instance.lastValue) {
            instance.lastValue.resize = true;
            this.renderValue(el, instance.lastValue, instance);
        }
    },
    //function called everytime input parameters change
    renderValue: function(el, x, instance) {

        instance.lastValue = x;
        //the render function behaves differently when called repeatedly
        //the first run is actually still a part of the initialization step
        if(x.new_sample){
            console.log("new sample");
            x.new_sample = false;
  		    var intens = x["intens"];
            var intens_rev = "";
            var rev = 0;  //offset on labels in case we have alternative reference
            if(x["intens_rev"] !== null){
                intens_rev = x["intens_rev"];
                rev = 1;
            }
            var single_rev = x["single_rev"];
            //console.log("single_rev:"+single_rev);
            if(instance.instanceCounter>=1){ //cleanup after previous sample
                instance.focus.selectAll(".area").remove();
                instance.focus.selectAll(".line_f").remove();
                instance.focus.selectAll(".line_r").remove();
                instance.focus.selectAll(".peak_label").remove();
                //instance.focus.selectAll(".zoom").remove();
                instance.context.selectAll(".context").remove();
//                 //instance.focus.selectAll(".excl").remove();
                instance.updateFilt(0,0,single_rev);
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
//
  			widthScale.domain([0,domain_x]);
  			width2Scale.domain([0,domain_x]);
  			heightScale.domain([0,domain_y]);
  	        instance.heightScale_fwd_split.domain([0,domain_y]);
            instance.heightScale_rev_split.domain([0,domain_y]);
//
//   			//visualise fullseq width
  			instance.full_line
  				.append("line").attr("class","context fullSeqW")
  				.attr("x1",0)
  				.attr("y1",25)
  				.attr("x2",widthScale(domain_x))
  				.attr("y2",25)
  				.attr("stroke-width",3).attr("stroke","rgba(180,180,180,1.0)");

//             //intron/exon boxes
            instance.setIntrexBoxes(intrex);
  		    //brush.x(width2Scale);
//             var brush_fw_height = 140;
//             if(intens_rev == ""){
//                 brush_fw_height = 280;
//             }
//             if(single_rev==true){
//                 brush_fw_height = 0;
//             }
//             instance.finish_fwBrushInit(x["brush_fw"],rev,brush_fw_height);
//             //console.log("brush_rv: ", x["brush_rv"])
//             if(rev!=0){
//                 instance.finish_rvBrushInit(x["brush_rv"],domain_x,rev);
//             }else if(single_rev){
//                 instance.finish_rvBrushInit(x["brush_fw"],domain_x,rev);
//
//             }
//             else{
//                 instance.finish_rvBrushInit(0,0);
//             }
            context.append("g")
//				.attr("class", "x brush")
                .attr("class","brush context brush_context")
                .call(brush)
                //.call(brush.move)
      	 		.selectAll(".selection")
                //.attr("class","main_brush selection")
                //.attr("fill","rgba(70,130,180,0.05)") //steal blue
                .attr("fill","none")
                .attr("stroke-width",3).attr("stroke","rgb(70,130,180)")//.attr("stroke","#ff4444")//.attr("stroke-dasharray","3,6")
  		 		.attr("y", -16)
                .attr("height", 82) //height2 + 10)
                .attr("rx",3)
  		 		.attr("ry",3)
                .attr("shape-rendering","geometricPrecision")
                .style("filter", "url(#drop-shadow)")
  		 		.attr("opacity",1);

            instance.updateLine([intens["A"]],"A",false);
            instance.updateLine([intens["C"]],"C",false);
            instance.updateLine([intens["G"]],"G",false);
            instance.updateLine([intens["T"]],"T",false);
            if(single_rev){
                instance.updateFilt(0,x["brush_fw"],single_rev);
            }else{
                instance.updateFilt(x["brush_fw"]);
            }

            //noise indicator
            var a_noise_fwd = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
            //reverse strand
            if(intens_rev != ""){
                var a_noise_rev = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_rev"]]);
                instance.updateLine([intens_rev["A"]],"A",true);
                instance.updateLine([intens_rev["C"]],"C",true);
                instance.updateLine([intens_rev["G"]],"G",true);
                instance.updateLine([intens_rev["T"]],"T",true);
                instance.updateFilt(x["brush_fw"],x["brush_rv"]);
            }
            instance.setNoiseArea(a_noise_fwd,a_noise_rev);
            //on single strand always show "join view"
            if(intens_rev != ""){
                instance.joinView();
            }else{
                instance.joinView("TRUE");
            }
//
//             instance.scope_g.selectAll("scope").data(calls).enter()        //scope (position indicator)
//                  .append("rect").attr("class",function(d){return "scope ".concat("scope").concat(d["trace_peak"]);})
//                  .attr("x",function(d){return (widthScale(d["trace_peak"])-12)})
//                  .attr("y",20).attr("rx",2).attr("ry",2)
//                  .attr("width",24).attr("height",instance.height-90)
//                  .attr("fill","rgba(155, 155, 255, 0.12)").attr("opacity",0);
            focus.append("g").selectAll("text.seq.codon").data(calls).enter() //codon stuff
                .append("text").attr("class",function(d){if(d["coding_seq"]%10==0){return "peak_label coding_ten"}else{return "peak_label";}})
                .text(function(d){
//                    if   (d["coding_seq"] > 0){return d["coding_seq"]+" : "+d["codon"]+"."+d["ord_in_cod"];}
                    if   (d["coding_seq"] > 0){return "c." + d["coding_seq"];}
                    else {                     return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .attr("y",(instance.label_pos["codon"]+rev))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
                .append("text").attr("class",function(d){if((d["coding_seq"]%10==0)&&(d["coding_seq"]>0)){return "peak_label coding_ten"}else{return "peak_label"}})
                .text(function(d){return "g." + d["gen_coord"];})
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
//             //console.log(x);
//             //console.log(x["qual_present"]);
            if(x["qual_present"]){
                instance.setQualityLabels(calls,rev);
            }
            var show_qual  = x["show_qual"];
            if(show_qual){
                instance.show_qual = "visible";
                d3.selectAll(".q").attr("visibility","visible");
            }else{
                instance.show_qual = "hidden";
                d3.selectAll(".q").attr("visibility","hidden");
            }
//             //default
            focus.selectAll(".call").attr("opacity",instance.call_opacity);
            instance.setPeakLabel(calls,"user_sample");
            instance.setPeakLabel(calls,"user_mut");
            instance.setCodingLabel(calls);
            //d3.selectAll(".peak_label").filter(function(d){return d["trace_peak"] < x["brush_fw"]}).attr("stroke","black");
// /*
//             //horizontal line on top of the chrom
//                 focus
//       				.append("line")
//       				.attr("x1",0)
//       				.attr("y1",intens_guide_line)
//       				.attr("x2",1200)
//       				.attr("y2",intens_guide_line)
//       				.attr("stroke-width",1).attr("stroke","rgba(0,0,0,0.6)").attr("stroke-dasharray",2);
// */
//
            instance.showVarInMinimap(choices);
//             instance.showNoisyNbr(noisy_neighbors);
// 	        //zooming in so that the first view is not ugly dense graph
//
            if (typeof choices[0] !== 'undefined') {
                from = choices[0]["trace_peak"]-200;
                to   = choices[0]["trace_peak"]+220;
                if(from < 0) {from = 0;to = 420}
                instance.setBrush(from,to);
            }else{
                instance.setBrush(200,1200);
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
            //     var noisy_neighbors = HTMLWidgets.dataframeToD3(x["noisy_neighbors"]);
                var show_calls  = x["show_calls"];
                if(show_calls){ instance.call_opacity = 0.8; }
                else{ instance.call_opacity = 0; }
                instance.focus.selectAll(".call").attr("opacity",instance.call_opacity);
                instance.setIntrexBoxes(HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"]));
                instance.showVarInMinimap(choices);
                instance.choices = x.choices;
            //     instance.showNoisyNbr(noisy_neighbors);
                instance.setPeakLabel(calls,"user_sample");
                instance.setPeakLabel(calls,"user_mut");
                instance.setPeakLabel(calls,"reference")
                instance.setPeakLabel(calls,"mut_call_fwd",instance.call_opacity);
                if(rev){
                    instance.setPeakLabel(calls,"mut_call_rev",instance.call_opacity);
                }
                instance.setCodingLabel(calls);
                var show_qual  = x["show_qual"];
                if(show_qual){
                    instance.show_qual = "visible";
                    d3.selectAll(".q").attr("visibility","visible");
                }else{
                    instance.show_qual = "hidden";
                    d3.selectAll(".q").attr("visibility","hidden");
                }
            //
  			}else if(x.resize == true){
                instance.lastValue.resize = false;
                instance.redraw();
                instance.setIntrexBoxes(HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"]));
                instance.showVarInMinimap(HTMLWidgets.dataframeToD3(x["choices"]));
            //
            }else { console.log(x) }
        }
    }
});
