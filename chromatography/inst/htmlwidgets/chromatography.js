HTMLWidgets.widget({

    name: 'chromatography',
    type: 'output',

    initialize: function(el, w, h) {

        console.log("height:");
        console.log(h);

        var instanceCounter = 0;
        var intrex = "";
        var max_x = 0;
        var max_y = 0;
        var margin  = {top: 65, right: 10,bottom: 20, left:10},
            margin2 = {top: 20, right: 10,bottom: 420,  left:10},
            width   = w - margin.left - margin.right,
            height  = h - margin.top  - margin.bottom,
            half_height = height/1.6;
        var height2 = h - margin2.top - margin2.bottom;
        var widthScale   = d3.scaleLinear().range([0,width]),
            width2Scale  = d3.scaleLinear().range([0,width]),  //remains constant, to be used with context
            heightScale  = d3.scaleLinear().range([height,0]),
            height2Scale = d3.scaleLinear().range([height2,0]),
            heightScale_fwd_split = d3.scaleLinear().range([half_height,(2*half_height -  height)]),
            heightScale_rev_split = d3.scaleLinear().range([height,half_height]),
            heightScale_fwd = heightScale_fwd_split,
            heightScale_rev = heightScale_rev_split;


        var call_opacity = 0
        var mult = 2;
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

        // should be 'drag' not 'zoom'
        var zoom = d3.zoom()
                     .scaleExtent([1, Infinity])
                     .translateExtent([[0, 0], [width, zoom_height]])
                     .extent([[0, 0], [width, zoom_height]])
                     .on("start",zoom_start)
                     .on("zoom", zoomed)
                     .on("end",zoom_end);

        var filt_fwd_start = 0;
        var filt_rev_start = 0;
        var filt_fwd_end = 0;
        var filt_rev_end = 0;
        var filt_has_rev = false;
        var filt_fwd_start_mini = svg.append("line")
                            .attr("x1", 10).attr("y1", 18)
                            .attr("x2", 250).attr("y2", 18)
                            .attr("clip-path","url(#clip)")
                            .attr("stroke-width", 2)
                            .style("stroke-dasharray", ("3, 3"))
                            .attr("stroke", "red").attr("opacity",1);
        var filt_fwd_end_mini = svg.append("line")
                            .attr("y1", 18).attr("y2", 18)
                            //.attr("x1", 10).attr("x2", 250) // don't know before init
                            .attr("clip-path","url(#clip)")
                            .attr("stroke-width", 2)
                            .style("stroke-dasharray", ("3, 3"))
                            .attr("stroke", "red").attr("opacity",1);

        var filt_rev_start_mini = svg.append("line")
                            .attr("x1", 10).attr("y1", 72)
                            .attr("x2", 250).attr("y2", 72)
                            .attr("clip-path","url(#clip)")
                            .attr("stroke-width", 2)
                            .style("stroke-dasharray", ("3, 3"))
                            .attr("stroke", "red").attr("opacity",1);

        var filt_rev_end_mini = svg.append("line")
                            .attr("y1", 72).attr("y2", 72)
                            //.attr("x1", 10).attr("x2", 250) // don't know before init
                            .attr("clip-path","url(#clip)")
                            .attr("stroke-width", 2)
                            .style("stroke-dasharray", ("3, 3"))
                            .attr("stroke", "red").attr("opacity",1);

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
        var selected_pos = focus.append("rect")
                                        .attr("class","selected_pos")
                                        .attr("x",10)
                                        .attr("y", 25 + 113)
                                        .attr("width",20)
                                        .attr("height",280)
                                        .attr("stroke-width", 0)
                                        .attr("fill","lightblue")
                                        .attr("opacity",0.5);
        var selected_pos_x = -30;

        var quals_g    = focus.append("g");
        var quals_g_r  = focus.append("g");
        var quals_txt  = focus.append("g");
        var gLines     = focus.append("g");

        const lines    = new LineSet(gLines,height,widthScale,h);


        var aa_ref_line       = focus.append("g");
        var aa_sample_line    = focus.append("g");
        var aa_mut_line       = focus.append("g");

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

        var aa_ref            = focus.append("g");
        var aa_sample         = focus.append("g");
        var aa_mut            = focus.append("g");

        focus.append("rect")
                 .attr("class", "zoom")
                 .attr("y",25 +113)
                 .attr("width", width)
                 .attr("height", zoom_height - 108)
                 .call(zoom).on("mousewheel.zoom", null)
             .on("DOMMouseScroll.zoom", null) // disables older versions of Firefox
             .on("wheel.zoom", null)

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
           lines.redrawLines();

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
           lines.updateFilt();
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
           var s = widthScale.range().map(t.invertX,t)
           context.select(".selection").attr("x",s[0])
                                       .attr("width",s[1]-s[0])
           lines.redrawLines();
       }
       function zoom_end() {
           if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return;
           if(d3.event.sourceEvent==null) return;
           var t = d3.event.transform;
           widthScale.domain(t.rescaleX(width2Scale).domain());
           context.select(".brush").call(brush.move, widthScale.range().map(t.invertX, t));
       }

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
            focus.selectAll(".selected_pos").attr("width",widthScale(10)-widthScale(0));
            focus.selectAll(".selected_pos").attr("x",widthScale(selected_pos_x));


            focus.selectAll(".var_noise_indic").attr("stroke-width",widthScale(12)-widthScale(0));
            focus.selectAll(".var_noise_indic").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                    .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            lines.updateFilt();
        }

        function redraw()  {

            var wsd = widthScale.domain();
            var wsr = widthScale.range();
            var bp = (wsd[1]-wsd[0])/(wsr[1]*11);  //base per
            var bpt = bp*1000; //base per 1000

            //focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
            //focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
            //focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
            //focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");
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
            focus.selectAll(".qual_fwd").attr("x",function(d){return widthScale(d["trace_peak"])-2;});
            focus.selectAll(".qual_rev").attr("x",function(d){return widthScale(d["trace_peak"])+2;});
            focus.selectAll(".line").attr("x1",function(d){return widthScale(d["trace_peak"]);})
                                    .attr("x2",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".var_noise_indic").attr("stroke-width",widthScale(12)-widthScale(0));
            //conditional visibility
            //console.log(w);

            if(bpt<22){
                var this_show_qual = show_qual;
                focus.selectAll(".peak_label").attr("visibility","visible");
                if(bp===0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                this_show_qual = "hidden";
                focus.selectAll(".peak_label").attr("visibility","hidden");
                focus.selectAll(".ref").attr("visibility","visible");
                focus.selectAll(".user").attr("visibility","visible");
                focus.selectAll(".var_noise_indic").attr("visibility","visible");

                if(bpt<70){
                  this_show_qual = show_qual;
                  focus.selectAll(".qual_fwd").attr("visibility","visible");
                  focus.selectAll(".qual_rev").attr("visibility","visible");
                  focus.selectAll(".aa").attr("visibility","visible");
                  focus.selectAll(".coding_ten").attr("visibility","visible");
                  //focus.selectAll(".aa").attr("visibility","visible");
                }
                //if(bpt<25){focus.selectAll(".short").attr("visibility","visible");}
            }
            focus.selectAll(".q").attr("visibility",this_show_qual);
            focus.selectAll(".call").attr("visibility",'visible');

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
                //focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
                //focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
                //focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
                //focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");

            }else if(j==="TRUE"){
                heightScale_fwd = heightScale;
                heightScale_rev = heightScale;
                //focus.selectAll("g").selectAll(".line_f").attr("d",line_fwd);
                //focus.selectAll("g").selectAll(".line_r").attr("d",line_rev);
                //focus.selectAll("g").selectAll(".area_fwd").attr("d",noise_area_fwd).attr("visibility","visible");
                //focus.selectAll("g").selectAll(".area_rev").attr("d",noise_area_rev).attr("visibility","visible");
            }

        }


        function updateFilt_mini(fwd_start,fwd_end,rev_start,rev_end,has_rev){

            if(fwd_start != undefined){filt_fwd_start = fwd_start;}
            if(fwd_end != undefined){filt_fwd_end = fwd_end;}
            if(rev_start != undefined){filt_rev_start = rev_start;}
            if(rev_end != undefined){filt_rev_end = rev_end;}
            if(has_rev!=undefined){filt_has_rev = has_rev;}

            filt_fwd_start_mini.attr("x2",width2Scale(filt_fwd_start));
            filt_fwd_end_mini.attr("x1",width2Scale(filt_fwd_end)).attr("x2",width);
            filt_rev_start_mini.attr("x2",width2Scale(filt_rev_start));
            filt_rev_end_mini.attr("x1",width2Scale(filt_rev_end)).attr("x2",width);


            if(filt_rev_end == 0  ){
                //clip_rev.attr("x",0).attr("width",0);
                filt_rev_end_mini.attr("x1",0).attr("x2",0);
                //filt_line_rev.attr("x1",0).attr("x2",0)
                //             .attr("y1",0).attr("y2",0);
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
            text.enter().append("text")                         //Enter
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
                .on("click",function(d,i){posClick(d["id"],d["trace_peak"]);})
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

        }

        function setCodingLabel(calls){
            //console.log(calls);
            var aa_r = aa_ref.selectAll("text").data(calls);
            aa_r.exit().remove();
            aa_r.enter().append("text")
                .merge(aa_r).attr("class","peak_label short aa")
                .text(function(d){
                    //if   (d["ord_in_cod"] == 1) {return d["aa_ref"].toUpperCase()+""+d["codon"];}
                    //if   (d["ord_in_cod"] == 1) {return d["aa_ref"]+""+d["codon"];}
                    if   (d["ord_in_cod"] == 2) {return d["aa_ref"]+d["codon"];}
                    else {                       return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .on("click",function(d,i){posClick(d["id"],d["trace_peak"]);})
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
                .on("click",function(d,i){posClick(d["id"],d["trace_peak"]);})
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
                .on("click",function(d,i){posClick(d["id"],d["trace_peak"]);})
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
                    .on("click",function(d,i){posClick(d["id"],d["trace_peak"]);})
                    .attr("y",label_pos["qual"] -2)
                    .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "10px");
            qt.exit().remove();
            var qg = quals_g.selectAll("rect").data(calls);
            qg.enter().append("rect")
            .merge(qg).attr("class","peak_label qual_fwd q")
                    .attr("x",function(d){return (widthScale(d["trace_peak"])-2);})
                    .attr("y",function(d){return label_pos["qual"] - 11 -d[q]/3}).attr("rx",1).attr("ry",1)
                    .attr("width",2+(!rev)*2)
                    .attr("height",function(d){return d[q]/3;})
                    .attr("fill", "rgba(50,50,50,0.8)");
            qg.exit().remove();
            if(rev!=0){
                var qg_r = quals_g_r.selectAll("rect").data(calls);
                qg_r.enter().append("rect")
                .merge(qg_r).attr("class","peak_label qual_rev q")
                    .attr("x",function(d){return (widthScale(d["trace_peak"])+2);})
                    .attr("y",function(d){return label_pos["qual"] - 11 -d["quality_rev"]/3}).attr("rx",1).attr("ry",1)
                    .attr("width",2)
                    .attr("height",function(d){return  d["quality_rev"]/3;})
                    .attr("fill", "rgba(50,50,50,0.8)");
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
                .text(function(d){if(d["coding_seq"]==null){return d["id"]}else{return d["id"]+" (c."+d["coding_seq"]+")";}})
                .attr("fill","black");
        }
        function setQual(value){
            show_qual = value;
        }

        Shiny.addCustomMessageHandler("goto",
            function(message) {
                setBrush(Number(message)-180,Number(message)+200);
            }
        );

        Shiny.addCustomMessageHandler("s2n_min",
            function(message) {
                mult = Number(message);
                redraw();
            }
        );

        Shiny.addCustomMessageHandler("join",
            function(message){
                join = message;
                joinView(message);
                }
        );

        function posClick(id,tp){
            var zoom = 400
            var scope = widthScale.domain()
            var pos = (id-1)*12 - 2;
            if(pos < scope[0]){
                var from = Math.max(0,pos - 200)
                setBrush(from,from+400);
            }else if(pos > scope[1]){
                var to = Math.min(pos + 200,width2Scale.domain()[1])
                setBrush(to-400,to);
            }
            focus.selectAll(".selected_pos").attr("x",widthScale(tp-5));
            selected_pos_x = tp-5;
            Shiny.onInputChange("pos_click", {id: id});
        };
        //Expose the posClick function
        window.posClick = posClick;

        function get_color(col){
            if      (col === "A"){ return "#00A100"; }
            else if (col === "C"){ return "#2985EA"; }
            else if (col === "G"){ return "#6C6A6C"; }
            else if (col === "T"){ return "#EA2929"; }
            else if (col === "-"){ return "white"; }
            return "yellow";
        }

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
            lines: lines,
            updateFilt_mini: updateFilt_mini,
            rescaleWidth: rescaleWidth,

            brush:   brush,
            redraw:  redraw,
            join:    join,
            setQual: setQual,
            joinView:joinView,

            reHeight: reHeight,
            setBrush: setBrush,
            setPeakLabel: setPeakLabel,
            showVarInMinimap: showVarInMinimap,
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

        console.log(x);
        instance.lastValue = x;

        //the render function behaves differently when called repeatedly
        //the first run is actually still a part of the initialization step
        if(x.new_sample){
            x.new_sample = false;

            if(instance.instanceCounter>=1){ //cleanup after previous sample
                instance.lines.destroy();
                instance.focus.selectAll(".peak_label").remove();
                instance.context.selectAll(".context").remove();
            }

            instance.instanceCounter = instance.instanceCounter+1;
            var calls       = HTMLWidgets.dataframeToD3(x["samples"][0]["calls"]);
            var choices     = HTMLWidgets.dataframeToD3(x["samples"][0]["choices"]);
            //var noisy_neighbors     = HTMLWidgets.dataframeToD3(x["noisy_neighbors"]);
            var show_calls  = x["show_calls"];
            if(show_calls){
                instance.call_opacity = 0.8;
            }else{
                instance.call_opacity = 0;
            }

            var domain_y    = x["samples"][0]["intrexdat"]["max_y"];
            console.log(domain_y)
            instance.max_y  = domain_y;
            var domain_x    = x["samples"][0]["intrexdat"]["max_x"];
            instance.max_x  = domain_x;
            var intrex      = HTMLWidgets.dataframeToD3(x["samples"][0]["intrexdat"]["intrex"]);
            instance.intrex = intrex;

            var focus   = instance.focus;
            var context = instance.context;
            var brush   = instance.brush;
            var widthScale  = instance.widthScale;
            //var heightScale = instance.heightScale;
            var width2Scale = instance.width2Scale;

            widthScale.domain([0,domain_x]);
            width2Scale.domain([0,domain_x]);

            instance.lines.create(x["num_samples"],domain_x,domain_y);
            //heightScale.domain([0,domain_y]);
            //instance.heightScale_fwd_split.domain([0,domain_y]);
            //instance.heightScale_rev_split.domain([0,domain_y]);

            //visualise fullseq width
            instance.full_line
                .append("line").attr("class","context fullSeqW")
                .attr("x1",0)
                .attr("y1",25)
                .attr("x2",widthScale(domain_x))
                .attr("y2",25)
                .attr("stroke-width",3).attr("stroke","rgba(180,180,180,1.0)");

            //intron/exon boxes
            instance.setIntrexBoxes(intrex);

            context.append("g")
                .attr("class","brush context brush_context")
                .call(brush)
                .selectAll(".selection")
                .attr("fill","none")
                .attr("stroke-width",3).attr("stroke","rgb(70,130,180)")//.attr("stroke","#ff4444")//.attr("stroke-dasharray","3,6")
                .attr("y", -16)
                .attr("height", 82) //height2 + 10)
                .attr("rx",3)
                .attr("ry",3)
                .attr("shape-rendering","geometricPrecision")
                .style("filter", "url(#drop-shadow)")
                .attr("opacity",1);

            instance.lines.updateLine(x.samples);

            //noise indicator
            //var a_noise_fwd = HTMLWidgets.dataframeToD3([x["calls"]["trace_peak"],x["calls"]["noise_abs_fwd"]]);
            instance.lines.setNoiseArea(x.samples);

            //todo solve filters in mini when multiple samples
            if(x.num_samples == 1){
                instance.updateFilt_mini(x['samples'][0]['brush_fwd_start'],x['samples'][0]['brush_fwd_end'],x['samples'][0]['brush_rev_start'],x['samples'][0]['brush_rev_end'],x['samples'][0]['intens_rev'] != "");
                instance.lines.updateFilt(x['samples'][0]['brush_fwd_start'],x['samples'][0]['brush_fwd_end'],x['samples'][0]['brush_rev_start'],x['samples'][0]['brush_rev_end'],x['samples'][0]['intens_rev'] != "");
            }

            //on single strand always show "join view"
            // todo move this to chromat_lines.js
            //if(intens_rev != ""){
            //    instance.joinView();
            //}else{
            //    instance.joinView("TRUE");
            //}

            //console.log("before labels");
            focus.append("g").selectAll("text.seq.codon").data(calls).enter() //coding coord
                .append("text").attr("class",function(d){if(d["coding_seq"]%10==0){return "peak_label coding_ten"}else{return "peak_label";}})
                .text(function(d){
                    if   (d["coding_seq"] > 0){return "c." + d["coding_seq"];}
                    else {                     return "";}})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .attr("y",(instance.label_pos["codon"]))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            //console.log("inbetween labels");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
                .append("text").attr("class",function(d){if((d["coding_seq"]%10==0)&&(d["coding_seq"]>0)){return "peak_label coding_ten"}else{return "peak_label"}})
                .text(function(d){return "g." + d["gen_coord"];})
                .attr("text-anchor", "middle")
                .attr("x",function(d){return widthScale(d["trace_peak"]);})
                .attr("y",(instance.label_pos["gen_coord"]))
                .attr("fill", "black").attr("opacity", 0.8).attr("font-family", "sans-serif").attr("font-size", "11px");
            //console.log("after of labels");
            instance.setPeakLabel(calls,"reference");
            instance.setPeakLabel(calls,"call");

            //todo move sample specific labels to chromat_lines.js
            instance.setPeakLabel(calls,"mut_call_fwd");
            if(x['samples'][0]['intens_rev'] != 'undefined'){
                instance.setPeakLabel(calls,"call_rev");
                instance.setPeakLabel(calls,"mut_call_rev");
            }
            if(x.num_samples == 1){
                if(x['samples'][0]["qual_present"]){
                    console.log("set qual labels");
                    instance.setQualityLabels(calls,x['samples'][0]['intens_rev'] != 'undefined');
                }
            }
            var show_qual  = x["show_qual"];
            if(show_qual){
                instance.setQual("visible");
                d3.selectAll(".q").attr("visibility","visible");
            }else{
                instance.setQual("hidden");
                d3.selectAll(".q").attr("visibility","hidden");
            }
            //default
            focus.selectAll(".call").attr("opacity",instance.call_opacity);
            instance.setPeakLabel(calls,"user_sample");
            instance.setPeakLabel(calls,"user_mut");
            instance.setCodingLabel(calls);

            instance.showVarInMinimap(choices);
            instance.lines.updateVariants(choices,1);

            if (typeof choices[0] !== 'undefined') {
                var from = choices[0]["trace_peak"]-200;
                var to   = choices[0]["trace_peak"]+220;
                if(from < 0) {from = 0;to = 420}
                instance.setBrush(from,to);
            }else{
                instance.setBrush(200,1200);
            }

        }else{
            //console.log("render");

            var calls   = HTMLWidgets.dataframeToD3(x['samples'][0]["calls"]);

            instance.lines.setNoiseArea(x.samples)

            //todo solve filters in mini when multiple samples
            if(x.num_samples == 1){
                instance.updateFilt_mini(x['samples'][0]['brush_fwd_start'],x['samples'][0]['brush_fwd_end'],x['samples'][0]['brush_rev_start'],x['samples'][0]['brush_rev_end'],x['samples'][0]['intens_rev'] != "");
                instance.lines.updateFilt(x['samples'][0]['brush_fwd_start'],x['samples'][0]['brush_fwd_end'],x['samples'][0]['brush_rev_start'],x['samples'][0]['brush_rev_end'],x['samples'][0]['intens_rev'] != "");
            }

            if(x['samples'][0]["qual_present"]){
                console.log("set qual labels");
                instance.setQualityLabels(calls,x['samples'][0]['intens_rev'] != 'undefined');
            }
            if(instance.max_x != x['samples'][0]["intrexdat"]["max_x"]){
                instance.max_x = x['samples'][0]["intrexdat"]["max_x"];
                var domain_x    = x['samples'][0]["intrexdat"]["max_x"];
                instance.width2Scale.domain([0,domain_x]);
                var intrex      = HTMLWidgets.dataframeToD3(x['samples'][0]["intrexdat"]["intrex"]);
                instance.setIntrexBoxes(intrex);
            }
            if(x['samples'][0]["intrexdat"]["max_y"]!= instance.max_y){
                instance.reHeight(x["intrexdat"]["max_y"]);
                instance.max_y = x["intrexdat"]["max_y"];
            }else if(x['samples'][0].choices != instance.choices){
                var choices = HTMLWidgets.dataframeToD3(x['samples'][0]["choices"]);
                var show_calls  = x["show_calls"];
                if(show_calls){
                    instance.call_opacity = 0.8; }
                else{ instance.call_opacity = 0; }
                instance.focus.selectAll(".call").attr("opacity",instance.call_opacity);
                instance.setIntrexBoxes(HTMLWidgets.dataframeToD3(x['samples'][0]["intrexdat"]["intrex"]));
                instance.showVarInMinimap(choices);
                instance.choices = x['samples'][0].choices;
                instance.setPeakLabel(calls,"user_sample");
                instance.setPeakLabel(calls,"user_mut");
                instance.setPeakLabel(calls,"reference")
                instance.setPeakLabel(calls,"mut_call_fwd",instance.call_opacity);
                if(x['samples'][0]['intens_rev'] != 'undefined'){
                    instance.setPeakLabel(calls,"mut_call_rev",instance.call_opacity);
                }
                instance.setCodingLabel(calls);
                var show_qual  = x["show_qual"];
                if(show_qual){
                    //instance.show_qual = "visible";
                    instance.setQual("visible")
                    d3.selectAll(".q").attr("visibility","visible");
                }else{
                    //instance.show_qual = "hidden";
                    instance.setQual("hidden");
                    d3.selectAll(".q").attr("visibility","hidden");
                }

            }else if(x.resize == true){
                instance.lastValue.resize = false;
                instance.redraw();
                instance.setIntrexBoxes(HTMLWidgets.dataframeToD3(x["intrexdat"]["intrex"]));
                instance.showVarInMinimap(HTMLWidgets.dataframeToD3(x["choices"]));
            }else { console.log(x) }
        }
    }
});
