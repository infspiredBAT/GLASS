(function(window){

    class LineSet {

        constructor(line_el,height,widthScale) {
            this.line_el = line_el ;
            this.count = 0;
            this.height = height;
            this.bases = ["A","C","G","T"];

            this.filt_fwd_start = 0;
            this.filt_fwd_end =0;
            this.filt_rev_start = 0;
            this.filt_rev_end = 0;
            this.filt_has_rev = false;

            let half_height = height/1.6;

            this.heightScale  = d3.scaleLinear().range([height,0]),
            this.heightScale_fwd_split = d3.scaleLinear().range([half_height,(2*half_height -  height)]),
            this.heightScale_rev_split = d3.scaleLinear().range([height,half_height]);
            this.heightScale_fwd = this.heightScale_fwd_split,
            this.heightScale_rev = this.heightScale_rev_split;

            this.widthScale = widthScale;

            //const widthScale = widthScale;
            const heightScale_fwd = this.heightScale_fwd;
            const heightScale_rev = this.heightScale_rev;

            this.line_fwd = d3.line()
                .x(function(d,i){return widthScale(i)})
                .y(function(d){return heightScale_fwd(d)});
            this.line_rev = d3.line()
                .x(function(d,i){return widthScale(i)})
                .y(function(d){return heightScale_rev(d)});

            const mult = 2;

            this.noise_area_fwd = d3.area()
                        .x(function(d){return widthScale(d[0]);})
                        .y0(function(d){return (heightScale_fwd(0)+2);})
                        .y1(function(d){return heightScale_fwd(d[1]*mult);});

            this.noise_area_rev = d3.area()
                        .x(function(d){return widthScale(d[0]);})
                        .y0(function(d){return heightScale_rev(0)+2;})
                        .y1(function(d){return heightScale_rev(d[1]*mult);});
        }
        create(n,domain_x,domain_y){
            this.count = n;
            let i;
            let svg = d3.select("svg");
            //filter for Grayscale to be used on filtered part of traces
            let filter_gs   = svg.append("filter")
                               .attr("id","monochrome");
            let colormatrix = filter_gs.append("feColorMatrix")
                                .attr("type","matrix")
                                .attr("values","2 0.5 0.5 0 0 0.5 2 0.5 0 0 0.5 0.5 2 0 0 0 0 0 1 0");
            for(i=1;i<=this.count;i++){
                this.line_el.append("g").attr("class","trace_line line_a_" + i);
                this.line_el.append("g").attr("class","trace_line line_c_" + i);
                this.line_el.append("g").attr("class","trace_line line_g_" + i);
                this.line_el.append("g").attr("class","trace_line line_t_" + i);

                this.line_el.append("g").attr("class","trace_line line_a_r_" + i);
                this.line_el.append("g").attr("class","trace_line line_c_r_" + i);
                this.line_el.append("g").attr("class","trace_line line_g_r_" + i);
                this.line_el.append("g").attr("class","trace_line line_t_r_" + i);

                this.line_el.append("g").attr("class","trace_line line_a_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_c_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_g_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_t_" + i + "_gray");

                this.line_el.append("g").attr("class","trace_line line_a_r_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_c_r_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_g_r_" + i + "_gray");
                this.line_el.append("g").attr("class","trace_line line_t_r_" + i + "_gray");

                this.line_el.append("g").attr("class","gNoise_fwd_" + i);
                this.line_el.append("g").attr("class","gNoise_rev_" + i);

                this.line_el.append("g").attr("class","var_ind_focus_" + i);

                this.clip_fwd = svg.append("clipPath")
                                .attr("id", "rect_fwd_clip_" + i)
                                .append("rect")
                                .attr("x", 0).attr("y", 0)
                                .attr("width", 0)
                                .attr("height",this.height)
                                .attr("clip-path","url(#clip)");

                this.clip_rev = svg.append("clipPath")
                                .attr("id", "rect_rev_clip_" + i)
                                .append("rect")
                                .attr("x", 0).attr("y", 0)
                                .attr("width", 0)
                                .attr("height",this.height);

                this.filt_line_fwd_start = this.line_el.append("line")
                                       .attr("x1", 0).attr("y1", 0)
                                       .attr("x2", 0).attr("y2", 0)
                                       .attr("stroke-width", 2)
                                       .attr("stroke", "red")
                                       .attr("opacity",1);
                this.filt_line_fwd_end = this.line_el.append("line")
                                       .attr("x1", 0).attr("y1", 0)
                                       .attr("x2", 0).attr("y2", 0)
                                       .attr("stroke-width", 2)
                                       .attr("stroke", "red")
                                       .attr("opacity",1);
                this.filt_line_rev_start = this.line_el.append("line")
                                       .attr("x1", 0).attr("y1", 0)
                                       .attr("x2", 0).attr("y2", 0)
                                       .attr("stroke-width", 2)
                                       .attr("stroke", "red")
                                       .attr("opacity",1);
                this.filt_line_rev_end = this.line_el.append("line")
                                       .attr("x1", 0).attr("y1", 0)
                                       .attr("x2", 0).attr("y2", 0)
                                       .attr("stroke-width", 2)
                                       .attr("stroke", "red")
                                       .attr("opacity",1);


            }
            this.widthScale.domain([0,domain_x])
            this.heightScale.domain([0,domain_y]);
            this.heightScale_fwd_split.domain([0,domain_y]);
            this.heightScale_rev_split.domain([0,domain_y]);
        }

        updateLine(data_all,rev){

            const n = this.count;
            let base,data;

            for(const b in this.bases){
                base = this.bases[b];
                data = [data_all[base]]

                switch(base) {        //g2 must go first so its under g1 (svg has no depth, only order of elements)
                    case "A": if(rev){var g2 = this.line_el.selectAll(".line_a_r_" + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_a_r_" + n); }
                              else   {var g2 = this.line_el.selectAll(".line_a_"   + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_a_"   + n); }
                              var col = "#00A100"; break;
                    case "C": if(rev){var g2 = this.line_el.selectAll(".line_c_r_" + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_c_r_" + n); }
                              else   {var g2 = this.line_el.selectAll(".line_c_"   + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_c_"   + n); }
                              var col = "#2985EA"; break;
                    case "G": if(rev){var g2 = this.line_el.selectAll(".line_g_r_" + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_g_r_" + n); }
                              else   {var g2 = this.line_el.selectAll(".line_g_"   + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_g_"   + n); }
                              var col = "#6C6A6C"; break;
                    case "T": if(rev){var g2 = this.line_el.selectAll(".line_t_r_" + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_t_r_" + n); }
                              else   {var g2 = this.line_el.selectAll(".line_t_"   + n + "_gray");
                                      var g1 = this.line_el.selectAll(".line_t_"   + n); }
                              var col = "#EA2929"; break;
                }
                if(rev){var c = "line_r"; var l = this.line_rev; var cl = "url(#rect_rev_clip_" + n +")"}
                else   {var c = "line_f"; var l = this.line_fwd; var cl = "url(#rect_fwd_clip_" + n +")"}

                var line = g1.selectAll("path").data(data); //UPDATE
                line.exit().remove();                       //EXIT
                line.enter().append("path")                 //enter
                    .merge(line).attr("class", "path " + c)
                    .attr("d", l)
                    .attr("clip-path","url(#clip)")
                    .attr("filter","url(#monochrome)")
                    .attr("fill","none").attr("stroke", col)
                    .attr("stroke-width", 2);         // on reverse attr("stroke-dasharray","20,3,10,1,10,1");
                var line = g2.selectAll("path").data(data); //UPDATE
                line.exit().remove();                       //EXIT
                line.enter().append("path")                 //enter
                    .merge(line).attr("class","path "+c)
                    .attr("d",l)
                    .attr("clip-path", cl)
                    .attr("fill","none").attr("stroke",col)
                    .attr("stroke-width",2);         // on reverse attr("stroke-dasharray","20,3,10,1,10,1");
            }
        }

        updateFilt(fwd_start,fwd_end,rev_start,rev_end,has_rev){

            //should also cycle through this.count and adjust height accordingly
            if(fwd_start != undefined){this.filt_fwd_start = fwd_start;}
            if(fwd_end != undefined){this.filt_fwd_end = fwd_end;}
            if(rev_start != undefined){this.filt_rev_start = rev_start;}
            if(rev_end != undefined){this.filt_rev_end = rev_end;}
            if(has_rev!=undefined){this.filt_has_rev = has_rev;}


            let fwd_x = this.widthScale(this.filt_fwd_start);
            let fwd_width = this.widthScale(this.filt_fwd_end) - fwd_x;

            this.clip_fwd.attr("x",fwd_x).attr("width",fwd_width);
            this.filt_line_fwd_start.attr("x1",this.widthScale(this.filt_fwd_start))
                        .attr("x2",this.widthScale(this.filt_fwd_start))
                        .attr("y1",140).attr("y2",280);
            this.filt_line_fwd_end.attr("x1",this.widthScale(this.filt_fwd_end))
                        .attr("x2",this.widthScale(this.filt_fwd_end))
                        .attr("y1",140).attr("y2",280);

            if(this.filt_has_rev){
                let rev_x = this.widthScale(this.filt_rev_start);
                let rev_width = this.widthScale(this.filt_rev_end)- rev_x;
                this.clip_rev.attr("x",rev_x).attr("width",rev_width);
                this.filt_line_rev_start.attr("x1",rev_x)
                            .attr("x2",rev_x)
                            .attr("y1",280).attr("y2",420);
                this.filt_line_rev_end.attr("x1",(this.widthScale(this.filt_rev_end)-1))
                            .attr("x2",(this.widthScale(this.filt_rev_end)-1))
                            .attr("y1",280).attr("y2",420);
            }else{
                console.log("no filt on reverse");
            }

        }
        updateVariants(choices,n){
            //console.log(choices);
            let g_var = this.line_el.select(".var_ind_focus_" + n);
            let widthScale = this.widthScale;
            let v = g_var.selectAll("line").data(choices);
            v.exit().remove();
            
            //strandness in choices 0 = undetermined (should not occure); 1 = forward strand; 2 = reverse strand; 3 = both, 4 = undetermined indel, 5 = indel forward only, 6 = indel reverse only, 7 = indel both strands
            v.enter().append("line").attr("class","enter")
                .merge(v).attr("class","peak_label short line var_noise_indic")
                .attr("x1",function(d){return widthScale(d["trace_peak"]);})
                .attr("y1",function(d){if(d["strand"]==2){return 280;}else{return 140;}})
                .attr("x2",function(d){return widthScale(d["trace_peak"]);})
                .attr("y2",function(d){if(d["strand"]==1){return 270;}else{return 420;}})
                .attr("stroke-width",widthScale(12) - widthScale(0))
                .attr("stroke","rgba(255,0,0,0.12)").attr("clip-path","url(#clip)");
        }

        setNoiseArea(fwd,rev,n){

            var gnf = this.line_el.select(".gNoise_fwd_" + n).selectAll("path").data([fwd]);
            gnf.enter().append("path")
                .merge(gnf).attr("class","area area_fwd").attr("d",this.noise_area_fwd)
               .attr("fill","#000000").attr("stroke","none").attr("opacity",0.15).attr("clip-path","url(#clip)");
            gnf.exit().remove();
            if(typeof rev !== 'undefined'){
                var gnr = this.line_el.select(".gNoise_rev_" + n).selectAll("path").data([rev]);
                gnr.enter().append("path")
                   .merge(gnr).attr("class","area area_rev").attr("d",this.noise_area_rev)
                   .attr("fill","#440000").attr("stroke","none").attr("opacity",0.15).attr("clip-path","url(#clip)");
                gnr.exit().remove();
            }
        }

        redrawLines(){
            this.line_el.selectAll(".line_f").attr("d",this.line_fwd);
            this.line_el.selectAll(".line_r").attr("d",this.line_rev);

            this.line_el.select(".area_fwd").attr("d",this.noise_area_fwd).attr("visibility","visible");
            this.line_el.select(".area_rev").attr("d",this.noise_area_rev).attr("visibility","visible");
        }

        destroy(){
            console.log('removing lines');
            console.log(this.count);
            for(let i=1;i<=this.count;i++){
                console.log(this.line_el.selectAll(".trace_line"));
                this.line_el.selectAll(".trace_line").remove();
                this.line_el.selectAll(".area").remove();
                d3.select("svg").selectAll("clipPath").remove();
            }
            this.count = 0;
        }

    }

    if(typeof(window.LineSet) === 'undefined'){
         window.LineSet = LineSet;
     }
})(window);
