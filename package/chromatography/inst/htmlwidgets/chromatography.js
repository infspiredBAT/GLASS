HTMLWidgets.widget({

    name: 'chromatography',
    type: 'output',

    initialize: function(el, w, h) {

        var instanceCounter = 0;
        var intrex = "";
        var max_x = 0;
        var max_y = 0;
        var margin  = {top: 10,  right: 10, bottom: 100, left: 40},
            margin2 = {top: 230, right: 10, bottom: 20,  left: 40},
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

        function redraw()  {
            widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            var w = brush.extent()[1]-brush.extent()[0] ;
            focus.selectAll("g").selectAll("path").attr("d", line);
            focus.selectAll(".peak_label").attr("x",function(d){return widthScale(d["trace_peak"]);});
            focus.selectAll(".q").attr("x",function(d){return widthScale(d["trace_peak"])-9;});
            //conditional visibility
            if(w<410){
                focus.selectAll(".peak_label").attr("visibility","visible");
                if(w==0){focus.selectAll(".peak_label").attr("visibility","hidden");}
            }else{
                focus.selectAll(".peak_label").attr("visibility","hidden");
                if(w<800){focus.selectAll(".q").attr("visibility","visible")}
                if(w<2000){focus.selectAll(".short").attr("visibility","visible");}
            }

        }
        function brushed() {
            redraw();
    		}
        function reHeight(domain_y){
            heightScale.domain([0,domain_y]);
            //console.log("reseting.height",domain_y);
            focus.selectAll("g").selectAll("path").attr("d", line);
        }
        //setting brush programmatically
        function setBrush(start,end){
            context.call(brush.extent([start,end]));
            redraw();
        }
        function showVarInMap(choices){

            //genomic
            console.log("changing choices");
            context.selectAll("lines.choices").data(choices).enter()
        			.append("line")
                    .attr("class","varInMinimap")
      				.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      				.attr("y1",6)
      				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      				.attr("y2",24)
      				.attr("stroke-width",3)
      				.attr("stroke",function(d) {
      				    if      (d["reference"] === "A"){ return "#33CC33"; }
      				    else if (d["reference"] === "C"){ return "#0000FF"; }
      				    else if (d["reference"] === "G"){ return "#000000"; }
      				    else if (d["reference"] === "T"){ return "#FF0000"; }
      				    else    {                         return "white";  }});
      			// user
      			context.selectAll("lines.choices").data(choices).enter() //function(d){return d["trace_peak"]*5;})
      				.append("line")
                    .attr("class","varInMinimap")
      				.attr("x1",function(d){return width2Scale(d["trace_peak"]);})
      				.attr("y1",26)
      				.attr("x2",function(d){return width2Scale(d["trace_peak"]);})
      				.attr("y2",44)
      				.attr("stroke-width",3)
      				.attr("stroke",function(d) {
      				    if      (d["call"] === "A"){ return "#33CC33"; }
      				    else if (d["call"] === "C"){ return "#0000FF"; }
      				    else if (d["call"] === "G"){ return "#000000"; }
      				    else if (d["call"] === "T"){ return "#FF0000"; }
      				    else    {                    return "white";  }});
      			context.selectAll("text.choices.coord").data(choices).enter()
      				.append("text")
                    .attr("class","varInMinimap")
      				.attr("x",function(d){return width2Scale(d["trace_peak"])+4;})
      				.attr("y",30)
      				.attr("opacity",0.6)
      				.text(function(d){return d["id"];})
      				.attr("fill","black");
            focus.append("g").selectAll("variance_indicator").data(choices).enter()  //variance indicator
      				.append("line").attr("class","peak_label q")
      				.attr("x1",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y1",function(d){return d["quality"]+30;})
      				.attr("x2",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y2",190)
      				.attr("stroke-width",1).attr("stroke","rgba(255,0,0,0.3)").attr("stroke-dasharray",2);
        }
        //passing arguments
        //this enables to access vars and functions from the render function as instance.*
        return {
            svg:     svg,
            line:    line,
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
            showVarInMap:showVarInMap
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
    //function called everytime input paramters change
    renderValue: function(el, x, instance) {

        instance.lastValue = x;
        //the render function behaves differently when called repeatedly
        //the first run is actually still a part initialization step
        if(instance.instanceCounter === 0){

      			console.log(x)
      			instance.instanceCounter = instance.instanCounter+1;
      			var intens = x["intens"];
      			var intens_guide_line = x["intens_guide_line"];
      			var calls = HTMLWidgets.dataframeToD3(x["calls"]);
      			var choices = HTMLWidgets.dataframeToD3(x["choices"]);
      			var domain_y = x["helperdat"]["max_y"];
      			instance.max_y = domain_y;
      			var domain_x = x["helperdat"]["max_x"];
      			instance.max_x = domain_x;
      			var intrex = HTMLWidgets.dataframeToD3(x["helperdat"]["helper_intrex"])
      			instance.intrex = intrex;

      			var svg     = instance.svg;
      			var line    = instance.line;
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
      			//visualise introns/exons
            //lines
      			context.selectAll("lines.intrex").data(intrex).enter()
      				.append("line")
      				.attr("x1",function(d){return widthScale(d["start"]);})
      				.attr("y1",0)
      				.attr("x2",function(d){return widthScale(d["start"]);})
      				.attr("y2",50)
      				.attr("stroke-width",2).attr("stroke","rgba(20,20,20,0.6)").attr("stroke-dasharray",2);
      //				.on("mouseover", function(){d3.select(this).style("fill", "white");})
      //				.on("mouseout",  function(){d3.select(this).style("fill", "rgba(200,200,200,0.2)");});
            //intron/exon boxes
      			context.selectAll("rect").data(intrex).enter()
      				.append("rect")
      				.attr("x",function(d){return widthScale(d["start"]);})
      				.attr("y",0).attr("rx",5).attr("ry",5)
      				.attr("width",function(d){return widthScale(d["end"]-d["start"]);})
      				.attr("height",50)
      				.attr("fill",function(d) {
      				    if (/exon/.test(d)){ return "blue";
      				    } else {             return "rgba(200,200,200,0.3)"; }
      				});
      //				.attr("stroke-width",1).attr("stroke","rgba(20,20,20,0.8)").attr("stroke-dasharray",2)
      //				.on("mouseover", function(){d3.select(this).style("fill", "white");})
      //				.on("mouseout",  function(){d3.select(this).style("fill", "rgba(200,200,200,0.2)");});

      			context.selectAll("text.intrex.name").data(intrex).enter()
      				.append("text")
      				.attr("x",function(d){return widthScale(d["start"]);})
      				.attr("y",-2)
      				.attr("opacity",0.6)
      				.text(function(d){return d["attr"];})
      				.attr("fill","black");
      			context.selectAll("text.intrex.start").data(intrex).enter()
      				.append("text")
      				.attr("x",function(d){return widthScale(d["start"]);})
      				.attr("y",60)
      				.attr("opacity",0.6)
      				.text(function(d){return d["id"];}) // position labels !extract sequence coords
      				.attr("fill","black");

      			brush.x(width2Scale);

      			var group_a = focus.append("g");  //why do I need a group for each line?
      			var group_c = focus.append("g");
      			var group_g = focus.append("g");
      			var group_t = focus.append("g");

      			group_a.selectAll("path").data([intens["A"]]).enter()
      				.append("path").attr("class","path")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#33CC33").attr("stroke-width",0.75);
      			group_c.selectAll("path").data([intens["C"]]).enter()
      				.append("path")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#0000FF").attr("stroke-width",0.75);
      			group_g.selectAll("path").data([intens["G"]]).enter()
      				.append("path")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#000000").attr("stroke-width",0.75);
      			group_t.selectAll("path").data([intens["T"]]).enter()
      				.append("path")
      				.attr("d",line)
      				.attr("fill","none")
      				.attr("stroke","#FF0000").attr("stroke-width",0.75);


            //trace peak labels
            focus.append("g").selectAll("qualities").data(calls).enter()  //quality box
      				.append("rect").attr("class","peak_label q")
      				.attr("x",function(d){return (widthScale(d["trace_peak"])-9);})
      				.attr("y",0).attr("rx",2).attr("ry",2)
      				.attr("width",18)
      				.attr("height",function(d){return d["quality"];})
      				.attr("fill", "rgba(200,200,200,0.4)");
            focus.append("g").selectAll("text.qualities").data(calls).enter() // quality number !do we need it?
      				.append("text").attr("class","peak_label")
      				.text(function(d){return d["quality"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",-1)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.seq.genomic").data(calls).enter() //reference
        			.append("text").attr("class","peak_label short")
      				.text(function(d){return d["reference"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",10)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px")
                    .attr("stroke",function(d) {
        			    if      (d["reference"] === "A"){ return "#33CC33"; }
      				    else if (d["reference"] === "C"){ return "#0000FF"; }
      				    else if (d["reference"] === "G"){ return "#000000"; }
      				    else if (d["reference"] === "T"){ return "#FF0000"; }
      				    else    {                         return "#000000"; }});
            focus.append("g").selectAll("text.seq.user").data(calls).enter() //call
      				.append("text").attr("class","peak_label short")
      				.text(function(d){return d["call"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",22)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px")
                    .attr("stroke",function(d) {
          		        if      (d["call"] === "A"){ return "#33CC33"; }
      				    else if (d["call"] === "C"){ return "#0000FF"; }
      				    else if (d["call"] === "G"){ return "#000000"; }
      				    else if (d["call"] === "T"){ return "#FF0000"; }
      				    else    {                    return "#000000"; }});
            focus.append("g").selectAll("text.seq.user").data(calls).enter() //user_mod
      				.append("text").attr("class","peak_label short")
      				.text(function(d){return d["user_mod"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",34)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px")
                    .attr("stroke",function(d) {
          		        if      (d["user_mod"] === "A"){ return "#33CC33"; }
      				    else if (d["user_mod"] === "C"){ return "#0000FF"; }
      				    else if (d["user_mod"] === "G"){ return "#000000"; }
      				    else if (d["user_mod"] === "T"){ return "#FF0000"; }
      				    else    {                        return "#000000"; }});
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter() //gen coord
      				.append("text").attr("class","peak_label")
      				.text(function(d){return d["gen_coord"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",50)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.exon_intron").data(calls).enter() //intrex
      				.append("text").attr("class","peak_label")
      				.text(function(d){return d["exon_intron"];})
      				.attr("text-anchor", "middle")
      				.attr("x",function(d){return widthScale(d["trace_peak"]);})
      				.attr("y",62)
      				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");

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
      				.call(brush)
      				.selectAll("rect")
      				.attr("y", -14)
      				.attr("height", 78) //height2 + 10)
      				.attr("rx",3)
      				.attr("ry",3)
      				.attr("fill","rgba(255,255,255,0.999)")
      				.attr("stroke-width",2).attr("stroke","red").attr("stroke-dasharray","2,6")
      				.attr("opacity",0.5);
      			//context.selectAll("g").selectAll("path").attr("opacity",0.5);

            //In SVG, z-index is defined by the order the element appears in the document
            //http://stackoverflow.com/questions/17786618/how-to-use-z-index-in-svg-elements
            //the vars on the minimap are shown last sto that they stay on to
            instance.showVarInMap(choices);
	          //zooming in so that the first view is not ugly dense graph

            if (typeof choices[0] !== 'undefined') {
                instance.setBrush((choices[0]["trace_peak"]-300),(choices[0]["trace_peak"]+320));
            }else{
                instance.setBrush(200,1000);
            }

        }else{

      			if(x["helperdat"]["max_y"]!= instance.max_y){
      				instance.reHeight(x["helperdat"]["max_y"]);
      			}else if(x.choices != instance.choices){
              var choices = HTMLWidgets.dataframeToD3(x["choices"]);
              instance.context.selectAll(".varInMinimap").remove();
              instance.showVarInMap(choices);
              instance.choices = x.choices;
      			}
        }
    }
});
