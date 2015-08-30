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
            width2Scale  = d3.scale.linear().range([0,width]),
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

		function brushed() {
	        widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
	        focus.selectAll("g").selectAll("path").attr("d", line);
		}
        function reHeight(domain_y){
            heightScale.domain([0,domain_y]);
            console.log("reseting.height",domain_y);
            focus.selectAll("g").selectAll("path").attr("d", line);
        }
        function reWidth(domain_x){
            brush.extent([0,domain_x]);
            widthScale.domain(brush.empty() ? width2Scale.domain() : brush.extent());
            focus.selectAll("g").selectAll("path").attr("d", line);
        }

        //passing arguments
        //this enables to access vars and functions from the render function as instance.*
        return {
            svg: svg,
            line: line,
            linec: linec,
		    context: context,
	        brush: brush,
            focus: focus,
            widthScale: widthScale,
            heightScale: heightScale,
            width: w,
            height: h,
	        height2: height2,
            instanceCounter: instanceCounter,
            reHeight: reHeight,
            reWidth: reWidth
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
        //a somewhat nasty hack, the render function behaves differently when called repeatedly
        //the first run is actually still an initialization step
        if(instance.instanceCounter === 0){

//			console.log("first drawing")
			instance.instanceCounter = instance.instanCounter+1;
			var intens = x["intens"];
			var calls = HTMLWidgets.dataframeToD3(x["call"]);
			var choices = HTMLWidgets.dataframeToD3(x["choices"]);
			var domain_y = x["helperdat"]["max_y"];
			instance.max_y = domain_y;
			var domain_x = x["helperdat"]["max_x"];
			instance.max_x = domain_x;
			var intrex = HTMLWidgets.dataframeToD3(x["helperdat"]["helper_intrex"])
			instance.intrex = intrex;

			var svg = instance.svg;
			var line = instance.line;
			var linec = instance.linec;
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
			context.selectAll("lines.intrex").data(intrex).enter()
				.append("line")
				.attr("x1",function(d){return widthScale(d["start"]);})
				.attr("y1",0)
				.attr("x2",function(d){return widthScale(d["start"]);})
				.attr("y2",50)
				.attr("stroke-width",2).attr("stroke","rgba(20,20,20,0.6)").attr("stroke-dasharray",2);
//				.on("mouseover", function(){d3.select(this).style("fill", "white");})
//				.on("mouseout",  function(){d3.select(this).style("fill", "rgba(200,200,200,0.2)");});
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
				.text(function(d){return Math.ceil(widthScale(d["start"]));})
				.attr("fill","black");
			// genomic
			context.selectAll("lines.choices").data(choices).enter()
				.append("line")
				.attr("x1",function(d){return d["id"];})
				.attr("y1",6)
				.attr("x2",function(d){return d["id"];})
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
				.attr("x1",function(d){return d["id"];})
				.attr("y1",26)
				.attr("x2",function(d){return d["id"];})
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
				.attr("x",function(d){return d["id"]+4;})
				.attr("y",30)
				.attr("opacity",0.6)
				.text(function(d){return d["id"];})
				.attr("fill","black");

			brush.x(width2Scale);

			var group_a = focus.append("g");
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
            focus.append("g").selectAll("qualities").data(calls).enter()
				.append("rect")
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",0).attr("rx",2).attr("ry",2)
				.attr("width",5*5)
				.attr("height",function(d){return d["quality"];})
				.attr("fill", "rgba(200,200,200,0.4)");
            focus.append("g").selectAll("text.qualities").data(calls).enter()
				.append("text")
				.text(function(d){return d["quality"];})
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",4)
				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.seq.user").data(calls).enter()
				.append("text")
				.text(function(d){return d["call"];})
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",16)
				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.seq.genomic").data(calls).enter()
				.append("text")
				.text(function(d){return d["reference"];})
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",28)
				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.coord.genomic").data(calls).enter()
				.append("text")
				.text(function(d){return d["gen_coord"];})
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",40)
				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
            focus.append("g").selectAll("text.exon_intron").data(calls).enter()
				.append("text")
				.text(function(d){return d["exon_intron"];})
				.attr("x",function(d){return d["trace_peak"]*5;})
				.attr("y",52)
				.attr("fill", "black").attr("opacity", 0.7).attr("font-family", "sans-serif").attr("font-size", "10px");
			focus
				.append("line")
				.attr("x1",0)
				.attr("y1",100)
				.attr("x2",1200)
				.attr("y2",100)
				.attr("stroke-width",1).attr("stroke","rgba(0,0,0,0.6)").attr("stroke-dasharray",2);

/*
			var group_ac = context.append("g");
			var group_cc = context.append("g");
			var group_gc = context.append("g");
			var group_tc = context.append("g");

			group_ac.selectAll("path")
				.data([intens["A"]]).enter()
				.append("path").attr("d",linec)
				.attr("fill","none")
				.attr("stroke","#33CC33")
				.attr("stroke-width",0.5);
			group_cc.selectAll("path")
				.data([intens["C"]]).enter()
				.append("path").attr("d",linec)
				.attr("fill","none")
				.attr("stroke","#0000FF")
				.attr("stroke-width",0.5);
			group_gc.selectAll("path")
				.data([intens["G"]]).enter()
				.append("path").attr("d",linec)
				.attr("fill","none")
				.attr("stroke","#000000")
				.attr("stroke-width",0.5);
			group_tc.selectAll("path")
				.data([intens["T"]]).enter()
				.append("path").attr("d",linec)
				.attr("fill","none")
				.attr("stroke","#FF0000")
				.attr("stroke-width",0.5);
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

	          //zooming in so that the first view is not dense ugly graph
	          //works but does not show the brush tool
	          //instance.reWidth(2000);

        }else{
			console.log("redrawing");
			if(x["helperdat"]["max_y"]!= instance.max_y){
				instance.reHeight(x["helperdat"]["max_y"]);
			}
        }
    }
});
