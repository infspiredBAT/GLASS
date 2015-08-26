// This recieves messages of type "testmessage" from the server.
  Shiny.addCustomMessageHandler("testmessage",
    function(message) {
      if(message == "ready") {
//        alert(JSON.stringify(message));
//        var bc = document.getElementById('go_load').style.backgroundColor;  //#DDD
        var bc = "#DDD";
//        var tc = document.getElementById('go_load').style.color;  //#DA230F
        var tc = "#DA230F";
        document.getElementById('go_load').style.backgroundColor = tc;
        document.getElementById('go_load').style.color = bc;
      } else {
        var bc = document.getElementById('go_load').style.backgroundColor;
        var tc = document.getElementById('go_load').style.color;
        document.getElementById('go_load').style.backgroundColor = tc;
        document.getElementById('go_load').style.color = bc;
      }
    }
  );