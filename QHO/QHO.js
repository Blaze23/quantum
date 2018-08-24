"use strict";
//Global Initial Parameters:
var currentHref = window.location.href;

var m = 1;
var w = 1;
var hbar = 1;
var N = 2000;
var dx;
var xRange = 20;
var coeffs = [1, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0];
var userCoeffs = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
var ENUM = 11;
var observed = false;
var userShow = false;

var hermCoeffs = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [-2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, -12, 0, 8, 0, 0, 0, 0, 0, 0, 0],
                [12, 0, -48, 0, 16, 0, 0, 0, 0, 0, 0, 0],
                [0, 120, 0, -160, 0, 32, 0, 0, 0, 0, 0],
                [-120, 0, 720, 0, -480, 0, 64, 0, 0, 0, 0],
                [0, -1680, 0, 3360, 0, -1344, 0, 128, 0, 0, 0],
                [1680, 0, -13440, 0, 13440, 0, -3584, 0, 256, 0, 0],
                [0, 30240, 0, -80640, 0, 48384, 0, -9216, 0, 512, 0],
                [-30240, 0, 302400, 0, -403200, 0, 161280, 0, -23040, 0, 1024]];

/*var hermCoeffs = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [3, 0, -6, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 15, 0, -10, 0, 1, 0, 0, 0, 0, 0],
                [-15, 0, 45, 0, -15, 0, 1, 0, 0, 0, 0],
                [0, -105, 0, 105, 0, -21, 0, 1, 0, 0, 0],
                [105, 0, -420, 0, 210, 0, -28, 0, 1, 0, 0],
                [0, 945, 0, -1260, 0, 378, 0, -36, 0, 1, 0],
                [-945, 0, 4725, 0, -3150, 0, 630, 0, -45, 0, 1]];*/

var layout = {
    width: 400, height: 400,
    margin: {l:50, r:50, t:50, b:50},
    hovermode: "closest",
    showlegend: false,
    xaxis: {label: 'x', range: [-10,10]},
    yaxis: {lavel: 'y', range: [-1,1]},
    aspectratio: {x:1, y:1}
};

function factorialRt(num){
    if(num == 1 || num == 0){
        return 1;
    }

    return Math.sqrt(num)*factorialRt(num-1);
}

function hermite(n, x){
    var hermVal = 0;
    for(var i = 0; i<=n; i++){
        hermVal += hermCoeffs[n][i]*Math.pow(x, i);
    }

    return hermVal;
}

function eigenstate(n, c, xArr, hermArr){
    var qhoVal;
    var f = factorialRt(n);

    for(var i=0; i<N; i++){
        var x = xArr[i]

        var c1 = 1/(Math.sqrt(Math.pow(2, n))*f);
        var c2 = Math.pow(m*w/Math.PI*hbar, 1/4);
        var c3 = Math.exp(-m*w*Math.pow(x, 2)/(2*hbar));
        var c4 = hermite(n, Math.sqrt(m*w/hbar)*x);

        qhoVal = Math.sqrt(c)*c1*c2*c3*c4;

        hermArr[i] += qhoVal;

    }

    return hermArr;

}

function QHOstate(coeffArray){
    var sum = 0;
    var data = setArr();

    for(var i=0; i<ENUM; i++){
        if(!isNaN(coeffArray[i])){
            sum += coeffArray[i];
        }
    }

    if(sum != 0 && !isNaN(sum)){
        for(var i=0; i<ENUM; i++){
            var c = coeffArray[i]/sum;
            data[1] = eigenstate(i, c, data[0], data[1]);
        }
    }

    return data;
}

function setArr(){
    var xArr = [];
    var hermArr = [];
    var vArr = [];
    
    for(var i=0; i<N; i++){
        var x = i*dx-(xRange/2);
        xArr.push(x);
        hermArr.push(0);
        vArr.push(0.01*(x+5)*(x-5));
    }

    return [xArr, hermArr, vArr];
}

//Plots
function initPlot() {
    Plotly.purge("graph");
    var graphData = QHOstate(coeffs);

    var data = [{
        x: graphData[0],
        y: graphData[1],
        type: 'line'
      }];
    
    data.push({
        x: graphData[0],
        y: graphData[2],
        type: 'line'
    })

    Plotly.newPlot("graph", data, layout);
    userShow = false;

    return;
}

function select(){
    var rnd = Math.random();
    var sum = 0;
    var runningSum = 0;
    var eigen = -1;

    for(var i = 0; i<ENUM; i++){
        sum += coeffs[i];
    }

    for(var i = 0; i<ENUM; i++){
        runningSum += coeffs[i]/sum;

        if(rnd < runningSum){
            if(eigen == -1){
                eigen = i;
            }
        }
    }

    return eigen;

}

function collapse() {
    if(!observed){
        observed = true;
        Plotly.purge("graph");

        var graphData = setArr();

        var eigenvalue = select();

        graphData[1] = eigenstate(eigenvalue, 1, graphData[0], graphData[1]);

        var data = [{
            x: graphData[0],
            y: graphData[1],
            type: 'line'
        }];

        data.push({
            x: graphData[0],
            y: graphData[2],
            type: 'line'
        })

        Plotly.newPlot("graph", data, layout);
        Plotly.plot;

        var energyUnits = eigenvalue + 0.5;

        var eigenString = energyUnits.toString();

        document.getElementById("eigenvalue").innerHTML = eigenString;
    }
    return;
}

function reset(){
    observed = false;
    document.getElementById("eigenvalue").innerHTML = "";
    initPlot();
}

function randomise(){
    var rnd;

    for(var i=0; i<ENUM; i++){
        rnd = Math.round(3*Math.random());
        coeffs[i] = rnd;
    }

    reset();
}


function main() {
    dx = xRange/N;

    initPlot();
}

function tabChange(idName, btnName){
    var hideTab = document.getElementsByClassName("tab-pane active");
    var showTab = document.getElementById(idName);
    var activeBtn = document.getElementById(btnName);
    var hideBtn = document.getElementsByClassName("button active");

    for(var i = 0; i<hideTab.length; i++){
        hideTab[i].className = "tab-pane";
    }

    for(var i = 0; i<hideBtn.length; i++){
        hideBtn[i].className = "button";
    }

    showTab.className = "tab-pane active";
    activeBtn.className = "button active";
}

function sliderChange(){
    userCoeffs[0] = parseInt($('#zero').val());
    userCoeffs[1] = parseInt($('#one').val());
    userCoeffs[2] = parseInt($('#two').val());
    userCoeffs[3] = parseInt($('#three').val());
    userCoeffs[4] = parseInt($('#four').val());
    userCoeffs[5] = parseInt($('#five').val());
    userCoeffs[6] = parseInt($('#six').val());
    userCoeffs[7] = parseInt($('#seven').val());
    userCoeffs[8] = parseInt($('#eight').val());
    userCoeffs[9] = parseInt($('#nine').val());

    plotUserWave();
}

function plotUserWave(){
    var graphData = QHOstate(userCoeffs);

    var data = [{
        x: graphData[0],
        y: graphData[1],
        type: 'line'
      }];

    if(userShow){
        Plotly.deleteTraces("graph", 2);
    }

    Plotly.plot("graph", data);  
    userShow = true; 
}

$(document).ready(main);