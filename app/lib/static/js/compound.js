document.addEventListener('DOMContentLoaded', function(){
curr_url = window.location.href;
url_name = curr_url.split("/");
process_id = url_name[url_name.length-1];
$.getJSON("{{url_for('network.net_result',process_id='ADDSHARE1') }}".replace("ADDSHARE1", process_id),
function(data) {
    edge_data = data;
});
$.getJSON("{{url_for('network.net_style',process_id='ADDSHARE1') }}".replace("ADDSHARE1", process_id),
function(data) {
    module_style = data;
});
var cy = window.cy = cytoscape({
	container: document.getElementById('cy'),
	ready: function() {
		let layoutUtilities = this.layoutUtilities({
			desiredAspectRatio: this.width() / this.height(),
			componentSpacing: 30
		});
		this.layout({
			name: 'fcose',
			animationEasing: 'ease-out'
		}).run();
	},
	layout: {
		name: 'fcose'
	},
	style: module_style,
	elements: edge_data
	
	
});
cy.on('tap', 'node',
function(evt) {
	var node = evt.target;
	console.log('tapped ' + node.id());
});
var node_list = [];
var clust_id = []
cy.ready(function() {
	cy.nodes().forEach(function(n) {
		if(!n.data('id').startsWith('clust')){
			node_list.push(n.data('id'));
			var $links = [{
				name: 'GeneCard',
				url: 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=' + n.data('id')
			},
			{
				name: 'UniProt search',
				url: 'http://www.uniprot.org/uniprot/?query=' + n.data('id') + '&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score'
			}].map(function(link) {
				return h('a', {
					target: '_blank',
					href: link.url,
					'class': 'tip-link'
				},
				[t(link.name)]);
			});

			var tippy = makeTippy(n, h('div', {},
			$links));

			n.data('tippy', tippy);

			n.on('click',
			function(e) {
				tippy.show();

				cy.nodes().not(n).forEach(hideTippy);
			});
		}else{	
			clust_id.push(n.data('id'));
		}
		
	});

});



var t = function(text){
  var el = document.createTextNode(text);
  return el;
};

var h = function(tag, attrs, children){
  var el = document.createElement(tag);

  Object.keys(attrs).forEach(function(key){
    var val = attrs[key];

    el.setAttribute(key, val);
  });

  children.forEach(function(child){
    el.appendChild(child);
  });

  return el;
};

var makeTippy = function(node, html){
  return tippy( node.popperRef(), {
    html: html,
    trigger: 'manual',
    arrow: true,
    placement: 'bottom',
    hideOnClick: false,
    interactive: true
  } ).tooltips[0];
};
//poper tip
var hideTippy = function(node){
  var tippy = node.data('tippy');

  if(tippy != null){
    tippy.hide();
  }
};

var hideAllTippies = function(){
  cy.nodes().forEach(hideTippy);
};

cy.on('tap',
function(e) {
	var evtTarget = e.target;
	if (evtTarget === cy) {
		hideAllTippies();
		console.log('tap on background');
	} else {
		console.log('tap on some element');
	}
});

cy.on('tap', 'edge', function(e){
  hideAllTippies();
});

cy.on('zoom pan', function(e){
  hideAllTippies();
});

autocomplete(document.getElementById("Search"), node_list);
document.getElementById("Search").addEventListener("keyup",
function(e) {
	if (e.keyCode == 13) {
		var input_protein = document.getElementById(this.id).value;
		cy.getElementById(input_protein).select();
	}
});
document.getElementById("fcoseButton").addEventListener("click",
function() {
	var layout = cy.layout({
		name: 'fcose',
		quality: 'proof',
		randomize: true,
		animate: true,
		animationEasing: 'ease-out'
	});
	layout.run();
});

document.getElementById("separate").addEventListener("click", function(){
	if(!(document.getElementById("separate").checked) & clust_id.length!=0){
		cy.nodes().children().move({parent: null});
		for(clust_num=0;clust_num<clust_id.length;clust_num++){
			cy.remove(cy.getElementById(clust_id[clust_num]));
		}	
	}
	
});



document.getElementById("a_dtn").addEventListener('click',
function() {
	var jpg64 = cy.jpg();
	document.querySelector('#a_dtn').setAttribute('href', jpg64);
	//document.querySelector('#jpg-eg').setAttribute('src', jpg64);
});
	//node manipulate
  
function autocomplete(inp, arr) {
	/*the autocomplete function takes two arguments,
  the text field element and an array of possible autocompleted values:*/
	var currentFocus;
	/*execute a function when someone writes in the text field:*/
	inp.addEventListener("input",
	function(e) {
		var a, b, i, val = this.value;
		/*close any already open lists of autocompleted values*/
		closeAllLists();
		if (!val) {
			return false;
		}
		currentFocus = -1;
		/*create a DIV element that will contain the items (values):*/
		a = document.createElement("DIV");
		a.setAttribute("id", this.id + "autocomplete-list");
		a.setAttribute("class", "autocomplete-items");
		a.style.opacity = 0.5;
		a.style.position = "fixed";
		a.style.zIndex = "99999";
		/*append the DIV element as a child of the autocomplete container:*/
		this.parentNode.appendChild(a);
		/*for each item in the array...*/
		for (i = 0; i < arr.length; i++) {
			/*check if the item starts with the same letters as the text field value:*/
			if (arr[i].substr(0, val.length).toUpperCase() == val.toUpperCase()) {
				/*create a DIV element for each matching element:*/
				b = document.createElement("DIV");
				/*make the matching letters bold:*/
				b.innerHTML = "<strong>" + arr[i].substr(0, val.length) + "</strong>";
				b.innerHTML += arr[i].substr(val.length);
				/*insert a input field that will hold the current array item's value:*/
				b.innerHTML += "<input type='hidden' value='" + arr[i] + "'>";
				/*execute a function when someone clicks on the item value (DIV element):*/
				b.addEventListener("click",
				function(e) {
					/*insert the value for the autocomplete text field:*/
					inp.value = this.getElementsByTagName("input")[0].value;
					/*close the list of autocompleted values,
			  (or any other open lists of autocompleted values:*/
					closeAllLists();
				});
				a.appendChild(b);
			}
		}
	});
	/*execute a function presses a key on the keyboard:*/
	inp.addEventListener("keydown",
	function(e) {
		var x = document.getElementById(this.id + "autocomplete-list");
		if (x) x = x.getElementsByTagName("div");
		if (e.keyCode == 40) {
			/*If the arrow DOWN key is pressed,
		increase the currentFocus variable:*/
			currentFocus++;
			/*and and make the current item more visible:*/
			addActive(x);
		} else if (e.keyCode == 38) { //up
			/*If the arrow UP key is pressed,
		decrease the currentFocus variable:*/
			currentFocus--;
			/*and and make the current item more visible:*/
			addActive(x);
		} else if (e.keyCode == 13) {
			/*If the ENTER key is pressed, prevent the form from being submitted,*/
			/*e.preventDefault();*/
			if (currentFocus > -1) {
				/*and simulate a click on the "active" item:*/
				if (x) x[currentFocus].click();
			}
		}
	});
	function addActive(x) {
		/*a function to classify an item as "active":*/
		if (!x) return false;
		/*start by removing the "active" class on all items:*/
		removeActive(x);
		if (currentFocus >= x.length) currentFocus = 0;
		if (currentFocus < 0) currentFocus = (x.length - 1);
		/*add class "autocomplete-active":*/
		x[currentFocus].classList.add("autocomplete-active");
	}
	function removeActive(x) {
		/*a function to remove the "active" class from all autocomplete items:*/
		for (var i = 0; i < x.length; i++) {
			x[i].classList.remove("autocomplete-active");
		}
	}
	function closeAllLists(elmnt) {
		/*close all autocomplete lists in the document,
	except the one passed as an argument:*/
		var x = document.getElementsByClassName("autocomplete-items");
		for (var i = 0; i < x.length; i++) {
			if (elmnt != x[i] && elmnt != inp) {
				x[i].parentNode.removeChild(x[i]);
			}
		}
	}
	/*execute a function when someone clicks in the document:*/
	document.addEventListener("click",
	function(e) {
		closeAllLists(e.target);
	});
}//end of autocomplete





});