class Atom {
  constructor(symbol, id, connections, connectionTypes) {
    let atom = "";
    for (let j = 0; j < elements.length; j++) {
      let json = elements[j];
      if (json['symbol']==symbol) {
        atom = json;
      }
    }
    this.mass = atom['atomic_mass'];
    this.EN = atom['electronegativity_pauling']==null?2:atom['electronegativity_pauling'];
    this.id = id;
    let shells = atom['shells'];
    this.connectionCount = shells.length==1?2-shells[0]:8-shells[shells.length-1];
    this.connections = connections;
    this.connectionTypes = connectionTypes;
    //this.position = p5.Vector.random3D();
    this.pos = p5.Vector.random3D().normalize();
    this.force = createVector(0, 0);
    let col = unhex("FF"+(atom['cpk-hex']==null?"FFFFFF":atom['cpk-hex']))
      this.c = col;
    let hex = (atom['cpk-hex']==null?"FFFFFF":atom['cpk-hex']),
      r = unhex(hex.substring(0, 2)),
      g = unhex(hex.substring(2, 4)),
      b = unhex(hex.substring(4, 6));
    this.symbol = symbol;
    this.pg = createGraphics(200, 200, RGB);
    this.pg.background(r, g, b);
    this.pg.textAlign(CENTER, CENTER);
    this.pg.textSize(60);
    this.pg.textFont('Arial');
    this.pg.noStroke();
    this.pg.fill(int((0.299 * r + 0.587 * g + 0.114 * b) < 128) * 255);
    this.pg.text(symbol, 50, 100);
    this.pg.text(symbol, 150, 100);
  }
  update() {
    if (!this.locked) {
      this.pos.add(p5.Vector.div(this.force.limit(100), this.mass));
    }
    this.connections.length = this.connectionCount;
    this.connectionTypes.length = this.connectionCount;
    push();
    translate(this.pos.x, this.pos.y, this.pos.z);
    texture(this.pg);
    textureMode(IMAGE);
    sphere(20);
    pop();
    this.force = createVector(0, 0);
  }
  rotateAX(sin, cos) {
    let pos = this.position;
    this.pos = createVector(pos.x, pos.y * cos + pos.z * sin, pos.y * -sin + pos.z * cos);
  }
  rotateAY(sin, cos) {
    let pos = this.position;
    this.pos = createVector(pos.x*cos+pos.z*-sin, pos.y, pos.x * sin + pos.z * cos);
  }
}
let distMoved = 0;
let atoms = [];
let elements = [];

// Put any asynchronous data loading in preload to complete before "setup" is run
function preload() {
  data = loadJSON("data/periodic-table-lookup.json");
}
let data = {}, slider;
function windowResized() {
  resizeCanvas(windowWidth+2, windowHeight+2);
}

function setup() {
  let cnv = createCanvas(windowWidth, windowHeight, WEBGL);
  cnv.position(0, 0);
  document.getElementById("picker").addEventListener("click", function() {
    document.getElementById("c").click();
  }
  );
  document.getElementById("upload").addEventListener("click", function() {
    document.getElementById("hiddenUpload").click();
  }
  );
  show();
  var ignoreClickOnMeElement = document.getElementById('settings');
  document.addEventListener('click', function(event) {
    var isClickInsideElement = ignoreClickOnMeElement.contains(event.target);
    if (!isClickInsideElement) {
      show();
    }
  }
  );

  perspective(PI/3.0, width/height, 1, 100000);
  setupSim();
}
let recovery = false;
function setupSim() {
  let json = data;
  let order = json["order"];
  elements = [];
  for (let i = 0; i < order.length; i++) {
    elements[i] = json[order[i]];
  }
  let input = "2,3dimethyl3,4diethylhexan";
  if (document.getElementById('tbiupac').value!="" && !recovery) {
    input=document.getElementById('tbiupac').value;
  }
  recovery = false;
  let alkane = ["methan",
    "ethan",
    "propan",
    "butan",
    "pentan",
    "hexan",
    "heptan",
    "octan",
    "nonan",
    "decan"];
  let rest = ["methyl",
    "ethyl",
    "propyl",
    "butyl",
    "pentyl",
    "hexyl",
    "heptyl",
    "octyl",
    "nonyl",
    "decyl"];
  let kind = -1;

  for (let i = 0; i < 10; i++) {
    if (input.includes(alkane[i])) {
      print(alkane[i], i);
      kind = i;
      break;
    }
  }
  if (kind==-1) {
    alert("Name of longest chain is incorrect. Correct it to continue");
    recovery = true;
    setupSim();
    return;
  }
  input = input.replaceAll(alkane[kind], "");
  let Hauptkette = [[]];
  for (let i = 0; i < kind; i++) {
    Hauptkette.push([]);
  }
  let loops = input.split("yl").length-1;
  for (let n = 0; n < loops; n++) {
    let i1 = 0, i2 = input.length;
    for (let i = 0; i<input.length; i++) {
      if (match(str(input.charAt(i)), "[0-9,]") == null) {
        i1 = i;
        break;
      }
    }
    for (let i = i1; i<input.length; i++) {
      if (match(str(input.charAt(i)), "[a-z]") == null) {
        i2 = i;
        break;
      }
    }
    let subKind = 0;
    let places = input.substring(0, i1).split(",");
    let part = input.substring(i1, i2);
    input = input.substring(i2);
    for (let i = 0; i < 10; i++) {
      if (part.includes(rest[i])) {
        print(rest[i], i);
        subKind = i+1;
        break;
      }
    }
    if (subKind == 0 | part.length==0) {
      alert("Fault in: "+ input.substring(0, i2));
      recovery = true;
      setupSim();
      return;
    }
    for (let h = 0; h < places.length&&match(places[h], "[0-9]")>0; h++) {
      Hauptkette[int(places[h])] = Hauptkette[int(places[h])].concat([subKind]);
    }
  }
  atoms=[];
  let len1 = 0, len2 = 0;
  for (let i = 0; i < Hauptkette.length; i++) {
    len1 = max(len1, Hauptkette[i]);
    for (let n = 0; n < Hauptkette[i].length; n++) {
      len2 = max(len2, Hauptkette[i][n]);
    }
  }
  let ketten = [];
  for (let i = 0; i < Hauptkette.length; i++) {
    for (let val = 0; val < Hauptkette[i].length; val++) {
      ketten.push([Hauptkette[i][val], i]);
    }
  }
  for (let j = 0; j < ketten.length; j++) {
    print(ketten[j][0].length);
    //for (let l = 0; l < ketten[j][0].length; l++) {
    let lange = ketten[j][0];
    if (lange==0) {
      continue;
    }
    let subkette = [[]];
    for (let i = 0; i < lange; i++) {
      subkette.push([0, 0, 0, 0]);
    }
    subkette[0][0] = -2;
    for (let i = 0; i < lange; i++) {
      if (i>0) {
        subkette[i][0] = i;
        subkette[i-1][3] = i+1;
      }
      subkette[i] = [subkette[i][0], -1, -1, -1];
    }
    let lid = -1;
    for (let i = 0; i < lange; i++) {
      let id = atoms.length;
      if (lid!=-1) {
        for (let n = 0; n < 4; n++) {
          if (atoms[lid].connections[n]==-1) {
            atoms[lid].connections[n]=id;
            atoms[lid].connectionTypes[n]=1;
            break;
          }
        }
      }
      if (lid==-1) {
        ketten[j] = [ketten[j][0], ketten[j][1], id];
      }
      atoms.push(new Atom("C", id, [lid, -1, -1, -1], [1, 1, 1, 1]));
      for (let n = 0; n < 4; n++) {
        if (subkette[i][n] == -1) {
          let idh = atoms.length;
          atoms.push(new Atom("H", idh, [id], [0]));
          atoms[id].connections[n]=idh;
        }
      }
      lid=id;
    }
  }
  let lid = -1;
  for (let i = 0; i < Hauptkette.length; i++) {
    let id = atoms.length;
    atoms.push(new Atom("C", id, [lid, -1, -1, -1], [1, 0, 0, 0]));
    if (lid!=-1) {
      for (let n = 0; n < 4; n++) {
        if (atoms[lid].connections[n] == -1) {
          atoms[lid].connections[n] = id;
          atoms[lid].connectionTypes[n] = 1;
          break;
        }
      }
    }
    for (let a = 0; a < ketten.length; a++) {
      arr = ketten[a];
      if (arr.length > 2) {
        let ar = arr[1];
        if (ar == i+1) {
          for (let n = 0; n < 4; n++) {
            if (atoms[i].connections[n]==-1) {
              atoms[id].connections[n]=arr[2];
              atoms[id].connectionTypes[n]=1;
              break;
            }
          }
          for (let n = 0; n < 4; n++) {
            if (atoms[arr[2]].connections[n]==-1) {
              atoms[arr[2]].connections[n]=id;
              atoms[arr[2]].connectionTypes[n]=1;
              break;
            }
          }
        }
      }
    }
    lid = id;
  }
  for (let k = 0; k < atoms.length; k++) {
    for (let n = 0; n < atoms[k].connectionCount; n++) {
      if (atoms[k].connections[n] == -1) {
        let idh = atoms.length;
        atoms.push(new Atom("H", idh, [atoms[k].id], [1]));
        atoms[k].connections[n] = idh;
        atoms[k].connectionTypes[n] = 1;
      }
    }
  }
}
function parse(file) {
  print(file);
  let atomList = JSON.parse(file)['atoms'], atomsTemp = [];
  atoms = [];
  for (let i = 0; i < atomList.length; i++) {
    atom = atomList[i];
    print(atom);
    let a = new Atom(atom['symbol'], atom['id'], atom['connections'], atom['connectionTypes']);
    for (let n = 0; n < a.connectionTypes.length; n++) {
      if (a.connectionTypes[n] == 0) {
        let id = atomList.length+atomsTemp.length;
        a.connections[n] = id;
        a.connectionTypes[n] = 1;
        atomsTemp.push(new Atom("H", id, [a.id], [1]));
      }
    }
    atoms.push(a);
  }
  atoms = concat(atoms, atomsTemp);
  print(atoms);
}
let Gain = 0.03, theta = 0, mp;
function touchStarted() {
  mp = true;
  mx = mouseX;
  my = mouseY;
}
function touchEnded() {
  mp =false;
}
let mx = 0, my = 0, alpha = 0, beta = 0, avg;
function draw() {
  background(0);
  avg = createVector(0, 0);
  for (let i = 0; i < atoms.length; i++) {
    avg.add(atoms[i].pos);
  }
  avg.div(atoms.length);
  z += zoom*20;
  camera(width/2.0, height/2.0, (height/2.0+z) / tan(PI*30.0 / 180.0), width/2.0, height/2.0, z, 0, 1, 0);
  let connections = [];
  if (((mouseX < 180 && mouseY < 250 && shown) | !mp) && document.getElementById("rotate").checked) {
    theta += pow(10, (int(document.getElementById("speed").value)/100)-3);
  }
  if (mp && !(mouseX < 180 && mouseY < 250)) {
    let dx = map(mouseX - mx, -width, width, -PI, PI)*2, dy = map(mouseY - my, -height, height, -PI, PI);
    beta += dx;
    alpha += dy;
    mx = mouseX;
    my = mouseY;
  }

  translate(width/2-avg.x, height/2-avg.y, -avg.z);
  rotateX(-alpha);
  rotateY(beta + theta);

  let col = document.getElementById('c').value.replace("#", "");
  let red = parseInt(col.substring(0, 2), 16), green = parseInt(col.substring(2, 4), 16), blue = parseInt(col.substring(4, 6), 16);

  for (let ai1 = 0; ai1 < atoms.length; ai1++) {
    let a1 = atoms[ai1];
    for (let ai2 = 0; ai2 < atoms.length; ai2++) {
      let a2 = atoms[ai2];
      if (a1.id==a2.id) {
        continue;
      }
      let attract = a2.connections.includes(a1.id) | a1.connections.includes(a2.id);
      let f = int(attract) * attraction(a1, a2) + coulomb(a1, a2);
      stroke(red, green, blue);
      strokeWeight(8);
      if (attract & connections.indexOf(int(str(min(a1.id, a2.id))+str(max(a1.id, a2.id))))==-1) {
        connections.push(int(str(min(a1.id, a2.id))+str(max(a1.id, a2.id))));
        let normal = p5.Vector.sub(a1.pos, a2.pos).cross(p5.Vector.sub(a2.pos, createVector(width/2, height/2, -2000))).setMag(8);
        let index = a1.connections.indexOf(a2.id);
        switch(index != -1 ? a1.connectionTypes[a1.connections.indexOf(a2.id)] : -1) {
        case 1:
          lineVec(a1.pos, a2.pos);
          break;
        case 2:
          let offsetdi = createVector(0, 8, 0);
          lineVec(p5.Vector.add(a1.pos, offsetdi), p5.Vector.add(a2.pos, offsetdi));
          lineVec(p5.Vector.sub(a1.pos, offsetdi), p5.Vector.sub(a2.pos, offsetdi));
          break;
        case 3:
          let offsettri2 = mat(p5.Vector.sub(a1.pos, a2.pos).copy().normalize(), normal, 0.8660254, -0.5),
            offsettri3 = mat(p5.Vector.sub(a1.pos, a2.pos).copy().normalize(), normal, -0.8660254, -0.5);
          lineVec(p5.Vector.add(a1.pos, normal), p5.Vector.add(a2.pos, normal));
          lineVec(p5.Vector.add(a1.pos, offsettri2), p5.Vector.add(a2.pos, offsettri2));
          lineVec(p5.Vector.add(a1.pos, offsettri3), p5.Vector.add(a2.pos, offsettri3));
          break;
        }
      }
      a1.force.add(p5.Vector.sub(a2.pos, a1.pos).setMag(f*5));
      a2.force.add(p5.Vector.sub(a1.pos, a2.pos).setMag(f*5));
    }
    noStroke();
    a1.update();
  }
}
function reset() {
  atoms = [];
  setupSim();
}
let z = 0;
let zoom = 0;
function keyPressed() {
  if (keyCode==38) {
    zoom = 1;
  }
  if (keyCode==40) {
    zoom = -1;
  }
  if (key == 'r') {
    reset();
  }
}
function keyReleased() {
  zoom = 0;
}
function coulomb(a1, a2) {
  let r = max(p5.Vector.dist(a1.pos, a2.pos), 0.01), K = -2000;
  return K*((a1.EN*a2.EN)/(r*r));
}
function attraction(a1, a2) {
  let r = max(p5.Vector.dist(a1.pos, a2.pos), 0.01), K = 0.008;
  return K*a1.EN*a2.EN*r;
}
function lineVec(p1, p2) {
  line(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
}
function mat(v, p, sin, cos) {
  cos = 1 - cos;
  let W = [[1, -v.z*sin, v.y*sin], [v.z*sin, 1, -v.x*sin], [-v.y*sin, v.x*sin, 1]], w2 = [[0, -v.z, v.y], [v.z, 0, -v.x], [-v.y, v.x, 0]], nmat = [[]];
  for (let i = 0; i < 3; i++) {
    for (let n = 0; n < 3; n++) {
      nmat[i][n] = createVector(w2[i][0], w2[i][1], w2[i][2]).dot(createVector(w2[0][n], w2[1][n], w2[2][n]))*cos+W[i][n];
    }
  }
  ret = [];
  for (let i = 0; i < 3; i++) {
    for (let n = 0; n < 3; n++) {
      ret[n]=ret[n]+nmat[i][n]*p.array()[i];
    }
  }
  return createVector(ret[0], ret[1], ret[2]);
}
function tbclick() {
  setupSim();
}
let shown = false;
function show() {
  shown = false;
  document.getElementById("picker").style.display = "none";
  document.getElementById("c").style.display = "none";
  document.getElementById("LabBind").style.display = "none";
  document.getElementById("rotate").style.display = "none";
  document.getElementById("LabRot").style.display = "none";
  document.getElementById("speLab").style.display = "none";
  document.getElementById("speed").style.display = "none";
  document.getElementById("tbiupac").style.display = "none";
  document.getElementById("upload").style.display = "none";
  document.getElementById("set").style.display = "initial";
  document.getElementById("rect").style.height = "50";
  document.getElementById("rect").style.width = "50";
}

function hide() {
  shown = true;
  document.getElementById("picker").style.display = "initial";
  document.getElementById("c").style.display = "initial";
  document.getElementById("LabBind").style.display = "initial";
  document.getElementById("rotate").style.display = "initial";
  document.getElementById("LabRot").style.display = "initial";
  document.getElementById("speLab").style.display = "initial";
  document.getElementById("speed").style.display = "initial";
  document.getElementById("tbiupac").style.display = "initial";
  document.getElementById("set").style.display = "none";
  document.getElementById("rect").className = "shadow2";
  document.getElementById("rect").style.height = "250";
  document.getElementById("rect").style.width = "180";
  document.getElementById("upload").style.display = "initial";
}
let dataRaw = "";
function clickUpload() {
  print("File submitted");
  document.getElementById('hiddenUpload').files[0].text().then(data=>parse(data));
}
function dropHandler(ev) {
  console.log('File(s) dropped');
  // Prevent default behavior (Prevent file from being opened)
  ev.preventDefault();
  if (ev.dataTransfer.items) {
    // Use DataTransferItemList interface to access the file(s)
    for (var i = 0; i < ev.dataTransfer.items.length; i++) {
      // If dropped items aren't files, reject them
      if (ev.dataTransfer.items[i].kind === 'file') {
        var file = ev.dataTransfer.items[i].getAsFile().text().then(data=>parse(data));
      }
    }
  } else {
    // Use DataTransfer interface to access the file(s)
    for (var f = 0; f < ev.dataTransfer.files.length; f++) {
      let file = ev.dataTransfer.items[f].getAsFile();
      file.text().then(data=>parse(data));
    }
  }
}
function dragOverHandler(ev) {
  console.log('File(s) in drop zone');

  // Prevent default behavior (Prevent file from being opened)
  ev.preventDefault();
}

//function matRot(v, sin, cos) {
//  cos = 1 - cos;
//  print(sin, cos);
//  let W = [[1, -v.z*sin, v.y*sin], [v.z*sin, 1, -v.x*sin], [-v.y*sin, v.x*sin, 1]], w2 = [[0, -v.z, v.y], [v.z, 0, -v.x], [-v.y, v.x, 0]], nmat = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
//  for (let i = 0; i < 3; i++) {
//    for (let n = 0; n < 3; n++) {
//      nmat[i][n] = createVector(w2[i][0], w2[i][1], w2[i][2]).dot(createVector(w2[0][n], w2[1][n], w2[2][n]))*cos+W[i][n];
//    }
//  }
//  return nmat;
//}
//function matWorld(n1, n2) {
//  let n3 = p5.Vector.cross(n1, n2);
//  return [concat(n1.array(), [0]), concat(n3.array(), [0]), concat(p5.Vector.cross(n1, n3).array(), [0]), [0, 0, 0, 1]];
//}
//function QuickInverse(m) { // Only for Rotation/Translation Matrices
//  let mat = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
//  mat[0][0] = m[0][0];
//  mat[0][1] = m[1][0];
//  mat[0][2] = m[2][0];
//  mat[0][3] = 0.0;
//  mat[1][0] = m[0][1];
//  mat[1][1] = m[1][1];
//  mat[1][2] = m[2][1];
//  mat[1][3] = 0.0;
//  mat[2][0] = m[0][2];
//  mat[2][1] = m[1][2];
//  mat[2][2] = m[2][2];
//  mat[2][3] = 0.0;
//  mat[3][0] = -(m[3][0] * mat[0][0] + m[3][1] * mat[1][0] + m[3][2] * mat[2][0]);
//  mat[3][1] = -(m[3][0] * mat[0][1] + m[3][1] * mat[1][1] + m[3][2] * mat[2][1]);
//  mat[3][2] = -(m[3][0] * mat[0][2] + m[3][1] * mat[1][2] + m[3][2] * mat[2][2]);
//  mat[3][3] = 1.0;
//  return mat;
//}
//function matMult(m1, m2) {
//  let mat = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
//  for (let i1 = 0; i1 < 4; i1++) {
//    for (let i2 = 0; i2 < 4; i2++) {
//      for (let i3 = 0; i3 < 4; i3++) {
//        mat[i1][i2] += m2[i3][i2] * m1[i1][i3];
//      }
//    }
//  }
//  return mat;
//}
