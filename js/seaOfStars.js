/**
 * Uniform random number in the interval [0, 1].
 */
 function random(){
	return Math.random();
}

/**
 * Uniform random number in the interval [-1, 1].
 */
function rand(){
	return Math.random() * 2.0 - 1.0;
}

/**
 * Samples disk at origin with radius.
 */
function sampleDisk(radius){
	let r = radius * Math.sqrt(random());
	let theta = random() * 2 * Math.PI;
	return [r * Math.cos(theta),  r * Math.sin(theta), 0];
}

//3D texture size.
let width = 64;
let height = 64;
let depth = 64;

let data = new Float32Array( width * height * depth);

let i = 0;
for ( let z = 0; z < depth; z ++ ) {
	for ( let y = 0; y < height; y ++ ) {
		for ( let x = 0; x < width; x ++ ) {
			data[i] = 0.9;
			i++;
		}
	}
}

/*=======================
	N-body simulation
	a_i = G * \sum_{j = 0 && j!=i}^{nBodies} m_j * (r_j * r_i) / |r_j - r_i|³
	where r is the position, a the acceleration and m the mass
	Adapted from: https://medium.com/swlh/create-your-own-n-body-simulation-with-python-f417234885e9
	Remember that 2^64 ~= 1 * 10^19
  =======================*/
var body = function(mass, position, velocity, acceleration){
	this.mass = mass; // in kg
	this.pos = position;// in m
	this.v = velocity; // in m/s
	this.acc = acceleration; // in m/s²
	return this;
};

//G = 6.67×10^-11 m³/kg²/s
const G = 6.67 * Math.pow(10, -11);
const epsilon = 0.0001;

const targetUpdate = 100000;
const dt = 1.0 / targetUpdate;

const maxMass = 1 * Math.pow(10, 6);
const minStarMass = Math.pow(10, 2);
const minMass = Math.pow(10, 1);
const maxAcc = [0.001, 0.001, 0.001];
const lightMaxIntensity = 10;
let radius = 0.5;

const nbodies = 40;
const maxLights = nbodies;
let bodies = [];

//Body with a high mass in the center simulating the sun.
mass = 1 * Math.pow(10, 15); 
pos = [0, 0, 0];
v = [0, 0, 0];
acc = [0,0,0];
bodies.push(new body(mass, pos, v, acc));

//Another big body
mass = 1 * Math.pow(10, 2); 
pos = [radius * 0.3, radius * 0.3, 0];
v = [0, 0, 0];
acc = [0,0,0];
bodies.push(new body(mass, pos, v, acc));

var light = function(position, transformWCS, Li, bodyIndex){
	this.position = position; 
	this.transformWCS = transformWCS;
	this.Li = Li;
	return this;
};

let lights = [];
lights[0] = new light(new THREE.Vector3(0,0,0), new THREE.Matrix4(), new THREE.Vector3(100, 10, 10), 0);
lights[1] = new light(new THREE.Vector3(0,0,0), new THREE.Matrix4(), new THREE.Vector3(1, 10, 100), 1);

//Initializes all bodies with random positions and masses in a disk with the sun at its center.
for(i = 2; i < nbodies; i++){
	mass = Math.max(minMass, Math.random() * maxMass); 
	lightIntensityFactor = (mass/maxMass) * lightMaxIntensity;
	pos = sampleDisk(radius);
	v = [0, 0, 0];
	acc = [0, 0, 0];

	lights[i] = new light(new THREE.Vector3(0,0,0), new THREE.Matrix4(), new THREE.Vector3(random() * lightIntensityFactor, random() * lightIntensityFactor, random() * lightIntensityFactor), i);

	bodies.push(new body(mass, pos, v, acc));
}

//Newton law of universal gravitation
//Calculates the acceleration of body_i
function calcAcc(i){
	acc = [0, 0, 0];
	
	for(j = 0; j < bodies.length; j++){
		if( i == j )
			continue;
		
		dx = bodies[j].pos[0] - bodies[i].pos[0];
		dy = bodies[j].pos[1] - bodies[i].pos[1];
		dz = bodies[j].pos[2] - bodies[i].pos[2];
		
		dx = dx / (Math.abs(dx * dx * dx) + epsilon);
		dy = dy / (Math.abs(dy * dy * dy) + epsilon);
		dz = dz / (Math.abs(dz * dz * dz) + epsilon);
		
		dx = bodies[j].mass * dx;
		dy = bodies[j].mass * dy;
		dz = bodies[j].mass * dz;
		
		acc[0] += dx;
		acc[1] += dy;
		acc[2] += dz;
		
	};
	
	//Ad hoc for visual aesthetics, not in the original formula
	bodies[i].acc[0] = G * acc[0];
	bodies[i].acc[1] = G * acc[1];
	bodies[i].acc[2] = G * acc[2];
};

for(i = 0; i < bodies.length; i++){
	calcAcc(i);
};

function updateBodies(){
	//Uses leapfrog integration
	for(i = 0; i < bodies.length; i++){
		// 1/2 kick
		bodies[i].v[0] += bodies[i].acc[0] * dt/2.0;
		bodies[i].v[1] += bodies[i].acc[1] * dt/2.0;
		bodies[i].v[2] += bodies[i].acc[2] * dt/2.0;
		
		//drift
		bodies[i].pos[0] += bodies[i].v[0] * dt;
		bodies[i].pos[1] += bodies[i].v[1] * dt;
		bodies[i].pos[2] += bodies[i].v[2] * dt;

		lights[i].position.x = 	bodies[i].pos[0];
		lights[i].position.y = 	bodies[i].pos[1];
		lights[i].position.z = 	bodies[i].pos[2];
		
		calcAcc(i);
		
		// 1/2 kick
		bodies[i].v[0] += bodies[i].acc[0] * dt/2.0;
		bodies[i].v[1] += bodies[i].acc[1] * dt/2.0;
		bodies[i].v[2] += bodies[i].acc[2] * dt/2.0;
	};
};

const scene = new THREE.Scene();
const renderer = new THREE.WebGLRenderer({ antialias: false });

/**
 * Volumetric definitions
 */
const texture = new THREE.DataTexture3D( data, width, height, depth );
texture.format = THREE.RedFormat;
texture.minFilter = THREE.LinearFilter;
texture.magFilter = THREE.LinearFilter;
texture.type = THREE.FloatType;
texture.wrapS = THREE.RepeatWrapping;
texture.wrapT = THREE.RepeatWrapping;
texture.wrapR = THREE.RepeatWrapping;

const scattering = new THREE.Vector3(0.01, 0.1, 1.1);
const absorption = new THREE.Vector3(0.1, 0.09, 0.09);

const volumeMaterial = new THREE.RawShaderMaterial( {
	glslVersion: THREE.GLSL3,
	uniforms: {
		volume: { value: texture },
		invmodel: { value:  new THREE.Matrix4()},
		cameraPosition: { value: new THREE.Vector3() },
		scattering: { value: scattering },
		absorption: { value: absorption },
		densityCoef: { value: 1 },
		g: { value: 0 },
		numberOfSteps: { value: 50 },
		numberOfShadowSteps: { value: 5 },
		lights: { value: lights }
	},
	defines: {
		NLIGHTS : lights.length 
	},
	vertexShader,
	fragmentShader,
	side: THREE.BackSide,
} );

const volumeAABB = new THREE.BoxGeometry(1, 1, 1);
volumeMesh = new THREE.Mesh( volumeAABB, volumeMaterial );
volumeMesh.scale.set(3, 1, 1);
volumeMesh.updateMatrix();
scene.add( volumeMesh );

const dummy = new THREE.Object3D();

function update(){

	volumeMesh.material.uniforms.cameraPosition.value.copy( camera.position );
	volumeMesh.material.uniforms.lights.value = lights ;
	volumeMesh.material.needsUpdate = true;

	dummy.position.set(0, 0, 0);
	dummy.scale.set(2, 1, 1);
	dummy.updateMatrix();
	
	volumeMesh.material.uniforms.invmodel.value.copy( volumeMesh.matrix.invert() );

	updateBodies();
}

function render(){
	renderer.render( scene, camera );
}

function coreLoop(){
	requestAnimationFrame( coreLoop );
	
	update();
	render();
}
		