const vertexShader = `
	in vec3 position;

	out vec3 rayDirection;
	out vec3 vertCoord;
	flat out vec3 translation;
	flat out vec3 scale;

	uniform mat4 modelMatrix;
	uniform mat4 modelViewMatrix;
	uniform mat4 projectionMatrix;
	uniform vec3 cameraPosition;

	void main(){

		scale.x = modelMatrix[0][0];
		scale.y = modelMatrix[1][1];
		scale.z = modelMatrix[2][2];

		translation.x = modelMatrix[3][0];
		translation.y = modelMatrix[3][1];
		translation.z = modelMatrix[3][2];
		
		vertCoord = position;
		vec3 posWCS = vec3(modelMatrix * vec4(position, 1.0));
		rayDirection = posWCS - cameraPosition;
		gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
	}
`;

const fragmentShader = `
	precision highp float;
	precision highp sampler3D;
	#define PI 3.1415926535897932384626433832795
	#define INV4PI 0.07957747154594766788

	out vec4 color; 
	in vec3 rayDirection;
	in vec3 vertCoord;
	flat in vec3 translation;
	flat in vec3 scale;

	uniform mat4 modelMatrix;
	uniform mat4 modelViewMatrix;
	uniform mat4 projectionMatrix;
	uniform mat4 invmodel;
	uniform vec3 cameraPosition;

	uniform sampler3D volume;

	/* Volume properties */
	uniform vec3 scattering;
	uniform vec3 absorption;
	#define extinction (absorption + scattering)
	#define albedo (scattering/extinction)
	uniform float densityCoef;

	/* Phase function properties */
	uniform float g;

	/* Method Params */
	uniform int numberOfSteps;
	uniform int numberOfShadowSteps;
	const float transmittanceThreshold = 0.01;

	/* Light properties */
	struct Light {    
		vec3 position;
		mat4 transformWCS;
		vec3 Li;
	};

	uniform Light lights[NLIGHTS];

	struct Ray{
		vec3 origin;
		vec3 direction;
	};

	//Boundaries of a unit cube centered at origin
	vec3 boxMin = ( vec4(-1.5f, -0.5f, -0.5f,1)).xyz; // TODO removed modelMatrix multiplication
	vec3 boxMax = ( vec4(1.5, 0.5, 0.5, 1)).xyz;

	vec3 getPointAt(Ray r, float t){
		return r.origin + t * r.direction;
	}

	/*
		Ray must be in WCS
	*/
	vec2 intersectBox(vec3 orig, vec3 dir, vec3 bmin, vec3 bmax) {
		//Line's/Ray's equation
		// o + t*d = y
		// t = (y - o)/d
		//when t is negative, the box is behind the ray origin
		vec3 tMinTemp = (bmin - orig) / dir;
		vec3 tmaxTemp = (bmax - orig) / dir;

		vec3 tMin = min(tMinTemp, tmaxTemp);
		vec3 tMax = max(tMinTemp, tmaxTemp);

		float t0 = max(tMin.x, max(tMin.y, tMin.z));
		float t1 = min(tMax.x, min(tMax.y, tMax.z));

		return vec2(t0, t1);
	}

	/*
		if(t.x == 0)
			origin inside AABB
		if(t.x > t.y)
			miss
	*/
	vec2 intersectBox(vec3 orig, vec3 dir) {
		return intersectBox(orig, dir, boxMin, boxMax);
	}

	float density(ivec3 gridPoint){
		return texelFetch( volume, gridPoint, 0).x;
	}

	float interpolatedDensity(vec3 gridPoint) {
		vec3 pSamples = vec3(gridPoint.x - .5f, gridPoint.y - .5f, gridPoint.z - .5f);
		ivec3 pi = ivec3(floor(pSamples.x), floor(pSamples.y), floor(pSamples.z));
		vec3 d = pSamples - vec3(pi);

		// Trilinearly interpolate density values to compute local density
		float d00 = mix(density(pi), density(pi + ivec3(1, 0, 0)), d.x);
		float d10 = mix(density(pi + ivec3(0, 1, 0)), density(pi + ivec3(1, 1, 0)), d.x);
		float d01 = mix(density(pi + ivec3(0, 0, 1)), density(pi + ivec3(1, 0, 1)), d.x);
		float d11 = mix(density(pi + ivec3(0, 1, 1)), density(pi + ivec3(1, 1, 1)), d.x);
		float d0 = mix(d00, d10, d.y);
		float d1 = mix(d01, d11, d.y);
		return mix(d0, d1, d.z);
	}

	float sampleVolume(vec3 pos){
		//tex must be in Texture Coordinate System [0,0,0]->[1,1,1]
		//and this cube model must have width = 1 and be centered at origin.
		vec3 tex = ( invmodel * vec4(pos,1)).xyz + 0.5;
		vec3 texGridPoint = tex * vec3(textureSize(volume, 0));

		return interpolatedDensity(texGridPoint);
	}

	vec3 evalHG(vec3 incoming, vec3 scattered) {
		float cosTheta = dot(normalize(incoming), normalize(scattered));
		float denom = 1.0f + g * g + 2.0f * g * cosTheta;
		
		float res = INV4PI * (1.0f - g * g) / (denom * sqrt(denom));
		return vec3(res, res, res);
	}

	//geometric term
	vec3 evaluateLight(Light light, in vec3 pos){
		float d = length(light.position - pos);
		return (light.Li/100.0) / (d*d);
	}

	vec3 volumetricShadow(in vec3 cubePos, in vec3 lightPos){
		vec3 transmittance = vec3(1.0);
		float distance = length(lightPos - cubePos) / float(numberOfShadowSteps);
		vec3 lightDir = normalize(lightPos - cubePos);
		float stepSizeShadow = 1.0 / float(numberOfShadowSteps);

		for(float tshadow = stepSizeShadow; tshadow < distance; tshadow += stepSizeShadow){
			vec3 cubeShadowPos = cubePos + tshadow * lightDir;
			

			float density = sampleVolume(cubeShadowPos);
			if(density == 0.0)
				continue;
			

			// e ^(-ext * density[i] + -ext * density[i+1] + ...)
			transmittance *= exp(-extinction * density);
			if(transmittance.x < transmittanceThreshold)
				break;
		}
		return transmittance;
	}

	bool integrate(Ray incomingRay, Light light, vec2 tHit, out vec3 transmittance, out vec3 inScattering){
		if (tHit.x > tHit.y) 
			return false;

		tHit.x = max(tHit.x, 0.0f);
		vec3 absDir = abs(incomingRay.direction);
		float dt = 1.0f / ( float(numberOfSteps) * max(absDir.x, max(absDir.y, absDir.z)));
		incomingRay.origin = getPointAt(incomingRay, tHit.x);
		vec3 totalTr = vec3(1.0f);
		vec3 sum = vec3(0.0f);
		vec3 currentTr = vec3(0);

		for(float t = tHit.x; t < tHit.y; t += dt){
			float density = sampleVolume(incomingRay.origin);
			if(density == 0.0){
				incomingRay.origin = getPointAt(incomingRay, dt);
				continue;
			}
			
			//If accumulated transmittance is near zero, there's no point in continue calculating transmittance
			if(totalTr.x > transmittanceThreshold){
				currentTr = exp(-extinction * density * densityCoef);
				totalTr *= currentTr;
			}
			
			vec3 Ls = evaluateLight(light, incomingRay.origin) 
						* volumetricShadow(incomingRay.origin, light.position) 
						* evalHG(incomingRay.direction, normalize(light.position - incomingRay.origin));
						
			//Integrate Ls from 0 to d
			Ls =  (Ls - Ls * currentTr) / (extinction * densityCoef); 

			sum += totalTr * (scattering * densityCoef) * Ls;
			incomingRay.origin = getPointAt(incomingRay, dt);
		}

		transmittance = totalTr;
		inScattering = sum;
		return true;
	}

	void main(){
		vec4 thisColor;
		vec3 transmittance;
		vec3 inScattering;
		vec4 bg = vec4(0,0,0,1);
		Ray incomingRay = Ray(cameraPosition, normalize(rayDirection));
		vec2 tHit = intersectBox(incomingRay.origin, incomingRay.direction); 
		
		for(int i = 0; i < NLIGHTS; i++){
			if(integrate(incomingRay, lights[i], tHit, transmittance, inScattering))
				thisColor += vec4(bg.xyz * transmittance + inScattering, 1.0f);
		}
		color = thisColor;
		//color = vec4(lights[1].position.xyz,1);
	}
`;