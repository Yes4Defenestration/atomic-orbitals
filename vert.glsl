#version 300 es
layout(location=0) in vec2 aPos; //unit cube [0,1]

uniform vec3 uCam;
uniform vec3 uTarget;
uniform vec3 uUp;
uniform float uFov;
uniform float uAspect;

out vec3 vRd;
flat out vec3 vRo;

void main(){
	gl_Position = vec4(aPos,0.0,1.0);
	vec3 forward = normalize(uTarget-uCam);
	vec3 right = normalize(cross(forward,uUp));
	vec3 up = cross(right,forward);
	
	vec2 screenPos = aPos;

        float tanHalfFov = tan(uFov * 0.5);

	// Perspective ray direction
	vRd = forward + right * screenPos.x * uAspect * tanHalfFov
                      + up    * screenPos.y * tanHalfFov;

	vRo = uCam;

	//vec3 vol_translation = vec3(0.5) - uScale * 0.5; //min corner
	//	gl_Position = uMVP * vec4(pos * uScale + vol_translation, 1);
	//vRo = (uCam - vol_translation) / uScale;
	//vRd = pos - vRo;
}
