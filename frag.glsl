#version 300 es
precision lowp float;
precision lowp int;

uniform float uC;
uniform float uN; //normalisation constants
uniform int u_n;
uniform int u_l;
uniform int u_m;
uniform float ua0;
uniform vec3 uDimensions;
uniform float uStepSf;
uniform float uMaxR;
uniform int uVolumePlaneSurface;
uniform int uAxis;
uniform float uPlanePos;

in vec3 vRd;
flat in vec3 vRo;

out vec4 fragColor;

float gen_laguerre(int n, float alpha, float x){
	if (n == 0) return 1.0;
	if (n == 1) return 1.0 + alpha - x;

	float prev2 = 1.0;
	float prev1 = 1.0 + alpha -x;
	float L = 0.0;

	for (int i = 2; i <= n; i++){
		float fi = float(i);
		L = ((2.0*fi-1.0 + alpha - x)*prev1 - (fi-1.0 + alpha)*prev2)/fi;
		prev2 = prev1;
		prev1 = L;
	}
	
	return L;
}

float assoc_legendre_recur(int m, int l, float x)
{
    float P_mm = 1.0; //P_0,0
    float fm = float(m);
    float fl = float(l);
    if (m > 0)
    {
        float termx2 = sqrt(1.0 - x*x);
        float coeff = 1.0;

        for (int i = 1; i <= m; i++)
        {
            P_mm *= -coeff * termx2;
            coeff += 2.0;
        }
    }

    if (l == m) return P_mm;

    float P_mmp1 = x * (2.0*fm + 1.0) * P_mm;

    if (l == m + 1) return P_mmp1;

    float P = 0.0;

    for (int n = m+2; n <=l; n++){
	float fn = float(n);
    	P = (x*(2.0*fn-1.0)*P_mmp1 - (fn+fm-1.0)*P_mm)/(fn-fm);
    	P_mm = P_mmp1;
    	P_mmp1 = P;
    }
    return P;
}


float psi(vec3 pos) {
	float f_un = float(u_n);
	float f_ul = float(u_l);
	float f_um = float(u_m);
	float r = length(pos);
	if (r < 1e-6) return 0.0;
	float theta = acos(pos.z/r); float phi = atan(pos.y,pos.x);

	//radial
	float rho = (2.0*r)/(f_un*ua0);
	float L = gen_laguerre(u_n - u_l -1, 2.0*f_ul+1.0, rho);
	float R = uC * exp(-rho * 0.5) * pow(rho,f_ul) * L;

	//angular
	int absM = abs(u_m);
	float fabsM = float(absM);
	float P = assoc_legendre_recur(absM,u_l,cos(theta));
	float Y;
	if (u_m < 0) {Y = sqrt(2.0) * uN * P * sin(fabsM*phi);} // check these defs lmao
	else if (u_m > 0) {Y = sqrt(2.0) * uN * P * cos(f_um*phi); }
	else {Y = uN * P;}

	return R*Y;
}

vec2 boxIntersection(vec3 ro, vec3 rd) {
	const vec3 box_min = vec3(0);
	const vec3 box_max = vec3(1);
	vec3 inv_rd = 1.0/rd;
	vec3 tmin_tmp = (box_min - ro) *inv_rd;
	vec3 tmax_tmp = (box_max - ro) *inv_rd;
	vec3 tmin = min(tmin_tmp, tmax_tmp);
	vec3 tmax = max(tmin_tmp, tmax_tmp);
	float t0 = max(tmin.x, max(tmin.y, tmin.z));
	float t1 = min(tmax.x, min(tmax.y, tmax.z));
	return vec2(t0,t1);
}
vec4 transfer(float val) {
	float opacity = clamp(abs(val)*5.0,0.0,1.0);
	vec3 color = val>= 0.0 ? vec3(1.0,0.2,0.1) : vec3(0.1,0.2,1.0); //change- this is just orange/blue for now
	return vec4(color,opacity);
}
vec4 raymarch(in vec3 ro, in vec3 rd, float t0, float t1){
	if (t0 > t1) return vec4(0.0);  
	t0 = max(t0,0.0);

	vec3 vec3dt = 1.0/(vec3(uDimensions) * abs(rd)); 
	//can add zoom factor yap
	float dt = min(vec3dt.x,min(vec3dt.y,vec3dt.z)) * uStepSf / uMaxR;
	
	vec3 p = ro + (t0 + dt*0.5) * rd; //start 1 step ahead to prevent issues with box clipping
	vec4 color = vec4(0.0);

	for (float t = t0 + dt*0.5; t < t1; t+=dt) {
		if (any(lessThan(p, vec3(0.0))) || any(greaterThan(p, vec3(1.0)))) break; // throw out if not in [0,1] cube
		vec3 phys_p = (p -0.5) *2.0 * uMaxR;
		float v = psi(phys_p); //psi wf is for an orbital centred at (0,0,0)
		if (abs(v) < 1e-4){
			p += rd * 4.0*dt;
		}
		vec4 val_color = transfer(v);
		val_color.rgb *= val_color.a;
		color += val_color * (1.0 - color.a);
		//color.rgb += (1.0 - color.a) * val_color.a * val_color.rgb;
		//color.a += (1.0 - color.a) * val_color.a;

		if (color.a >= 0.9) break;
		p += rd * dt;
	}
	return color;
}
vec4 plane(in vec3 ro, in vec3 rd, float t0, float t1) {
	float d = uPlanePos; //at origin rn
	vec3 N;
	//duplicated; prolly move into main()
	if (uAxis == 0) N = vec3(0.0,0.0,1.0); //xy 
	else if (uAxis == 1) N = vec3(0.0,1.0,0.0);  //xz
	else N = vec3(1.0,0.0,0.0);  //yz

	if (abs(dot(rd,N)) < 1e-6) return vec4 (0.0); //parallel to normal - avoid division by 0

	float tplane = (d-dot(ro,N)) /dot(rd,N);

	if (tplane <t0 || tplane > t1) return vec4(0.0);

	vec3 pplane = ro + tplane *rd;

	if (any(lessThan(pplane, vec3(0.0))) || any(greaterThan(pplane, vec3(1.0)))) return vec4(0.0); //unit cube

	vec3 phys_p = (pplane-0.5) *2.0 *uMaxR;
	float v = psi(phys_p);
	return transfer(v);
}
// can make vec3 isosurface  & cut planes
void main() {
	vec3 rd = normalize(vRd), ro = vRo; 

	vec2 hitray = boxIntersection(ro,rd);
	float t0 = hitray.x;
	float t1 = hitray.y;

	if (uVolumePlaneSurface == 0 || uVolumePlaneSurface == 2) fragColor = raymarch(ro,rd, t0,t1);
	else if (uVolumePlaneSurface == 1) fragColor = plane(ro,rd, t0,t1);
}
