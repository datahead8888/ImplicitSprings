//Modified from http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/cubeXform.frag
//and http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf
//and based on the code presented by Soumya Dutta in class

//This shader implements the equations for local illumination

uniform vec4 ambient_coef;
uniform vec4 diffuse_coef;
uniform vec4 specular_coef;
uniform float mat_shininess;
uniform vec4 lightAmbient;
uniform vec4 lightDiffuse;
uniform vec4 lightSpecular;
uniform vec4 lightPosition;
uniform vec4 eyePosition;

varying vec4 pcolor; 
varying vec3 v_normal;
varying vec4 pos_in_eye; 

void main() {
	//Prepare vectors
	vec3 lightVector = normalize(vec3(lightPosition - pos_in_eye));
	vec3 eyeVector = normalize(vec3(eyePosition - pos_in_eye));
	vec3 normal = normalize(vec3(v_normal));

	//Ambient
	vec4 ambient = ambient_coef * lightAmbient;

	//Diffuse
	float ndotl = max(dot(normal,lightVector),0);
	vec4 diffuse = diffuse_coef * lightDiffuse * ndotl;

	//Specular
	vec3 reflectVector = normalize(reflect(lightVector, normal));
	float rdote = max(dot(reflectVector, eyeVector),0);
	vec4 specular = specular_coef * lightSpecular * pow(rdote, mat_shininess);

	//Set the color
	vec4 finalColor = (ambient + diffuse + specular) * pcolor;

    gl_FragColor = finalColor; 
	//gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);

 } 